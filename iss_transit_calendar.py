#!/usr/bin/env python3
"""
iss_transit_calendar.py - Keep a Google Calendar in sync with upcoming ISS
(and, optionally, Tiangong) transits across the Sun and the Moon.

Transits are computed locally with https://github.com/beschne/transit-finder
(transit_finder.py must sit next to this script). Meant to run once a day.

Each run:
  1. computes the transits for the next N days (--days, default 30) using a
     fresh TLE from Celestrak. Beyond ~1 week forecasts are only indicative:
     they are flagged with "~" and get refined as the date approaches;
  2. renders a visibility corridor map (PNG) for each transit;
  3. creates / updates / deletes Google Calendar events without duplicates:
       - deterministic event ID at creation: <iss|css>transit<body><YYYYMMDD>
       - private extended properties: source, satellite, body, transit time
         (UTC), minimum separation, Drive file ID;
       - existing predictions are matched by body + time, because predictions
         drift a little from one day to the next;
  4. uploads the map to Google Drive (folder "Satellite Transits") and attaches it
     to the event (Drive attachment, plus a link in the description).

See README.md for installation, Google Cloud setup and scheduling.
"""
from __future__ import annotations

import argparse
import fcntl
import logging
import math
import os
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

from google.auth.exceptions import RefreshError
from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
from googleapiclient.http import MediaFileUpload
import numpy as np
from skyfield.api import EarthSatellite, Loader, wgs84

import transit_finder as tf  # https://github.com/beschne/transit-finder

SOURCE_TAG = "iss-transit-finder"  # historical name, kept so existing events stay tracked
DEFAULT_DRIVE_FOLDER = "Satellite Transits"

# Supported satellites. ISS is the default; Tiangong (CSS) is opt-in via --satellites.
SATELLITES = {
    "iss": {"label": "ISS", "norad": 25544, "id_prefix": "isstransit"},
    "css": {"label": "Tiangong", "norad": 48274, "id_prefix": "csstransit"},  # Tianhe core module
}
SATELLITE_ALIASES = {"tiangong": "css"}
SCOPES = [
    "https://www.googleapis.com/auth/calendar.events",
    "https://www.googleapis.com/auth/drive.file",  # only files created by this script
]
MATCH_TOL_S = 10800         # tolerance (3 h) when matching an existing prediction:
                            # long-range forecasts can drift by several minutes
FAR_FORECAST_DAYS = 7       # beyond this, the forecast is flagged as indicative
MAX_TLE_AGE_FOR_DELETE = 5  # days; with an older TLE nothing is ever deleted

BODY_NAME = {"sun": "Sun", "moon": "Moon"}
BODY_WITH_ARTICLE = {"sun": "the Sun", "moon": "the Moon"}
BODY_EMOJI = {"sun": "☀️", "moon": "🌙"}
BODY_COLOR = {"sun": "5", "moon": "8"}  # Google Calendar colors: yellow / graphite

log = logging.getLogger("iss-cal")


# -- Context -------------------------------------------------------------------
@dataclass
class Ctx:
    args: argparse.Namespace
    tz: ZoneInfo
    ts: object
    eph: object
    sats: dict          # satellite key -> skyfield EarthSatellite
    ages: dict          # satellite key -> TLE age in days
    cal: object
    drive: object
    folder_id: str | None
    tmpdir: str


# -- Google authentication -----------------------------------------------------
def get_credentials(workdir: Path) -> Credentials:
    token_file = workdir / "token.json"
    cred_file = workdir / "credentials.json"
    creds = None
    if token_file.exists():
        creds = Credentials.from_authorized_user_file(str(token_file), SCOPES)
    if creds and creds.expired and creds.refresh_token:
        try:
            creds.refresh(Request())
        except RefreshError:
            log.warning("Refresh token rejected, a new authorization is required.")
            creds = None
    if not creds or not creds.valid:
        if not sys.stdin.isatty():
            raise RuntimeError(
                "Google authorization missing or expired: run the script once "
                "manually to regenerate token.json."
            )
        flow = InstalledAppFlow.from_client_secrets_file(str(cred_file), SCOPES)
        creds = flow.run_local_server(port=0)
    token_file.write_text(creds.to_json())
    try:
        token_file.chmod(0o600)
    except OSError:
        pass
    return creds


# -- Search: observer grid + radius-aware pre-scan -----------------------------
# Three limitations of the original transit_finder.py, worked around here
# without modifying that file:
#  1. The observers it tests are only the centre plus 8 points on the circle's
#     circumference. The interior of the disc is not sampled, although the
#     visibility corridor is only a few km wide
#     -> use a grid of points spaced --grid-step-km apart.
#  2. Its pre-scan keeps a candidate only if the ISS passes within 2 degrees of
#     the Sun/Moon as seen from the CENTRE. Because of parallax (the satellite is
#     400-2000 km away), an observer 50 km off sees it shifted by several
#     degrees -> the threshold now depends on the search radius.
#  3. A 10 s time step lets the ISS travel up to ~10 degrees between samples
#     near the zenith, so a fast pass could fall between two samples
#     -> the step is reduced (--coarse-step).
MU_EARTH_KM3_S2 = 398600.4418
EARTH_RADIUS_KM = 6371.0


def min_range_km(sat) -> float:
    """Approximate observer-to-satellite distance when the satellite is at the
    zenith, i.e. its perigee altitude (minus a 10 km margin), from the TLE."""
    try:
        n = sat.model.no_kozai / 60.0  # rad/min -> rad/s
        a = (MU_EARTH_KM3_S2 / n ** 2) ** (1.0 / 3.0)
        return max(150.0, a * (1.0 - sat.model.ecco) - EARTH_RADIUS_KM - 10.0)
    except Exception:
        return 350.0  # conservative fallback for a low-Earth-orbit station


def make_observer_grid(step_km: float):
    def grid(lat: float, lon: float, radius_km: float) -> list[tuple[float, float]]:
        pts = [(lat, lon)]
        if radius_km <= 0.1:
            return pts
        n_rings = max(1, round(radius_km / step_km))
        for k in range(1, n_rings + 1):
            r = radius_km * k / n_rings
            n = max(6, math.ceil(2 * math.pi * r / step_km))
            offset = (k % 2) * 180.0 / n  # stagger every other ring
            for j in range(n):
                pts.append(tf.offset_point(lat, lon, 360.0 * j / n + offset, r))
        return pts
    return grid


def coarse_scan_v2(lat, lon, days, ts, eph, iss):
    """Like tf.coarse_scan, but processed in chunks (memory) and returning, for
    each close approach, the instant of MINIMUM separation rather than the
    first sample (with a wide threshold the first sample would fall too early
    for the window of the fine scan)."""
    earth = eph["earth"]
    observer = wgs84.latlon(lat, lon)
    step, thr = tf.COARSE_STEP_S, tf.APPROACH_DEG
    t0_tt = ts.from_datetime(datetime.now(timezone.utc)).tt
    n = int(days * 86400 / step) + 1
    hits: dict[str, list[tuple[float, float]]] = {"sun": [], "moon": []}
    chunk = 100_000
    for i0 in range(0, n, chunk):
        tt = t0_tt + np.arange(i0, min(n, i0 + chunk)) * (step / 86400.0)
        t_arr = ts.tt_jd(tt)
        iss_alt, iss_az, _ = (iss - observer).at(t_arr).altaz()
        earth_at = (earth + observer).at(t_arr)
        for body in ("sun", "moon"):
            b_alt, b_az, _ = earth_at.observe(eph[body]).apparent().altaz()
            sep = tf.sep_vec(iss_alt.degrees, iss_az.degrees, b_alt.degrees, b_az.degrees)
            ok = (iss_alt.degrees > 0) & (b_alt.degrees > -5) & (sep < thr)
            hits[body] += list(zip(tt[ok].tolist(), sep[ok].tolist()))

    clustered: list[tuple[float, str]] = []
    for body, lst in hits.items():
        group: list[tuple[float, float]] = []
        for t, sp in lst:
            if group and (t - group[-1][0]) * 86400 > step * 2.5:
                clustered.append((min(group, key=lambda x: x[1])[0], body))
                group = []
            group.append((t, sp))
        if group:
            clustered.append((min(group, key=lambda x: x[1])[0], body))
    clustered.sort(key=lambda x: x[0])
    log.info("Pre-scan: %d close approach(es) within %.1f deg (step %.0f s)",
             len(clustered), thr, step)
    return clustered, t0_tt


def configure_search(args, sat) -> None:
    tf.COARSE_STEP_S = args.coarse_step
    auto = math.degrees(math.atan(args.radius / min_range_km(sat))) + 0.55 * args.coarse_step
    tf.APPROACH_DEG = args.approach_deg or min(25.0, max(2.0, auto))
    tf.FINE_WINDOW_S = 40  # the candidate is now the instant of minimum separation
    tf.observer_ring = make_observer_grid(args.grid_step_km)
    tf.coarse_scan = coarse_scan_v2
    n_obs = len(tf.observer_ring(args.lat, args.lon, args.radius))
    log.info("Search: %d observer points (%.1f km grid), pre-scan threshold %.1f deg",
             n_obs, args.grid_step_km, tf.APPROACH_DEG)


# -- Transit computation -------------------------------------------------------
def compute_transits(args, workdir: Path):
    loader = Loader(str(workdir / "skyfield-data"))
    ts = loader.timescale()
    eph = loader("de421.bsp")
    transits: list[dict] = []
    sats: dict = {}
    ages: dict = {}
    for key in args.sat_keys:
        label = SATELLITES[key]["label"]
        try:
            name, l1, l2 = tf.fetch_tle(SATELLITES[key]["norad"])  # calls sys.exit on failure
        except SystemExit as e:
            raise RuntimeError(f"Could not fetch the {label} TLE: {e}") from None
        sat = EarthSatellite(l1, l2, name, ts)
        age = (datetime.now(timezone.utc) - sat.epoch.utc_datetime()).total_seconds() / 86400
        sats[key], ages[key] = sat, age
        log.info("%s: TLE %s, age %.1f d", label, name.strip(), age)
        if age > 3:
            log.warning("%s: TLE is %.1f days old: reduced accuracy", label, age)
        configure_search(args, sat)
        t_start = time.monotonic()
        found = tf.find_all_transits(args.lat, args.lon, args.radius, args.days, ts, eph, sat)
        log.info("%s: transit computation took %.0f s", label, time.monotonic() - t_start)
        for tr in found:
            if tr["body_alt"] >= args.min_alt:
                tr["satellite"] = key
                transits.append(tr)
    transits.sort(key=lambda t: t["time"])
    return transits, ts, eph, sats, ages


# -- Google Drive --------------------------------------------------------------
def ensure_folder(drive, name: str) -> str:
    q = (f"name='{name}' and mimeType='application/vnd.google-apps.folder' "
         "and trashed=false")
    res = drive.files().list(q=q, fields="files(id)", spaces="drive").execute()
    if res.get("files"):
        return res["files"][0]["id"]
    f = drive.files().create(
        body={"name": name, "mimeType": "application/vnd.google-apps.folder"},
        fields="id",
    ).execute()
    return f["id"]


def upload_map(ctx: Ctx, path: str, title: str) -> dict:
    media = MediaFileUpload(path, mimetype="image/png")
    meta = {"name": title, "parents": [ctx.folder_id]}
    f = ctx.drive.files().create(
        body=meta, media_body=media, fields="id,webViewLink"
    ).execute()
    if ctx.args.share_images:
        ctx.drive.permissions().create(
            fileId=f["id"], body={"type": "anyone", "role": "reader"}
        ).execute()
    return f


def delete_drive_file(ctx: Ctx, file_id: str | None) -> None:
    if not file_id:
        return
    try:
        ctx.drive.files().delete(fileId=file_id).execute()
    except HttpError as e:
        if e.resp.status != 404:
            log.warning("Could not delete Drive file %s: %s", file_id, e)


def render_map(ctx: Ctx, tr: dict, corridor: dict) -> str | None:
    """Call tf.generate_map. Its map title is hard-coded to "ISS ...", so for
    other satellites the title is rewritten on the fly."""
    key = tr["satellite"]
    label = SATELLITES[key]["label"]
    out_dir = os.path.join(ctx.tmpdir, key)
    orig_set_title = axes_mod = None
    if label != "ISS":
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.axes as axes_mod
            orig_set_title = axes_mod.Axes.set_title

            def set_title(self, text, *a, **kw):
                if isinstance(text, str) and text.startswith("ISS "):
                    text = label + text[3:]
                return orig_set_title(self, text, *a, **kw)

            axes_mod.Axes.set_title = set_title
        except ImportError:
            orig_set_title = None  # tf.generate_map will report the missing library
    try:
        return tf.generate_map(tr, corridor, out_dir, ctx.args.lat, ctx.args.lon, ctx.args.name)
    finally:
        if orig_set_title is not None:
            axes_mod.Axes.set_title = orig_set_title


def make_and_upload(ctx: Ctx, tr: dict) -> tuple[dict | None, str | None]:
    """Render the corridor map and upload it. Returns (attachment, Drive file ID)."""
    if ctx.args.no_maps:
        return None, None
    try:
        sat = ctx.sats[tr["satellite"]]
        corridor = tf.compute_corridor(tr, ctx.ts, ctx.eph, sat)
        path = render_map(ctx, tr, corridor)
        if not path:
            log.warning("Map not generated (missing libraries?)")
            return None, None
        local = tr["time"].astimezone(ctx.tz)
        label = SATELLITES[tr["satellite"]]["label"]
        title = f"{label} {BODY_NAME[tr['body']]} {local:%Y-%m-%d %H:%M:%S}.png"
        f = upload_map(ctx, path, title)
        att = {"fileUrl": f["webViewLink"], "title": title,
               "mimeType": "image/png", "fileId": f["id"]}
        return att, f["id"]
    except Exception as e:  # the map must never block the sync
        log.warning("Map/Drive failure, event will have no image: %s", e)
        return None, None


# -- Event construction --------------------------------------------------------
def event_id(tr: dict) -> str:
    # Google only allows lowercase a-v and 0-9 in event IDs (base32hex)
    prefix = SATELLITES[tr["satellite"]]["id_prefix"]
    return f"{prefix}{tr['body']}{tr['time']:%Y%m%d}"


def build_event(ctx: Ctx, tr: dict, att: dict | None, drive_id: str | None) -> dict:
    a = ctx.args
    body = tr["body"]
    sat_key = tr["satellite"]
    label = SATELLITES[sat_key]["label"]
    t_utc = tr["time"]
    local = t_utc.astimezone(ctx.tz)
    days_ahead = (t_utc - datetime.now(timezone.utc)).total_seconds() / 86400
    far = days_ahead > FAR_FORECAST_DAYS
    start = local - timedelta(minutes=a.lead_min)
    end = local + timedelta(minutes=2)

    desc = [
        f"<b>{label} transit across {BODY_WITH_ARTICLE[body]}</b>",
        f"Exact time: <b>{local:%H:%M:%S}</b> ({ctx.tz.key}) - {t_utc:%H:%M:%S} UTC",
        f"Duration: {tr['duration_s']:.1f} s",
        f"Minimum separation: {tr['min_sep_am']:.0f}\" (disc radius: {tr['r_am']:.0f}\")",
        f"Chord coverage: {tr['coverage']:.0f} % of the disc diameter",
        f"{BODY_NAME[body]}: altitude {tr['body_alt']:.1f}°, azimuth {tr['body_az']:.0f}° "
        f"({tf.azimuth_label(tr['body_az'])})",
        f"{label}: altitude {tr['iss_alt']:.1f}°, azimuth {tr['iss_az']:.0f}° "
        f"({tf.azimuth_label(tr['iss_az'])})",
        f"Best observing spot: "
        f"<a href=\"https://www.google.com/maps?q={tr['obs_lat']:.5f},{tr['obs_lon']:.5f}\">"
        f"{tr['obs_lat']:.4f}, {tr['obs_lon']:.4f}</a> "
        f"(~{tr['obs_dist_km']:.0f} km from the search centre)",
        f"TLE age: {ctx.ages[sat_key]:.1f} d - prediction made "
        f"{datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC",
    ]
    if far:
        desc.append(
            f"<b>Indicative forecast (D+{days_ahead:.0f})</b>: the uncertainty on the "
            f"satellite position grows with forecast time (about 1-2 km per day for the ISS), "
            f"while the corridor is only 5-15 km wide. Time and location will be refined "
            f"every day.")
    if att:
        desc.append(f"Corridor map: <a href=\"{att['fileUrl']}\">open image</a> "
                    "(also attached)")
    if body == "sun":
        desc.append("<b>⚠ A certified solar filter is mandatory.</b>")

    private = {
        "source": SOURCE_TAG,
        "satellite": sat_key,
        "body": body,
        "transitUtc": t_utc.isoformat(),
        "minSepAm": f"{tr['min_sep_am']:.1f}",
    }
    if drive_id:
        private["driveFileId"] = drive_id

    ev = {
        "summary": (f"{BODY_EMOJI[body]} {label} transit - {BODY_NAME[body]} "
                    f"{'~' if far else ''}{local:%H:%M:%S}"),
        "description": "<br>".join(desc),
        "location": f"{tr['obs_lat']:.5f}, {tr['obs_lon']:.5f}",
        "start": {"dateTime": start.isoformat(), "timeZone": ctx.tz.key},
        "end": {"dateTime": end.isoformat(), "timeZone": ctx.tz.key},
        "colorId": BODY_COLOR[body],
        "reminders": {"useDefault": False, "overrides": [
            {"method": "popup", "minutes": 60},
        ]},
        "extendedProperties": {"private": private},
    }
    if att:
        ev["attachments"] = [att]
    return ev


# -- Synchronisation -----------------------------------------------------------
def list_tracked(ctx: Ctx, t_min: datetime, t_max: datetime) -> list[dict]:
    events, token = [], None
    while True:
        r = ctx.cal.events().list(
            calendarId=ctx.args.calendar_id,
            timeMin=t_min.isoformat(), timeMax=t_max.isoformat(),
            privateExtendedProperty=f"source={SOURCE_TAG}",
            singleEvents=True, showDeleted=False, maxResults=250, pageToken=token,
        ).execute()
        events += r.get("items", [])
        token = r.get("nextPageToken")
        if not token:
            return events


def priv(ev: dict) -> dict:
    return ev.get("extendedProperties", {}).get("private", {})


def sat_of(ev: dict) -> str:
    # Events created before multi-satellite support carry no "satellite" property: ISS
    return priv(ev).get("satellite", "iss")


def tracked_time(ev: dict) -> datetime | None:
    try:
        return datetime.fromisoformat(priv(ev)["transitUtc"])
    except (KeyError, ValueError):
        return None


def match_existing(existing: list[dict], tr: dict, used: set[str]) -> dict | None:
    best, best_dt = None, MATCH_TOL_S + 1
    for ev in existing:
        if (ev["id"] in used or priv(ev).get("body") != tr["body"]
                or sat_of(ev) != tr["satellite"]):
            continue
        t_old = tracked_time(ev)
        if t_old is None:
            continue
        dt = abs((t_old - tr["time"]).total_seconds())
        if dt < best_dt:
            best, best_dt = ev, dt
    return best if best_dt <= MATCH_TOL_S else None


def needs_update(ctx: Ctx, ev: dict, tr: dict) -> bool:
    p = priv(ev)
    t_old = tracked_time(ev)
    if t_old is None:
        return True
    if abs((t_old - tr["time"]).total_seconds()) > ctx.args.time_tol:
        return True
    if abs(float(p.get("minSepAm", 0)) - tr["min_sep_am"]) > ctx.args.sep_tol:
        return True
    if not ctx.args.no_maps and not p.get("driveFileId"):
        return True  # retry the map
    return False


def create_event(ctx: Ctx, tr: dict) -> bool:
    att, drive_id = make_and_upload(ctx, tr)
    body = build_event(ctx, tr, att, drive_id)
    body["id"] = event_id(tr)
    try:
        ctx.cal.events().insert(
            calendarId=ctx.args.calendar_id, body=body, supportsAttachments=True
        ).execute()
    except HttpError as e:
        delete_drive_file(ctx, drive_id)
        if e.resp.status == 409:
            log.info("ID %s already in use (event deleted by hand?): skipped", body["id"])
            return False
        raise
    return True


def update_event(ctx: Ctx, ev: dict, tr: dict) -> None:
    old_drive = priv(ev).get("driveFileId")
    att, drive_id = make_and_upload(ctx, tr)
    body = build_event(ctx, tr, att, drive_id)
    if not att and old_drive:  # map not regenerated: keep the previous one
        body["extendedProperties"]["private"]["driveFileId"] = old_drive
    ctx.cal.events().patch(
        calendarId=ctx.args.calendar_id, eventId=ev["id"], body=body,
        supportsAttachments=True,
    ).execute()
    if drive_id and old_drive:
        delete_drive_file(ctx, old_drive)


def delete_event(ctx: Ctx, ev: dict) -> None:
    ctx.cal.events().delete(calendarId=ctx.args.calendar_id, eventId=ev["id"]).execute()
    delete_drive_file(ctx, priv(ev).get("driveFileId"))


def sync(ctx: Ctx, transits: list[dict], now: datetime, horizon: datetime) -> None:
    existing = list_tracked(ctx, now, horizon + timedelta(hours=1))
    log.info("%d tracked event(s) already in the calendar", len(existing))
    used: set[str] = set()
    created = updated = unchanged = deleted = 0

    for tr in transits:
        ev = match_existing(existing, tr, used)
        if ev:
            used.add(ev["id"])
            if needs_update(ctx, ev, tr):
                update_event(ctx, ev, tr)
                updated += 1
                log.info("UPDATED  %s %s %s", tr["satellite"], tr["body"], tr["time"])
            else:
                unchanged += 1
        elif create_event(ctx, tr):
            created += 1
            log.info("CREATED  %s %s %s", tr["satellite"], tr["body"], tr["time"])

    # Delete predictions that disappeared: only inside the reliable window, only for
    # satellites computed in this run, and never when that satellite's TLE is stale
    warned: set[str] = set()
    for ev in existing:
        sat_key = sat_of(ev)
        t_old = tracked_time(ev)
        if ev["id"] in used or not t_old or sat_key not in ctx.sats:
            continue
        if ctx.ages[sat_key] > MAX_TLE_AGE_FOR_DELETE:
            if sat_key not in warned:
                log.warning("%s TLE too old: no deletion performed", SATELLITES[sat_key]["label"])
                warned.add(sat_key)
            continue
        if now + timedelta(minutes=10) < t_old < horizon - timedelta(hours=2):
            delete_event(ctx, ev)
            deleted += 1
            log.info("DELETED  %s %s %s (no longer predicted)",
                     sat_key, priv(ev).get("body"), t_old)

    log.info("Summary: %d created, %d updated, %d unchanged, %d deleted",
             created, updated, unchanged, deleted)


# -- Entry point ---------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="ISS / Tiangong transits -> Google Calendar")
    p.add_argument("--lat", type=float, required=True, help="latitude (degrees, north positive)")
    p.add_argument("--lon", type=float, required=True, help="longitude (degrees, east positive)")
    p.add_argument("--radius", type=float, default=50.0, help="search radius (km)")
    p.add_argument("--days", type=int, default=30,
                   help="forecast horizon in days (default 30; beyond ~7 days "
                        "forecasts are only indicative)")
    p.add_argument("--name", default="", help="location name (shown on the map)")
    p.add_argument("--satellites", default="iss",
                   help="comma-separated satellites to track: iss, css (Tiangong). "
                        "Default: iss. Example: --satellites iss,css")
    p.add_argument("--drive-folder", default=DEFAULT_DRIVE_FOLDER,
                   help="Google Drive folder that receives the maps")
    p.add_argument("--grid-step-km", type=float, default=3.0,
                   help="spacing of the observer grid inside the radius (default 3 km, "
                        "about the minimum corridor width; smaller = more reliable "
                        "but slower)")
    p.add_argument("--coarse-step", type=float, default=2.0,
                   help="pre-scan time step in seconds (default 2)")
    p.add_argument("--approach-deg", type=float, default=0.0,
                   help="pre-scan threshold in degrees (0 = automatic, from the radius)")
    p.add_argument("--calendar-id", default="primary", help="target Google Calendar ID")
    p.add_argument("--tz", default="UTC",
                   help="time zone used in event titles/descriptions, e.g. Europe/Paris")
    p.add_argument("--min-alt", type=float, default=5.0,
                   help="minimum altitude (deg) of the Sun/Moon to keep a transit")
    p.add_argument("--lead-min", type=int, default=5,
                   help="the event starts N minutes before the transit")
    p.add_argument("--time-tol", type=float, default=3.0,
                   help="time shift (s) above which an event is updated")
    p.add_argument("--sep-tol", type=float, default=3.0,
                   help="separation change (arcsec) above which an event is updated")
    p.add_argument("--share-images", action="store_true",
                   help="make the Drive maps readable by anyone with the link")
    p.add_argument("--no-maps", action="store_true", help="skip map generation and upload")
    p.add_argument("--dry-run", action="store_true",
                   help="compute and print transits without touching Google")
    p.add_argument("--workdir", default=str(Path(__file__).resolve().parent),
                   help="folder for credentials.json / token.json / Skyfield data")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()
    args.sat_keys = parse_satellites(p, args.satellites)
    return args


def parse_satellites(parser: argparse.ArgumentParser, value: str) -> list[str]:
    keys: list[str] = []
    for item in value.lower().split(","):
        item = SATELLITE_ALIASES.get(item.strip(), item.strip())
        if not item:
            continue
        if item not in SATELLITES:
            parser.error(f"unknown satellite '{item}' (choose from: {', '.join(SATELLITES)})")
        if item not in keys:
            keys.append(item)
    if not keys:
        parser.error("--satellites must list at least one satellite")
    return keys


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    fmt = "%Y-%m-%d %H:%M:%S %Z"
    t_begin = time.monotonic()
    log.info("=== RUN STARTED: %s ===", datetime.now().astimezone().strftime(fmt))
    rc = 1  # stays 1 if an unexpected exception escapes
    try:
        rc = run(args)
    finally:
        d = int(time.monotonic() - t_begin)
        log.info("=== RUN FINISHED: %s - duration %dh%02dm%02ds - exit code %d ===",
                 datetime.now().astimezone().strftime(fmt), d // 3600, d % 3600 // 60, d % 60, rc)
    return rc


def run(args: argparse.Namespace) -> int:
    workdir = Path(args.workdir)
    os.chdir(workdir)  # transit_finder writes its maps to a relative path

    # Lock: a fine-grained search can take a long time, avoid overlapping runs
    lock_file = open(workdir / ".lock", "w")
    try:
        fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        log.error("Previous run still in progress: aborting")
        return 0
    now = datetime.now(timezone.utc)
    horizon = now + timedelta(days=args.days)
    if args.days > 60:
        log.warning("Horizon of %d days: forecasts are very uncertain beyond a few weeks "
                    "(TLE drift, ISS reboost maneuvers)", args.days)

    try:
        transits, ts, eph, sats, ages = compute_transits(args, workdir)
    except Exception as e:
        # Computation failed: never touch the calendar
        log.error("Computation failed, calendar left untouched: %s", e)
        return 1
    log.info("%d transit(s) predicted over %d days", len(transits), args.days)

    if args.dry_run:
        for tr in transits:
            print(f"{tr['time']:%Y-%m-%d %H:%M:%S} UTC  {tr['satellite']:3s}  {tr['body']:5s}  "
                  f"duration {tr['duration_s']:.1f}s  sep. {tr['min_sep_am']:.0f}\"  "
                  f"alt {tr['body_alt']:.0f}°  ({tr['obs_lat']:.4f}, {tr['obs_lon']:.4f})")
        return 0

    try:
        creds = get_credentials(workdir)
        cal = build("calendar", "v3", credentials=creds, cache_discovery=False)
        drive = build("drive", "v3", credentials=creds, cache_discovery=False)
        folder_id = None if args.no_maps else ensure_folder(drive, args.drive_folder)
        with tempfile.TemporaryDirectory() as tmp:
            ctx = Ctx(args, ZoneInfo(args.tz), ts, eph, sats, ages, cal, drive, folder_id, tmp)
            sync(ctx, transits, now, horizon)
    except Exception as e:
        log.exception("Synchronisation failed: %s", e)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
