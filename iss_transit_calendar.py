#!/usr/bin/env python3
"""
iss_transit_calendar.py - Keep a Google Calendar in sync with upcoming ISS
transits across the Sun and the Moon.

Transits are computed locally with https://github.com/beschne/transit-finder
(transit_finder.py must sit next to this script). Meant to run once a day.

Each run:
  1. computes the transits for the next N days (--days, default 30) using a
     fresh TLE from Celestrak. Beyond ~1 week forecasts are only indicative:
     they are flagged with "~" and get refined as the date approaches;
  2. renders a visibility corridor map (PNG) for each transit;
  3. creates / updates / deletes Google Calendar events without duplicates:
       - deterministic event ID at creation: isstransit<body><YYYYMMDD>
       - private extended properties: source, body, transit time (UTC),
         minimum separation, Drive file ID;
       - existing predictions are matched by body + time, because predictions
         drift a little from one day to the next;
  4. uploads the map to Google Drive (folder "ISS Transits") and attaches it
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

SOURCE_TAG = "iss-transit-finder"
DRIVE_FOLDER = "ISS Transits"
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
    iss: object
    tle_age_days: float
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
#     the Sun/Moon as seen from the CENTRE. Because of parallax (the ISS is
#     400-2000 km away), an observer 50 km off sees the ISS shifted by several
#     degrees -> the threshold now depends on the search radius.
#  3. A 10 s time step lets the ISS travel up to ~10 degrees between samples
#     near the zenith, so a fast pass could fall between two samples
#     -> the step is reduced (--coarse-step).
ISS_MIN_RANGE_KM = 408.0  # observer-to-ISS distance when the ISS is at the zenith


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


def configure_search(args) -> None:
    tf.COARSE_STEP_S = args.coarse_step
    auto = math.degrees(math.atan(args.radius / ISS_MIN_RANGE_KM)) + 0.55 * args.coarse_step
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
    try:
        name, l1, l2 = tf.fetch_tle()  # calls sys.exit on failure
    except SystemExit as e:
        raise RuntimeError(f"Could not fetch the TLE: {e}") from None
    iss = EarthSatellite(l1, l2, name, ts)
    age = (datetime.now(timezone.utc) - iss.epoch.utc_datetime()).total_seconds() / 86400
    log.info("TLE %s, age %.1f d", name.strip(), age)
    if age > 3:
        log.warning("TLE is %.1f days old: reduced accuracy", age)
    configure_search(args)
    t_start = time.monotonic()
    transits = tf.find_all_transits(args.lat, args.lon, args.radius, args.days, ts, eph, iss)
    log.info("Transit computation took %.0f s", time.monotonic() - t_start)
    transits = [t for t in transits if t["body_alt"] >= args.min_alt]
    return transits, ts, eph, iss, age


# -- Google Drive --------------------------------------------------------------
def ensure_folder(drive) -> str:
    q = (f"name='{DRIVE_FOLDER}' and mimeType='application/vnd.google-apps.folder' "
         "and trashed=false")
    res = drive.files().list(q=q, fields="files(id)", spaces="drive").execute()
    if res.get("files"):
        return res["files"][0]["id"]
    f = drive.files().create(
        body={"name": DRIVE_FOLDER, "mimeType": "application/vnd.google-apps.folder"},
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


def make_and_upload(ctx: Ctx, tr: dict) -> tuple[dict | None, str | None]:
    """Render the corridor map and upload it. Returns (attachment, Drive file ID)."""
    if ctx.args.no_maps:
        return None, None
    try:
        corridor = tf.compute_corridor(tr, ctx.ts, ctx.eph, ctx.iss)
        path = tf.generate_map(tr, corridor, ctx.tmpdir, ctx.args.lat, ctx.args.lon,
                               ctx.args.name)
        if not path:
            log.warning("Map not generated (missing libraries?)")
            return None, None
        local = tr["time"].astimezone(ctx.tz)
        title = f"ISS {BODY_NAME[tr['body']]} {local:%Y-%m-%d %H:%M:%S}.png"
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
    return f"isstransit{tr['body']}{tr['time']:%Y%m%d}"


def build_event(ctx: Ctx, tr: dict, att: dict | None, drive_id: str | None) -> dict:
    a = ctx.args
    body = tr["body"]
    t_utc = tr["time"]
    local = t_utc.astimezone(ctx.tz)
    days_ahead = (t_utc - datetime.now(timezone.utc)).total_seconds() / 86400
    far = days_ahead > FAR_FORECAST_DAYS
    start = local - timedelta(minutes=a.lead_min)
    end = local + timedelta(minutes=2)

    desc = [
        f"<b>ISS transit across {BODY_WITH_ARTICLE[body]}</b>",
        f"Exact time: <b>{local:%H:%M:%S}</b> ({ctx.tz.key}) - {t_utc:%H:%M:%S} UTC",
        f"Duration: {tr['duration_s']:.1f} s",
        f"Minimum separation: {tr['min_sep_am']:.0f}\" (disc radius: {tr['r_am']:.0f}\")",
        f"Chord coverage: {tr['coverage']:.0f} % of the disc diameter",
        f"{BODY_NAME[body]}: altitude {tr['body_alt']:.1f}°, azimuth {tr['body_az']:.0f}° "
        f"({tf.azimuth_label(tr['body_az'])})",
        f"ISS: altitude {tr['iss_alt']:.1f}°, azimuth {tr['iss_az']:.0f}° "
        f"({tf.azimuth_label(tr['iss_az'])})",
        f"Best observing spot: "
        f"<a href=\"https://www.google.com/maps?q={tr['obs_lat']:.5f},{tr['obs_lon']:.5f}\">"
        f"{tr['obs_lat']:.4f}, {tr['obs_lon']:.4f}</a> "
        f"(~{tr['obs_dist_km']:.0f} km from the search centre)",
        f"TLE age: {ctx.tle_age_days:.1f} d - prediction made "
        f"{datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC",
    ]
    if far:
        desc.append(
            f"<b>Indicative forecast (D+{days_ahead:.0f})</b>: the uncertainty on the ISS "
            f"position grows by roughly 1-2 km per day of forecast, while the corridor is "
            f"only 5-15 km wide. Time and location will be refined every day.")
    if att:
        desc.append(f"Corridor map: <a href=\"{att['fileUrl']}\">open image</a> "
                    "(also attached)")
    if body == "sun":
        desc.append("<b>⚠ A certified solar filter is mandatory.</b>")

    private = {
        "source": SOURCE_TAG,
        "body": body,
        "transitUtc": t_utc.isoformat(),
        "minSepAm": f"{tr['min_sep_am']:.1f}",
    }
    if drive_id:
        private["driveFileId"] = drive_id

    ev = {
        "summary": (f"{BODY_EMOJI[body]} ISS transit - {BODY_NAME[body]} "
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


def tracked_time(ev: dict) -> datetime | None:
    try:
        return datetime.fromisoformat(priv(ev)["transitUtc"])
    except (KeyError, ValueError):
        return None


def match_existing(existing: list[dict], tr: dict, used: set[str]) -> dict | None:
    best, best_dt = None, MATCH_TOL_S + 1
    for ev in existing:
        if ev["id"] in used or priv(ev).get("body") != tr["body"]:
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
                log.info("UPDATED  %s %s", tr["body"], tr["time"])
            else:
                unchanged += 1
        elif create_event(ctx, tr):
            created += 1
            log.info("CREATED  %s %s", tr["body"], tr["time"])

    # Delete predictions that disappeared (only inside the reliable window)
    if ctx.tle_age_days <= MAX_TLE_AGE_FOR_DELETE:
        for ev in existing:
            t_old = tracked_time(ev)
            if (ev["id"] not in used and t_old
                    and now + timedelta(minutes=10) < t_old < horizon - timedelta(hours=2)):
                delete_event(ctx, ev)
                deleted += 1
                log.info("DELETED  %s %s (no longer predicted)", priv(ev).get("body"), t_old)
    else:
        log.warning("TLE too old: no deletion performed")

    log.info("Summary: %d created, %d updated, %d unchanged, %d deleted",
             created, updated, unchanged, deleted)


# -- Entry point ---------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="ISS transits -> Google Calendar")
    p.add_argument("--lat", type=float, required=True, help="latitude (degrees, north positive)")
    p.add_argument("--lon", type=float, required=True, help="longitude (degrees, east positive)")
    p.add_argument("--radius", type=float, default=50.0, help="search radius (km)")
    p.add_argument("--days", type=int, default=30,
                   help="forecast horizon in days (default 30; beyond ~7 days "
                        "forecasts are only indicative)")
    p.add_argument("--name", default="", help="location name (shown on the map)")
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
    return p.parse_args()


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
        transits, ts, eph, iss, age = compute_transits(args, workdir)
    except Exception as e:
        # Computation failed: never touch the calendar
        log.error("Computation failed, calendar left untouched: %s", e)
        return 1
    log.info("%d transit(s) predicted over %d days", len(transits), args.days)

    if args.dry_run:
        for tr in transits:
            print(f"{tr['time']:%Y-%m-%d %H:%M:%S} UTC  {tr['body']:5s}  "
                  f"duration {tr['duration_s']:.1f}s  sep. {tr['min_sep_am']:.0f}\"  "
                  f"alt {tr['body_alt']:.0f}°  ({tr['obs_lat']:.4f}, {tr['obs_lon']:.4f})")
        return 0

    try:
        creds = get_credentials(workdir)
        cal = build("calendar", "v3", credentials=creds, cache_discovery=False)
        drive = build("drive", "v3", credentials=creds, cache_discovery=False)
        folder_id = None if args.no_maps else ensure_folder(drive)
        with tempfile.TemporaryDirectory() as tmp:
            ctx = Ctx(args, ZoneInfo(args.tz), ts, eph, iss, age, cal, drive, folder_id, tmp)
            sync(ctx, transits, now, horizon)
    except Exception as e:
        log.exception("Synchronisation failed: %s", e)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
