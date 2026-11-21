# ISS Transit Calendar

Keeps a Google Calendar in sync with upcoming **ISS transits across the Sun and the Moon** (and, optionally, **Tiangong** transits), seen from a location and search radius of your choice. Designed to run once a day (cron, NAS task scheduler, ...).

Transits are computed locally by [beschne/transit-finder](https://github.com/beschne/transit-finder) (Skyfield + a live TLE from Celestrak). This project adds the calendar integration, a denser search grid, and a more robust pre-scan on top of it.

> **Safety:** never look at the Sun without a certified solar filter. Solar transit events carry a reminder in their description.

## Features

- Computes transits for the next N days (`--days`, default 30).
- Tracks the **ISS** by default; **Tiangong** (CSS, NORAD 48274) is opt-in with `--satellites iss,css`. Each satellite is searched separately and tracked independently in the calendar.
- Creates, updates and deletes Google Calendar events **without duplicates**:
  - deterministic event ID at creation: `<iss|css>transit<sun|moon><YYYYMMDD>`;
  - private extended properties store the satellite, the body, the transit time (UTC), the minimum separation and the Drive file ID;
  - existing events are matched by satellite + body + time (±3 h), because forecasts drift from day to day;
  - events are only patched when the prediction moved by more than a threshold (`--time-tol`, `--sep-tol`);
  - predictions that disappear are deleted, but only inside the reliable part of the horizon, only for satellites computed in that run, and never if that satellite's TLE is stale or the computation failed.
- Renders a **visibility corridor map** for each transit, uploads it to Google Drive (folder `Satellite Transits`, see `--drive-folder`) and attaches it to the event.
- Event description: exact time, duration, minimum separation, chord coverage, Sun/Moon and ISS altitude/azimuth, best observing spot (Google Maps link), TLE age.
- Forecasts beyond 7 days are flagged with `~` in the title and labelled "indicative".
- One popup reminder, 60 minutes before the event.
- **Weather outlook and interest score** for transits in the next 10 days (`--weather-days`): a cloud forecast from [Open-Meteo](https://open-meteo.com/) is turned into a 🟢 / 🟡 / 🔴 rating shown in the event title, the event colour and the description (see below).
- A lock file prevents overlapping runs; start/end timestamps, duration and exit code are logged.

## Requirements

- Linux or macOS (the lock uses `fcntl`), Python **3.10+**.
- A Google account and a Google Cloud project (free) with the Calendar and Drive APIs enabled.
- Internet access (TLE from Celestrak, ephemeris file on first run, map tiles from OpenStreetMap/CARTO).

## Installation

```bash
git clone https://github.com/leirn/iss-transit-calendar.git
cd iss-transit-calendar
python3 -m venv venv && . venv/bin/activate
pip install -r requirements.txt

# Fetch the transit computation module (MIT licensed, by beschne)
curl -O https://raw.githubusercontent.com/beschne/transit-finder/main/transit_finder.py
```

`transit_finder.py` is intentionally **not** included in this repository. This script overrides some of its internals (`observer_ring`, `coarse_scan` and a few constants), so if the upstream file changes significantly, check that everything still works with `--dry-run`.

## Google Cloud setup (once)

1. In the [Google Cloud console](https://console.cloud.google.com/), create a project and enable **Google Calendar API** and **Google Drive API**.
2. Configure the OAuth consent screen (type *External*) and **publish it ("In production")**. In *Testing* mode the refresh token expires after 7 days and the scheduled job will stop working.
3. Create an **OAuth client ID** of type *Desktop app* and download the JSON as `credentials.json` into the project folder.
4. Run the script once **by hand** on a machine with a browser to authorize access and create `token.json`:

   ```bash
   python iss_transit_calendar.py --lat 48.85 --lon 2.35 --radius 50 --tz Europe/Paris
   ```

   You can then copy `credentials.json` and `token.json` to the machine that runs the scheduled job.

Scopes requested: `calendar.events` and `drive.file` (the script can only see Drive files it created itself).

## Usage

```bash
# Test run: print the predicted transits without touching Google (no credentials needed)
python iss_transit_calendar.py --lat 48.85 --lon 2.35 --radius 50 --days 7 --dry-run

# Same, and also render the corridor maps locally
python iss_transit_calendar.py --lat 48.85 --lon 2.35 --radius 50 --dry-run --maps-dir ./test-maps

# Real run
python iss_transit_calendar.py --lat 48.85 --lon 2.35 --radius 50 --days 30 \
    --name "Home" --tz Europe/Paris

# Also track Tiangong
python iss_transit_calendar.py --lat 48.85 --lon 2.35 --radius 50 --satellites iss,css
```

| Option | Default | Description |
| --- | --- | --- |
| `--lat`, `--lon` | required | Centre of the search area (degrees) |
| `--radius` | 50 | Search radius (km) |
| `--days` | 30 | Forecast horizon (days) |
| `--name` | | Location name shown on the maps |
| `--satellites` | `iss` | Comma-separated list: `iss`, `css` (Tiangong). Example: `iss,css` |
| `--drive-folder` | `Satellite Transits` | Drive folder that receives the maps |
| `--grid-step-km` | 3 | Spacing of the observer grid inside the radius |
| `--coarse-step` | 2 | Pre-scan time step (seconds) |
| `--approach-deg` | 0 (auto) | Pre-scan threshold (degrees); auto depends on the radius |
| `--calendar-id` | `primary` | Target calendar |
| `--tz` | `UTC` | Time zone used in event titles and descriptions |
| `--min-alt` | 5 | Minimum Sun/Moon altitude (degrees) |
| `--lead-min` | 5 | The event starts N minutes before the transit |
| `--time-tol` | 3 | Time shift (s) that triggers an event update |
| `--sep-tol` | 3 | Separation change (arcsec) that triggers an update |
| `--share-images` | off | Make Drive maps readable by anyone with the link |
| `--no-maps` | off | Skip map generation and upload |
| `--weather-days` | 10 | Rate transits in the next N days with a cloud forecast (max 15, 0 disables) |
| `--weather-model` | | Open-Meteo model, e.g. `meteofrance_seamless`, `ecmwf_ifs025` (default: Open-Meteo's best match) |
| `--no-weather` | off | Skip the weather evaluation |
| `--dry-run` | off | Compute and print only: no Google credentials, calendar untouched |
| `--maps-dir` | | With `--dry-run`: also render the corridor maps into this folder |
| `--workdir` | script folder | Where `credentials.json`, `token.json` and Skyfield data live |

Exit codes: `0` success, `1` computation failed (calendar untouched), `2` Google synchronisation failed.

## Scheduling

```cron
15 3 * * * cd /opt/iss-transit-calendar && ./venv/bin/python iss_transit_calendar.py \
    --lat 48.85 --lon 2.35 --radius 50 --days 30 --name "Home" --tz Europe/Paris >> run.log 2>&1
```
<<<<<<< HEAD

Synology Task Scheduler (or any other scheduler) works the same way. The search can take a while with a fine grid and a large radius; a lock file prevents a new run from starting while the previous one is still going. Time a `--days 7 --dry-run` run first and extrapolate.

## Weather outlook and interest score

For transits within `--weather-days`, the script fetches the hourly forecast at the best observing spot (only those coordinates are sent to Open-Meteo), interpolates it at the transit time and rates the event:

- **Sky cover**: low and mid clouds count fully, high clouds count for half (a thin cirrus veil often lets the Sun through with a filter), layers are assumed independent. A rain probability of 70 % or more also pushes the rating down.
- **Weather level**: below 25 % sky cover is *good*, below 60 % *fair*, otherwise *poor*.
- **Geometry level**: *good* if the chord covers at least 60 % of the disc and the Sun/Moon is at least 15° high; *fair* from 30 % and 8°; otherwise *poor*.
- **Interest** is the lower of the two levels.
- **Trend**: beyond D+5 a cloud forecast is only a trend. The title shows `~` before the sky cover, and a *poor* weather rating is softened to *fair* so that a good transit is not written off on an uncertain forecast.

How it shows up in the calendar:

| Where | What |
| --- | --- |
| Title | `🟢 ☀️ ISS transit - Sun 07:43:46 · ☁ 10 %` (level first, so it is visible in month, agenda and widget views) |
| Event colour | Green / yellow / red for good / fair / poor; graphite when there is no weather data (outside the horizon, or a fetch failure) |
| Description | Interest line, low/mid/high cloud cover, rain probability, visibility, fetch time and the Open-Meteo credit |

The rating is refreshed on every run. When only the weather changed (level change, or sky cover moved by 15 points or more), the title, colour and description are patched and the map is left alone. If the weather service is unreachable, events are still created and updated, just without a rating.

=======
 
>>>>>>> 9f499b3ab0dee56a1716e1c1a0b29a975639420f
## How the search differs from the original transit-finder

The original tool tests the centre plus 8 points on the circle's circumference. For this use case that is not enough:

1. **Sampling.** The interior of the disc is never tested, and with a 100 km radius the 8 points are ~80 km apart while the corridor is only ~4-19 km wide. This script tests a grid of points every `--grid-step-km` (≈950 points for a 50 km radius at 3 km spacing).
2. **Pre-scan threshold.** The original keeps a candidate only if the ISS passes within 2° of the Sun/Moon seen from the centre. Because of parallax (the satellite is 400-2000 km away), an observer 50 km away sees it shifted by several degrees, so valid transits were discarded. The threshold now grows with the radius.
3. **Time step.** A 10 s step lets the ISS travel up to ~10° between samples near the zenith, so a fast pass could be missed. The default step is now 2 s, and the candidate time is the instant of minimum separation.

## Limitations

- ISS position accuracy degrades by roughly 1-2 km per day of forecast, and reboost maneuvers can shift things further. The corridor is only a few km wide: treat anything beyond ~7 days as a heads-up, and confirm a few days before travelling.
- Tiangong uses the Celestrak TLE for NORAD 48274 (Tianhe core module). `transit_finder.py` was written for the ISS: its own console messages say "ISS", the computation itself works from the TLE, and this script fixes the map titles for other satellites.
- Computation time is per satellite: tracking both roughly doubles the run time.
- Cloud forecasts become unreliable after about 5 days, and a transit lasts under a second: a 🟢 rating at D+8 is a good sign, not a promise. Re-check the event a day or two before travelling.
- Open-Meteo is free for non-commercial use and needs no API key; attribution is included in the event description.
- Events deleted by hand are not recreated (their ID stays reserved by Google).
- The map is a Drive attachment/link; Google Calendar cannot display images inline. Maps are only visible to you unless you use `--share-images`.
- Map generation downloads tiles; without network access events are created without a basemap or without an image.
- Not tested on Windows (`fcntl`).

## Credits

- [beschne/transit-finder](https://github.com/beschne/transit-finder) (MIT) - transit and corridor computation.
- [Skyfield](https://rhodesmill.org/skyfield/) and the JPL DE421 ephemeris.
- [Celestrak](https://celestrak.org/) for the satellite TLEs.
- [Open-Meteo](https://open-meteo.com/) for the weather forecasts (CC BY 4.0 data).

## License
<<<<<<< HEAD

[MIT](LICENSE)
=======
 
[MIT](LICENSE)
 
>>>>>>> 9f499b3ab0dee56a1716e1c1a0b29a975639420f
