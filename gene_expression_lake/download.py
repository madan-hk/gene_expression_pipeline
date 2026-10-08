"""Download a file safely: to a temporary name first, renamed only when complete."""
from __future__ import annotations

import time
from pathlib import Path

import requests


class DownloadError(RuntimeError):
    pass


def download(url: str, target: Path, attempts: int = 3) -> bool:
    """Fetch `url` into `target`. Returns True if downloaded, False if it was already there.

    The data goes to `<target>.part` and is renamed only when the whole file has
    arrived and its size matches what the server announced, so a cut-off download
    can never be mistaken for a finished one."""
    if target.exists():
        return False
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".part")
    last = "unknown error"
    for attempt in range(1, attempts + 1):
        try:
            with requests.get(url, stream=True, timeout=60) as r:
                if r.status_code == 404:
                    raise DownloadError(f"Not found on recount3: {url}")
                r.raise_for_status()
                expected = int(r.headers.get("Content-Length", -1))
                written = 0
                with open(part, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
                        written += len(chunk)
            if expected >= 0 and written != expected:
                raise DownloadError(f"Incomplete download: got {written} of {expected} bytes")
            part.replace(target)
            return True
        except DownloadError:
            part.unlink(missing_ok=True)
            raise
        except (requests.RequestException, OSError) as e:
            last = f"{type(e).__name__}: {e}"
            part.unlink(missing_ok=True)
            if attempt < attempts:
                time.sleep(2 * attempt)
    raise DownloadError(f"Download failed after {attempts} attempts ({last}): {url}")


def remote_size(url: str) -> int | None:
    """Size in bytes of the file at `url` (asked of the server, nothing is downloaded), or None."""
    try:
        r = requests.head(url, timeout=30, allow_redirects=True)
        return int(r.headers["Content-Length"]) if r.status_code == 200 and "Content-Length" in r.headers else None
    except (requests.RequestException, ValueError):
        return None
