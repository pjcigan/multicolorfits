"""Capture light/dark screenshots of the browser GUI for docs.

Starts a local ``mcf-web`` session (NGC 602 by default), waits for panels
and Fast preview to settle, then grabs PNGs with headless Chrome.

Usage (repo root)::

    python docs/capture_web_gui_screenshots.py
"""
from __future__ import annotations

import os
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from shutil import which

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = Path(__file__).resolve().parent / '_static' / 'gui'
NGC602 = ROOT / 'hidden' / 'testing' / 'ngc602'
CHROME = 'google-chrome'


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(('127.0.0.1', 0))
        return int(s.getsockname()[1])


def _wait_http(url: str, timeout: float = 60.0) -> None:
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            with urllib.request.urlopen(url, timeout=2) as r:
                if r.status == 200:
                    return
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            time.sleep(0.25)
    raise RuntimeError('server did not become ready: %s' % url)


def _chrome_screenshot(url: str, dest: Path, width: int = 1600, height: int = 1000,
                       settle_ms: int = 20000) -> None:
    """Headless Chrome screenshot after virtual-time settle."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Write to a temp path in CWD-friendly location; chrome --screenshot
    # defaults to cwd/screenshot.png when given a bare filename.
    tmp = dest.with_suffix('.tmp.png')
    if tmp.exists():
        tmp.unlink()
    cmd = [
        CHROME,
        '--headless=new',
        '--disable-gpu',
        '--no-sandbox',
        '--disable-dev-shm-usage',
        '--hide-scrollbars',
        '--force-device-scale-factor=1',
        f'--window-size={width},{height}',
        f'--virtual-time-budget={settle_ms}',
        f'--screenshot={tmp}',
        url,
    ]
    subprocess.run(cmd, check=True, cwd=str(dest.parent),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not tmp.is_file():
        # Older chrome may ignore absolute --screenshot paths; fall back.
        fallback = dest.parent / 'screenshot.png'
        if fallback.is_file():
            fallback.rename(tmp)
    if not tmp.is_file():
        raise RuntimeError('chrome did not write screenshot for %s' % url)
    tmp.replace(dest)
    print('wrote', dest, '(%d bytes)' % dest.stat().st_size)


def _build_session_json(path: Path, data_dir: Path) -> None:
    os.environ.setdefault('MPLBACKEND', 'Agg')
    import multicolorfits as mcf

    s = mcf.McfSession(n_panels=3)
    files = [
        data_dir / 'ngc602_ir.fits',
        data_dir / 'ngc602_optical_R.fits',
        data_dir / 'ngc602_optical_B.fits',
    ]
    colors = ['#BE599E', '#DEA215', '#77C0F9']
    labels = ['IR', 'R', 'B']
    for i, (fp, color, label) in enumerate(zip(files, colors, labels)):
        s.panels[i].load_fits(str(fp))
        s.panels[i].stretch = 'linear'
        s.panels[i].set_percentiles(1.0, 99.5)
        s.panels[i].color = color
        s.panels[i].label = label
    s.compose.combine_mode = 'lab'
    s.compose.combine_blend = 'screen'
    s.compose.gamma = 2.2
    s.compose.show_legend = True
    s.compose.title = 'NGC 602'
    s.save_state(str(path))


def _start_server(port: int, state: Path) -> subprocess.Popen:
    env = os.environ.copy()
    env['MPLBACKEND'] = 'Agg'
    env.pop('QT_QPA_PLATFORM', None)
    cmd = [
        sys.executable, '-c',
        (
            'import multicolorfits as mcf;'
            'mcf.gui(host="127.0.0.1", port=%d, open_browser=False, state=%r)'
        ) % (port, str(state)),
    ]
    return subprocess.Popen(
        cmd, cwd=str(ROOT), env=env,
        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
    )


def main() -> int:
    if not NGC602.is_dir():
        print('NGC 602 data missing at', NGC602, file=sys.stderr)
        return 1
    if which(CHROME) is None:
        print('google-chrome not found on PATH', file=sys.stderr)
        return 1

    # Short paths in the panel filepath fields (docs screenshots).
    data_dir = Path(tempfile.mkdtemp(prefix='mcf-ngc602-'))
    for name in ('ngc602_ir.fits', 'ngc602_optical_R.fits', 'ngc602_optical_B.fits'):
        os.symlink(NGC602 / name, data_dir / name)

    state = data_dir / 'ngc602_docs.json'
    _build_session_json(state, data_dir)

    port = _free_port()
    server = _start_server(port, state)
    try:
        _wait_http(f'http://127.0.0.1:{port}/', timeout=120)
        # Let the server finish reading FITS into panel state.
        time.sleep(4.0)
        for theme, name in (('light', 'gui_web_walkthrough_light.png'),
                            ('dark', 'gui_web_walkthrough_dark.png')):
            url = f'http://127.0.0.1:{port}/?theme={theme}'
            _chrome_screenshot(url, OUT_DIR / name, settle_ms=25000)
    finally:
        server.terminate()
        try:
            server.wait(timeout=8)
        except subprocess.TimeoutExpired:
            server.kill()
        err = ''
        if server.stderr is not None:
            err = server.stderr.read().decode(errors='replace')
        if err and 'Error' in err:
            print(err[-2500:], file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
