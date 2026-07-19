"""
Launcher for the multicolorfits browser GUI: starts a local uvicorn server
and opens the page in a web browser.
"""

import argparse
import socket
import threading
import webbrowser


def _port_available(host, port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex((host, port)) != 0


def _build_session(session=None, state=None, files=None, colors=None, labels=None):
    """Return a McfSession, optionally pre-loaded from state or FITS paths."""
    from ..session import McfSession

    if session is not None:
        s = session
    else:
        s = McfSession()
    if state:
        s.load_state(state)
    elif files:
        s.load_files(files, colors=colors, labels=labels)
    return s


def resolve_port(host, port, max_tries=50):
    """Return the first port in ``port, port+1, ...`` that is free on *host*."""
    tries = 0
    while not _port_available(host, port) and tries < max_tries:
        port += 1
        tries += 1
    return port


def open_gui_url(url, browser=None):
    """
    Open *url* in a browser.

    Parameters
    ----------
    url : str
    browser : str or None
        Browser name understood by :func:`webbrowser.get` (e.g. ``'firefox'``,
        ``'google-chrome'``, ``'chromium'``, ``'safari'``), or a full path to a
        browser executable.  ``None`` uses the system default (on Unix, the
        ``BROWSER`` environment variable is also consulted by the stdlib).
    """
    try:
        if browser:
            webbrowser.get(browser).open(url)
        else:
            webbrowser.open(url)
    except webbrowser.Error as exc:
        print('Could not open browser %r (%s). Open manually: %s' % (browser, exc, url))


def run(host='127.0.0.1', port=8321, open_browser=True, browser=None,
        session=None, state=None, files=None, colors=None, labels=None):
    """Start the multicolorfits web GUI (blocks until the server is stopped).

    Parameters
    ----------
    host, port : server bind address
    open_browser : bool
        If True, open a tab after the server starts.
    browser : str or None
        Which browser to open when ``open_browser`` is True (see
        :func:`open_gui_url`).  Ignored when ``open_browser`` is False.
    session, state, files, colors, labels
        Optional pre-load of an :class:`~multicolorfits.McfSession`.
    """
    import uvicorn
    from .server import create_app

    port = resolve_port(host, port)
    app = create_app(session=_build_session(session=session, state=state,
                                              files=files, colors=colors,
                                              labels=labels))
    url = 'http://%s:%d' % (host, port)
    print('multicolorfits web GUI:  %s   (Ctrl+C to quit)' % url)
    if open_browser:
        threading.Timer(0.8, lambda: open_gui_url(url, browser=browser)).start()
    uvicorn.run(app, host=host, port=port, log_level='warning')


def main():
    parser = argparse.ArgumentParser(description='multicolorfits browser GUI')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', type=int, default=8321)
    browser_group = parser.add_mutually_exclusive_group()
    browser_group.add_argument(
        '--no-browser', action='store_true',
        help="don't auto-open a browser tab",
    )
    browser_group.add_argument(
        '--browser', default=None, metavar='NAME',
        help="browser to open (e.g. firefox, google-chrome, chromium, safari); "
             "default is the system browser (Unix: also respects $BROWSER)",
    )
    parser.add_argument('--state', default=None,
                        help='JSON session file to restore on startup')
    parser.add_argument('--files', nargs='*', default=None,
                        help='FITS files to load into panels on startup')
    parser.add_argument('--colors', nargs='*', default=None,
                        help='Hex colors for --files (one per file)')
    parser.add_argument('--labels', nargs='*', default=None,
                        help='Channel labels for --files')
    args = parser.parse_args()
    run(host=args.host, port=args.port,
        open_browser=not args.no_browser, browser=args.browser,
        state=args.state, files=args.files, colors=args.colors, labels=args.labels)


if __name__ == '__main__':
    main()
