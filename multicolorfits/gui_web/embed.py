"""
Embed the multicolorfits browser GUI inside Jupyter, Colab, VS Code notebooks,
and similar environments.

Typical usage::

    import multicolorfits as mcf

    s = mcf.McfSession()
    s.load_files(['a.fits', 'b.fits'], colors=['#f00', '#0ff'])
    mcf.gui_embed(session=s)          # auto-detect notebook environment
    mcf.gui_embed(session=s, height=950, display='colab')
"""

from __future__ import annotations

import os
import socket
import threading
import time
from dataclasses import dataclass
from typing import Optional

from .launcher import _build_session, resolve_port


@dataclass
class GuiServer:
  """Handle for a background multicolorfits web GUI server."""

  host: str
  port: int
  url: str
  thread: threading.Thread
  _server: object = None

  def stop(self):
    """Request a graceful server shutdown (best-effort)."""
    srv = self._server
    if srv is not None:
      srv.should_exit = True


def _server_listening(host: str, port: int) -> bool:
  with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    return s.connect_ex((host, port)) == 0


def _wait_until_listening(host: str, port: int, timeout: float = 10.0) -> bool:
  deadline = time.time() + timeout
  while time.time() < deadline:
    if _server_listening(host, port):
      return True
    time.sleep(0.05)
  return False


def start_gui_server(host='127.0.0.1', port=8321, session=None,
                     state=None, files=None, colors=None, labels=None,
                     *, wait: bool = True, timeout: float = 10.0) -> GuiServer:
  """
  Start the web GUI in a background daemon thread.

  Returns a :class:`GuiServer` with ``url``, ``port``, and ``host``.  The
  thread keeps running until the kernel exits or :meth:`GuiServer.stop` is
  called.

  Parameters
  ----------
  host, port
      Bind address.  ``port`` is advanced automatically when busy (same as
      :func:`~multicolorfits.gui`).
  session, state, files, colors, labels
      Same pre-load options as :func:`~multicolorfits.gui`.
  wait : bool
      When True (default), block until the server accepts connections.
  timeout : float
      Seconds to wait when ``wait`` is True.
  """
  import uvicorn
  from .server import create_app

  port = resolve_port(host, port)
  app = create_app(session=_build_session(session=session, state=state,
                                            files=files, colors=colors,
                                            labels=labels))
  config = uvicorn.Config(app, host=host, port=port, log_level='warning')
  server = uvicorn.Server(config)
  thread = threading.Thread(target=server.run, daemon=True, name='mcf-web-gui')
  thread.start()
  if wait and not _wait_until_listening(host, port, timeout=timeout):
    server.should_exit = True
    raise RuntimeError('multicolorfits web GUI did not start on %s:%d within %.1fs'
                       % (host, port, timeout))
  url = 'http://%s:%d' % (host, port)
  return GuiServer(host=host, port=port, url=url, thread=thread, _server=server)


def detect_notebook_environment():
  """
  Guess the notebook host environment, or return ``None`` outside IPython.

  Returns
  -------
  str or None
      One of ``'colab'``, ``'jupyter'``, ``'vscode'``, ``'binder'``, or
      ``None`` when not in an IPython kernel.
  """
  try:
    ip = get_ipython()  # noqa: F821
  except NameError:
    return None
  shell = getattr(ip, '__class__', type('x', (), {})).__name__
  if shell == 'GoogleColabShell':
    return 'colab'
  if os.environ.get('BINDER_SERVICE_URL') or os.environ.get('BINDER_REPO_URL'):
    return 'binder'
  if os.environ.get('VSCODE_PID') or os.environ.get('VSCODE_INJECTION'):
    return 'vscode'
  if shell in ('ZMQInteractiveShell', 'TerminalInteractiveShell'):
    return 'jupyter'
  return 'jupyter'


def jupyter_proxy_url(port: int, path: str = '') -> Optional[str]:
  """
  Build a JupyterHub / Binder ``/proxy/{port}/`` URL when available.

  Returns ``None`` when no proxy prefix is detected (local Jupyter can use the
  direct ``http://127.0.0.1:{port}`` URL instead).
  """
  prefix = os.environ.get('JUPYTERHUB_SERVICE_PREFIX', '')
  if not prefix and os.environ.get('BINDER_SERVICE_URL'):
    prefix = '/'
  if not prefix:
    return None
  if not prefix.endswith('/'):
    prefix += '/'
  rel = 'proxy/%d/' % int(port)
  if path:
    rel += path.lstrip('/')
  return prefix + rel


def embed_url(server: GuiServer, *, proxy_url: Optional[str] = None) -> str:
  """Pick the URL to load in an iframe for *server*."""
  if proxy_url:
    return proxy_url
  proxied = jupyter_proxy_url(server.port)
  if proxied:
    return proxied
  return server.url


def gui_embed(host='127.0.0.1', port=8321, session=None,
              state=None, files=None, colors=None, labels=None,
              *, display='auto', height=900, width='100%',
              proxy_url: Optional[str] = None, return_url: bool = False,
              print_url: bool = False) -> GuiServer | str:
  """
  Start the web GUI and embed it in the current notebook output.

  This is the recommended entry point for demos in Jupyter, Google Colab,
  VS Code notebooks, and Binder.  It starts uvicorn in a background thread and
  renders an iframe (or Colab's port proxy) when IPython display is available.

  Parameters
  ----------
  host, port, session, state, files, colors, labels
      Same as :func:`~multicolorfits.gui`.
  display : {'auto', 'iframe', 'colab', 'url', 'none'}
      * ``auto`` — Colab proxy when in Colab, iframe in other IPython kernels,
        otherwise print the URL.
      * ``iframe`` — ``IPython.display.IFrame`` (local Jupyter / VS Code).
      * ``colab`` — ``google.colab.output.serve_kernel_port_as_iframe``.
      * ``url`` — print the URL only (open manually).
      * ``none`` — start the server but do not display anything.
  height, width
      iframe dimensions (ignored for Colab except ``height``).
  proxy_url : str or None
      Override the iframe ``src`` (advanced: custom JupyterHub proxy paths).
  return_url : bool
      When True, return the server URL string instead of displaying.
  print_url : bool
      When True, also ``print`` the direct local URL (useful with ``display='colab'``).

  Returns
  -------
  GuiServer or str
      Server handle, unless ``return_url=True`` (then the URL string).

  Examples
  --------
  Preload a demo session and embed inline::

      import multicolorfits as mcf
      s = mcf.McfSession()
      s.load_files(['halpha.fits', 'oiii.fits'], colors=['#f00', '#0ff'])
      mcf.gui_embed(session=s, height=950)
  """
  server = start_gui_server(host=host, port=port, session=session, state=state,
                            files=files, colors=colors, labels=labels)
  url = embed_url(server, proxy_url=proxy_url)

  if print_url or display == 'url':
    print('multicolorfits web GUI:  %s' % server.url)
    if proxy_url or jupyter_proxy_url(server.port):
      print('Notebook iframe URL:     %s' % url)

  if return_url:
    return url

  if display == 'auto':
    env = detect_notebook_environment()
    if env == 'colab':
      display = 'colab'
    elif env in ('jupyter', 'vscode', 'binder'):
      display = 'iframe'
    else:
      display = 'url'

  if display == 'none':
    return server

  if display == 'colab':
    try:
      from google.colab import output
    except ImportError as exc:
      raise ImportError(
          "display='colab' requires google.colab (run inside Google Colab)"
      ) from exc
    output.serve_kernel_port_as_iframe(server.port, height=int(height))
    return server

  if display == 'iframe':
    try:
      from IPython.display import IFrame, display as ipy_display
    except ImportError as exc:
      raise ImportError(
          "display='iframe' requires IPython (install ipython or run in Jupyter)"
      ) from exc
    ipy_display(IFrame(src=url, width=width, height=int(height)))
    return server

  # display == 'url' already printed above
  return server
