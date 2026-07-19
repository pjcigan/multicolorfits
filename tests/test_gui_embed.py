"""Tests for notebook / Colab GUI embedding helpers."""

import socket
import threading

import pytest


def _free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


class TestResolvePort:
    def test_resolve_port_returns_first_free(self):
        from multicolorfits.gui_web.launcher import resolve_port, _port_available

        port = _free_port()
        assert resolve_port('127.0.0.1', port) == port
        assert _port_available('127.0.0.1', port)


class TestNotebookDetect:
    def test_outside_ipython_is_none(self):
        from multicolorfits.gui_web.embed import detect_notebook_environment

        assert detect_notebook_environment() is None

    def test_jupyter_proxy_url_with_prefix(self, monkeypatch):
        from multicolorfits.gui_web.embed import jupyter_proxy_url

        monkeypatch.setenv('JUPYTERHUB_SERVICE_PREFIX', '/user/foo/')
        assert jupyter_proxy_url(8321) == '/user/foo/proxy/8321/'

    def test_jupyter_proxy_url_none_locally(self, monkeypatch):
        from multicolorfits.gui_web.embed import jupyter_proxy_url

        monkeypatch.delenv('JUPYTERHUB_SERVICE_PREFIX', raising=False)
        monkeypatch.delenv('BINDER_SERVICE_URL', raising=False)
        assert jupyter_proxy_url(8321) is None


@pytest.fixture
def web_extra():
    pytest.importorskip('fastapi')
    pytest.importorskip('uvicorn')


class TestStartGuiServer:
    def test_start_and_fetch_state(self, web_extra):
        from multicolorfits.gui_web.embed import start_gui_server

        port = _free_port()
        server = start_gui_server(port=port)
        try:
            import urllib.request
            with urllib.request.urlopen(server.url + '/api/state', timeout=5) as resp:
                assert resp.status == 200
            assert server.port == port
            assert server.url == 'http://127.0.0.1:%d' % port
        finally:
            server.stop()

    def test_preloaded_session(self, web_extra, test_fits_file):
        import json
        import urllib.request
        from multicolorfits.gui_web.embed import start_gui_server
        from multicolorfits.session import McfSession

        s = McfSession()
        s.load_files([test_fits_file], colors=['#ff0000'])
        port = _free_port()
        server = start_gui_server(session=s, port=port)
        try:
            with urllib.request.urlopen(server.url + '/api/state', timeout=5) as resp:
                data = json.loads(resp.read().decode())
            assert data['panels'][0]['loaded'] is True
        finally:
            server.stop()


class TestGuiEmbedDisplay:
    def test_return_url(self, web_extra, monkeypatch):
        from multicolorfits.gui_web import embed as embed_mod

        port = _free_port()
        fake = embed_mod.GuiServer('127.0.0.1', port, 'http://127.0.0.1:%d' % port,
                                   threading.Thread())

        monkeypatch.setattr(embed_mod, 'start_gui_server', lambda **kw: fake)
        monkeypatch.setattr(embed_mod, 'detect_notebook_environment', lambda: None)

        import multicolorfits as mcf
        url = mcf.gui_embed(port=port, return_url=True)
        assert url == fake.url

    def test_display_none_returns_server(self, web_extra, monkeypatch):
        from multicolorfits.gui_web import embed as embed_mod

        port = _free_port()
        fake = embed_mod.GuiServer('127.0.0.1', port, 'http://127.0.0.1:%d' % port,
                                   threading.Thread())

        monkeypatch.setattr(embed_mod, 'start_gui_server', lambda **kw: fake)

        import multicolorfits as mcf
        srv = mcf.gui_embed(port=port, display='none')
        assert srv is fake
