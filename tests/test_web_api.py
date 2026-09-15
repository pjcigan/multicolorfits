"""Web GUI backend: endpoint behavior via FastAPI's TestClient."""

import os

import numpy as np
import pytest

fastapi = pytest.importorskip('fastapi')
from fastapi.testclient import TestClient

from multicolorfits.gui_web.server import create_app


@pytest.fixture
def client():
    return TestClient(create_app())


@pytest.fixture
def loaded_client(client, test_fits_file):
    r = client.post('/api/panel/0/load', json={'path': test_fits_file})
    assert r.status_code == 200, r.text
    return client


class TestState:
    def test_initial_state(self, client):
        r = client.get('/api/state')
        assert r.status_code == 200
        data = r.json()
        assert len(data['panels']) == 4
        assert data['compose']['gamma'] == 2.2
        assert 'asinh' in data['stretches']
        assert data['max_panels'] == 16
        assert data['min_panels'] == 1

    def test_index_served(self, client):
        r = client.get('/')
        assert r.status_code == 200
        assert 'MultiColorFits' in r.text
        assert 'btn-add-panel' in r.text


class TestDynamicPanels:
    def test_add_and_remove_panel(self, client):
        r = client.post('/api/panels/add', json={})
        assert r.status_code == 200, r.text
        body = r.json()
        assert body['index'] == 4
        assert len(body['state']['panels']) == 5

        r2 = client.post('/api/panels/4/remove')
        assert r2.status_code == 200, r2.text
        assert len(r2.json()['state']['panels']) == 4

    def test_cannot_remove_last_panel(self, client):
        # Shrink to one via repeated remove from a fresh 4-panel session.
        for idx in (3, 2, 1):
            assert client.post('/api/panels/%d/remove' % idx).status_code == 200
        assert len(client.get('/api/state').json()['panels']) == 1
        r = client.post('/api/panels/0/remove')
        assert r.status_code == 400


class TestPanelLoad:
    def test_load_from_path(self, client, test_fits_file):
        r = client.post('/api/panel/0/load', json={'path': test_fits_file})
        assert r.status_code == 200
        p = r.json()
        assert p['loaded'] is True
        assert p['shape'] == [64, 64]

    def test_load_missing_file(self, client):
        r = client.post('/api/panel/0/load', json={'path': '/nonexistent/file.fits'})
        assert r.status_code == 400

    def test_bad_panel_index(self, client):
        r = client.post('/api/panel/9/load', json={'path': 'x'})
        assert r.status_code == 404

    def test_upload(self, client, test_fits_file, tmp_path, monkeypatch):
        # Browser uploads land in ~/.cache/.../uploads; point that at a
        # writable temp dir so a root-owned home cache cannot fail this test.
        import multicolorfits.gui_web.server as web_server
        monkeypatch.setattr(web_server, 'UPLOAD_CACHE_DIR', str(tmp_path / 'uploads'))
        with open(test_fits_file, 'rb') as f:
            r = client.post('/api/panel/1/upload', files={'file': ('synthetic.fits', f, 'application/fits')})
        assert r.status_code == 200, r.text
        assert r.json()['loaded'] is True
        assert (tmp_path / 'uploads' / 'synthetic.fits').is_file()

    def test_clear(self, loaded_client):
        r = loaded_client.post('/api/panel/0/clear')
        assert r.status_code == 200
        assert loaded_client.get('/api/state').json()['panels'][0]['loaded'] is False


class TestPanelParams:
    def test_set_stretch_color(self, loaded_client):
        r = loaded_client.post('/api/panel/0/params', json={'stretch': 'asinh', 'color': '#C11B17'})
        p = r.json()
        assert p['stretch'] == 'asinh'
        assert p['color'] == '#C11B17'

    def test_bad_stretch(self, loaded_client):
        r = loaded_client.post('/api/panel/0/params', json={'stretch': 'bogus'})
        assert r.status_code == 400

    def test_vmin_updates_percentile(self, loaded_client):
        state = loaded_client.get('/api/state').json()['panels'][0]
        mid = (state['data_min'] + state['data_max']) / 2.
        p = loaded_client.post('/api/panel/0/params', json={'vmin': mid}).json()
        assert p['vmin'] == pytest.approx(mid)
        assert p['percent_min'] > 0

    def test_percentiles_update_limits(self, loaded_client):
        p = loaded_client.post('/api/panel/0/params', json={'percent_min': 5, 'percent_max': 95}).json()
        assert p['vmin'] > p['data_min']
        assert p['vmax'] < p['data_max']

    def test_zscale_and_minmax(self, loaded_client):
        z = loaded_client.post('/api/panel/0/zscale').json()
        assert z['vmin'] < z['vmax']
        m = loaded_client.post('/api/panel/0/minmax').json()
        assert m['vmin'] == pytest.approx(m['data_min'])
        assert m['vmax'] == pytest.approx(m['data_max'])


class TestViews:
    def test_preview_png(self, loaded_client):
        r = loaded_client.get('/api/panel/0/preview.png')
        assert r.status_code == 200
        assert r.headers['content-type'] == 'image/png'
        assert r.content[:8] == b'\x89PNG\r\n\x1a\n'

    def test_preview_unloaded_404(self, client):
        assert client.get('/api/panel/2/preview.png').status_code == 404

    def test_preview_buffer(self, loaded_client):
        r = loaded_client.get('/api/panel/0/preview_buffer?max_size=64')
        assert r.status_code == 200
        data = r.json()
        assert len(data['shape']) == 2
        assert data['shape'][0] > 0 and data['shape'][1] > 0
        import base64
        raw = base64.b64decode(data['bytes'])
        assert len(raw) == data['shape'][0] * data['shape'][1] * 4

    def test_preview_buffer_unloaded_404(self, client):
        assert client.get('/api/panel/1/preview_buffer').status_code == 404

    def test_preview_buffer_bin(self, loaded_client):
        r = loaded_client.get('/api/panel/0/preview_buffer.bin?max_size=64')
        assert r.status_code == 200
        raw = r.content
        assert len(raw) >= 8
        ny = int.from_bytes(raw[0:4], 'little')
        nx = int.from_bytes(raw[4:8], 'little')
        assert len(raw) == 8 + ny * nx * 4

    def test_histogram(self, loaded_client):
        r = loaded_client.get('/api/panel/0/histogram?bins=32')
        data = r.json()
        assert len(data['counts']) == 32
        assert len(data['edges']) == 33

    def test_header_roundtrip(self, loaded_client):
        text = loaded_client.get('/api/panel/0/header').text
        assert 'CRVAL1' in text
        r = loaded_client.post('/api/panel/0/header', json={'text': text})
        assert r.status_code == 200

    def test_combined_png(self, loaded_client):
        r = loaded_client.get('/api/combined.png')
        assert r.status_code == 200
        assert r.content[:4] == b'\x89PNG'

    def test_combined_png_preview(self, loaded_client):
        r = loaded_client.get('/api/combined.png?preview=true&max_size=128')
        assert r.status_code == 200
        assert r.content[:4] == b'\x89PNG'

    def test_combined_plot_png(self, loaded_client):
        r = loaded_client.get('/api/combined_plot.png')
        assert r.status_code == 200
        assert r.content[:4] == b'\x89PNG'

    def test_combined_no_data_400(self, client):
        assert client.get('/api/combined.png').status_code == 400


class TestComposeAndOutput:
    def test_set_compose(self, client):
        r = client.post('/api/compose', json={'gamma': 1.5, 'inverse': True, 'tickcolor': '#FF00FF'})
        c = r.json()
        assert c['gamma'] == 1.5
        assert c['inverse'] is True
        assert c['tickcolor'] == '#FF00FF'

    def test_set_combine_mode_blend(self, client):
        r = client.post('/api/compose', json={'combine_mode': 'lab', 'combine_blend': 'max'})
        c = r.json()
        assert c['combine_mode'] == 'lab'
        assert c['combine_blend'] == 'max'

    def test_set_facecolor(self, client):
        c = client.post('/api/compose', json={'facecolor': '#204060'}).json()
        assert c['facecolor'] == '#204060'

    def test_facecolor_empty_becomes_none(self, client):
        c = client.post('/api/compose', json={'facecolor': ''}).json()
        assert c['facecolor'] == 'none'

    def test_save_image_transparent(self, loaded_client, tmp_path):
        loaded_client.post('/api/compose', json={'facecolor': 'none'})
        path = str(tmp_path / 'transparent.png')
        r = loaded_client.post('/api/save/image', json={'path': path, 'dpi': 50})
        assert r.status_code == 200, r.text
        from matplotlib import image as mpl_image
        arr = mpl_image.imread(path)
        assert arr.shape[-1] == 4 and (arr[..., 3] == 0).any()

    def test_export_cutout_download(self, loaded_client):
        r = loaded_client.post('/api/export/cutout', json={
            'size': 48, 'crop': 'auto', 'download': True,
            'alpha_lo': 40, 'alpha_hi': 99, 'alpha_gamma': 0.5,
        })
        assert r.status_code == 200, r.text
        assert r.headers['content-type'].startswith('image/png')
        assert len(r.content) > 100

    def test_export_cutout_save_path(self, loaded_client, tmp_path):
        path = str(tmp_path / 'cutout.png')
        r = loaded_client.post('/api/export/cutout', json={
            'size': 32, 'path': path, 'download': False, 'crop': 'none',
        })
        assert r.status_code == 200, r.text
        assert r.json()['ok'] is True
        assert os.path.isfile(path)

    def test_invalid_combine_mode(self, client):
        assert client.post('/api/compose', json={'combine_mode': 'bogus'}).status_code == 400

    def test_invalid_combine_blend(self, client):
        assert client.post('/api/compose', json={'combine_blend': 'bogus'}).status_code == 400

    def test_set_combine_mode_ryb(self, client):
        c = client.post('/api/compose', json={'combine_mode': 'ryb'}).json()
        assert c['combine_mode'] == 'ryb'

    def test_set_combine_mode_cmyk(self, client):
        c = client.post('/api/compose', json={'combine_mode': 'cmyk'}).json()
        assert c['combine_mode'] == 'cmyk'

    def test_session_save_load(self, loaded_client, tmp_path):
        path = str(tmp_path / 's.json')
        r = loaded_client.post('/api/session/save', json={'path': path})
        assert r.status_code == 200
        assert (tmp_path / 's.json').stat().st_size > 0
        c2 = TestClient(create_app())
        r2 = c2.post('/api/session/load', json={'path': path})
        assert r2.status_code == 200
        assert c2.get('/api/state').json()['compose']['gamma'] == 2.2

    def test_session_json_download(self, loaded_client):
        r = loaded_client.get('/api/session.json')
        assert r.status_code == 200
        assert r.headers['content-type'].startswith('application/json')
        body = r.text
        assert len(body) > 20
        assert 'multicolorfits_state_version' in body

    def test_session_reset(self, loaded_client):
        r = loaded_client.post('/api/session/reset')
        assert r.status_code == 200
        st = r.json()['state']
        assert not any(p.get('loaded') for p in st['panels'])
        assert st['compose']['gamma'] == 2.2
        assert loaded_client.get('/api/state').json()['panels'][0]['loaded'] is False

    def test_cursor_api(self, loaded_client):
        st = loaded_client.get('/api/state').json()
        ny, nx = st['panels'][0]['shape']
        r = loaded_client.get('/api/cursor', params={'x': nx // 2, 'y': ny // 2})
        assert r.status_code == 200
        assert r.json()['layers']

    def test_set_tick_direction(self, client):
        c = client.post('/api/compose', json={'tick_direction': 'out'}).json()
        assert c['tick_direction'] == 'out'

    def test_set_axis_labels(self, client):
        c = client.post('/api/compose', json={
            'xlabel': 'RA (J2000)', 'ylabel': 'Dec (J2000)',
        }).json()
        assert c['xlabel'] == 'RA (J2000)'
        assert c['ylabel'] == 'Dec (J2000)'

    def test_invalid_tick_direction(self, client):
        assert client.post('/api/compose', json={'tick_direction': 'sideways'}).status_code == 400

    def test_set_tick_direction_inout(self, client, test_fits_file):
        client.post('/api/panel/0/load', json={'path': test_fits_file})
        c = client.post('/api/compose', json={'tick_direction': 'inout'}).json()
        assert c['tick_direction'] == 'inout'
        r = client.get('/api/combined_plot.png')
        assert r.status_code == 200 and r.content[:4] == b'\x89PNG'

    def test_preview_swatch(self, client, test_fits_file):
        client.post('/api/panel/0/load', json={'path': test_fits_file})
        client.post('/api/panel/1/load', json={'path': test_fits_file})
        r = client.post('/api/preview_swatch', json={'delta_deg': 45})
        assert r.status_code == 200 and r.content[:4] == b'\x89PNG'

    def test_rotate_hues(self, client, test_fits_file):
        client.post('/api/panel/0/load', json={'path': test_fits_file})
        client.post('/api/panel/1/load', json={'path': test_fits_file})
        base = client.get('/api/state').json()['panels']
        colors = [p['color'] for p in base if p.get('loaded')]
        r = client.post('/api/rotate_hues', json={'delta_deg': 30, 'base_colors': colors})
        assert r.status_code == 200
        assert len(r.json()['result']['colors']) == 2

    def test_set_bare_plot(self, client):
        c = client.post('/api/compose', json={'bare_plot': True}).json()
        assert c['bare_plot'] is True

    def test_set_panel_preview_hd(self, client):
        c = client.post('/api/compose', json={'panel_preview_hd': True}).json()
        assert c['panel_preview_hd'] is True
        c2 = client.post('/api/compose', json={'panel_preview_hd': False}).json()
        assert c2['panel_preview_hd'] is False

    def test_set_combine_background(self, client):
        for bg in ('white', 'transparent', '#204060'):
            c = client.post('/api/compose', json={'combine_background': bg}).json()
            assert c['combine_background'] == bg

    def test_invalid_combine_background(self, client):
        assert client.post('/api/compose',
                           json={'combine_background': 'notacolor'}).status_code == 400

    def test_combined_preview_ryb(self, loaded_client):
        loaded_client.post('/api/compose', json={'combine_mode': 'ryb',
                                                 'combine_background': 'white'})
        r = loaded_client.get('/api/combined.png?preview=true&max_size=64')
        assert r.status_code == 200 and r.content[:4] == b'\x89PNG'

    def test_combined_transparent_has_alpha(self, loaded_client):
        from matplotlib import image as mpl_image
        import io
        loaded_client.post('/api/compose', json={'combine_background': 'transparent'})
        r = loaded_client.get('/api/combined.png')
        assert r.status_code == 200
        arr = mpl_image.imread(io.BytesIO(r.content))
        assert arr.shape[-1] == 4

    def test_export_script_ryb_and_background(self, loaded_client):
        loaded_client.post('/api/compose', json={'combine_mode': 'ryb',
                                                 'combine_background': 'white'})
        r = loaded_client.get('/api/export_script')
        assert 'combine_multicolor_alpha' in r.text and "mode='ryb'" in r.text
        compile(r.text, '<x>', 'exec')

    def test_export_script_rgb_background(self, loaded_client):
        loaded_client.post('/api/compose', json={'combine_background': 'transparent'})
        r = loaded_client.get('/api/export_script')
        assert 'composite_over_background' in r.text and 'background=None' in r.text
        compile(r.text, '<x>', 'exec')

    def test_combined_preview_with_lab(self, loaded_client):
        loaded_client.post('/api/compose', json={'combine_mode': 'lab', 'combine_blend': 'screen'})
        r = loaded_client.get('/api/combined.png?preview=true&max_size=64')
        assert r.status_code == 200
        assert r.content[:4] == b'\x89PNG'

    def test_export_script_colorspace(self, loaded_client):
        loaded_client.post('/api/compose', json={'combine_mode': 'hsv', 'combine_blend': 'mean'})
        r = loaded_client.get('/api/export_script')
        assert 'combine_multicolor_colorspace' in r.text
        assert "colorspace='hsv'" in r.text
        assert "blend='mean'" in r.text
        compile(r.text, '<x>', 'exec')

    def test_save_fits(self, loaded_client, tmp_path):
        import astropy.io.fits as pyfits
        path = str(tmp_path / 'out.fits')
        r = loaded_client.post('/api/save/fits', json={'path': path})
        assert r.status_code == 200
        hdr = pyfits.getheader(path)
        assert pyfits.getdata(path).shape == (3, 64, 64)
        assert hdr['MCFMODE'] == 'rgb'
        assert hdr['MCFTYPE'] == 'combined'

    def test_save_fits_lab_mode_annotated(self, loaded_client, tmp_path):
        import astropy.io.fits as pyfits
        loaded_client.post('/api/compose', json={'combine_mode': 'lab', 'combine_blend': 'screen'})
        path = str(tmp_path / 'lab.fits')
        r = loaded_client.post('/api/save/fits', json={'path': path, 'kind': 'combined'})
        assert r.status_code == 200
        hdr = pyfits.getheader(path)
        assert hdr['MCFMODE'] == 'lab'
        assert hdr['MCFBLEND'] == 'screen'

    def test_fits_targets_lists_layers(self, loaded_client):
        r = loaded_client.get('/api/fits_targets')
        assert r.status_code == 200
        targets = r.json()['targets']
        kinds = [t['kind'] for t in targets]
        assert kinds[0] == 'combined'
        assert 'data' in kinds and 'color' in kinds

    def test_save_component_data_fits(self, loaded_client, tmp_path):
        import astropy.io.fits as pyfits
        path = str(tmp_path / 'layer_data.fits')
        r = loaded_client.post('/api/save/fits',
                               json={'path': path, 'kind': 'data', 'index': 0})
        assert r.status_code == 200
        assert pyfits.getdata(path).ndim == 2
        assert pyfits.getheader(path)['MCFTYPE'] == 'component-data'

    def test_save_component_color_fits(self, loaded_client, tmp_path):
        import astropy.io.fits as pyfits
        path = str(tmp_path / 'layer_color.fits')
        r = loaded_client.post('/api/save/fits',
                               json={'path': path, 'kind': 'color', 'index': 0})
        assert r.status_code == 200
        assert pyfits.getdata(path).shape[0] == 3
        assert pyfits.getheader(path)['MCFTYPE'] == 'component-color'

    def test_save_component_unloaded_panel_errors(self, loaded_client, tmp_path):
        path = str(tmp_path / 'nope.fits')
        r = loaded_client.post('/api/save/fits',
                               json={'path': path, 'kind': 'data', 'index': 3})
        assert r.status_code == 400

    def test_save_image(self, loaded_client, tmp_path):
        import os
        path = str(tmp_path / 'out.png')
        r = loaded_client.post('/api/save/image', json={'path': path, 'dpi': 80})
        assert r.status_code == 200
        assert os.path.getsize(path) > 0

    def test_save_image_format_from_extension(self, loaded_client, tmp_path):
        import os
        path = str(tmp_path / 'out.jpg')
        r = loaded_client.post('/api/save/image', json={'path': path, 'dpi': 80})
        assert r.status_code == 200
        assert os.path.isfile(path)
        assert os.path.getsize(path) > 0

    def test_params_text(self, loaded_client):
        r = loaded_client.get('/api/params')
        assert 'image1' in r.text

    def test_export_script(self, loaded_client):
        r = loaded_client.get('/api/export_script')
        assert 'mcf.combine_multicolor' in r.text
        compile(r.text, '<x>', 'exec')

    def test_browse(self, client, tmp_path, test_fits_file):
        import os
        r = client.get('/api/browse?path=' + os.path.dirname(test_fits_file))
        data = r.json()
        assert 'synthetic.fits' in data['files']


@pytest.fixture
def mismatched_fits(tmp_path):
    """A second FITS file with a different grid than test_fits_file (64x64)."""
    import astropy.io.fits as pyfits
    from conftest import make_test_header, make_test_data
    hdr = make_test_header(100, 80)
    hdr['CRVAL1'] = 150.02
    path = tmp_path / 'mismatch.fits'
    pyfits.writeto(str(path), make_test_data(100, 80), hdr)
    return str(path)


class TestAlignment:
    def test_grid_status_single_aligned(self, loaded_client):
        r = loaded_client.get('/api/grid_status')
        assert r.status_code == 200
        s = r.json()
        assert s['aligned'] is True
        assert 'reproject_available' in s
        assert 'reference' in s['align_targets'] or 'reference' in s['align_targets']

    def test_grid_status_detects_mismatch(self, loaded_client, mismatched_fits):
        loaded_client.post('/api/panel/1/load', json={'path': mismatched_fits})
        s = loaded_client.get('/api/grid_status').json()
        assert s['aligned'] is False
        assert 1 in s['mismatched']

    def test_invalid_align_target(self, loaded_client):
        r = loaded_client.post('/api/align', json={'target': 'bogus'})
        assert r.status_code == 400

    @pytest.mark.skipif(
        __import__('multicolorfits').reproject_available() is False,
        reason='reproject not installed')
    def test_align_to_reference(self, loaded_client, mismatched_fits):
        loaded_client.post('/api/panel/1/load', json={'path': mismatched_fits})
        r = loaded_client.post('/api/align', json={'target': 'reference'})
        assert r.status_code == 200, r.text
        body = r.json()
        assert body['result']['ok'] is True
        # after aligning, grids match
        assert loaded_client.get('/api/grid_status').json()['aligned'] is True

    @pytest.mark.skipif(
        __import__('multicolorfits').reproject_available() is False,
        reason='reproject not installed')
    def test_align_to_icrs(self, loaded_client, mismatched_fits):
        loaded_client.post('/api/panel/1/load', json={'path': mismatched_fits})
        r = loaded_client.post('/api/align', json={'target': 'icrs'})
        assert r.status_code == 200, r.text
        assert loaded_client.get('/api/grid_status').json()['aligned'] is True


class TestPalettes:
    def test_list_palettes(self, client):
        r = client.get('/api/palettes')
        assert r.status_code == 200
        names = [p['name'] for p in r.json()['palettes']]
        assert 'pob' in names and 'rgb' in names
        assert 'tol' in names
        menu = r.json().get('menu', [])
        assert menu and menu[0]['group']
        for p in r.json()['palettes']:
            assert all(c.startswith('#') for c in p['colors'])

    def test_apply_hue_pattern(self, client, test_fits_file):
        client.post('/api/panel/0/load', json={'path': test_fits_file})
        client.post('/api/panel/1/load', json={'path': test_fits_file})
        client.post('/api/panel/2/load', json={'path': test_fits_file})
        r = client.post('/api/apply_palette', json={'name': 'triad'})
        assert r.status_code == 200
        assert len(r.json()['result']['colors']) == 3

    def test_apply_palette(self, client, test_fits_file):
        client.post('/api/panel/0/load', json={'path': test_fits_file})
        client.post('/api/panel/1/load', json={'path': test_fits_file})
        r = client.post('/api/apply_palette', json={'name': 'ryb'})
        assert r.status_code == 200, r.text
        body = r.json()
        assert len(body['result']['colors']) == 2
        assert body['state']['panels'][0]['color'] == '#C11B17'
        assert 'colorblind' in body

    def test_apply_perceptual(self, client, test_fits_file):
        client.post('/api/panel/0/load', json={'path': test_fits_file})
        r = client.post('/api/apply_palette', json={'name': 'perceptual'})
        assert r.status_code == 200
        assert len(r.json()['result']['colors']) == 1

    def test_invalid_palette(self, loaded_client):
        r = loaded_client.post('/api/apply_palette', json={'name': 'nope'})
        assert r.status_code == 400

    def test_colorblind_endpoint(self, loaded_client):
        r = loaded_client.get('/api/colorblind')
        assert r.status_code == 200
        body = r.json()
        assert 'ok' in body and 'kinds' in body

    def test_colorblind_flags_confusable(self, client, test_fits_file):
        client.post('/api/panel/0/load', json={'path': test_fits_file})
        client.post('/api/panel/1/load', json={'path': test_fits_file})
        client.post('/api/panel/0/params', json={'color': '#FF0000'})
        client.post('/api/panel/1/params', json={'color': '#FA0505'})
        body = client.get('/api/colorblind').json()
        assert body['ok'] is False


class TestLegend:
    def test_set_show_legend(self, loaded_client):
        r = loaded_client.post('/api/compose', json={'show_legend': True,
                                                     'legend_loc': 'lower left'})
        assert r.status_code == 200
        body = r.json()
        assert body['show_legend'] is True
        assert body['legend_loc'] == 'lower left'

    def test_set_panel_label(self, loaded_client):
        r = loaded_client.post('/api/panel/0/params', json={'label': 'H-alpha'})
        assert r.status_code == 200
        assert r.json()['label'] == 'H-alpha'

    def test_label_defaults_to_filename(self, client, test_fits_file):
        r = client.post('/api/panel/0/load', json={'path': test_fits_file})
        assert r.json()['label'] == 'synthetic'

    def test_legend_in_state(self, client):
        s = client.get('/api/state').json()
        assert 'show_legend' in s['compose']
        assert 'legend_loc' in s['compose']


class TestSwatchAndBandLabels:
    def test_set_combo_swatch(self, loaded_client):
        r = loaded_client.post('/api/compose', json={
            'show_combo_swatch': True, 'combo_swatch_loc': 'upper left',
            'combo_swatch_labels': True,
            'combo_swatch_inset_scale': 0.32, 'combo_swatch_size': 256})
        assert r.status_code == 200
        body = r.json()
        assert body['show_combo_swatch'] is True
        assert body['combo_swatch_loc'] == 'upper left'
        assert body['combo_swatch_labels'] is True
        assert body['combo_swatch_inset_scale'] == pytest.approx(0.32)
        assert body['combo_swatch_size'] == 256

    def test_set_band_labels(self, loaded_client):
        r = loaded_client.post('/api/compose', json={
            'show_band_labels': True, 'band_labels_loc': 'lower right'})
        assert r.status_code == 200
        body = r.json()
        assert body['show_band_labels'] is True
        assert body['band_labels_loc'] == 'lower right'

    def test_swatch_fields_in_state(self, client):
        c = client.get('/api/state').json()['compose']
        for k in ('show_combo_swatch', 'combo_swatch_loc', 'combo_swatch_labels',
                  'combo_swatch_inset_scale', 'combo_swatch_size',
                  'show_band_labels', 'band_labels_loc'):
            assert k in c

    def test_combined_figure_with_annotations(self, loaded_client):
        loaded_client.post('/api/compose', json={
            'show_combo_swatch': True, 'show_band_labels': True})
        r = loaded_client.get('/api/combined_plot.png')
        assert r.status_code == 200
        assert r.headers['content-type'] == 'image/png'


class TestOverlays:
    def test_set_overlay_flags(self, loaded_client):
        r = loaded_client.post('/api/compose', json={
            'show_compass': True, 'show_beam': True, 'show_scale_bar': True})
        assert r.status_code == 200
        body = r.json()
        assert body['show_compass'] is True
        assert body['show_beam'] is True
        assert body['show_scale_bar'] is True

    def test_overlay_fields_in_state(self, client):
        state = client.get('/api/state').json()
        for k in ('show_compass', 'show_beam', 'show_scale_bar'):
            assert k in state['compose']
        assert 'overlays_available' in state
        assert isinstance(state['overlays_available'], bool)

    def test_combined_plot_with_overlays(self, loaded_client):
        loaded_client.post('/api/compose', json={'show_compass': True, 'show_scale_bar': True})
        r = loaded_client.get('/api/combined_plot.png')
        assert r.status_code == 200
        assert r.headers['content-type'] == 'image/png'
