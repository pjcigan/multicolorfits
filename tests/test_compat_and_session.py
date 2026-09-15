"""Backward-compat API surface and the shared McfSession controller."""

import numpy as np
import pytest

import multicolorfits as mcf
from multicolorfits import PanelState, ComposeState, McfSession
from conftest import make_test_header, make_test_data

# Every public name from the v2.1.3 flat module that scripts may rely on
V2_PUBLIC_API = [
    'scaling_fns',
    'adjust_gamma', 'drawProgressBar', 'nanpercofscore',
    'force_hdr_to_2D', 'force_hdr_to_3D', 'force_hdr_floats',
    'deg2dms', 'dms2deg', 'deg2hour', 'hour2deg',
    'dec2sex', 'sex2dec', 'angulardistance',
    'getcdelts', 'convsky2pix', 'convpix2sky', 'makesimpleheader',
    'getcdmatrix', 'getdegperpix', 'getasecperpix', 'getsteradperpix',
    'squeeze_image', 'header_coord_grids',
    'McfHeader', 'as_mcfheader', 'save_header', 'load_header',
    'beampars_asec_fromhdr', 'pixperbeam_from_hdr',
    'reproject2D', 'reproject3D',
    'cropfits2D', 'cropfits3D', 'cropfits2D_coords', 'cropfits3D_coords',
    'smooth_image',
    'hex_to_rgb', 'rgb_to_hex', 'hex_to_hsv',
    'rgb_to_hsv_vectorized', 'hsv_to_rgb_vectorized',
    'hexinv', 'to_hex',
    'greyRGBize_image', 'colorize_image_direct_rgb', 'colorize_image',
    'combine_multicolor',
    'plotsinglemulticolorRGB', 'comparemulticolorRGB_pureRGB', 'saveRGBfits',
    'gui', 'mcf_gui',
    'overlays_available',
]

V3_OPTIONAL_API = ['overlays']


class TestCompat:
    @pytest.mark.parametrize('name', V2_PUBLIC_API)
    def test_name_present(self, name):
        assert hasattr(mcf, name), 'v2.x public API name %r missing from package' % name

    @pytest.mark.parametrize('name', V3_OPTIONAL_API)
    def test_optional_submodule_present(self, name):
        assert hasattr(mcf, name), 'optional submodule %r missing' % name

    def test_no_gui_toolkit_imported(self):
        """Importing multicolorfits must not drag in traits/pyqt/pyside/wx.

        Checked in a fresh subprocess so the result is independent of any GUI
        toolkit that other test modules (e.g. the Qt GUI tests) may import.
        """
        import subprocess
        import sys
        code = (
            "import multicolorfits, sys\n"
            "bad = [m for m in ('traits','traitsui','pyface','PyQt5','PyQt6','PySide6','wx')\n"
            "       if m in sys.modules]\n"
            "print(','.join(bad))\n"
            "sys.exit(1 if bad else 0)\n"
        )
        r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
        assert r.returncode == 0, \
            'GUI toolkit imported at package import time: %s' % (r.stdout.strip() or r.stderr.strip())

    def test_scaling_fns_keys(self):
        assert set(mcf.scaling_fns.keys()) == {'linear', 'sqrt', 'squared', 'log', 'power', 'sinh', 'asinh'}

    def test_version(self):
        assert mcf.__version__ == '3.1.0'


class TestLegacyAliases:
    """The central old->new rename map in multicolorfits.compat."""

    def test_every_alias_resolves(self):
        for old, new in mcf.LEGACY_ALIASES.items():
            assert hasattr(mcf, old), 'legacy alias %r missing' % old
            assert hasattr(mcf, new), 'new name %r missing' % new

    def test_alias_docstrings_flag_deprecation(self):
        for old in mcf.LEGACY_ALIASES:
            obj = getattr(mcf, old)
            if callable(obj):
                assert 'Deprecated alias' in (obj.__doc__ or ''), old

    def test_non_callable_aliases_are_same_object(self):
        assert mcf.scaling_fns is mcf.stretch_functions

    def test_alias_results_match_new_names(self):
        data = make_test_data()
        assert np.allclose(mcf.greyRGBize_image(data, rescalefn='asinh'),
                           mcf.to_grey_rgb(data, rescalefn='asinh'))
        assert mcf.hexinv('#C11B17') == mcf.hex_complement('#C11B17')
        hdr = make_test_header()
        assert mcf.getdegperpix(hdr) == mcf.deg_per_pixel(hdr)
        assert np.allclose(mcf.convsky2pix(hdr, 150.0, 2.0),
                           mcf.sky_to_pixel(hdr, 150.0, 2.0))

    def test_v2_quickstart_runs_on_old_names(self):
        """The classic v2.x three-step recipe, verbatim on legacy names."""
        data = make_test_data()
        grey = mcf.greyRGBize_image(data, rescalefn='asinh', scaletype='perc', min_max=[0., 100.])
        colr = mcf.colorize_image(grey, '#C11B17', colorintype='hex')
        combined = mcf.combine_multicolor([colr], gamma=2.2)
        assert combined.shape == data.shape + (3,)
        assert combined.min() >= 0 and combined.max() <= 1


class TestPanelState:
    def test_load_fits(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        assert p.in_use
        assert p.data.shape == (64, 64)
        assert p.vmin == p.data_min
        assert p.vmax == p.data_max

    def test_set_data_drops_extra_axes(self, test_header):
        p = PanelState()
        cube = make_test_data()[np.newaxis, :, :]  # [1, ny, nx]
        hdr = test_header.copy()
        hdr['NAXIS'] = 3; hdr['NAXIS3'] = 1
        p.set_data(cube, hdr)
        assert p.data.ndim == 2
        assert p.header['NAXIS'] == 2

    def test_set_data_cleans_ghost_higher_axis_cards_from_2d_header(self, test_header):
        p = PanelState()
        hdr = test_header.copy()
        hdr['NAXIS3'] = 1
        hdr['CTYPE3'] = 'STOKES'
        hdr['CRVAL3'] = 1.0
        hdr['WCSAXES'] = 3
        p.set_data(make_test_data(), hdr)
        assert p.data.ndim == 2
        assert p.header['NAXIS'] == 2
        assert p.header['WCSAXES'] == 2
        assert 'NAXIS3' not in p.header
        assert 'CTYPE3' not in p.header

    def test_set_limits_updates_percentiles(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        med = float(np.nanmedian(p.data))
        p.set_limits(vmin=med)
        assert 40 < p.percent_min < 60  # median should be ~50th percentile

    def test_set_percentiles_updates_limits(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        p.set_percentiles(pmin=5, pmax=95)
        assert p.vmin > p.data_min
        assert p.vmax < p.data_max

    def test_zscale(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        p.apply_zscale()
        assert p.vmin < p.vmax

    def test_reset_minmax(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        p.set_percentiles(pmin=10, pmax=90)
        p.reset_minmax()
        assert p.vmin == p.data_min
        assert p.vmax == p.data_max
        assert p.percent_max >= 99.0

    def test_render_display_range(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        p.color = '#C11B17'
        disp = p.render_display(gamma=2.2)
        assert disp.shape == (64, 64, 3)
        assert disp.min() >= 0 and disp.max() <= 1

    def test_render_inverted(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        disp = p.render_display(gamma=2.2, inverse=True)
        # Inverted single image: background (low values) should be bright
        assert disp.mean() > 0.5

    def test_render_with_smoothing(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        p.smooth = True
        p.smooth_sigma = 2.0
        disp = p.render_display()
        assert disp.shape == (64, 64, 3)

    def test_header_string_roundtrip(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        s = p.header_string()
        assert 'CRVAL1' in s
        p.apply_header_string(s)
        assert np.isclose(p.header['CRVAL1'], 150.0)

    def test_clear(self, test_fits_file):
        p = PanelState()
        p.load_fits(test_fits_file)
        p.clear()
        assert not p.in_use
        assert p.data is None

    def test_render_without_data_raises(self):
        with pytest.raises(ValueError):
            PanelState().render_color_rgb()


class TestMcfSession:
    def _loaded_session(self, test_fits_file):
        s = McfSession()
        s.panels[0].load_fits(test_fits_file)
        s.panels[0].color = '#FF0000'
        s.panels[1].load_fits(test_fits_file)
        s.panels[1].color = '#0000FF'
        return s

    def test_render_combined(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        combined = s.render_combined()
        assert combined.shape == (64, 64, 3)
        assert combined.min() >= 0 and combined.max() <= 1

    def test_render_combined_inverse(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        combined = s.render_combined(inverse=True)
        # Inverse: background should be bright (white-ish)
        assert combined.mean() > 0.5

    def test_no_panels_raises(self):
        with pytest.raises(ValueError):
            McfSession().render_combined()

    def test_render_combined_white_background(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        s.compose.combine_background = 'white'
        combined = s.render_combined()
        assert combined.shape == (64, 64, 3)
        # Zero-signal corners should be (near) white.
        assert combined.mean() > 0.5

    def test_render_combined_transparent_is_rgba(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        s.compose.combine_background = 'transparent'
        combined = s.render_combined()
        assert combined.shape == (64, 64, 4)

    def test_render_combined_ryb(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        s.compose.combine_mode = 'ryb'
        s.compose.combine_background = 'white'
        combined = s.render_combined()
        assert combined.shape == (64, 64, 3)
        assert combined.min() >= 0 and combined.max() <= 1

    def test_save_rgb_fits_transparent_flattened(self, tmp_path, test_fits_file):
        import astropy.io.fits as pyfits
        s = self._loaded_session(test_fits_file)
        s.compose.combine_background = 'transparent'
        path = str(tmp_path / 'combined_tp.fits')
        s.save_rgb_fits(path)
        assert pyfits.getdata(path).shape == (3, 64, 64)   # alpha flattened out
        assert pyfits.getheader(path)['MCFBKG'].lower() == 'transparent'

    def test_export_script_background_and_ryb(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        s.compose.combine_mode = 'ryb'
        s.compose.combine_background = 'white'
        src = s.export_script()
        assert 'combine_multicolor_alpha' in src and "mode='ryb'" in src
        compile(src, '<x>', 'exec')
        s.compose.combine_mode = 'rgb'
        s.compose.combine_background = 'transparent'
        src = s.export_script()
        assert 'composite_over_background' in src and 'background=None' in src
        compile(src, '<x>', 'exec')

    def test_save_rgb_fits(self, tmp_path, test_fits_file):
        import astropy.io.fits as pyfits
        s = self._loaded_session(test_fits_file)
        path = str(tmp_path / 'combined.fits')
        s.save_rgb_fits(path)
        assert pyfits.getdata(path).shape == (3, 64, 64)

    def test_save_rgb_fits_records_mode(self, tmp_path, test_fits_file):
        import astropy.io.fits as pyfits
        s = self._loaded_session(test_fits_file)
        s.compose.combine_mode = 'lab'
        s.compose.combine_blend = 'max'
        path = str(tmp_path / 'combined_lab.fits')
        s.save_rgb_fits(path)
        hdr = pyfits.getheader(path)
        assert hdr['MCFMODE'] == 'lab'
        assert hdr['MCFBLEND'] == 'max'
        assert hdr['MCFTYPE'] == 'combined'

    def test_save_rgb_fits_has_provenance(self, tmp_path, test_fits_file):
        import astropy.io.fits as pyfits
        s = self._loaded_session(test_fits_file)
        path = str(tmp_path / 'combined.fits')
        s.save_rgb_fits(path)
        hdr = pyfits.getheader(path)
        blob = ' '.join(str(c) for c in hdr.cards).lower()
        assert 'multicolorfits' in blob

    def test_load_state_restores_overlay_and_swatch_offset(self, tmp_path, test_fits_file):
        s = self._loaded_session(test_fits_file)
        s.compose.show_compass = True
        s.compose.show_scale_bar = True
        s.compose.scale_bar_asec = 10.0
        s.compose.combo_swatch_label_offset = 0.45
        path = str(tmp_path / 'state.json')
        s.save_state(path)
        s2 = McfSession()
        s2.load_state(path)
        assert s2.compose.show_compass is True
        assert s2.compose.show_scale_bar is True
        assert s2.compose.scale_bar_asec == pytest.approx(10.0)
        assert s2.compose.combo_swatch_label_offset == pytest.approx(0.45)

    def test_load_state_restores_colors_and_percentiles(self, tmp_path, test_fits_file):
        s = McfSession()
        s.panels[0].load_fits(test_fits_file)
        s.panels[0].color = '#BE599E'
        s.panels[0].set_percentiles(pmin=5, pmax=95)
        s.panels[1].load_fits(test_fits_file)
        s.panels[1].color = '#DEA215'
        s.panels[1].set_percentiles(pmin=10, pmax=90)
        path = str(tmp_path / 'state.json')
        s.save_state(path)
        s2 = McfSession()
        s2.load_state(path)
        assert s2.panels[0].color == '#BE599E'
        assert s2.panels[1].color == '#DEA215'
        assert s2.panels[0].percent_min == pytest.approx(5.0)
        assert s2.panels[1].percent_max == pytest.approx(90.0)

    def test_load_state_resolves_basenames_via_sibling_dir(self, tmp_path, test_fits_file):
        import shutil
        d = tmp_path / 'data'
        d.mkdir()
        f0 = d / 'layer_a.fits'
        f1 = d / 'layer_b.fits'
        shutil.copy(test_fits_file, f0)
        shutil.copy(test_fits_file, f1)
        state = {
            'multicolorfits_state_version': 1,
            'panels': [
                {'loaded': True, 'filepath': str(f0), 'stretch': 'linear',
                 'color': '#FF0000', 'percent_min': 0, 'percent_max': 100},
                {'loaded': True, 'filepath': 'layer_b.fits', 'stretch': 'linear',
                 'color': '#0000FF', 'percent_min': 0, 'percent_max': 100},
                {'loaded': False, 'filepath': ''},
                {'loaded': False, 'filepath': ''},
            ],
            'compose': ComposeState().to_dict(),
        }
        s = McfSession()
        warnings = s.load_state(state)
        assert warnings == []
        assert s.panels[1].in_use
        assert s.panels[1].color == '#0000FF'

    def test_save_state_uses_relative_paths_for_colocated_fits(self, tmp_path, test_fits_file):
        import json
        import shutil
        d = tmp_path / 'share'
        d.mkdir()
        f0 = d / 'a.fits'
        f1 = d / 'b.fits'
        shutil.copy(test_fits_file, f0)
        shutil.copy(test_fits_file, f1)
        s = McfSession()
        s.panels[0].load_fits(str(f0))
        s.panels[1].load_fits(str(f1))
        json_path = d / 'session.json'
        s.save_state(str(json_path))
        with open(json_path) as fh:
            state = json.load(fh)
        assert state['panels'][0]['filepath'] == './a.fits'
        assert state['panels'][1]['filepath'] == './b.fits'
        # Relocate the share folder and still load via relative paths.
        moved = tmp_path / 'elsewhere'
        shutil.move(str(d), str(moved))
        s2 = McfSession()
        warnings = s2.load_state(str(moved / 'session.json'))
        assert warnings == []
        assert s2.panels[0].in_use and s2.panels[1].in_use

    def test_make_portable_path_keeps_external_absolute(self, tmp_path, test_fits_file):
        import os
        outside = os.path.abspath(test_fits_file)
        session_dir = str(tmp_path / 'share')
        os.makedirs(session_dir, exist_ok=True)
        assert McfSession._make_portable_path(outside, session_dir) == outside
        local = os.path.join(session_dir, 'local.fits')
        assert McfSession._make_portable_path(local, session_dir) == './local.fits'

    def test_to_dict_uses_absolute_filepath(self, test_fits_file):
        import os
        p = PanelState()
        p.load_fits(test_fits_file)
        d = p.to_dict()
        assert os.path.isabs(d['filepath'])
        assert d['filepath'] == os.path.abspath(test_fits_file)


    def test_resolve_fits_path(self, tmp_path, test_fits_file):
        import shutil
        d = tmp_path / 'fits'
        d.mkdir()
        target = d / 'target.fits'
        shutil.copy(test_fits_file, target)
        assert McfSession._resolve_fits_path(str(target)) == str(target.resolve())
        assert McfSession._resolve_fits_path('target.fits', search_dirs=[str(d)]) == str(target.resolve())

    def test_save_component_data(self, tmp_path, test_fits_file):
        import astropy.io.fits as pyfits
        s = self._loaded_session(test_fits_file)
        path = str(tmp_path / 'layer1.fits')
        s.save_component_fits(0, path, kind='data')
        data = pyfits.getdata(path)
        assert data.ndim == 2 and data.shape == (64, 64)
        assert pyfits.getheader(path)['MCFTYPE'] == 'component-data'

    def test_save_component_color(self, tmp_path, test_fits_file):
        import astropy.io.fits as pyfits
        s = self._loaded_session(test_fits_file)
        path = str(tmp_path / 'layer1_color.fits')
        s.save_component_fits(0, path, kind='color')
        assert pyfits.getdata(path).shape == (3, 64, 64)
        hdr = pyfits.getheader(path)
        assert hdr['MCFTYPE'] == 'component-color'
        assert hdr['MCFCOLOR'] == '#FF0000'

    def test_save_component_bad_kind(self, tmp_path, test_fits_file):
        s = self._loaded_session(test_fits_file)
        with pytest.raises(ValueError):
            s.save_component_fits(0, str(tmp_path / 'x.fits'), kind='bogus')

    def test_save_component_unloaded_panel(self, tmp_path, test_fits_file):
        s = self._loaded_session(test_fits_file)
        with pytest.raises(ValueError):
            s.save_component_fits(3, str(tmp_path / 'x.fits'), kind='data')

    def test_fits_save_targets(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        targets = s.fits_save_targets()
        assert targets[0]['kind'] == 'combined'
        # two loaded layers => 1 combined + 2*(data,color) = 5 entries
        assert len(targets) == 5
        assert {t['kind'] for t in targets} == {'combined', 'data', 'color'}

    def test_params_text(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        text = s.params_text()
        assert 'image1' in text
        assert 'image2' in text
        assert "'#FF0000'" in text
        assert 'gamma = 2.2' in text

    def test_export_script_is_valid_python(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        script = s.export_script()
        compile(script, '<exported>', 'exec')  # Must be syntactically valid
        assert 'mcf.to_grey_rgb' in script
        assert 'mcf.combine_multicolor' in script
        assert "'#FF0000'" in script

    def test_exported_script_reproduces_combined(self, tmp_path, test_fits_file):
        """Running the exported script must reproduce render_combined() output."""
        s = self._loaded_session(test_fits_file)
        expected = s.render_combined()
        script = s.export_script()
        # Strip the final plotting call (matplotlib not needed for equivalence check)
        lines = [ln for ln in script.split('\n') if not ln.startswith('mcf.plot_combined_rgb')]
        namespace = {}
        exec('\n'.join(lines), namespace)
        assert np.allclose(namespace['combined'], expected, atol=1e-6)

    def test_to_dict(self, test_fits_file):
        s = self._loaded_session(test_fits_file)
        d = s.to_dict()
        assert len(d['panels']) == 4
        assert d['panels'][0]['loaded'] is True
        assert d['panels'][2]['loaded'] is False
        assert d['compose']['gamma'] == 2.2
