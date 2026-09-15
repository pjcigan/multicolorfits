"""Tests for palettes, pipeline helpers, session colorspace, and provenance."""

import json
import os

import numpy as np
import pytest
import astropy.io.fits as pyfits

import multicolorfits as mcf
from multicolorfits import McfSession


class TestPalettes:
    def test_list_and_get(self):
        names = mcf.list_palettes()
        assert 'pob' in names
        assert len(mcf.get_palette('rgb')) == 3

    def test_suggest_colors_count(self):
        colors = mcf.suggest_colors(4)
        assert len(colors) == 4
        assert all(c.startswith('#') for c in colors)

    def test_pob_palette_passes_cvd_check(self):
        ok, failures = mcf.palettes.check_palette_colorblind(mcf.get_palette('pob'))
        assert ok is True
        assert failures == []

    def test_cvd_check_flags_similar_colors(self):
        ok, failures = mcf.palettes.check_palette_colorblind(['#FF0000', '#FF1010'])
        assert ok is False
        assert len(failures) == 1

    def test_get_palette_n_truncates(self):
        assert len(mcf.get_palette('rgb', n=2)) == 2


class TestPipeline:
    def test_combine_layers_rgb(self, test_data):
        d2 = test_data * 0.7
        out = mcf.combine_layers([test_data, d2], ['#FF0000', '#0000FF'])
        assert out.shape == test_data.shape + (3,)
        assert 0 <= out.min() and out.max() <= 1

    def test_combine_layers_lab(self, test_data):
        out = mcf.combine_layers([test_data, test_data * 0.5],
                                 ['#C11B17', '#4CC417'], colorspace='lab')
        assert out.shape == test_data.shape + (3,)

    def test_combine_from_files(self, test_fits_file, tmp_path):
        pytest.importorskip('reproject')
        path2 = str(tmp_path / 'copy.fits')
        data, hdr = pyfits.getdata(test_fits_file, header=True)
        pyfits.writeto(path2, data, hdr, overwrite=True)
        combined, hdr_out = mcf.combine_from_files(
            [test_fits_file, path2], ['#FF0000', '#0000FF'], align_reference=0)
        assert combined.shape[2] == 3
        assert hdr_out is not None


class TestFloat32Preview:
    def test_preview_returns_float32(self, test_data):
        out = mcf.combine_layers([test_data, test_data * 0.6],
                                 ['#FF0000', '#0000FF'], preview=True)
        assert out.dtype == np.float32

    def test_preview_default_is_float64(self, test_data):
        out = mcf.combine_layers([test_data], ['#00FF00'])
        assert out.dtype == np.float64

    def test_preview_matches_full_precision(self, test_data):
        colors = ['#BE599E', '#DEA215', '#77C0F9']
        layers = [test_data, test_data * 0.7, test_data * 0.4]
        f64 = mcf.combine_layers(layers, colors)
        f32 = mcf.combine_layers(layers, colors, preview=True)
        assert np.allclose(f64, f32.astype(np.float64), atol=1e-4)

    def test_preview_lab_colorspace(self, test_data):
        out = mcf.combine_layers([test_data, test_data * 0.5],
                                 ['#C11B17', '#4CC417'], colorspace='lab', preview=True)
        assert out.dtype == np.float32
        assert np.isfinite(out).all()

    def test_explicit_dtype_overrides_preview(self, test_data):
        out = mcf.combine_layers([test_data], ['#00FF00'], preview=True,
                                 dtype=np.float64)
        assert out.dtype == np.float64

    def test_combine_multicolor_dtype(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear', dtype=np.float32)
        assert grey.dtype == np.float32
        colr = mcf.colorize_image(grey, '#FF0000', colorintype='hex', dtype=np.float32)
        combined = mcf.combine_multicolor([colr], gamma=2.2, dtype=np.float32)
        assert combined.dtype == np.float32
        assert combined.min() >= 0.0 and combined.max() <= 1.0


class TestSessionColorspace:
    def test_render_combined_lab(self, test_fits_file):
        s = McfSession()
        s.panels[0].load_fits(test_fits_file)
        s.panels[0].color = '#E8A04C'
        s.compose.combine_mode = 'lab'
        s.compose.combine_blend = 'screen'
        out = s.render_combined()
        assert out.shape[2] == 3
        assert out.max() <= 1.0


class TestSessionPreview:
    def test_render_combined_preview_downsamples_and_float32(self, test_fits_file):
        s = McfSession()
        s.panels[0].load_fits(test_fits_file)
        s.panels[0].color = '#FF0000'
        full = s.render_combined()
        prev = s.render_combined(preview=True, max_size=16)
        assert prev.dtype == np.float32
        assert max(prev.shape[:2]) <= 16
        assert full.shape[2] == prev.shape[2] == 3

    def test_render_combined_preview_close_to_full_at_same_size(self, test_fits_file):
        # With a max_size >= image size, no downsampling: preview ~ full precision
        s = McfSession()
        s.panels[0].load_fits(test_fits_file)
        s.panels[0].color = '#4CC417'
        full = s.render_combined()
        prev = s.render_combined(preview=True, max_size=4096)
        assert prev.shape == full.shape
        assert np.allclose(full, prev.astype(np.float64), atol=1e-4)


class TestStateSaveRestore:
    def test_roundtrip(self, test_fits_file, tmp_path):
        s = McfSession()
        s.panels[0].load_fits(test_fits_file)
        s.panels[0].color = '#C11B17'
        s.panels[0].stretch = 'asinh'
        s.compose.gamma = 1.8
        path = str(tmp_path / 'state.json')
        s.save_state(path)
        s2 = McfSession()
        warnings = s2.load_state(path)
        assert warnings == []
        assert s2.panels[0].in_use
        assert s2.panels[0].color == '#C11B17'
        assert s2.compose.gamma == 1.8


class TestProvenance:
    def test_annotate_adds_history(self, test_header):
        hdr = test_header.copy()
        mcf.annotate_provenance_header(hdr, {'panels': [{'loaded': True, 'stretch': 'linear',
                                                          'vmin': 0, 'vmax': 1, 'color': '#FF0000'}],
                                           'compose': {'gamma': 2.2, 'combine_mode': 'lab'}})
        hist = '\n'.join(str(c) for c in hdr.cards if c[0] == 'HISTORY')
        assert 'multicolorfits' in hist
        assert 'layer 1' in hist

    def test_annotate_accepts_note_string(self, test_header):
        hdr = test_header.copy()
        mcf.annotate_provenance_header(hdr, 'f560w north-up crop')
        hist = '\n'.join(str(c[1]) for c in hdr.cards if c[0] == 'HISTORY')
        assert 'f560w north-up crop' in hist

    def test_save_combined(self, tmp_path, test_data, test_header):
        combined = mcf.combine_layers([test_data], ['#00FF00'])
        path = str(tmp_path / 'out.fits')
        mcf.save_combined(path, combined, test_header)
        assert os.path.isfile(path)
        assert pyfits.getdata(path).shape[0] == 3
