"""Tests for the stateful McfHeader wrapper."""

import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import pytest

import multicolorfits as mcf
from multicolorfits.core.mcfheader import McfHeader, as_mcfheader

from conftest import make_test_header


class TestConstruction:
    def test_from_astropy_header(self, test_header):
        h = McfHeader(test_header)
        assert h['CTYPE1'] == 'RA---TAN'
        assert h.shape == (test_header['NAXIS2'], test_header['NAXIS1'])

    def test_none_yields_placeholder(self):
        h = McfHeader(None)
        assert h.shape is None
        assert h.is_celestial is False

    def test_copy_is_independent(self, test_header):
        h = McfHeader(test_header)
        h2 = h.copy()
        h2['CRVAL1'] = 99.0
        assert h['CRVAL1'] != 99.0

    def test_does_not_mutate_source(self, test_header):
        original = test_header['CRVAL1']
        h = McfHeader(test_header)
        h['CRVAL1'] = 12.3
        assert test_header['CRVAL1'] == original

    def test_exported_top_level(self):
        assert hasattr(mcf, 'McfHeader')
        assert hasattr(mcf, 'as_mcfheader')

    def test_as_mcfheader_passthrough(self, test_header):
        h = McfHeader(test_header)
        assert as_mcfheader(h) is h
        assert isinstance(as_mcfheader(test_header), McfHeader)


class TestNormalization:
    def test_string_wcs_cards_coerced_to_float(self):
        hdr = make_test_header()
        hdr['CDELT1'] = '-0.0002777778'  # stringly typed, as some writers emit
        hdr['CRVAL1'] = '150.0'
        h = McfHeader(hdr)
        assert isinstance(h['CDELT1'], float)
        assert isinstance(h['CRVAL1'], float)

    def test_forces_2d(self):
        hdr = make_test_header()
        hdr['NAXIS'] = 3
        hdr['NAXIS3'] = 1
        hdr['CTYPE3'] = 'FREQ'
        hdr['CRVAL3'] = 1.4e9
        h = McfHeader(hdr)
        assert int(h['NAXIS']) == 2
        assert 'NAXIS3' not in h
        assert 'CTYPE3' not in h

    def test_forces_2d_for_nominally_2d_header_with_ghost_axis_cards(self):
        hdr = make_test_header()
        hdr['NAXIS3'] = 1
        hdr['CTYPE3'] = 'STOKES'
        hdr['CRVAL3'] = 1.0
        hdr['WCSAXES'] = 3
        h = McfHeader(hdr)
        assert int(h['NAXIS']) == 2
        assert int(h['WCSAXES']) == 2
        assert 'NAXIS3' not in h
        assert 'CTYPE3' not in h

    def test_float_wcsaxes_coerced_to_int_without_warning(self):
        import warnings
        from astropy.wcs.wcs import FITSFixedWarning

        hdr = make_test_header()
        hdr['WCSAXES'] = 2.0
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always', FITSFixedWarning)
            h = McfHeader(hdr)
            _ = h.wcs
        assert h['WCSAXES'] == 2
        assert isinstance(h['WCSAXES'], int)
        assert not [w for w in caught if issubclass(w.category, FITSFixedWarning)]

    def test_normalize_false_preserves_strings(self):
        hdr = make_test_header()
        hdr['CRVAL1'] = '150.0'
        h = McfHeader(hdr, normalize=False, force_2d=False)
        assert h['CRVAL1'] == '150.0'


class TestDerivedState:
    def test_wcs_cached(self, test_header):
        h = McfHeader(test_header)
        w1 = h.wcs
        w2 = h.wcs
        assert w1 is w2
        assert isinstance(w1, pywcs.WCS)

    def test_wcs_matches_astropy(self, test_header):
        h = McfHeader(test_header)
        ref = pywcs.WCS(test_header).celestial
        assert np.allclose(h.wcs.wcs.crval, ref.wcs.crval)

    def test_frame(self, test_header):
        assert McfHeader(test_header).frame == 'fk5'

    def test_pixscale_asec(self, test_header):
        h = McfHeader(test_header)
        assert h.pixscale_asec == pytest.approx(1.0, rel=1e-3)

    def test_is_celestial(self, test_header):
        assert McfHeader(test_header).is_celestial is True

    def test_cache_invalidated_on_setitem(self, test_header):
        h = McfHeader(test_header)
        _ = h.wcs
        h['CRVAL1'] = 200.0
        assert h.wcs.wcs.crval[0] == pytest.approx(200.0)


class TestSerialization:
    def test_roundtrip_tostring_fromstring(self, test_header):
        h = McfHeader(test_header)
        s = h.tostring()
        h2 = McfHeader.fromstring(s, normalize=False, force_2d=False)
        assert h2['CTYPE1'] == 'RA---TAN'

    def test_astropy_property_is_header(self, test_header):
        h = McfHeader(test_header)
        assert isinstance(h.astropy, pyfits.Header)


class TestMinimalWcsHeader:
    def test_strips_non_wcs_cards(self, test_header):
        hdr = make_test_header()
        hdr['HISTORY'] = 'processed by some pipeline'
        hdr['COMMENT'] = 'a chatty comment'
        hdr['BUNIT'] = 'Jy/beam'
        mini = McfHeader(hdr).minimal_wcs_header()
        assert 'HISTORY' not in mini
        assert 'COMMENT' not in mini
        assert 'BUNIT' not in mini

    def test_keeps_grid_defining_cards(self, test_header):
        mini = McfHeader(test_header).minimal_wcs_header()
        assert int(mini['NAXIS']) == 2
        assert mini['NAXIS1'] == test_header['NAXIS1']
        assert mini['NAXIS2'] == test_header['NAXIS2']
        assert mini['CTYPE1'].startswith('RA')

    def test_defines_same_grid(self, test_header):
        h = McfHeader(test_header)
        assert h.matches_grid(McfHeader(h.minimal_wcs_header()))


class TestHeaderPortability:
    def test_save_load_roundtrip(self, test_header, tmp_path):
        p = tmp_path / 'grid.hdr'
        h = McfHeader(test_header)
        mcf.save_header(h, str(p))
        h2 = mcf.load_header(str(p))
        assert isinstance(h2, McfHeader)
        assert h.matches_grid(h2)

    def test_save_minimal_is_compact(self, test_header, tmp_path):
        hdr = make_test_header()
        for i in range(20):
            hdr['HISTORY'] = 'line %d of provenance' % i
        full_path = tmp_path / 'full.hdr'
        mini_path = tmp_path / 'mini.hdr'
        h = McfHeader(hdr)
        mcf.save_header(h, str(full_path), minimal=False)
        mcf.save_header(h, str(mini_path), minimal=True)
        assert mini_path.read_text().count('\n') < full_path.read_text().count('\n')
        # minimal file still reloads to the same grid
        assert h.matches_grid(mcf.load_header(str(mini_path)))

    def test_totextfile_fromtextfile(self, test_header, tmp_path):
        p = tmp_path / 'grid.hdr'
        h = McfHeader(test_header)
        h.totextfile(str(p))
        h2 = McfHeader.fromtextfile(str(p))
        assert h.matches_grid(h2)

    def test_exported_top_level(self):
        assert hasattr(mcf, 'save_header')
        assert hasattr(mcf, 'load_header')


class TestGridMatch:
    def test_same_grid_matches(self, test_header):
        assert McfHeader(test_header).matches_grid(McfHeader(test_header))

    def test_different_shape_does_not_match(self, test_header):
        assert not McfHeader(test_header).matches_grid(McfHeader(make_test_header(nx=128, ny=128)))

    def test_different_crval_does_not_match(self, test_header):
        other = make_test_header()
        other['CRVAL1'] = 42.0
        assert not McfHeader(test_header).matches_grid(McfHeader(other))


class TestToFrame:
    def test_to_galactic(self, test_header):
        reproject = pytest.importorskip('reproject')  # noqa: F841
        h = McfHeader(test_header)
        gal = h.to_frame('galactic')
        assert isinstance(gal, McfHeader)
        assert gal.frame == 'galactic'


class TestSessionIntegration:
    def test_panel_header_is_mcfheader(self, test_fits_file):
        sess = mcf.McfSession()
        sess.panels[0].load_fits(test_fits_file)
        assert isinstance(sess.panels[0].header, McfHeader)

    def test_common_header_is_mcfheader(self, test_fits_file):
        sess = mcf.McfSession()
        sess.panels[0].load_fits(test_fits_file)
        assert isinstance(sess.common_header, McfHeader)
        assert sess.common_header.is_celestial

    def test_header_string_roundtrip(self, test_fits_file):
        sess = mcf.McfSession()
        p = sess.panels[0]
        p.load_fits(test_fits_file)
        s = p.header_string()
        assert 'CTYPE1' in s
        p.apply_header_string(s)
        assert isinstance(p.header, McfHeader)
        assert p.header['CTYPE1'] == 'RA---TAN'

    def test_save_rgb_fits_still_works(self, test_fits_file, tmp_path):
        sess = mcf.McfSession()
        sess.panels[0].load_fits(test_fits_file)
        sess.panels[0].color = '#FF0000'
        out = tmp_path / 'combined.fits'
        sess.save_rgb_fits(str(out))
        assert out.exists()
