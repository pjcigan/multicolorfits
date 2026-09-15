"""Tests for header tidying, north-up rotation, overlap crop, and beam matching."""

import warnings

import numpy as np
import pytest
import astropy.io.fits as pyfits

import multicolorfits as mcf
from conftest import make_test_header, make_test_data


def _pc_header(theta_deg=40.0, flipped=False):
    hdr = make_test_header()
    th = np.deg2rad(theta_deg)
    c, s = np.cos(th), np.sin(th)
    # Keep the CDELT scales; put the rotation in PC.
    scale = abs(float(hdr['CDELT1']))
    hdr['CDELT1'] = -scale
    hdr['CDELT2'] = scale
    hdr['PC1_1'] = c
    hdr['PC1_2'] = s
    hdr['PC2_1'] = -s if not flipped else s
    hdr['PC2_2'] = c if not flipped else -c
    return hdr


class TestTidyHeader:
    def test_drops_conflicting_cd_and_crota(self, test_header):
        hdr = test_header.copy()
        hdr['PC1_1'] = 1.0
        hdr['PC2_2'] = 1.0
        hdr['CD1_1'] = hdr['CDELT1']
        hdr['CD2_2'] = hdr['CDELT2']
        hdr['CROTA2'] = 12.0
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            out = mcf.tidy_header(hdr, warn_flip=False)
        assert 'CD1_1' not in out
        assert 'CROTA2' not in out
        assert 'PC1_1' in out
        assert any('CD' in str(w.message) for w in caught)
        assert mcf.wcs_is_flipped(test_header) is False

    def test_flip_warning(self):
        hdr = _pc_header(flipped=True)
        if not mcf.wcs_is_flipped(hdr):
            pytest.skip('constructed header was not mirrored')
        with pytest.warns(UserWarning, match='mirrored'):
            mcf.tidy_header(hdr)


class TestRotatedHeader:
    def test_north_up_keeps_equatorial_frame(self, test_header):
        out = mcf.make_north_up_header(test_header)
        assert str(out['CTYPE1']).startswith('RA')
        assert str(out['CTYPE2']).startswith('DEC')
        assert float(out['CDELT1']) < 0
        assert float(out['CDELT2']) > 0

    def test_galactic_stays_galactic(self):
        hdr = make_test_header()
        hdr['CTYPE1'] = 'GLON-TAN'
        hdr['CTYPE2'] = 'GLAT-TAN'
        hdr['CRVAL1'] = 120.0
        hdr['CRVAL2'] = -10.0
        out = mcf.make_north_up_header(hdr)
        assert str(out['CTYPE1']).startswith('GLON')
        assert str(out['CTYPE2']).startswith('GLAT')
        assert str(out['CTYPE1']).startswith('RA') is False

    def test_oversample_shrinks_pixels(self, test_header):
        out = mcf.make_rotated_header(test_header, oversample=2)
        assert abs(float(out['CDELT1'])) == pytest.approx(abs(float(test_header['CDELT1'])) / 2)

    def test_reproject_keeps_peak_near_center(self, test_data, test_header):
        pytest.importorskip('reproject')
        hdr = _pc_header(theta_deg=35.0)
        data, hdr2 = mcf.reproject_north_up(test_data, hdr, order=1)
        assert data.shape == (hdr2['NAXIS2'], hdr2['NAXIS1'])
        assert str(hdr2['CTYPE1']).startswith('RA')
        finite = np.isfinite(data)
        assert finite.any()
        peak = np.unravel_index(np.nanargmax(data), data.shape)
        cy, cx = (data.shape[0] - 1) / 2.0, (data.shape[1] - 1) / 2.0
        assert abs(peak[0] - cy) < data.shape[0] * 0.25
        assert abs(peak[1] - cx) < data.shape[1] * 0.25


class TestCropToOverlap:
    def test_trims_nan_border(self, test_header):
        a = np.ones((20, 30))
        b = np.ones((20, 30))
        a[:4, :] = np.nan
        b[:, -3:] = np.nan
        cropped, hdr = mcf.crop_to_overlap([a, b], test_header)
        assert cropped[0].shape[0] == 16
        assert cropped[0].shape[1] == 27
        assert np.isfinite(cropped[0]).all()
        assert np.isfinite(cropped[1]).all()
        assert int(hdr['NAXIS1']) == 27


class TestBeams:
    def test_identity_kernel_is_zero(self):
        maj, min_, pa = mcf.convolution_beam(5, 3, 20, 5, 3, 20)
        assert maj == 0 and min_ == 0

    def test_larger_target_has_positive_kernel(self):
        maj, min_, _pa = mcf.convolution_beam(2, 2, 0, 4, 4, 0)
        assert maj > 0 and min_ > 0

    def test_smaller_target_raises(self):
        with pytest.raises(ValueError):
            mcf.convolution_beam(4, 4, 0, 2, 2, 0)

    def test_match_updates_header(self, test_data, test_header):
        hdr = test_header.copy()
        hdr['BMAJ'] = 4.0 / 3600.0
        hdr['BMIN'] = 4.0 / 3600.0
        hdr['BPA'] = 0.0
        out, hdrout = mcf.match_beam(test_data, hdr, 8.0, 8.0, 0.0)
        assert out.shape == test_data.shape
        assert float(hdrout['BMAJ']) == pytest.approx(8.0 / 3600.0)
        assert np.isfinite(out).all()
