"""WCS/header helpers, cropping, FITS round-trip, and lazy reproject imports."""

import numpy as np
import astropy.io.fits as pyfits
import pytest

import multicolorfits as mcf
from conftest import make_test_header, make_test_data


class TestCoordConversions:
    def test_sex2dec_roundtrip(self):
        ra, dec = 150.123456, -30.654321
        sexstrings = mcf.dec2sex(ra, dec, as_string=True, decimal_places=4)
        ra2, dec2 = mcf.sex2dec(*sexstrings)
        assert np.isclose(ra, ra2, atol=1e-5)
        assert np.isclose(dec, dec2, atol=1e-5)

    def test_dec2sex_negative_zero_deg(self):
        # DEC between 0 and -1 must keep the minus sign
        out = mcf.dec2sex(10.0, -0.5, as_string=True)
        assert out[1].startswith('-00')

    def test_pix_sky_roundtrip(self, test_header):
        ra, dec = mcf.pixel_to_sky(test_header, 20, 40)
        x, y = mcf.sky_to_pixel(test_header, ra, dec, precise=True)
        assert np.isclose(x, 20, atol=0.01)
        assert np.isclose(y, 40, atol=0.01)

    def test_dec2sex_rollover_no_60_seconds(self):
        # Ported edge-case from skyplothelper: -12.7 should not print 60 sec.
        out = mcf.dec2sex(10.0, -12.7, as_string=True, decimal_places=2)
        assert out[1] == '-12:42:00.00'

    def test_sex2dec_parses_hms_dms_delimiters(self):
        ra, dec = mcf.sex2dec('05h34m31.9s', '-22d00m52.2s')
        assert np.isclose(ra, 83.63291666666667)
        assert np.isclose(dec, -22.0145)

    def test_deg2dms_dms2deg_roundtrip_negative_subdegree(self):
        v = -0.001
        back = mcf.dms2deg(mcf.deg2dms(v))
        assert back == pytest.approx(v, abs=1e-10)

    def test_angulardistance_zero(self):
        sep = mcf.angular_distance([10.0, -3.0], [10.0, -3.0])
        assert sep == pytest.approx(0.0, abs=1e-15)


class TestHeaderHelpers:
    def test_getcdelts(self, test_header):
        cd1, cd2 = mcf.get_cdelts(test_header)
        assert np.isclose(cd1, -0.0002777778)
        assert np.isclose(cd2, 0.0002777778)

    def test_getcdelts_from_cd_matrix(self):
        hdr = make_test_header()
        cdelt = hdr['CDELT1']
        del hdr['CDELT1'], hdr['CDELT2']
        hdr['CD1_1'] = cdelt; hdr['CD1_2'] = 0.
        hdr['CD2_1'] = 0.; hdr['CD2_2'] = -cdelt
        cd1, cd2 = mcf.get_cdelts(hdr)
        assert np.isclose(cd1, cdelt)
        assert np.isclose(cd2, -cdelt)

    def test_deg_per_pixel_and_arcsec(self, test_header):
        assert np.isclose(mcf.deg_per_pixel(test_header), 0.0002777778)
        assert np.isclose(mcf.arcsec_per_pixel(test_header), 1.0, atol=1e-4)

    def test_sterad_per_pixel(self, test_header):
        expected = (0.0002777778 * np.pi / 180.)**2
        assert np.isclose(mcf.sterad_per_pixel(test_header), expected)

    def test_force_header_2d(self):
        hdr = make_test_header()
        hdr['NAXIS'] = 3
        hdr['WCSAXES'] = 3
        hdr['NAXIS3'] = 10
        hdr['CRVAL3'] = 1.0e9
        hdr2 = mcf.force_header_2d(hdr)
        assert hdr2['NAXIS'] == 2
        assert hdr2['WCSAXES'] == 2
        assert 'NAXIS3' not in hdr2
        assert 'CRVAL3' not in hdr2

    def test_force_header_3d(self):
        hdr = make_test_header()
        hdr['NAXIS'] = 4
        hdr['WCSAXES'] = 4
        hdr['NAXIS3'] = 7
        hdr['NAXIS4'] = 2
        hdr['CRVAL4'] = 1.0
        hdr3 = mcf.force_header_3d(hdr)
        assert hdr3['NAXIS'] == 3
        assert hdr3['WCSAXES'] == 3
        assert 'NAXIS4' not in hdr3
        assert 'CRVAL4' not in hdr3

    def test_force_header_floats(self):
        hdr = make_test_header()
        hdr['CRPIX1'] = '32.0'  # string that should become float
        hdr['PC1_1'] = '1.0'    # ensure PC cards are coerced too
        mcf.force_header_floats(hdr)
        assert isinstance(hdr['CRPIX1'], float)
        assert isinstance(hdr['PC1_1'], float)

    def test_make_simple_header(self, test_header):
        simple = mcf.make_simple_header(test_header)
        assert simple['NAXIS'] == 2
        assert np.isclose(simple['CRVAL1'], test_header['CRVAL1'])
        assert np.isclose(simple['CRVAL2'], test_header['CRVAL2'])

    def test_getcdmatrix(self, test_header):
        cd1_1, cd1_2, cd2_1, cd2_2 = mcf.get_cd_matrix(test_header)
        assert np.isclose(cd1_1, test_header['CDELT1'])
        assert np.isclose(cd2_2, test_header['CDELT2'])
        assert np.isclose(cd1_2, 0.)

    def test_deg_per_pixel_unequal_returns_mean(self):
        hdr = make_test_header()
        hdr['CDELT1'] = -0.001
        hdr['CDELT2'] = 0.0015
        with pytest.warns(UserWarning):
            out = mcf.deg_per_pixel(hdr)
        assert out == pytest.approx(np.mean([0.001, 0.0015]))

    def test_beam_helpers(self, test_header):
        test_header['BMAJ'] = 2. / 3600.  # 2 arcsec
        test_header['BMIN'] = 1. / 3600.
        test_header['BPA'] = 45.
        bmaj, bmin, bpa = mcf.beam_params_arcsec(test_header)
        assert np.isclose(bmaj, 2.0)
        assert np.isclose(bmin, 1.0)
        assert bpa == 45.
        ppb = mcf.pixels_per_beam(test_header)
        assert ppb > 0


class TestCrop:
    def test_crop_image_shape(self, test_data, test_header):
        crop, crophdr = mcf.crop_image(test_data, test_header, [10, 29], [20, 39])
        assert crop.shape == (20, 20)
        assert crophdr['NAXIS1'] == 20
        assert crophdr['NAXIS2'] == 20

    def test_crop_image_wcs_preserved(self, test_data, test_header):
        """Sky coordinate of a given feature must be identical before and after crop."""
        crop, crophdr = mcf.crop_image(test_data, test_header, [10, 49], [10, 49])
        ra_orig, dec_orig = mcf.pixel_to_sky(test_header, 30, 30)
        ra_crop, dec_crop = mcf.pixel_to_sky(crophdr, 20, 20)
        assert np.isclose(ra_orig, ra_crop, atol=1e-8)
        assert np.isclose(dec_orig, dec_crop, atol=1e-8)

    def test_crop_image_negative_bounds_raise(self, test_data, test_header):
        with pytest.raises(Exception):
            mcf.crop_image(test_data, test_header, [-5, 20], [0, 20])

    def test_crop_image_sky(self, test_data, test_header):
        ra, dec = mcf.pixel_to_sky(test_header, 32, 32)
        crop, crophdr = mcf.crop_image_sky(test_data, test_header, [ra, dec], 10.)
        # 10 asec radius at 1 asec/pix -> ~20 pixels wide
        assert abs(crop.shape[0] - 21) <= 1
        assert abs(crop.shape[1] - 21) <= 1


class TestUtilityHelpers:
    def test_squeeze_image_2d_passthrough(self, test_data, test_header):
        out, hdr = mcf.squeeze_image(test_data, test_header, verbose=False)
        assert out.shape == test_data.shape
        assert hdr['NAXIS'] == 2

    def test_squeeze_image_cleans_ghost_higher_axis_cards_from_2d_header(self, test_data, test_header):
        hdr = test_header.copy()
        hdr['NAXIS3'] = 1
        hdr['CTYPE3'] = 'STOKES'
        hdr['CRVAL3'] = 1.0
        hdr['WCSAXES'] = 3
        out, hdr2 = mcf.squeeze_image(test_data, hdr, verbose=False)
        assert out.shape == test_data.shape
        assert hdr2['NAXIS'] == 2
        assert hdr2['WCSAXES'] == 2
        assert 'NAXIS3' not in hdr2
        assert 'CTYPE3' not in hdr2

    def test_squeeze_image_degenerate_4d(self, test_header):
        data = np.ones((1, 1, 8, 9))
        out, hdr = mcf.squeeze_image(data, test_header, verbose=False)
        assert out.shape == (8, 9)
        assert hdr['NAXIS'] == 2

    def test_squeeze_image_raises_on_real_cube(self, test_header):
        data = np.ones((4, 8, 9))
        with pytest.raises(ValueError):
            mcf.squeeze_image(data, test_header, verbose=False)

    def test_header_coord_grids_shape(self, test_header):
        lon, lat = mcf.header_coord_grids(test_header)
        assert lon.shape == (test_header['NAXIS2'], test_header['NAXIS1'])
        assert lat.shape == (test_header['NAXIS2'], test_header['NAXIS1'])


class TestFitsRoundtrip:
    def test_save_rgb_fits_roundtrip(self, tmp_path, test_data, test_header):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        colr = mcf.colorize_image(grey, '#00FF00', colorintype='hex')
        combined = mcf.combine_multicolor([colr], gamma=2.2)
        path = str(tmp_path / 'rgb.fits')
        mcf.save_rgb_fits(path, combined, test_header)
        loaded = pyfits.getdata(path)
        # Written as [3, ny, nx] cube
        assert loaded.shape == (3, test_data.shape[0], test_data.shape[1])
        restored = np.moveaxis(loaded, 0, -1)
        assert np.allclose(restored, combined, atol=1e-6)

    def test_save_rgb_fits_adds_provenance_comment_on_bare_header(self, tmp_path, test_data, test_header):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        colr = mcf.colorize_image(grey, '#00FF00', colorintype='hex')
        combined = mcf.combine_multicolor([colr], gamma=2.2)
        path = str(tmp_path / 'rgb.fits')
        mcf.save_rgb_fits(path, combined, test_header)
        blob = ' '.join(str(c) for c in pyfits.getheader(path).cards).lower()
        assert 'created with multicolorfits' in blob

    def test_plot_combined_rgb(self, tmp_path, test_data, test_header):
        import matplotlib
        matplotlib.use('Agg')
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        colr = mcf.colorize_image(grey, '#00FF00', colorintype='hex')
        combined = mcf.combine_multicolor([colr], gamma=2.2)
        path = str(tmp_path / 'plot.png')
        mcf.plot_combined_rgb(combined, test_header, 'test', path)
        import os
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0


class TestLazyReproject:
    def test_helpful_error_when_missing(self, test_data, test_header, monkeypatch):
        """If reproject isn't installed, a clear ImportError should be raised (only on use)."""
        import sys
        # Setting a module's sys.modules entry to None makes its import raise ImportError
        monkeypatch.setitem(sys.modules, 'reproject', None)
        with pytest.raises(ImportError, match='reproject'):
            mcf.reproject_image(test_data, test_header, make_test_header(32, 32), method='interp')

    def test_reproject_if_available(self, test_data, test_header):
        pytest.importorskip('reproject')
        hdrto = make_test_header(32, 32)
        out = mcf.reproject_image(test_data, test_header, hdrto, method='interp')
        assert out.shape == (32, 32)
