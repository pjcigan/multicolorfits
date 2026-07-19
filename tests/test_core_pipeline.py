"""Image pipeline: scaling, greyRGBize, colorize, combine, NaN handling."""

import numpy as np
import pytest

import multicolorfits as mcf


class TestScaling:
    def test_rescale_linear_range(self, test_data):
        scaled = mcf.rescale_image(test_data, 'linear')
        assert np.nanmin(scaled) >= 0.0
        assert np.nanmax(scaled) <= 1.0

    @pytest.mark.parametrize('stretch', ['linear', 'sqrt', 'squared', 'log', 'power', 'sinh', 'asinh'])
    def test_all_stretches_run(self, test_data, stretch):
        scaled = mcf.rescale_image(test_data, stretch)
        assert scaled.shape == test_data.shape
        assert np.nanmax(scaled) <= 1.0 + 1e-6

    def test_rescale_matches_astropy_directly(self, test_data):
        from astropy.visualization import SqrtStretch, ManualInterval
        vmin, vmax = 0.0, 1.0
        expected = (SqrtStretch() + ManualInterval(vmin=vmin, vmax=vmax))(test_data)
        result = mcf.rescale_image(test_data, 'sqrt', vmin=vmin, vmax=vmax)
        assert np.allclose(result, expected, equal_nan=True)

    def test_adjust_gamma_nan_safe(self):
        arr = np.array([0.25, np.nan, 1.0])
        out = mcf.adjust_gamma(arr, 2.0)
        assert np.isnan(out[1])
        assert np.isclose(out[0], 0.25**2.0)

    def test_zscale_limits(self, test_data):
        vmin, vmax = mcf.zscale_limits(test_data)
        assert vmin < vmax

    def test_nan_percentile_of_score(self, test_data_nans):
        # Should not raise or return nan despite NaNs in data
        score = mcf.nan_percentile_of_score(test_data_nans.ravel(), np.nanmedian(test_data_nans))
        assert 0 <= score <= 100
        assert not np.isnan(score)


class TestGreyRGBize:
    def test_shape_and_range(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        assert grey.shape == test_data.shape + (3,)
        assert np.nanmin(grey) >= 0.0
        assert np.nanmax(grey) <= 1.0

    def test_channels_equal(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='asinh')
        assert np.allclose(grey[..., 0], grey[..., 1], equal_nan=True)
        assert np.allclose(grey[..., 0], grey[..., 2], equal_nan=True)

    def test_percentile_scaletype(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear', scaletype='perc', min_max=[1., 99.])
        assert grey.shape == test_data.shape + (3,)

    def test_nan_input(self, test_data_nans):
        grey = mcf.to_grey_rgb(test_data_nans, rescalefn='sqrt')
        assert grey.shape == test_data_nans.shape + (3,)
        # Non-NaN region should be finite
        assert np.isfinite(grey[10:, 10:]).all()


class TestColorize:
    def test_white_preserves_grey(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        colr = mcf.colorize_image(grey, '#FFFFFF', colorintype='hex')
        assert np.allclose(np.nan_to_num(colr), np.nan_to_num(grey), atol=1e-3)

    def test_pure_red(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        colr = mcf.colorize_image(grey, '#FF0000', colorintype='hex')
        # Green/blue channels should be zero for pure red tint
        assert np.nanmax(colr[..., 1]) < 1e-6
        assert np.nanmax(colr[..., 2]) < 1e-6

    def test_direct_rgb_matches_hsv_path(self, test_data):
        """The fast direct-RGB path (gammacorr_color=1) should closely match the HSV path."""
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        fast = mcf.colorize_image_direct_rgb(grey, '#C11B17')
        # Force the HSV path by using colorintype='hsv' with the same color
        hsv_color = mcf.hex_to_hsv('#C11B17')
        slow = mcf.colorize_image(grey, hsv_color, colorintype='hsv')
        assert np.allclose(np.nan_to_num(fast), np.nan_to_num(slow), atol=5e-3)

    def test_rgb_colorintype(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        colr = mcf.colorize_image(grey, (193, 27, 23), colorintype='rgb')
        assert colr.shape == grey.shape

    def test_invalid_colorintype_raises(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        with pytest.raises(Exception):
            mcf.colorize_image(grey, '#FFFFFF', colorintype='bogus')

    def test_gammacorr_color_path(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        colr = mcf.colorize_image(grey, '#C11B17', colorintype='hex', gammacorr_color=2.2)
        assert colr.shape == grey.shape
        assert np.nanmax(colr) <= 1.0 + 1e-6


class TestCombine:
    def test_single_image(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        colr = mcf.colorize_image(grey, '#FF0000', colorintype='hex')
        combined = mcf.combine_multicolor([colr], gamma=2.2)
        assert combined.shape == colr.shape
        assert combined.min() >= 0.0
        assert combined.max() <= 1.0

    def test_two_images(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        c1 = mcf.colorize_image(grey, '#FF0000', colorintype='hex')
        c2 = mcf.colorize_image(grey, '#0000FF', colorintype='hex')
        combined = mcf.combine_multicolor([c1, c2], gamma=2.2)
        assert combined.shape == c1.shape
        # Red + blue should give purple-ish: both R and B present
        assert combined[..., 0].max() > 0.5
        assert combined[..., 2].max() > 0.5

    def test_inverse_background_white(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear', gamma=1. / 2.2)
        colr = mcf.colorize_image(grey, mcf.hex_complement('#FF0000'), colorintype='hex')
        combined = mcf.combine_multicolor([colr], gamma=2.2, inverse=True)
        assert combined.min() >= 0.0
        assert combined.max() <= 1.0


class TestSmooth:
    def test_2d(self, test_data):
        sm = mcf.smooth_image(test_data, sigma=2)
        assert sm.shape == test_data.shape
        assert sm.std() < test_data.std()

    def test_3d_multichannel(self, test_data):
        grey = mcf.to_grey_rgb(test_data, rescalefn='linear')
        sm = mcf.smooth_image(np.nan_to_num(grey), sigma=2)
        assert sm.shape == grey.shape
