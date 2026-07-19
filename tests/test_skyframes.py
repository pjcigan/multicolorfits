"""Tests for experimental sky-frame and stack-alignment helpers."""

import numpy as np
import pytest
import astropy.io.fits as pyfits

import multicolorfits as mcf
from conftest import make_test_header, make_test_data


class TestHeaderFrameName:
    def test_equatorial_default(self, test_header):
        assert mcf.header_frame_name(test_header) in ('fk5', 'icrs', 'FK5', 'ICRS')

    def test_galactic_ctype(self):
        hdr = make_test_header()
        hdr['CTYPE1'] = 'GLON-TAN'
        hdr['CTYPE2'] = 'GLAT-TAN'
        assert mcf.header_frame_name(hdr) == 'galactic'


class TestConvertHeaderFrame:
    def test_galactic_header_shape(self, test_header):
        hdr_gal = mcf.convert_header_frame(test_header, frame='galactic')
        assert hdr_gal['NAXIS'] == 2
        assert str(hdr_gal['CTYPE1']).startswith('GLON')
        assert str(hdr_gal['CTYPE2']).startswith('GLAT')
        assert int(hdr_gal['NAXIS1']) >= int(test_header['NAXIS1'])

    def test_invalid_frame_raises(self, test_header):
        with pytest.raises(ValueError):
            mcf.convert_header_frame(test_header, frame='not_a_frame')


class TestReprojectToFrame:
    def test_reproject_if_available(self, test_data, test_header):
        pytest.importorskip('reproject')
        data_gal, hdr_gal = mcf.reproject_to_galactic(test_data, test_header)
        assert data_gal.shape == (hdr_gal['NAXIS2'], hdr_gal['NAXIS1'])
        assert str(hdr_gal['CTYPE1']).startswith('GLON')

    def test_alias_matches_explicit(self, test_data, test_header):
        pytest.importorskip('reproject')
        a, ha = mcf.reproject_to_galactic(test_data, test_header)
        b, hb = mcf.reproject_to_frame(test_data, test_header, frame='galactic')
        assert a.shape == b.shape
        assert ha['CTYPE1'] == hb['CTYPE1']


class TestStackTools:
    def test_downsample_unchanged_when_small(self):
        small = np.zeros((32, 32))
        assert mcf.downsample_for_preview(small, max_size=512) is small

    def test_downsample_reduces_large(self):
        big = np.zeros((1000, 800))
        out = mcf.downsample_for_preview(big, max_size=256)
        assert max(out.shape[:2]) <= 256

    def test_downsample_rgb(self):
        rgb = np.zeros((600, 400, 3))
        out = mcf.downsample_for_preview(rgb, max_size=200)
        assert out.shape[2] == 3
        assert max(out.shape[:2]) <= 200

    def test_reproject_stack_to_reference(self, test_data, test_header):
        pytest.importorskip('reproject')
        # Slightly shifted header as second "layer"
        hdr2 = test_header.copy()
        hdr2['CRVAL1'] = float(hdr2['CRVAL1']) + 0.01
        images = [(test_data, test_header), (test_data * 0.5, hdr2)]
        data_list, hdr_out = mcf.reproject_stack_to_reference(images, reference=0)
        assert len(data_list) == 2
        assert data_list[0].shape == (hdr_out['NAXIS2'], hdr_out['NAXIS1'])
        assert data_list[1].shape == data_list[0].shape

    def test_align_stack(self, test_data, test_header):
        pytest.importorskip('reproject')
        hdr2 = test_header.copy()
        hdr2['CRVAL2'] = float(hdr2['CRVAL2']) + 0.005
        images = [(test_data, test_header), (test_data * 0.8, hdr2)]
        aligned = mcf.align_stack(images, reference=0)
        assert len(aligned) == 2
        assert aligned[0][0].shape == aligned[1][0].shape
        assert aligned[0][1] is aligned[1][1]  # same header object


class TestOptimalHeader:
    def test_optimal_if_available(self, test_data, test_header):
        pytest.importorskip('reproject')
        hdr2 = test_header.copy()
        hdr2['CRVAL1'] = float(hdr2['CRVAL1']) + 0.02
        images = [(test_data, test_header), (test_data, hdr2)]
        hdr = mcf.optimal_common_header(images)
        assert hdr['NAXIS'] == 2
        assert int(hdr['NAXIS1']) > 0
