import numpy as np
import astropy.io.fits as pyfits
import pytest


def make_test_header(nx=64, ny=64):
    """Simple valid celestial WCS header."""
    hdr = pyfits.Header()
    hdr['NAXIS'] = 2
    hdr['NAXIS1'] = nx
    hdr['NAXIS2'] = ny
    hdr['CTYPE1'] = 'RA---TAN'
    hdr['CTYPE2'] = 'DEC--TAN'
    hdr['CRPIX1'] = nx / 2.
    hdr['CRPIX2'] = ny / 2.
    hdr['CRVAL1'] = 150.0
    hdr['CRVAL2'] = -30.0
    hdr['CDELT1'] = -0.0002777778  # 1 arcsec/pix
    hdr['CDELT2'] = 0.0002777778
    hdr['RADESYS'] = 'FK5'
    hdr['EQUINOX'] = 2000.0
    return hdr


def make_test_data(nx=64, ny=64, with_nans=False, seed=42):
    """Synthetic image: 2D gaussian blob plus noise."""
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:ny, 0:nx]
    blob = np.exp(-(((x - nx / 2.) / (nx / 8.))**2 + ((y - ny / 2.) / (ny / 8.))**2))
    data = blob + rng.normal(0, 0.05, size=(ny, nx))
    if with_nans:
        data[0:3, 0:3] = np.nan
    return data


@pytest.fixture
def test_header():
    return make_test_header()


@pytest.fixture
def test_data():
    return make_test_data()


@pytest.fixture
def test_data_nans():
    return make_test_data(with_nans=True)


@pytest.fixture
def test_fits_file(tmp_path, test_data, test_header):
    """A synthetic FITS file on disk."""
    path = tmp_path / 'synthetic.fits'
    pyfits.writeto(str(path), make_test_data(), make_test_header())
    return str(path)
