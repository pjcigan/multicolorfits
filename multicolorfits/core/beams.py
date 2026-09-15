"""
Elliptical-Gaussian beam matching.

Ported from Phil Cigan's ``tasks.convolve2Dgaus`` (and the
``*_matchhdr`` wrappers). The kernel is the one implied by J.P. Wild,
Australian Journal of Physics, vol. 23, pp. 113–115 (1970): if the desired
beam F0 is the convolution of a kernel F1 with the current beam F2,

    F0 = F1 ⊛ F2

then with D = maj² − min²,

    D1² = D0² + D2² − 2 D0 D2 cos(2 (PA0 − PA2))
    maj1² = ½ [(maj0² + min0²) − (maj2² + min2²) + D1]
    min1² = ½ [(maj0² + min0²) − (maj2² + min2²) − D1]
    tan(2 PA1) = [D0 sin(2 PA0) − D2 sin(2 PA2)] / [D0 cos(2 PA0) − D2 cos(2 PA2)]

Convolution can only degrade resolution (target beam must be at least as
large as the input). This does not reproject — match beams on a common
grid, then reproject, or the reverse, depending on which pixel scale you
want the kernel measured in. For Jy/beam images set ``per_beam=True`` so
the output is scaled by the ratio of beam areas (same pixel scale both
sides, because the array is not resampled here).
"""

import numpy as np

from .wcs_tools import get_cdelts

__all__ = [
    'convolution_beam',
    'match_beam',
    'match_beam_to_header',
    'convolve2Dgaus',
    'convolve2Dgaus_matchhdr',
]


def convolution_beam(bmaj_from_asec, bmin_from_asec, bpa_from_deg,
                     bmaj_to_asec, bmin_to_asec, bpa_to_deg):
    """
    Kernel FWHM (arcsec, arcsec, deg) that takes the input beam to the target.

    Raises ValueError if the target is smaller than the input (cannot
    deconvolve by convolution).
    """
    maj0 = float(bmaj_to_asec) / 3600.0
    min0 = float(bmin_to_asec) / 3600.0
    maj2 = float(bmaj_from_asec) / 3600.0
    min2 = float(bmin_from_asec) / 3600.0
    pa0 = np.deg2rad(float(bpa_to_deg))
    pa2 = np.deg2rad(float(bpa_from_deg))
    if abs(maj0 - maj2) < 1e-15 and abs(min0 - min2) < 1e-15 and abs(pa0 - pa2) < 1e-12:
        return 0.0, 0.0, float(bpa_to_deg)

    D0 = maj0 ** 2 - min0 ** 2
    D2 = maj2 ** 2 - min2 ** 2
    D1sq = D0 ** 2 + D2 ** 2 - 2.0 * D0 * D2 * np.cos(2.0 * (pa0 - pa2))
    D1 = np.sqrt(max(0.0, D1sq))
    maj1_sq = 0.5 * ((maj0 ** 2 + min0 ** 2) - (maj2 ** 2 + min2 ** 2) + D1)
    min1_sq = 0.5 * ((maj0 ** 2 + min0 ** 2) - (maj2 ** 2 + min2 ** 2) - D1)
    if maj1_sq < -1e-18 or min1_sq < -1e-18:
        raise ValueError(
            'target beam is smaller than the input beam along at least one axis; '
            'convolution can only degrade resolution')
    maj1 = np.sqrt(max(0.0, maj1_sq))
    min1 = np.sqrt(max(0.0, min1_sq))
    if abs(D0) < 1e-30 and abs(D2) < 1e-30:
        pa1 = 0.0
    else:
        # Published tan(2 PA1) formula, evaluated with arctan2 so the quadrant
        # is preserved (the handwritten tasks.py draft used sin(PA2) in one
        # branch; this follows the comment / Wild formula).
        pa1 = 0.5 * np.arctan2(
            D0 * np.sin(2.0 * pa0) - D2 * np.sin(2.0 * pa2),
            D0 * np.cos(2.0 * pa0) - D2 * np.cos(2.0 * pa2))
    return float(maj1 * 3600.0), float(min1 * 3600.0), float(np.rad2deg(pa1))


def _beam_area_ratio(bmaj_to, bmin_to, bmaj_from, bmin_from):
    return (float(bmaj_to) * float(bmin_to)) / (float(bmaj_from) * float(bmin_from))


def match_beam(data, hdr, bmaj_to_asec, bmin_to_asec, bpa_to_deg,
               bmaj_from_asec=None, bmin_from_asec=None, bpa_from_deg=None,
               per_beam=False, pad_pix=20, update_header=True):
    """
    Convolve a 2D map so its beam matches the requested target.

    Beam major/minor are FWHM in arcsec; PA is degrees east of north.
    Missing ``*_from`` values are read from ``BMAJ`` / ``BMIN`` / ``BPA``
    (degrees in the header).

    per_beam : bool
        Scale by (target beam area) / (input beam area). Use this for
        Jy/beam maps. Both areas use the header FWHM values — this function
        does not reproject, so the pixel scale is unchanged.
    """
    if bmaj_from_asec is None:
        bmaj_from_asec = float(hdr['BMAJ']) * 3600.0
        bmin_from_asec = float(hdr['BMIN']) * 3600.0
        bpa_from_deg = float(hdr['BPA']) if 'BPA' in hdr else 0.0
    maj1, min1, pa1 = convolution_beam(
        bmaj_from_asec, bmin_from_asec, bpa_from_deg,
        bmaj_to_asec, bmin_to_asec, bpa_to_deg)
    ratio = 1.0
    if per_beam:
        ratio = _beam_area_ratio(bmaj_to_asec, bmin_to_asec, bmaj_from_asec, bmin_from_asec)
    arr = np.asarray(data, dtype=float)
    if arr.ndim != 2:
        raise ValueError('match_beam expects a 2D image (got shape %s)' % (arr.shape,))
    if maj1 < 1e-6 and min1 < 1e-6:
        out = arr * ratio
    else:
        out = _convolve_gaussian(arr, hdr, maj1, min1, pa1, pad_pix=pad_pix) * ratio
    if not update_header:
        return out
    hdrout = hdr.copy()
    hdrout['BMAJ'] = float(bmaj_to_asec) / 3600.0
    hdrout['BMIN'] = float(bmin_to_asec) / 3600.0
    hdrout['BPA'] = float(bpa_to_deg)
    return out, hdrout


def match_beam_to_header(data, hdr, hdr_target, per_beam=False, pad_pix=20,
                         update_header=True):
    """Match ``data`` to the BMAJ/BMIN/BPA of ``hdr_target``."""
    return match_beam(
        data, hdr,
        float(hdr_target['BMAJ']) * 3600.0,
        float(hdr_target['BMIN']) * 3600.0,
        float(hdr_target['BPA']) if 'BPA' in hdr_target else 0.0,
        per_beam=per_beam, pad_pix=pad_pix, update_header=update_header)


def _convolve_gaussian(data, hdr, maj_asec, min_asec, pa_deg, pad_pix=20):
    from astropy.convolution import convolve_fft
    from astropy.convolution import Model2DKernel
    from astropy.modeling.models import Gaussian2D

    pixsizex, pixsizey = np.abs(get_cdelts(hdr)[:2])
    maj = float(maj_asec) / 3600.0
    min_ = float(min_asec) / 3600.0
    pa = np.deg2rad(float(pa_deg))
    sigmax = maj / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    sigmay = min_ / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    kernel_size = int(np.ceil(0.5 * np.sqrt(np.pi / np.log(2.0)) * maj / pixsizex)) + 2 * int(pad_pix)
    if kernel_size % 2 == 0:
        kernel_size += 1
    kernel_size = max(kernel_size, 3)
    model = Gaussian2D(1.0, 0.0, 0.0, sigmax / pixsizex, sigmay / pixsizey, pa)
    kernel = Model2DKernel(model, x_size=kernel_size, y_size=kernel_size)
    filled = np.nan_to_num(data)
    return convolve_fft(filled, kernel, normalize_kernel=True, allow_huge=True)


def convolve2Dgaus(datain, headerin, beammaj_from_asec, beammin_from_asec, beampa_from_deg,
                   beammaj_to_asec, beammin_to_asec, beampa_to_deg,
                   pad_pix=20, perbeamunitratio=1.0, use_astropy=True,
                   just_plot_ellipses=False):
    """
    tasks.py-compatible entry point. Always uses astropy.convolution
    (the old pixel-loop kernel is not reproduced). ``just_plot_ellipses``
    is accepted and ignored.
    """
    del use_astropy, just_plot_ellipses
    out, _hdr = match_beam(
        datain, headerin, beammaj_to_asec, beammin_to_asec, beampa_to_deg,
        bmaj_from_asec=beammaj_from_asec, bmin_from_asec=beammin_from_asec,
        bpa_from_deg=beampa_from_deg, per_beam=False, pad_pix=pad_pix,
        update_header=True)
    return out * float(perbeamunitratio)


def convolve2Dgaus_matchhdr(datain, headerin, header_to_match, pad_pix=20,
                            perbeamunits=False, use_astropy=True,
                            just_plot_ellipses=False):
    """tasks.py-compatible wrapper around match_beam_to_header."""
    del use_astropy, just_plot_ellipses
    out, _hdr = match_beam_to_header(
        datain, headerin, header_to_match, per_beam=perbeamunits,
        pad_pix=pad_pix, update_header=True)
    return out
