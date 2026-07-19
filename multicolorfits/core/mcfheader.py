"""
Stateful FITS-header wrapper for multicolorfits.

Historically (v2.x) header treatment was purely functional: every consumer
re-parsed a raw ``astropy.io.fits.Header`` and rebuilt a ``WCS`` on demand
(``force_header_floats``, ``force_header_2d``, ``get_cdelts``, ``header_frame_name``
...).  ``McfHeader`` keeps all of those functions working unchanged, but wraps a
header in a small object that:

* normalizes once on construction (float-coerces WCS cards, drops >2D axes),
* lazily builds and **caches** the celestial ``WCS``, frame name, and pixel
  scale, and
* still behaves like a ``Header`` (indexing / ``in`` / ``get`` / iteration) and
  exposes the underlying object via :attr:`astropy` for any code (or the classic
  functional API, ``save_rgb_fits``, ``reproject_image`` ...) that wants the raw
  ``astropy.io.fits.Header``.

This makes it cheap for each GUI panel to retain its own header-derived state.
"""

import astropy.io.fits as pyfits
import astropy.wcs as pywcs

from .wcs_tools import (force_header_2d, force_header_floats, normalize_wcs_meta,
                        _has_higher_axis_cards)

__all__ = ['McfHeader', 'as_mcfheader', 'save_header', 'load_header']

_EMPTY_HEADER_STRING = 'COMMENT  No header'


class McfHeader:
    """
    A stateful wrapper around an ``astropy.io.fits.Header``.

    Parameters
    ----------
    header : astropy.io.fits.Header, McfHeader, str, dict, or None
        Source header.  ``None`` yields a minimal placeholder header.
    normalize : bool
        Float-coerce WCS cards (``force_header_floats``) and integer-coerce
        ``NAXIS`` / ``WCSAXES`` (``normalize_wcs_meta``) on construction.
    force_2d : bool
        Strip >2D structure cards (``force_header_2d``) when ``NAXIS > 2``.
    """

    def __init__(self, header=None, normalize=True, force_2d=True):
        if header is None:
            header = pyfits.Header.fromstring(_EMPTY_HEADER_STRING)
        elif isinstance(header, McfHeader):
            header = header.astropy
        elif isinstance(header, str):
            header = pyfits.Header.fromstring(header, sep='\n')
        elif not isinstance(header, pyfits.Header):
            header = pyfits.Header(header)
        header = header.copy()
        if normalize:
            force_header_floats(header)
        if force_2d and (int(header.get('NAXIS', 2)) > 2 or _has_higher_axis_cards(header)):
            header = force_header_2d(header)
        elif normalize:
            header = normalize_wcs_meta(header)
        self._hdr = header
        self._invalidate_cache()

    def _invalidate_cache(self):
        self._wcs = None
        self._frame = None
        self._pixscale = None
        self._pixscale_done = False

    # ------------------------------------------------------------------ #
    # interop with the raw astropy header / WCS
    # ------------------------------------------------------------------ #
    @property
    def astropy(self):
        """The underlying ``astropy.io.fits.Header`` (for the functional API)."""
        return self._hdr

    @property
    def wcs(self):
        """Cached celestial ``astropy.wcs.WCS`` built from this header."""
        if self._wcs is None:
            self._wcs = pywcs.WCS(self._hdr).celestial
        return self._wcs

    @property
    def frame(self):
        """Cached celestial frame name (e.g. ``'fk5'``, ``'galactic'``)."""
        if self._frame is None:
            from .skyframes import header_frame_name
            self._frame = header_frame_name(self._hdr)
        return self._frame

    @property
    def is_celestial(self):
        """True if the header defines a usable celestial WCS."""
        try:
            return bool(self.wcs.has_celestial)
        except Exception:
            return False

    @property
    def pixscale_asec(self):
        """Pixel scale in arcsec/pixel, or ``None`` if undefined/anisotropic."""
        if not self._pixscale_done:
            self._pixscale_done = True
            from .wcs_tools import arcsec_per_pixel
            try:
                self._pixscale = abs(float(arcsec_per_pixel(self._hdr)))
            except Exception:
                self._pixscale = None
        return self._pixscale

    @property
    def shape(self):
        """Image shape ``(ny, nx)`` from NAXIS2/NAXIS1, or ``None``."""
        try:
            return (int(self._hdr['NAXIS2']), int(self._hdr['NAXIS1']))
        except Exception:
            return None

    def matches_grid(self, other, rtol=1e-6):
        """
        True if ``other`` shares this header's pixel grid (shape + WCS), so the
        two images can be combined without reprojection.
        """
        other = as_mcfheader(other)
        if self.shape != other.shape:
            return False
        try:
            import numpy as np
            a, b = self.wcs.wcs, other.wcs.wcs
            return (np.allclose(a.crpix, b.crpix, rtol=rtol, atol=rtol) and
                    np.allclose(a.crval, b.crval, rtol=rtol, atol=rtol) and
                    np.allclose(self.wcs.pixel_scale_matrix,
                                other.wcs.pixel_scale_matrix, rtol=rtol, atol=1e-12))
        except Exception:
            return self.shape == other.shape

    # ------------------------------------------------------------------ #
    # frame conversion (delegates to skyframes)
    # ------------------------------------------------------------------ #
    def to_frame(self, frame='galactic', **kwargs):
        """
        Return a new :class:`McfHeader` describing the same sky area in a
        different celestial frame (wraps ``skyframes.convert_header_frame``).
        Resample the data separately with ``reproject_image`` / ``reproject_to_frame``.
        """
        from .skyframes import convert_header_frame
        return McfHeader(convert_header_frame(self._hdr, frame=frame, **kwargs))

    # ------------------------------------------------------------------ #
    # serialization
    # ------------------------------------------------------------------ #
    def minimal_wcs_header(self):
        """
        A stripped-down ``astropy.io.fits.Header`` containing only the cards
        needed to define the pixel grid -- the celestial WCS plus ``NAXIS``,
        ``NAXIS1`` and ``NAXIS2``.  HISTORY/COMMENT and other non-WCS cards are
        dropped, which keeps exported/portable common headers compact while
        still being sufficient for ``reproject_image`` and WCS axes.
        """
        hdr = self.wcs.to_header(relax=True)
        shape = self.shape
        hdr.set('NAXIS', 2, before=0)
        if shape is not None:
            hdr.set('NAXIS1', int(shape[1]), after='NAXIS')
            hdr.set('NAXIS2', int(shape[0]), after='NAXIS1')
        return hdr

    def tostring(self, sep='\n'):
        return self._hdr.tostring(sep=sep)

    @classmethod
    def fromstring(cls, text, sep='\n', **kwargs):
        return cls(pyfits.Header.fromstring(text, sep=sep), **kwargs)

    def totextfile(self, path, overwrite=True):
        """Write this header to a portable plain-text ``.hdr`` file."""
        self._hdr.totextfile(path, overwrite=overwrite)

    @classmethod
    def fromtextfile(cls, path, **kwargs):
        """Load an :class:`McfHeader` from a plain-text header file."""
        return cls(pyfits.Header.fromtextfile(path), **kwargs)

    def copy(self):
        return McfHeader(self._hdr)

    # ------------------------------------------------------------------ #
    # mapping delegation -- behave like a Header
    # ------------------------------------------------------------------ #
    def __getitem__(self, key):
        return self._hdr[key]

    def __setitem__(self, key, value):
        self._hdr[key] = value
        self._invalidate_cache()

    def __delitem__(self, key):
        del self._hdr[key]
        self._invalidate_cache()

    def __contains__(self, key):
        return key in self._hdr

    def __iter__(self):
        return iter(self._hdr)

    def __len__(self):
        return len(self._hdr)

    def get(self, key, default=None):
        return self._hdr.get(key, default)

    def keys(self):
        return self._hdr.keys()

    def set(self, *args, **kwargs):
        self._hdr.set(*args, **kwargs)
        self._invalidate_cache()

    def __repr__(self):
        return 'McfHeader(frame=%r, shape=%r, pixscale_asec=%r)' % (
            self.frame, self.shape, self.pixscale_asec)


def as_mcfheader(header, **kwargs):
    """Return ``header`` as an :class:`McfHeader` (pass-through if it already is)."""
    if isinstance(header, McfHeader):
        return header
    return McfHeader(header, **kwargs)


def save_header(header, path, overwrite=True, minimal=False):
    """
    Write a header to a portable plain-text ``.hdr`` file.

    Parameters
    ----------
    header : McfHeader or astropy.io.fits.Header
    path : str
    overwrite : bool
    minimal : bool
        If True and ``header`` is an :class:`McfHeader`, save only the WCS +
        NAXIS cards (see :meth:`McfHeader.minimal_wcs_header`).
    """
    if isinstance(header, McfHeader):
        hdr = header.minimal_wcs_header() if minimal else header.astropy
    else:
        hdr = header
    hdr.totextfile(path, overwrite=overwrite)


def load_header(path):
    """Load an :class:`McfHeader` from a plain-text header file written by
    :func:`save_header` / :meth:`McfHeader.totextfile`."""
    return McfHeader(pyfits.Header.fromtextfile(path))
