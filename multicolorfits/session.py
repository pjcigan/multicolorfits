"""
Toolkit-agnostic session model and render controller.

Encodes the same processing pipeline historically used by the TraitsUI GUI,
with no GUI toolkit dependencies, so the browser GUI, Qt GUI, and scripts
share one implementation::

    scaled   = stretch(data, vmin, vmax)               # rescale_image
    greyRGB  = intensity_to_grey_rgb(adjust_gamma(scaled, gamma))
    colorRGB = colorize_image(greyRGB, hexcolor, gammacorr_color=gamma)
    display  = colorRGB ** (1/gamma)                   # per-panel preview
    combined = combine_multicolor([colorRGB, ...], gamma)

Typical scripting flow::

    session = McfSession()
    session.load_files(['a.fits', 'b.fits'], colors=['#C11B17', '#4CC417'])
    session.align_panels(reference=0)   # if grids differ
    rgb = session.render_combined()
    session.save_rgb_fits('out.fits')
"""

import os
from dataclasses import dataclass
from typing import Optional, List

import numpy as np
import astropy.io.fits as pyfits

from .core import (
    stretch_functions,
    adjust_gamma,
    nan_percentile_of_score,
    rescale_image,
    zscale_limits,
    McfHeader,
    colorize_image,
    hex_complement,
    smooth_image,
    save_rgb_fits,
    combine_colorized_layers,
    combo_swatch,
    squeeze_image,
)
from .core.colorize import _intensity_to_grey_rgb
from .core.io_output import annotate_provenance_header


def _json_float(value, default=None):
    """Return a finite float for JSON, or *default* if missing / NaN / Inf.

    Starlette / FastAPI serialize with ``allow_nan=False``, so a single NaN in
    a panel or compose field will abort the response and stall the web GUI.
    """
    if value is None:
        return default
    try:
        v = float(value)
    except (TypeError, ValueError):
        return default
    if not np.isfinite(v):
        return default
    return v


from .core.stack_tools import downsample_for_preview
from .core.colormix import COLORSPACES
from .palettes import (
    PALETTES, list_palettes, list_palette_menu, get_palette, suggest_colors,
    colors_for_hue_pattern, rotate_colors_from_base, is_palette_name,
    check_palette_colorblind,
)

__all__ = ['PanelState', 'ComposeState', 'McfSession', 'reproject_available',
           'STRETCHES', 'COMBINE_MODES', 'ALIGN_FRAMES',
           'DEFAULT_N_PANELS', 'MIN_PANELS', 'MAX_PANELS']

STRETCHES = list(stretch_functions.keys())  # linear, sqrt, squared, log, power, sinh, asinh

# User-facing compositing modes (GUI + combine_colorized_layers).
COMBINE_MODES = ('rgb',) + tuple(COLORSPACES) + ('ryb', 'cmyk')

# Celestial frames offered for the "reproject / align" tools, in menu order.
ALIGN_FRAMES = ('icrs', 'galactic', 'fk5', 'fk4', 'ecliptic')

# CVD types checked by the palette colorblind report, in menu order.
CVD_KINDS = ('deuteranopia', 'protanopia', 'tritanopia')

# Panel-count bounds shared by McfSession and both GUIs.
DEFAULT_N_PANELS = 4
MIN_PANELS = 1
MAX_PANELS = 16


def reproject_available():
    """Return True if the optional ``reproject`` package is installed (needed to align layers)."""
    import importlib.util
    return importlib.util.find_spec('reproject') is not None


@dataclass
class PanelState:
    """
    State for one image layer (one GUI panel).

    Attributes
    ----------
    filepath : str
        Path of the loaded FITS file, if any.
    data, header : array / McfHeader
        Image data and WCS/header after optional squeeze / reproject.
    stretch, vmin, vmax : str / float
        Display stretch and absolute limits for :func:`rescale_image`.
    percent_min, percent_max : float
        Percentile limits (updated when the user adjusts percentile controls).
    color, label : str
        Layer color (hex) and optional legend label.
    smooth, smooth_sigma : bool / float
        Optional Gaussian smoothing before colorize.
    in_use : bool
        Whether this panel contributes to the combined image.
    reprojected : bool
        True after the layer has been aligned off its on-disk pixel grid.
    """
    filepath: str = ''
    data: Optional[np.ndarray] = None
    header: Optional[McfHeader] = None
    stretch: str = 'linear'
    vmin: float = 0.0
    vmax: float = 1.0
    percent_min: float = 0.0
    percent_max: float = 100.0
    data_min: float = 0.0    # full data range, set at load time
    data_max: float = 1.0
    color: str = '#FFFFFF'
    label: str = ''          # channel name for the combined-figure legend
    smooth: bool = False
    smooth_sigma: float = 3.0
    in_use: bool = False
    reprojected: bool = False  # True once aligned/reprojected off its on-disk grid

    # ---------- data loading ----------

    def load_fits(self, filepath):
        """
        Load a FITS image into this panel.

        Extra axes on >2D cubes are squeezed out via :func:`squeeze_image`.
        Sets ``vmin`` / ``vmax`` from the data range and marks the panel in use.
        """
        data, hdr = pyfits.getdata(filepath, header=True)
        self.set_data(np.asarray(data, dtype=float), hdr, filepath=filepath)

    def set_data(self, data, header=None, filepath=''):
        """
        Assign image data (and optional header) from arrays.

        Parameters
        ----------
        data : array
            2D intensity map (or squeezeable cube).
        header : fits header or McfHeader or None
            Accompanying WCS/header.
        filepath : str
            Optional source path recorded for provenance / export scripts.
        """
        data = np.asarray(data, dtype=float)
        if header is not None:
            data, header = squeeze_image(data, header, verbose=False)
        elif data.ndim > 2:
            # Preserve the legacy fallback for header-less arrays.
            data = data.reshape(data.shape[-2], data.shape[-1]) if all(
                s == 1 for s in data.shape[:-2]) else data[(0,) * (data.ndim - 2)]
        self.data = data
        # McfHeader normalizes (float-coerces WCS cards) and forces 2D on load,
        # and caches the derived WCS / frame / pixel scale for this panel.
        self.header = McfHeader(header)
        if filepath:
            self.filepath = os.path.abspath(os.path.expanduser(str(filepath)))
        else:
            self.filepath = ''
        # Default legend label = filename stem (user can override in the GUI).
        if filepath:
            self.label = os.path.splitext(os.path.basename(filepath))[0]
        self.data_min = float(np.nanmin(data))
        self.data_max = float(np.nanmax(data))
        self.vmin = self.data_min
        self.vmax = self.data_max
        self.percent_min = 0.0
        self.percent_max = 100.0
        self.in_use = True
        self.reprojected = False

    def set_reprojected(self, data, header):
        """
        Replace this panel's data and header with a reprojected version
        (from an align operation), preserving the user's display settings
        (stretch, color, smoothing, and absolute vmin/vmax scaling limits).
        """
        self.data = np.asarray(data, dtype=float)
        self.header = McfHeader(header)
        if np.isfinite(self.data).any():
            self.data_min = float(np.nanmin(self.data))
            self.data_max = float(np.nanmax(self.data))
        self.in_use = True
        self.reprojected = True

    def clear(self):
        """Reset the panel to its unloaded state."""
        self.__dict__.update(PanelState().__dict__)

    # ---------- limits ----------

    def set_limits(self, vmin=None, vmax=None):
        """Set absolute scaling limits, updating the matching percentiles."""
        if vmin is not None:
            self.vmin = float(vmin)
            if self.data is not None:
                self.percent_min = float(np.round(nan_percentile_of_score(self.data.ravel(), self.vmin, kind='strict'), 2))
        if vmax is not None:
            self.vmax = float(vmax)
            if self.data is not None:
                self.percent_max = float(np.round(nan_percentile_of_score(self.data.ravel(), self.vmax, kind='strict'), 2))

    def set_percentiles(self, pmin=None, pmax=None):
        """Set scaling limits from data percentiles."""
        if self.data is None:
            return
        if pmin is not None:
            self.percent_min = float(pmin)
            self.vmin = float(np.nanpercentile(self.data, self.percent_min))
        if pmax is not None:
            self.percent_max = float(pmax)
            self.vmax = float(np.nanpercentile(self.data, self.percent_max))

    def reset_minmax(self):
        """Reset limits to the full data min/max."""
        self.set_limits(self.data_min, self.data_max)

    def apply_zscale(self):
        """Set limits with the IRAF zscale algorithm."""
        if self.data is None:
            return
        vmin, vmax = zscale_limits(self.data)
        self.set_limits(vmin, vmax)

    # ---------- header editing ----------

    def header_string(self):
        """The panel header serialized as newline-separated cards."""
        if self.header is None:
            return 'COMMENT  No image loaded'
        return self.header.tostring(sep='\n')

    def apply_header_string(self, hdrstring):
        """Replace the panel header from newline-separated card text."""
        # Respect the user's edits verbatim (no float coercion / 2D forcing),
        # but still get a cached-WCS McfHeader.
        self.header = McfHeader.fromstring(hdrstring, sep='\n',
                                           normalize=False, force_2d=False)

    # ---------- rendering ----------

    def render_color_rgb(self, gamma=2.2, inverse=False, dtype=None, max_size=None):
        """
        Produce the colorized RGB array for this panel (the quantity that
        combine_multicolor consumes).  Not display-gamma-corrected.

        Parameters
        ----------
        gamma : float
        inverse : bool
        dtype : numpy dtype or None
            Working dtype.  Pass numpy.float32 for a faster/low-memory preview.
        max_size : int or None
            If set, downsample the source data so its longest side is at most
            max_size pixels *before* the compose (fast interactive preview).
            The full-resolution result uses max_size=None.
        """
        if self.data is None:
            raise ValueError('No data loaded in this panel')
        data = self.data
        if max_size is not None:
            data = downsample_for_preview(data, max_size=max_size)
        scaled = rescale_image(data, self.stretch, vmin=self.vmin, vmax=self.vmax)
        if dtype is not None:
            scaled = scaled.astype(dtype)
        # Same greyscale-RGB expansion as to_grey_rgb (broadcast, not stack/gray2rgb).
        grey_rgb = _intensity_to_grey_rgb(adjust_gamma(scaled, gamma), dtype=dtype)
        color = hex_complement(self.color) if inverse else self.color
        color_rgb = colorize_image(grey_rgb, color, colorintype='hex',
                                   gammacorr_color=gamma,
                                   dtype=dtype if dtype is not None else np.float64)
        if self.smooth:
            color_rgb = smooth_image(color_rgb, sigma=self.smooth_sigma)
        return color_rgb

    def render_display(self, gamma=2.2, inverse=False, dtype=None, max_size=None):
        """
        Display-ready [0..1] RGB preview of this single panel
        (gamma-corrected; inverted single images are white-background).

        Pass ``dtype=numpy.float32`` and ``max_size`` for fast GUI previews.
        """
        color_rgb = self.render_color_rgb(gamma=gamma, inverse=inverse,
                                          dtype=dtype, max_size=max_size)
        disp = color_rgb ** (1. / gamma)
        if inverse:
            disp = 1. - disp
        return np.clip(np.nan_to_num(disp), 0, 1)

    # ---------- serialization ----------

    def to_dict(self):
        fp = self.filepath
        if fp and self.in_use:
            exp = os.path.expanduser(fp)
            if os.path.isfile(exp):
                fp = os.path.abspath(exp)
            elif not os.path.isabs(exp):
                fp = os.path.abspath(exp)
        return {
            'filepath': fp,
            'loaded': self.in_use,
            'stretch': self.stretch,
            'vmin': _json_float(self.vmin),
            'vmax': _json_float(self.vmax),
            'percent_min': _json_float(self.percent_min, 0.0),
            'percent_max': _json_float(self.percent_max, 100.0),
            'data_min': _json_float(self.data_min),
            'data_max': _json_float(self.data_max),
            'color': self.color,
            'label': self.label,
            'smooth': self.smooth,
            'smooth_sigma': _json_float(self.smooth_sigma, 0.0),
            'shape': list(self.data.shape) if self.data is not None else None,
        }


@dataclass
class ComposeState:
    """
    Combined-image settings and plot styling for a session.

    Key compositing fields:

    * ``gamma`` — display gamma used when building and combining layers.
    * ``combine_mode`` — ``'rgb'``, ``'lab'``, ``'hsv'``, ``'hsl'``, ``'ryb'``,
      or ``'cmyk'``.
    * ``combine_blend`` — lightness blend for Lab/HSV/HSL (``'screen'``, …).
    * ``combine_background`` — ``'black'``, ``'white'``, ``'transparent'``, or
      a custom color.  A non-black background replaces the legacy
      ``inverse=True`` path (which complemented hues via ``1 - image``).
    * ``bare_plot`` — hide ticks, title, legend, and overlays for image-only
      export.
    * Overlay and legend toggles (`show_legend`, `show_compass`, …).
    """
    gamma: float = 2.2
    inverse: bool = False
    tickcolor: str = '0.9'
    facecolor: str = 'white'   # figure/canvas background; 'none' => transparent
    minorticks: bool = True
    tick_major_size: float = 8.0
    tick_minor_size: float = 4.0
    tick_major_width: float = 1.0
    tick_minor_width: float = 0.5   # plain-axes fallback; WCS minor width follows major
    tick_direction: str = 'in'      # 'in', 'out', or 'inout'
    coord_style: str = 'sexagesimal'   # 'sexagesimal' or 'decimal'
    x_format: str = 'hh:mm:ss.ss'
    y_format: str = 'dd:mm:ss.ss'
    title: str = ''
    xlabel: str = ''   # empty => auto from FITS CTYPE
    ylabel: str = ''
    show_legend: bool = False        # channel legend (which color = which image)
    legend_loc: str = 'upper right'  # matplotlib legend location
    # Overlapping color-combination swatch (composite-mode-aware circles).
    show_combo_swatch: bool = False
    combo_swatch_loc: str = 'lower right'
    combo_swatch_labels: bool = False   # short labels next to each circle
    combo_swatch_label_offset: float = 0.62  # radial label push (fraction of radius)
    combo_swatch_inset_scale: float = 0.24   # corner inset size (fraction of axes)
    combo_swatch_size: int = 320             # swatch render resolution (pixels/side)
    # Per-band colored text labels (skyplothelper-style add_bandlabels).
    show_band_labels: bool = False
    band_labels_loc: str = 'upper left'
    # Optional skyplothelper overlays (require multicolorfits[overlays]).
    show_compass: bool = False
    compass_loc: str = 'lower left'
    show_beam: bool = False
    beam_loc: str = 'lower left'
    beam_style: str = 'crosshair'
    show_scale_bar: bool = False
    scale_bar_asec: float = 0.0   # 0 => auto (~15% of field width)
    scale_bar_loc: int = 4
    scale_bar_color: str = 'white'
    scale_bar_stroke_color: str = 'black'
    scale_bar_stroke_lw: float = 1.75
    # Compositing beyond classic RGB addition.
    combine_mode: str = 'rgb'       # see COMBINE_MODES
    combine_blend: str = 'screen'   # lightness blend when combine_mode in lab/hsv/hsl
    # Compositing background at zero signal: 'black' (classic), 'white',
    # 'transparent', or any matplotlib color.  A non-black value composites the
    # image over that background (coverage = per-pixel brightness) and replaces
    # the legacy inverse=True path (which complemented hues).  RYB/CMYK modes
    # use the same background setting.
    combine_background: str = 'black'
    # Image-only export: hide ticks, labels, title, and annotations.
    bare_plot: bool = False
    # Panel thumbnail previews: fast downsampled float32 by default; HD uses full resolution.
    panel_preview_hd: bool = False
    panel_preview_max_size: int = 512

    def to_dict(self):
        return {
            'gamma': _json_float(self.gamma, 2.2),
            'inverse': self.inverse,
            'tickcolor': self.tickcolor,
            'facecolor': self.facecolor,
            'minorticks': self.minorticks,
            'tick_major_size': _json_float(self.tick_major_size, 8.0),
            'tick_minor_size': _json_float(self.tick_minor_size, 3.0),
            'tick_major_width': _json_float(self.tick_major_width, 1.0),
            'tick_minor_width': _json_float(self.tick_minor_width, 0.8),
            'tick_direction': self.tick_direction,
            'coord_style': self.coord_style,
            'x_format': self.x_format,
            'y_format': self.y_format,
            'title': self.title,
            'xlabel': self.xlabel,
            'ylabel': self.ylabel,
            'show_legend': self.show_legend,
            'legend_loc': self.legend_loc,
            'show_combo_swatch': self.show_combo_swatch,
            'combo_swatch_loc': self.combo_swatch_loc,
            'combo_swatch_labels': self.combo_swatch_labels,
            'combo_swatch_label_offset': _json_float(self.combo_swatch_label_offset, 0.62),
            'combo_swatch_inset_scale': _json_float(self.combo_swatch_inset_scale, 0.24),
            'combo_swatch_size': int(_json_float(self.combo_swatch_size, 320) or 320),
            'show_band_labels': self.show_band_labels,
            'band_labels_loc': self.band_labels_loc,
            'show_compass': self.show_compass,
            'compass_loc': self.compass_loc,
            'show_beam': self.show_beam,
            'beam_loc': self.beam_loc,
            'beam_style': self.beam_style,
            'show_scale_bar': self.show_scale_bar,
            'scale_bar_asec': _json_float(self.scale_bar_asec, 0.0),
            'scale_bar_loc': int(_json_float(self.scale_bar_loc, 4) or 4),
            'scale_bar_color': self.scale_bar_color,
            'scale_bar_stroke_color': self.scale_bar_stroke_color,
            'scale_bar_stroke_lw': _json_float(self.scale_bar_stroke_lw, 0.0),
            'combine_mode': self.combine_mode,
            'combine_blend': self.combine_blend,
            'combine_background': self.combine_background,
            'bare_plot': self.bare_plot,
            'panel_preview_hd': self.panel_preview_hd,
            'panel_preview_max_size': int(_json_float(self.panel_preview_max_size, 512) or 512),
        }


class McfSession:
    """
    Full editing session: N image panels plus compose settings.

    Methods cover loading and aligning FITS layers, rendering the combined
    image / previews, saving FITS or figures, exporting transparent cutouts,
    and writing a standalone recreation script
    (:meth:`export_script`).

    Examples
    --------
    >>> s = McfSession(n_panels=3)  # doctest: +SKIP
    >>> s.load_files(['ir.fits', 'r.fits', 'b.fits'],
    ...              colors=['#BE599E', '#DEA215', '#77C0F9'])
    >>> rgb = s.render_combined()
    """

    def __init__(self, n_panels=DEFAULT_N_PANELS):
        n = int(n_panels)
        if n < MIN_PANELS or n > MAX_PANELS:
            raise ValueError('n_panels must be between %d and %d' % (MIN_PANELS, MAX_PANELS))
        self.panels: List[PanelState] = [PanelState() for _ in range(n)]
        self.compose = ComposeState()

    def reset(self):
        """Unload all panels and restore default compose settings.

        Keeps the current number of panel slots (cleared empty); call
        :meth:`set_n_panels` afterward if you want to shrink back to the
        default count.
        """
        for p in self.panels:
            p.clear()
        self.compose = ComposeState()

    # ---------- panels ----------

    def active_panels(self):
        """Panels that currently have data loaded."""
        return [p for p in self.panels if p.in_use and p.data is not None]

    def active_indices(self):
        """Indices of panels that currently have data loaded."""
        return [i for i, p in enumerate(self.panels) if p.in_use and p.data is not None]

    def set_n_panels(self, n):
        """
        Grow or shrink the panel list to exactly *n* slots.

        Growing appends empty :class:`PanelState` instances.  Shrinking
        drops trailing panels (loaded data on those indices is discarded).

        Returns
        -------
        int
            The new panel count.
        """
        n = int(n)
        if n < MIN_PANELS or n > MAX_PANELS:
            raise ValueError('n_panels must be between %d and %d' % (MIN_PANELS, MAX_PANELS))
        while len(self.panels) < n:
            self.panels.append(PanelState())
        while len(self.panels) > n:
            self.panels.pop()
        return len(self.panels)

    def add_panel(self, color=None, label=''):
        """
        Append an empty panel slot (up to :data:`MAX_PANELS`).

        Parameters
        ----------
        color : str or None
            Optional hex color for the new layer (default white).
        label : str
            Optional legend label.

        Returns
        -------
        int
            Index of the new panel.
        """
        if len(self.panels) >= MAX_PANELS:
            raise ValueError('At most %d panels are supported' % MAX_PANELS)
        panel = PanelState()
        if color is not None:
            panel.color = str(color)
        if label:
            panel.label = str(label)
        self.panels.append(panel)
        return len(self.panels) - 1

    def remove_panel(self, index):
        """
        Remove the panel at *index* (indices after it shift down).

        At least :data:`MIN_PANELS` panel must remain.

        Returns
        -------
        int
            The new panel count.
        """
        if len(self.panels) <= MIN_PANELS:
            raise ValueError('Need at least %d panel' % MIN_PANELS)
        idx = int(index)
        if idx < 0 or idx >= len(self.panels):
            raise IndexError('No such panel: %s' % index)
        del self.panels[idx]
        return len(self.panels)

    # ---------- grid alignment ----------

    def grid_report(self):
        """
        Check whether all loaded panels share a common pixel grid (so they can
        be combined directly).  The first loaded panel is the reference.

        Returns
        -------
        dict with keys:
            aligned : bool        -- True if all layers match (or < 2 loaded)
            reference : int|None  -- index of the reference panel
            mismatched : list     -- indices of panels that don't match
            active : list         -- indices of all loaded panels
            shapes : dict         -- {index: [ny, nx]} for loaded panels
        """
        idx = self.active_indices()
        shapes = {i: list(self.panels[i].data.shape) for i in idx}
        if len(idx) < 2:
            return {'aligned': True, 'reference': (idx[0] if idx else None),
                    'mismatched': [], 'active': idx, 'shapes': shapes}
        ref = idx[0]
        ref_hdr = self.panels[ref].header
        mismatched = []
        for i in idx[1:]:
            p = self.panels[i]
            try:
                ok = bool(ref_hdr.matches_grid(p.header))
            except Exception:
                ok = (self.panels[ref].data.shape == p.data.shape)
            if not ok:
                mismatched.append(i)
        return {'aligned': not mismatched, 'reference': ref,
                'mismatched': mismatched, 'active': idx, 'shapes': shapes}

    def align_panels(self, target='reference', reference=None, method='interp', scale=False):
        """
        Reproject all loaded panels onto a single common grid so they can be
        combined.  Preserves each panel's display settings.

        Requires the optional 'reproject' package.

        Parameters
        ----------
        target : str
            'reference' -- reproject every layer onto the reference panel's grid.
            'icrs' / 'galactic' / 'fk5' / 'fk4' / 'ecliptic' -- build the common
            grid in that celestial frame (from the reference header) first, so
            e.g. the result has Galactic north up.
        reference : int or None
            Index of the reference panel (defaults to the first loaded panel).
        method : str
            reproject_image method: 'interp' (default), 'spi', or 'kapteyn'.
        scale : bool
            Flux-conserving reprojection (see reproject_image).

        Returns
        -------
        dict
            {'ok', 'changed' (indices), 'target', 'reference', 'shape', 'message'}
        """
        from .core.stack_tools import align_stack

        idx = self.active_indices()
        if len(idx) < 2:
            return {'ok': True, 'changed': [], 'target': target, 'reference': None,
                    'shape': None, 'message': 'Fewer than two images loaded; nothing to align.'}
        ref = idx[0] if reference is None else int(reference)
        if ref not in idx:
            raise ValueError('Reference panel %d has no image loaded' % (ref + 1))
        frame = None if target in (None, 'reference') else str(target).lower()

        images = [(self.panels[i].data, self.panels[i].header.astropy) for i in idx]
        aligned = align_stack(images, reference=idx.index(ref), frame=frame,
                              method=method, scale=scale)
        for i, (data, hdr) in zip(idx, aligned):
            self.panels[i].set_reprojected(data, hdr)

        tgt_desc = ('reference panel %d' % (ref + 1) if frame is None
                    else "the '%s' frame" % frame)
        return {'ok': True, 'changed': idx, 'target': target, 'reference': ref,
                'shape': list(self.panels[idx[0]].data.shape),
                'message': 'Aligned %d layers to %s.' % (len(idx), tgt_desc)}

    # ---------- colors & palettes ----------

    def panel_colors(self):
        """Current hex colors of the loaded panels, in panel order."""
        return [self.panels[i].color for i in self.active_indices()]

    def legend_entries(self):
        """
        (color, label) pairs for the loaded panels, for the combined-figure
        channel legend.  Labels fall back to 'Image N' when unset.
        """
        entries = []
        for i in self.active_indices():
            p = self.panels[i]
            entries.append((p.color, p.label or ('Image %d' % (i + 1))))
        return entries

    def render_combo_swatch(self, size=None):
        """
        Render the composite-mode-aware color-combination swatch: one circle per
        active panel (in its channel color), arranged so they overlap, with the
        overlap regions showing the *honest* combined colors.

        Each circle is a synthetic full-signal "band" image pushed through the
        exact same colorize + combine pipeline as the displayed image (RGB / RYB /
        Lab / HSV, the blend, gamma, inverse, and background), so the mixed colors
        match what the user sees.  The area outside the circles is transparent so
        the swatch floats over the figure.

        Parameters
        ----------
        size : int or None
            Pixel size of the (square) swatch image.  ``None`` uses
            ``compose.combo_swatch_size`` (default 320).

        Returns
        -------
        dict or None
            ``{'rgba': (size,size,4) float array in [0,1], 'circles': [...],
            'size': size}`` where each circle entry is
            ``{'x', 'y', 'r', 'color', 'label'}`` in swatch pixel coords
            (origin lower-left).  Returns None if no panels are loaded.
        """
        if size is None:
            size = int(getattr(self.compose, 'combo_swatch_size', 320) or 320)
        size = max(64, min(int(size), 1024))
        idx = self.active_indices()
        if not idx:
            return None
        colors = [self.panels[i].color for i in idx]
        labels = [self.panels[i].label or ('Image %d' % (i + 1)) for i in idx]
        return combo_swatch(
            colors,
            mode=getattr(self.compose, 'combine_mode', 'rgb') or 'rgb',
            blend=getattr(self.compose, 'combine_blend', 'screen') or 'screen',
            gamma=self.compose.gamma,
            background=getattr(self.compose, 'combine_background', 'black') or 'black',
            inverse=self.compose.inverse,
            size=size,
            labels=labels,
        )

    def render_combo_swatch_colors(self, colors, labels=None, size=320):
        """Like :meth:`render_combo_swatch` but for an arbitrary color list (preview)."""
        colors = list(colors)
        if not colors:
            return None
        if labels is None:
            labels = ['' for _ in colors]
        return combo_swatch(
            colors,
            mode=getattr(self.compose, 'combine_mode', 'rgb') or 'rgb',
            blend=getattr(self.compose, 'combine_blend', 'screen') or 'screen',
            gamma=self.compose.gamma,
            background=getattr(self.compose, 'combine_background', 'black') or 'black',
            inverse=self.compose.inverse,
            size=size,
            labels=list(labels),
        )

    def palette_colors_for(self, name, n):
        """
        Resolve ``n`` hex colors for a palette name.  ``'perceptual'`` (or
        ``'suggest'`` / ``'auto'``) generates perceptually even colors; any
        other name looks up a curated palette, padding with perceptual
        suggestions if the palette has fewer than ``n`` colors.
        """
        if n < 1:
            return []
        if name in ('perceptual', 'suggest', 'auto'):
            return suggest_colors(n)
        if name in ('complement', 'triad', 'split', 'square', 'analogous', 'even'):
            return colors_for_hue_pattern(name, n)
        avail = get_palette(name)
        if len(avail) >= n:
            return list(avail[:n])
        return list(avail) + suggest_colors(n)[len(avail):]

    def apply_palette(self, name):
        """
        Assign colors from a curated (or perceptual) palette to the loaded
        panels, in panel order.

        Returns
        -------
        dict : {'colors', 'indices', 'name'}
        """
        idx = self.active_indices()
        colors = self.palette_colors_for(name, len(idx))
        for gi, c in zip(idx, colors):
            self.panels[gi].color = c
        return {'colors': colors, 'indices': idx, 'name': name}

    def rotate_panel_hues(self, delta_deg, base_colors=None):
        """
        Rotate active panel colors around the hue wheel while preserving
        relative spacing.  Used by the interactive hue-wheel tuner.

        Parameters
        ----------
        delta_deg : float
            Rotation in degrees.
        base_colors : list of str or None
            Starting colors (one per active panel).  Defaults to current
            panel colors.

        Returns
        -------
        dict : {'colors', 'indices', 'delta_deg'}
        """
        idx = self.active_indices()
        if not idx:
            return {'colors': [], 'indices': [], 'delta_deg': float(delta_deg)}
        if base_colors is None:
            base = [self.panels[i].color for i in idx]
        else:
            base = list(base_colors)
        colors = rotate_colors_from_base(base, float(delta_deg))
        for gi, c in zip(idx, colors):
            self.panels[gi].color = c
        return {'colors': colors, 'indices': idx, 'delta_deg': float(delta_deg)}

    def colorblind_report(self, min_distance=25.):
        """
        Check whether the loaded panels' colors stay distinguishable under
        color-vision deficiency (protan/deutan/tritan simulation).

        Returns
        -------
        dict with keys:
            ok : bool             -- True if every CVD type passes (or < 2 colors)
            colors : list         -- active panel hex colors
            indices : list        -- global indices of the active panels
            kinds : dict          -- {kind: {'ok': bool, 'failures': [(a, b, dE)]}}
                                     where a, b are 1-based panel numbers
        """
        idx = self.active_indices()
        colors = [self.panels[i].color for i in idx]
        report = {'ok': True, 'colors': colors, 'indices': idx, 'kinds': {}}
        if len(colors) < 2:
            for kind in CVD_KINDS:
                report['kinds'][kind] = {'ok': True, 'failures': []}
            return report
        for kind in CVD_KINDS:
            ok, failures = check_palette_colorblind(colors, kind=kind,
                                                    min_distance=min_distance)
            # Map active-list positions back to 1-based panel numbers
            mapped = [(idx[a] + 1, idx[b] + 1, round(d, 1)) for a, b, d in failures]
            report['kinds'][kind] = {'ok': ok, 'failures': mapped}
            if not ok:
                report['ok'] = False
        return report

    @property
    def common_header(self):
        """
        McfHeader of the first active panel (all images share a pixel grid).
        Use ``.astropy`` for the raw astropy Header when needed.
        """
        active = self.active_panels()
        return active[0].header if active else None

    # ---------- rendering ----------

    def render_combined(self, inverse=None, preview=False, dtype=None, max_size=None):
        """
        Combine all active panels into the final multicolor RGB image.

        Parameters
        ----------
        inverse : bool or None
            Override the session's compose.inverse flag if not None.
        preview : bool
            Fast preview mode: run the compose in float32 and downsample the
            source data (see max_size).  Convenience for interactive updates on
            large images; the final export should use preview=False.
        dtype : numpy dtype or None
            Explicit working dtype (overrides the preview default of float32).
        max_size : int or None
            Downsample source data so its longest side is <= max_size before
            compositing.  Defaults to 1024 when preview=True and unset.

        Returns
        -------
        array
            Combined RGB image [ypixels, xpixels, 3] in [0..1]
        """
        active = self.active_panels()
        if not active:
            raise ValueError('No images loaded in any panel')
        inverse = self.compose.inverse if inverse is None else inverse
        gamma = self.compose.gamma
        if preview:
            if dtype is None:
                dtype = np.float32
            if max_size is None:
                max_size = 1024
        _, _, eff_inverse = self._resolve_background(inverse)
        colorized = [p.render_color_rgb(gamma=gamma, inverse=eff_inverse,
                                        dtype=dtype, max_size=max_size)
                     for p in active]
        return self._combine_colorized(colorized, dtype=dtype, inverse=inverse)

    def _resolve_background(self, inverse=None):
        """
        Resolve the effective compositing background from compose state.

        Returns
        -------
        (bg_arg, use_bg, eff_inverse)
            bg_arg : the background passed to the compositors (None => transparent,
                     a color string otherwise; 'black' collapses to no-op here).
            use_bg : True when a non-black background must be composited/flattened.
            eff_inverse : the inverse flag actually applied during colorize
                          (a non-black background replaces the legacy complementary inverse).
        """
        inverse = self.compose.inverse if inverse is None else inverse
        bg_setting = getattr(self.compose, 'combine_background', 'black') or 'black'
        bg_key = str(bg_setting).strip().lower()
        use_bg = bg_key not in ('black', '')
        bg_arg = None if bg_key in ('transparent', 'none') else bg_setting
        eff_inverse = False if use_bg else inverse
        return bg_arg, use_bg, eff_inverse

    def _combine_colorized(self, colorized, dtype=None, inverse=None):
        """
        Combine a list of already-colorized RGB layers into the final image,
        applying the current compose mode / blend / gamma / inverse / background.

        Thin wrapper over ``core.combine_colorized_layers`` (the single source of
        truth for the compositing math), shared by ``render_combined`` (real image
        data) and ``render_combo_swatch`` (the synthetic overlapping-circle
        legend), so the legend's overlap colors match the displayed image.
        """
        inverse = self.compose.inverse if inverse is None else inverse
        return combine_colorized_layers(
            colorized,
            mode=getattr(self.compose, 'combine_mode', 'rgb') or 'rgb',
            blend=getattr(self.compose, 'combine_blend', 'screen') or 'screen',
            gamma=self.compose.gamma,
            background=getattr(self.compose, 'combine_background', 'black') or 'black',
            inverse=inverse,
            dtype=dtype if dtype is not None else np.float64,
        )

    # ---------- output ----------

    def _combined_output_header(self):
        """Common header annotated with combine-mode provenance for FITS output."""
        hdr = self.common_header
        out = hdr.astropy.copy() if hdr is not None else pyfits.Header()
        mode = getattr(self.compose, 'combine_mode', 'rgb') or 'rgb'
        blend = getattr(self.compose, 'combine_blend', 'screen') or 'screen'
        bg = getattr(self.compose, 'combine_background', 'black') or 'black'
        out['MCFTYPE'] = ('combined', 'multicolorfits output type')
        out['MCFMODE'] = (mode, 'color-mix space (rgb/lab/hsv/hsl/ryb/cmyk)')
        if mode in ('lab', 'hsv', 'hsl'):
            out['MCFBLEND'] = (blend, 'lightness blend across layers')
        out['MCFGAMMA'] = (float(self.compose.gamma), 'gamma used for the color mix')
        out['MCFINVRT'] = (bool(self.compose.inverse), 'legacy inverse (complement hues) mode')
        out['MCFBKG'] = (str(bg)[:68], 'compositing background at zero signal')
        # The saved cube is ALWAYS display-ready RGB regardless of mix space, so
        # it opens correctly as an RGB cube; the mode is recorded for provenance.
        out.add_comment('multicolorfits combined image: 3-plane, display-ready RGB '
                        'cube (planes = R,G,B in [0,1]).')
        if mode in ('lab', 'hsv', 'hsl'):
            out.add_comment("Colors were mixed in '%s' space (blend='%s'); the stored "
                            'planes are the resulting RGB and display normally.'
                            % (mode, blend))
        elif mode == 'ryb':
            out.add_comment("Colors were mixed subtractively in RYB (paint) space; the "
                            'stored planes are the resulting RGB and display normally.')
        elif mode == 'cmyk':
            out.add_comment("Colors were mixed subtractively in CMYK (print-ink) space; the "
                            'stored planes are the resulting RGB and display normally.')
        if str(bg).strip().lower() not in ('black', ''):
            out.add_comment("Composited over a '%s' background (transparency flattened "
                            'to white for this RGB cube).' % bg)
        annotate_provenance_header(out, info=self.to_dict())
        return out

    def save_rgb_fits(self, savepath, inverse=None, overwrite=True):
        """
        Save the combined image as a display-ready RGB FITS cube.

        The header records the color-mix mode/blend/gamma (``MCFMODE`` etc.) and
        HISTORY provenance.  The cube is standard RGB regardless of mix space, so
        it opens normally as an RGB cube in DS9 etc.
        """
        combined = self.render_combined(inverse=inverse)
        combined = np.asarray(combined)
        if combined.ndim == 3 and combined.shape[-1] == 4:
            # Transparent output: flatten over white so the RGB cube is opaque.
            a = combined[..., 3:4]
            combined = combined[..., :3] * a + (1.0 - a)
        save_rgb_fits(savepath, np.clip(combined, 0, 1), self._combined_output_header(),
                    overwrite=overwrite)

    def plot_component_mosaic(self, **kwargs):
        """
        Combined image plus colorized component frames (shared WCS).

        See :func:`~multicolorfits.figures.make_component_mosaic`.
        """
        from .figures import make_component_mosaic
        return make_component_mosaic(self, **kwargs)

    def export_transparent_cutout(self, size=None, savepath=None, inverse=None,
                                  max_size=None, render_background='auto',
                                  **kwargs):
        """
        Render the combined image and build a transparent-sky RGBA cutout.

        By default the *color* render uses the session's current background
        when it is classic black RGB (so the cutout matches what you see in
        the GUI / a press-release figure).  Only when the session is already
        on a non-black background, or ``render_background='transparent'`` is
        passed, is the compose temporarily switched to transparent — that
        path uses the alpha compositor and can shift colors relative to a
        black-background ``combine_multicolor`` result.

        Parameters
        ----------
        size : int or None
            Max side length in pixels for the cutout (passed through).
        savepath : str or None
            If set, write PNG/TIFF via :func:`~multicolorfits.core.cutout.save_transparent_cutout`.
        inverse : bool or None
            Override compose.inverse for this render only.
        max_size : int or None
            Optional downsample of the *source* render before the cutout matte
            (faster for huge images when you only need a small stamp).
        render_background : {'auto', 'session', 'transparent', 'black'}
            How to render colors before building the sky matte:

            * ``'auto'`` (default) — keep the session background when it is
              black (classic RGB math); otherwise render transparent.
            * ``'session'`` — always use the session's stored background.
            * ``'black'`` / ``'transparent'`` — force that background for
              this render only.
        **kwargs
            Passed to :func:`~multicolorfits.core.cutout.make_transparent_cutout`
            (``alpha_lo``, ``alpha_hi``, ``alpha_gamma``, ``crop``, ``matte``, …).

        Returns
        -------
        array
            Float RGBA ``(ny, nx, 4)`` in [0, 1].
        """
        from .core.cutout import make_transparent_cutout, save_transparent_cutout

        saved_bg = getattr(self.compose, 'combine_background', 'black') or 'black'
        mode = str(render_background or 'auto').strip().lower()
        if mode == 'auto':
            bg_key = str(saved_bg).strip().lower()
            # Classic black RGB: preserve combine_multicolor color math.
            # Non-black / transparent sessions already use the alpha path.
            render_bg = saved_bg if bg_key in ('black', '') else 'transparent'
        elif mode == 'session':
            render_bg = saved_bg
        elif mode in ('black', 'transparent'):
            render_bg = mode
        else:
            raise ValueError(
                "render_background must be 'auto', 'session', 'black', or "
                "'transparent' (got %r)" % render_background)

        try:
            self.compose.combine_background = render_bg
            combined = self.render_combined(inverse=inverse, max_size=max_size)
        finally:
            self.compose.combine_background = saved_bg
        combined = np.asarray(combined)
        if combined.ndim == 3 and combined.shape[-1] == 4:
            # Use RGB + replace alpha with the cutout matte by default.
            kwargs.setdefault('existing_alpha', 'replace')
            rgb = combined
        else:
            rgb = combined
        rgba = make_transparent_cutout(rgb, size=size, **kwargs)
        if savepath:
            save_transparent_cutout(rgba, savepath)
        return rgba

    def save_component_fits(self, index, savepath, kind='data',
                            inverse=None, overwrite=True):
        """
        Save a single loaded panel to FITS.

        Parameters
        ----------
        index : int
            Panel index (0-based).
        savepath : str
        kind : str
            ``'data'``  -> the panel's current 2D data + header.  Useful to
            capture a layer *after* in-GUI reprojection/alignment (the on-disk
            original no longer matches the displayed grid).
            ``'color'`` -> the panel's display-ready colorized RGB cube + header
            (the single-layer color as it contributes to the composite).
        inverse : bool or None
            Override compose.inverse (only affects ``kind='color'``).
        overwrite : bool
        """
        if index < 0 or index >= len(self.panels):
            raise ValueError('Panel index %d out of range' % index)
        p = self.panels[index]
        if not p.in_use or p.data is None:
            raise ValueError('Panel %d has no image loaded' % (index + 1))
        hdr = p.header.astropy.copy() if p.header is not None else pyfits.Header()
        hdr['MCFLAYER'] = (int(index + 1), 'multicolorfits source layer number')
        if p.label:
            hdr['MCFLABEL'] = (str(p.label)[:68], 'layer label')
        if p.reprojected:
            hdr.add_comment('Data reprojected/aligned in multicolorfits (grid differs '
                            'from the original on-disk file).')
        if kind == 'data':
            hdr['MCFTYPE'] = ('component-data', 'multicolorfits output type')
            hdr.add_history('multicolorfits layer %d: raw 2D data (stretch/color '
                            'NOT applied)' % (index + 1))
            pyfits.writeto(savepath, np.asarray(p.data), hdr, overwrite=overwrite)
        elif kind == 'color':
            gamma = self.compose.gamma
            inv = self.compose.inverse if inverse is None else inverse
            disp = p.render_display(gamma=gamma, inverse=inv)
            hdr['MCFTYPE'] = ('component-color', 'multicolorfits output type')
            hdr['MCFCOLOR'] = (str(p.color), 'layer color (hex)')
            hdr['MCFGAMMA'] = (float(gamma), 'gamma used for colorize')
            hdr['MCFINVRT'] = (bool(inv), 'inverse (white-background) mode')
            hdr.add_comment('multicolorfits component image: display-ready RGB cube '
                            '(single colorized layer, planes = R,G,B in [0,1]).')
            hdr.add_history('multicolorfits layer %d: %s scaled=%s [%.4g, %.4g] color=%s'
                            % (index + 1, p.label or '', p.stretch, p.vmin, p.vmax, p.color))
            annotate_provenance_header(hdr, info=self.to_dict())
            save_rgb_fits(savepath, disp, hdr, overwrite=overwrite)
        else:
            raise ValueError("kind must be 'data' or 'color', got %r" % kind)

    def fits_save_targets(self):
        """
        Describe the FITS outputs the GUIs can offer, in menu order: the combined
        RGB cube, plus for each loaded layer a raw 2D-data option and a
        display-ready colorized-RGB option.  Each entry is a dict with keys
        ``kind`` ('combined'|'data'|'color'), ``index`` (-1 for combined), and
        ``label``.
        """
        targets = [{'kind': 'combined', 'index': -1, 'label': 'Combined RGB cube'}]
        for i in self.active_indices():
            name = self.panels[i].label or ('Image %d' % (i + 1))
            targets.append({'kind': 'data', 'index': i,
                            'label': 'Layer %d: %s -- 2D data' % (i + 1, name)})
            targets.append({'kind': 'color', 'index': i,
                            'label': 'Layer %d: %s -- colorized RGB' % (i + 1, name)})
        return targets

    def save_fits_target(self, kind, index=-1, savepath=None, overwrite=True):
        """Dispatch a FITS save described by :meth:`fits_save_targets`."""
        if not savepath:
            raise ValueError('No save path given')
        if kind == 'combined':
            self.save_rgb_fits(savepath, overwrite=overwrite)
        elif kind in ('data', 'color'):
            self.save_component_fits(index, savepath, kind=kind, overwrite=overwrite)
        else:
            raise ValueError('Unknown FITS save kind %r' % kind)

    def params_text(self):
        """Human-readable summary of the current plot parameters."""
        lines = ['RGB Image plot params:']
        for i, p in enumerate(self.panels, start=1):
            if p.in_use:
                lines.append('image%i:  (%s)' % (i, p.filepath or 'array input'))
                lines.append('    vmin = %.3e , vmax = %.3e, scale = %s' % (p.vmin, p.vmax, p.stretch))
                lines.append("    image color = '%s'" % p.color)
                if p.smooth:
                    lines.append('    smoothing sigma = %.2f pixels' % p.smooth_sigma)
        lines.append("gamma = %.1f , tick color = '%s'" % (self.compose.gamma, self.compose.tickcolor))
        return '\n'.join(lines)

    def export_script(self):
        """
        Generate a standalone Python script that reproduces the current
        combined image with the multicolorfits scripting API.

        Full-fidelity: if any panel was aligned/reprojected in the GUI, the
        script embeds the common WCS grid (stripped to the minimal cards) and
        replays the ``reproject_image`` step so the exported result matches the
        on-screen image exactly.  When no reprojection was needed the script
        stays lean and just uses the first image's on-disk header for the WCS.
        """
        gamma = self.compose.gamma
        active = [(i, p) for i, p in enumerate(self.panels, start=1) if p.in_use]
        needs_repro = any(p.reprojected for _, p in active)
        common = self.common_header

        lines = [
            '#!/usr/bin/env python',
            '# Generated by multicolorfits -- recreate the current GUI state in a script.',
            'import astropy.io.fits as pyfits',
            'import multicolorfits as mcf',
            '',
        ]

        hdr_var = None
        if needs_repro and common is not None:
            # Panels were aligned/reprojected onto a common grid in the GUI.
            # Embed that grid's WCS (only the cards needed to define the pixel
            # grid -- no HISTORY/COMMENT clutter) and replay the reprojection.
            cards = [c.rstrip() for c in common.minimal_wcs_header().tostring(sep='\n').split('\n')]
            cards = [c for c in cards if c and not c.startswith('END')]
            lines.append('# Common WCS grid the panels were aligned/reprojected onto.')
            lines.append("# (Reprojection needs the 'reproject' package: pip install reproject)")
            lines.append('_common_hdr = pyfits.Header.fromstring("\\n".join([')
            for c in cards:
                lines.append('    %r,' % c)
            lines.append(']), sep="\\n")')
            lines.append('')
            hdr_var = '_common_hdr'

        varnames = []
        for i, p in active:
            var = 'img%i' % i
            varnames.append(var)
            path = p.filepath if p.filepath else '<path_to_image_%i.fits>' % i
            lines.append("%s_data, %s_hdr = pyfits.getdata('%s', header=True)" % (var, var, path))
            if needs_repro and p.reprojected and hdr_var is not None:
                lines.append("%s_data = mcf.reproject_image(%s_data, %s_hdr, %s)  # align onto common grid"
                             % (var, var, var, hdr_var))
            lines.append("%s_grey = mcf.to_grey_rgb(%s_data, rescalefn='%s', scaletype='abs', "
                         "min_max=[%.6g, %.6g], gamma=%.3g)" % (var, var, p.stretch, p.vmin, p.vmax, gamma))
            lines.append("%s_color = mcf.colorize_image(%s_grey, '%s', colorintype='hex', gammacorr_color=%.3g)"
                         % (var, var, p.color, gamma))
            if p.smooth:
                lines.append('%s_color = mcf.smooth_image(%s_color, sigma=%.3g)' % (var, var, p.smooth_sigma))
            lines.append('')
        color_list = ', '.join('%s_color' % v for v in varnames)
        mode = getattr(self.compose, 'combine_mode', 'rgb') or 'rgb'
        bg_setting = getattr(self.compose, 'combine_background', 'black') or 'black'
        bg_key = str(bg_setting).strip().lower()
        use_bg = bg_key not in ('black', '')
        bg_repr = 'None' if bg_key in ('transparent', 'none') else repr(bg_setting)
        # A non-black background supersedes the legacy inverse flag.
        inverse_arg = ', inverse=True' if (self.compose.inverse and not use_bg) else ''
        if mode == 'ryb':
            combined_call = ("mcf.combine_multicolor_alpha([%s], background=%s, "
                             "mode='ryb', gamma=%.3g)" % (color_list, bg_repr, gamma))
        elif mode == 'cmyk':
            combined_call = ("mcf.combine_multicolor_alpha([%s], background=%s, "
                             "mode='cmyk', gamma=%.3g)" % (color_list, bg_repr, gamma))
        elif mode == 'rgb':
            combined_call = 'mcf.combine_multicolor([%s], gamma=%.3g%s)' % (
                color_list, gamma, inverse_arg)
        else:
            blend = getattr(self.compose, 'combine_blend', 'screen') or 'screen'
            combined_call = ("mcf.combine_multicolor_colorspace([%s], colorspace='%s', "
                             "blend='%s', gamma=%.3g%s)" % (
                                 color_list, mode, blend, gamma, inverse_arg))
        lines.append('combined = %s' % combined_call)
        if mode != 'ryb' and mode != 'cmyk' and use_bg:
            lines.append('combined = mcf.composite_over_background(combined, background=%s)'
                         % bg_repr)
        lines.append('')
        first = varnames[0] if varnames else 'img1'
        if hdr_var is None:
            hdr_var = '%s_hdr' % first
        from .figures import contrast_color
        facecolor = getattr(self.compose, 'facecolor', 'white') or 'white'
        labelcolor = contrast_color(facecolor)
        legend_arg = ''
        if getattr(self.compose, 'show_legend', False):
            entries = self.legend_entries()
            entry_src = ', '.join('(%r, %r)' % (c, lbl) for c, lbl in entries)
            legend_arg = (", legend=[%s], legend_loc=%r, legend_inverse=%r"
                          % (entry_src, self.compose.legend_loc, bool(self.compose.inverse)))
        band_arg = ''
        if getattr(self.compose, 'show_band_labels', False):
            entries = self.legend_entries()
            entry_src = ', '.join('(%r, %r)' % (c, lbl) for c, lbl in entries)
            band_arg = (", band_labels=[%s], band_labels_loc=%r"
                        % (entry_src, self.compose.band_labels_loc))
        swatch_arg = ''
        if getattr(self.compose, 'show_combo_swatch', False):
            entries = self.legend_entries()
            colors_src = ', '.join('%r' % c for c, _ in entries)
            labels_src = ', '.join('%r' % lbl for _, lbl in entries)
            blend = getattr(self.compose, 'combine_blend', 'screen') or 'screen'
            # Build the honest, composite-mode-aware swatch the same way the GUI does.
            lines.append("_swatch = mcf.combo_swatch([%s], mode=%r, blend=%r, gamma=%.3g, "
                         "background=%s, inverse=%r, size=%d, labels=[%s])"
                         % (colors_src, mode, blend, gamma, bg_repr,
                            bool(self.compose.inverse),
                            int(getattr(self.compose, 'combo_swatch_size', 320) or 320),
                            labels_src))
            swatch_arg = (", swatch=_swatch, swatch_loc=%r, swatch_labels=%r, "
                          "swatch_label_offset=%.4g, swatch_inset_scale=%.4g"
                          % (self.compose.combo_swatch_loc,
                             bool(self.compose.combo_swatch_labels),
                             float(getattr(self.compose, 'combo_swatch_label_offset', 0.62)),
                             float(getattr(self.compose, 'combo_swatch_inset_scale', 0.24))))
        overlay_arg = ''
        if (getattr(self.compose, 'show_compass', False)
                or getattr(self.compose, 'show_beam', False)
                or getattr(self.compose, 'show_scale_bar', False)):
            overlay_arg = (", show_compass=%r, compass_loc=%r, show_beam=%r, beam_loc=%r, "
                           "beam_style=%r, show_scale_bar=%r, scale_bar_asec=%.4g, scale_bar_loc=%d"
                           % (bool(self.compose.show_compass), self.compose.compass_loc,
                              bool(self.compose.show_beam), self.compose.beam_loc,
                              self.compose.beam_style, bool(self.compose.show_scale_bar),
                              float(self.compose.scale_bar_asec or 0.0),
                              int(self.compose.scale_bar_loc)))
            lines.append("# Overlays require: pip install multicolorfits[overlays]")
        labelcolor = contrast_color(facecolor)
        axis_arg = ''
        if (self.compose.xlabel or '').strip() or (self.compose.ylabel or '').strip():
            xl = (self.compose.xlabel or '').strip() or 'RA'
            yl = (self.compose.ylabel or '').strip() or 'DEC'
            axis_arg = ", xaxislabel=%r, yaxislabel=%r" % (xl, yl)
        bare_arg = ', bare_plot=True' if getattr(self.compose, 'bare_plot', False) else ''
        lines.append("mcf.plot_combined_rgb(combined, %s, '%s', './multicolor_image.png', "
                     "tickcolor='%s', facecolor='%s', labelcolor='%s', dpi=150%s%s%s%s%s%s)"
                     % (hdr_var, self.compose.title, self.compose.tickcolor,
                        facecolor, labelcolor, axis_arg, legend_arg, band_arg, swatch_arg,
                        overlay_arg, bare_arg))
        lines.append('')
        lines.append("# Save the combined RGB image as a FITS cube (portable common WCS):")
        lines.append("# mcf.save_rgb_fits('./multicolor_image.fits', combined, %s)" % hdr_var)
        lines.append("# Save just the common WCS grid to a portable text header:")
        lines.append("# mcf.save_header(%s, './multicolor_image.hdr')" % hdr_var)
        lines.append('')
        return '\n'.join(lines)

    def to_dict(self):
        return {
            'panels': [p.to_dict() for p in self.panels],
            'compose': self.compose.to_dict(),
        }

    def data_shape(self):
        """(ny, nx) of the first loaded panel, or None."""
        for p in self.panels:
            if p.in_use and p.data is not None:
                return p.data.shape
        return None

    def apply_bare_plot_style(self, transparent=False):
        """
        Configure this session for publication / image-only output.

        See :func:`~multicolorfits.figures.apply_bare_plot_style`.
        """
        from .figures import apply_bare_plot_style
        apply_bare_plot_style(self.compose, transparent=transparent)
        return self

    def load_files(self, paths, colors=None, labels=None):
        """
        Load FITS paths into successive panels (growing the session if needed,
        up to :data:`MAX_PANELS`).

        Parameters
        ----------
        paths : sequence of str
            FITS file paths.
        colors : sequence of str or None
            Optional hex colors, one per path.
        labels : sequence of str or None
            Optional channel labels.

        Returns
        -------
        list of str
            Warnings for paths that could not be loaded.
        """
        paths = list(paths)
        colors = list(colors) if colors is not None else []
        labels = list(labels) if labels is not None else []
        if len(paths) > len(self.panels):
            self.set_n_panels(min(len(paths), MAX_PANELS))
        warnings = []
        limit = min(len(paths), len(self.panels))
        if len(paths) > limit:
            warnings.append('Only the first %d of %d paths fit within MAX_PANELS=%d'
                            % (limit, len(paths), MAX_PANELS))
        for i, path in enumerate(paths[:limit]):
            try:
                self.panels[i].load_fits(path)
                if i < len(colors):
                    self.panels[i].color = colors[i]
                if i < len(labels):
                    self.panels[i].label = labels[i]
            except Exception as exc:
                warnings.append('Panel %d: could not load %s (%s)' % (i + 1, path, exc))
        return warnings

    def sample_at_pixel(self, x, y):
        """
        Sample raw data values at integer pixel ``(x, y)`` (origin lower).

        Returns per-layer fluxes and sky coordinates when WCS is available.
        """
        xi, yi = int(round(x)), int(round(y))
        layers = []
        ref_hdr = None
        for i, p in enumerate(self.panels):
            if not p.in_use or p.data is None:
                continue
            ny, nx = p.data.shape
            if xi < 0 or yi < 0 or xi >= nx or yi >= ny:
                continue
            label = p.label or ('Image %d' % (i + 1))
            layers.append({
                'index': i,
                'label': label,
                'color': p.color,
                'value': _json_float(p.data[yi, xi]),
            })
            if ref_hdr is None and p.header is not None:
                ref_hdr = p.header
        sky = None
        if ref_hdr is not None:
            try:
                from .core.wcs_tools import pixel_to_sky
                ra, dec = pixel_to_sky(ref_hdr.astropy, xi, yi)
                sky = {'ra': float(ra), 'dec': float(dec)}
            except Exception:
                pass
        return {'x': xi, 'y': yi, 'layers': layers, 'sky': sky}

    def display_to_data_pixel(self, ix, iy, img_w, img_h):
        """
        Map a cursor position on a displayed RGB image (origin upper-left,
        size ``img_w`` x ``img_h``) to data pixel coordinates (origin lower).
        """
        shape = self.data_shape()
        if shape is None or img_w < 1 or img_h < 1:
            return None, None
        ny, nx = shape
        fx = float(ix) / max(img_w - 1, 1)
        fy = float(iy) / max(img_h - 1, 1)
        x = int(round(fx * (nx - 1)))
        y = int(round((1.0 - fy) * (ny - 1)))
        return x, y

    # ---------- save / restore state (JSON, DS9 backup-style) ----------

    STATE_VERSION = 1

    @staticmethod
    def _resolve_fits_path(filepath, search_dirs=None):
        """
        Locate a FITS file on disk, tolerating relative paths and basenames.

        Tries the path as given, then each directory in *search_dirs* (session
        JSON folder, ``_session_dir``, sibling panel directories, etc.).
        """
        if not filepath:
            return None
        fp = os.path.expanduser(str(filepath).strip())
        if os.path.isfile(fp):
            return os.path.abspath(fp)
        name = os.path.basename(fp)
        dirs = []
        seen = set()
        for d in (search_dirs or []):
            if not d:
                continue
            d = os.path.abspath(os.path.expanduser(d))
            if d not in seen and os.path.isdir(d):
                seen.add(d)
                dirs.append(d)
        for d in dirs:
            for cand in (os.path.normpath(os.path.join(d, fp)),
                         os.path.join(d, name)):
                if os.path.isfile(cand):
                    return os.path.abspath(cand)
        if name != fp and os.path.isfile(name):
            return os.path.abspath(name)
        return None

    @staticmethod
    def _session_paths_root(panels):
        """Common parent directory of loaded panel filepaths, if any."""
        loaded_dirs = []
        for p in panels:
            if not p.get('loaded'):
                continue
            fp = p.get('filepath', '')
            if not fp:
                continue
            fp = os.path.expanduser(fp)
            if os.path.isfile(fp):
                loaded_dirs.append(os.path.dirname(os.path.abspath(fp)))
            elif os.path.isabs(fp):
                loaded_dirs.append(os.path.dirname(os.path.abspath(fp)))
        if not loaded_dirs:
            return None
        try:
            return os.path.commonpath(loaded_dirs)
        except ValueError:
            return loaded_dirs[0]

    @staticmethod
    def _make_portable_path(filepath, session_dir):
        """
        Prefer a ``./``-relative path when *filepath* lives under *session_dir*.

        Returns an absolute path when the file is outside that directory (or
        when *session_dir* is unset).  Uses forward slashes so a shared JSON
        stays readable across platforms.
        """
        if not filepath:
            return filepath
        abs_fp = os.path.abspath(os.path.expanduser(filepath))
        if not session_dir:
            return abs_fp
        abs_dir = os.path.abspath(os.path.expanduser(session_dir))
        try:
            rel = os.path.relpath(abs_fp, abs_dir)
        except ValueError:
            return abs_fp
        if rel.startswith('..') or os.path.isabs(rel):
            return abs_fp
        rel = rel.replace('\\', '/')
        return './' + rel if not rel.startswith('./') else rel

    @staticmethod
    def _initial_search_dirs(state, base_dir=None, extra_search_dirs=None):
        dirs = []
        if base_dir:
            dirs.append(os.path.abspath(os.path.expanduser(base_dir)))
        # Prefer the JSON / share folder before any machine-local absolute root.
        root = state.get('_session_dir')
        if root:
            dirs.append(os.path.abspath(os.path.expanduser(root)))
        for p in state.get('panels', []):
            fp = p.get('filepath', '')
            if not fp:
                continue
            fp = os.path.expanduser(fp)
            if os.path.isabs(fp):
                d = os.path.dirname(fp)
                if os.path.isdir(d):
                    dirs.append(d)
        for d in (extra_search_dirs or []):
            if d:
                dirs.append(os.path.abspath(os.path.expanduser(d)))
        seen = set()
        out = []
        for d in dirs:
            if d not in seen:
                seen.add(d)
                out.append(d)
        return out

    def to_state_dict(self, session_dir=None):
        """
        Full session state as a JSON-serializable dict.

        Parameters
        ----------
        session_dir : str or None
            Directory the JSON will live in (or the folder you will share).
            FITS files under this directory are stored as ``./...`` relative
            paths so a tarball of ``session.json`` + FITS is portable.  Files
            outside it keep absolute paths.  When *session_dir* is omitted,
            the common parent of all loaded FITS paths is used when possible.
        """
        panels = [p.to_dict() for p in self.panels]
        root = None
        if session_dir:
            root = os.path.abspath(os.path.expanduser(session_dir))
        else:
            root = self._session_paths_root(panels)
        if root:
            for p in panels:
                if p.get('loaded') and p.get('filepath'):
                    p['filepath'] = self._make_portable_path(p['filepath'], root)
        out = {
            'multicolorfits_state_version': self.STATE_VERSION,
            'panels': panels,
            'compose': self.compose.to_dict(),
        }
        # Same-machine fallback only; load prefers the JSON's own directory.
        if root:
            out['_session_dir'] = root
        return out

    def save_state(self, path):
        """
        Save the full session state (file paths, scaling, colors, compose
        settings) to a JSON file, so a session can be restored later with
        load_state().

        FITS files that live in the same directory as *path* (or under it)
        are stored as ``./``-relative paths for easy sharing.  Other files
        keep absolute paths.  Pixel data is not embedded -- the FITS files
        must still be present when the state is restored.
        """
        import json
        path = os.path.abspath(os.path.expanduser(path))
        session_dir = os.path.dirname(path)
        with open(path, 'w') as f:
            json.dump(self.to_state_dict(session_dir=session_dir), f, indent=2)

    def load_state(self, path_or_dict, base_dir=None, extra_search_dirs=None):
        """
        Restore a session saved with :meth:`save_state` or a state dict.

        Panels whose FITS files no longer exist are left empty and reported.

        Parameters
        ----------
        path_or_dict : str, path-like, or dict
            JSON file path or a state dict (as from :meth:`to_state_dict`).
        base_dir : str or None
            Optional directory to search for relative FITS paths (defaults to
            the JSON file's directory when loading from a path).

        Returns
        -------
        list of str
            Warnings for anything that could not be restored (empty list if
            everything loaded cleanly).
        """
        import json
        import os
        if isinstance(path_or_dict, dict):
            state = path_or_dict
        else:
            state_path = os.path.abspath(os.path.expanduser(path_or_dict))
            with open(state_path) as f:
                state = json.load(f)
            if base_dir is None:
                base_dir = os.path.dirname(state_path)
        warnings = []
        c = state.get('compose', {})
        for key, default in ComposeState().to_dict().items():
            if key in c:
                setattr(self.compose, key, c[key])
        search_dirs = self._initial_search_dirs(state, base_dir=base_dir,
                                                extra_search_dirs=extra_search_dirs)
        panel_states = list(state.get('panels', []))
        if panel_states:
            # Match the saved slot count (clamped to MAX_PANELS) so restored
            # sessions with more than the default four panels keep their tabs.
            want = min(len(panel_states), MAX_PANELS)
            if len(panel_states) > MAX_PANELS:
                warnings.append('State has more panels (%d) than MAX_PANELS=%d; extras ignored'
                                % (len(panel_states), MAX_PANELS))
            self.set_n_panels(want)
        for i, pstate in enumerate(panel_states[:len(self.panels)]):
            panel = self.panels[i]
            panel.clear()
            if not pstate.get('loaded'):
                continue
            fp = pstate.get('filepath', '')
            resolved = self._resolve_fits_path(fp, search_dirs=search_dirs)
            if not resolved:
                warnings.append('Panel %d: file not found: %s' % (i + 1, fp))
                continue
            try:
                panel.load_fits(resolved)
            except Exception as exc:
                warnings.append('Panel %d: could not load %s (%s)' % (i + 1, fp, exc))
                continue
            panel.stretch = pstate.get('stretch', panel.stretch)
            panel.color = pstate.get('color', panel.color)
            panel.label = pstate.get('label', panel.label)
            panel.smooth = bool(pstate.get('smooth', panel.smooth))
            panel.smooth_sigma = float(pstate.get('smooth_sigma', panel.smooth_sigma))
            if 'percent_min' in pstate or 'percent_max' in pstate:
                panel.set_percentiles(pstate.get('percent_min'), pstate.get('percent_max'))
            elif 'vmin' in pstate or 'vmax' in pstate:
                panel.set_limits(pstate.get('vmin'), pstate.get('vmax'))
            d = os.path.dirname(resolved)
            if d and d not in search_dirs:
                search_dirs.append(d)
        return warnings
