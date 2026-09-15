"""
MultiColorFits
==============

Colorize and combine FITS images to produce visually aesthetic scientific
plots -- with any number of image layers, in any colors.

Scripting quickstart (unchanged from v2.x)::

    import multicolorfits as mcf

    grey = mcf.to_grey_rgb(data, rescalefn='asinh', min_max=[0., 1.2])
    colr = mcf.colorize_image(grey, '#C11B17', colorintype='hex')
    combined = mcf.combine_multicolor([colr, ...], gamma=2.2)

GUIs (all optional):

* ``mcf.gui()`` or the ``mcf-web`` console script -- browser-based GUI
  (requires ``pip install multicolorfits[web]``)
* ``mcf.gui_qt()`` or the ``mcf-qt`` console script -- PySide6 desktop GUI
  (requires ``pip install multicolorfits[qt]``)

New here (human or AI agent)?  Run ``mcf.overview()`` for the mental model,
conventions, and a task index, then ``mcf.recipes('<keyword>')`` for
copy-paste code.  The same catalog is published as ``llms.txt`` /
``llms-full.txt`` for agent ingestion.
"""

__author__ = "Phil Cigan"
__version__ = "3.1.0"

# Agent- / newcomer-facing orientation (mcf.overview() / mcf.recipes()).
from ._overview import overview, recipes  # noqa: E402,F401

# Re-export the entire computational API under the historical flat namespace,
# so that `import multicolorfits as mcf` scripts from v2.x keep working.
from .core import (  # noqa: F401
    # scaling
    stretch_functions,
    SIGNED_STRETCHES,
    adjust_gamma,
    draw_progress_bar,
    nan_percentile_of_score,
    make_norm,
    rescale_image,
    zscale_limits,
    suggest_levels,
    describe_image,
    describe_images,
    # header / WCS / coordinates / cropping
    McfHeader,
    as_mcfheader,
    save_header,
    load_header,
    force_header_2d,
    force_header_3d,
    force_header_floats,
    deg2dms,
    dms2deg,
    deg2hour,
    hour2deg,
    dec2sex,
    sex2dec,
    angular_distance,
    get_cdelts,
    sky_to_pixel,
    pixel_to_sky,
    make_simple_header,
    get_cd_matrix,
    deg_per_pixel,
    arcsec_per_pixel,
    sterad_per_pixel,
    squeeze_image,
    header_coord_grids,
    beam_params_arcsec,
    pixels_per_beam,
    crop_image,
    crop_cube,
    crop_image_sky,
    crop_cube_sky,
    tidy_header,
    wcs_is_flipped,
    east_increases_right,
    describe_header,
    blank_missing,
    crop_to_overlap,
    # reprojection (lazy optional deps)
    reproject_image,
    reproject_cube,
    # Colorize / tint greyscale → colorized layers
    hex_to_rgb,
    rgb_to_hex,
    hex_to_hsv,
    rgb_to_hsv,
    hsv_to_rgb,
    hex_complement,
    to_hex,
    to_grey_rgb,
    colorize_image_direct_rgb,
    colorize_image,
    combine_multicolor,
    smooth_image,
    # Output helpers
    plot_combined_rgb,
    compare_multicolor_vs_rgb,
    save_rgb_fits,
    # Perceptual compositing (Lab supported; HSV/HSL experimental) + accessibility
    rgb_to_hsl,
    hsl_to_rgb,
    rgb_to_lab,
    lab_to_rgb,
    combine_multicolor_colorspace,
    mix_colors_hex,
    greyscale_image,
    simulate_colorblindness,
    # Subtractive RYB/CMYK + background / transparent compositing
    rgb_to_ryb,
    ryb_to_rgb,
    hex_to_ryb,
    ryb_to_hex,
    rgb_to_cmyk,
    cmyk_to_rgb,
    hex_to_cmyk,
    cmyk_to_hex,
    layer_coverage,
    composite_over_background,
    combine_multicolor_alpha,
    # Composite-mode swatch + shared combine dispatch
    combine_colorized_layers,
    swatch_layout,
    combo_swatch,
    # Transparent cutouts for slides / stamps
    luminance,
    alpha_from_intensity,
    make_transparent_cutout,
    save_transparent_cutout,
    make_checkerboard,
    preview_cutout_on_backgrounds,
    batch_transparent_cutouts,
    # Deblend a flattened composite back to RGBA (and re-flatten)
    flatten_rgba,
    deblend_background,
    # Celestial frame conversion (optional reproject)
    header_frame_name,
    convert_header_frame,
    reproject_to_frame,
    reproject_to_galactic,
    optimal_common_header,
    make_rotated_header,
    make_north_up_header,
    reproject_to_rotation,
    reproject_north_up,
    # Stack alignment and interactive preview downsampling
    reproject_stack_to_header,
    reproject_stack_to_reference,
    align_stack,
    downsample_for_preview,
    prep_layers,
    convolution_beam,
    match_beam,
    match_beam_to_header,
    convolve2Dgaus,
    convolve2Dgaus_matchhdr,
    annotate_provenance_header,
    read_fits,
    write_fits,
)

# Backward-compatible v2.x names.  Every legacy alias lives in compat.py
# (one place to corral for eventual removal); LEGACY_ALIASES maps old -> new.
from .compat import *  # noqa: F401,F403
from .compat import LEGACY_ALIASES  # noqa: F401

from .session import (  # noqa: F401
    PanelState, ComposeState, McfSession, reproject_available, ALIGN_FRAMES,
    STRETCHES, COMBINE_MODES, DEFAULT_N_PANELS, MIN_PANELS, MAX_PANELS,
)
from .core.colormix import COLORSPACES, BLENDS  # noqa: F401
from . import palettes  # noqa: F401
from .palettes import get_palette, list_palettes, suggest_colors  # noqa: F401
from .palettes import (  # noqa: F401
    resolve_palette_colors, palette_colorblind_report, CVD_KINDS,
    colors_from_hsv, colors_from_hue_angles, colors_for_hue_pattern,
)
from .pipeline import combine_layers, combine_from_files, save_combined  # noqa: F401
from .figures import (  # noqa: F401
    apply_bare_plot_style, apply_bare_axes, make_component_mosaic,
    make_combined_figure, setup_combined_axes,
)
from .palette_preview import preview_palette, PalettePreview  # noqa: F401
from . import overlays  # noqa: F401


def overlays_available():
    """True when the optional skyplothelper overlays extra is installed."""
    from .overlays import overlays_available as _fn
    return _fn()


def gui(host='127.0.0.1', port=8321, open_browser=True, browser=None,
        session=None, state=None, files=None, colors=None, labels=None):
    """
    Launch the (browser-based) multicolorfits GUI.

    Parameters
    ----------
    host, port
        Local uvicorn bind address (default ``127.0.0.1:8321``; port auto-
        increments if busy).
    open_browser : bool
        If True (default), open a browser tab when the server is ready.
    browser : str or None
        Which browser to open when ``open_browser`` is True.  Names understood
        by :func:`webbrowser.get` work (e.g. ``'firefox'``, ``'google-chrome'``,
        ``'chromium'``, ``'safari'``), as does a full path to an executable.
        ``None`` uses the system default; on Unix the ``BROWSER`` environment
        variable is also consulted by the stdlib.  Ignored when
        ``open_browser=False``.
    session : McfSession or None
        Optional pre-configured session (panels + compose state).
    state : str or dict or None
        JSON session file path or state dict to restore on startup
        (see :meth:`McfSession.save_state` / :meth:`McfSession.load_state`).
    files : sequence of str or None
        FITS paths to load into successive panels on startup.
    colors, labels : sequences or None
        Optional per-file colors and channel labels (with ``files``).

    Examples
    --------
    >>> mcf.gui()  # doctest: +SKIP
    >>> mcf.gui(browser='firefox')  # doctest: +SKIP
    >>> mcf.gui(open_browser=False)  # then open the printed URL yourself

    CLI equivalents: ``mcf-web --browser firefox``, ``mcf-web --no-browser``.

    Requires the web extra:  pip install multicolorfits[web]
    """
    try:
        from .gui_web.launcher import run
    except ImportError as exc:
        raise ImportError(
            "The browser GUI requires the 'web' extra dependencies. "
            "Install them with:  pip install multicolorfits[web]"
        ) from exc
    run(host=host, port=port, open_browser=open_browser, browser=browser,
        session=session, state=state, files=files, colors=colors, labels=labels)


def start_gui_server(host='127.0.0.1', port=8321, session=None,
                     state=None, files=None, colors=None, labels=None, **kwargs):
    """
    Start the web GUI in a background thread (non-blocking).

    Returns a :class:`~multicolorfits.gui_web.embed.GuiServer` with ``url``,
    ``port``, and ``host``.  See :func:`gui_embed` for notebook embedding.

    Requires the web extra:  pip install multicolorfits[web]
    """
    try:
        from .gui_web.embed import start_gui_server as _start
    except ImportError as exc:
        raise ImportError(
            "The browser GUI requires the 'web' extra dependencies. "
            "Install them with:  pip install multicolorfits[web]"
        ) from exc
    return _start(host=host, port=port, session=session, state=state,
                  files=files, colors=colors, labels=labels, **kwargs)


def gui_embed(host='127.0.0.1', port=8321, session=None,
              state=None, files=None, colors=None, labels=None, **kwargs):
    """
    Embed the browser GUI in Jupyter, Colab, VS Code notebooks, or Binder.

    Starts uvicorn in a background thread and displays an inline iframe (or
    Colab's port proxy).  Use ``display='url'`` to print the URL only, or
    ``return_url=True`` to get the URL string without displaying.

    Parameters
    ----------
    host, port, session, state, files, colors, labels
        Same as :func:`gui`.
    display : {'auto', 'iframe', 'colab', 'url', 'none'}
        Embedding mode; ``auto`` picks Colab vs Jupyter vs URL print.
    height, width
        iframe size (``height`` also used for Colab).
    proxy_url : str or None
        Override iframe ``src`` for custom JupyterHub proxy layouts.
    return_url, print_url
        URL-only helpers for scripting and debugging.

    Examples
    --------
    ::

        import multicolorfits as mcf
        s = mcf.McfSession()
        s.load_files(['a.fits', 'b.fits'], colors=['#f00', '#0ff'])
        mcf.gui_embed(session=s, height=950)

    Requires the web extra:  pip install multicolorfits[web]
    """
    try:
        from .gui_web.embed import gui_embed as _embed
    except ImportError as exc:
        raise ImportError(
            "The browser GUI requires the 'web' extra dependencies. "
            "Install them with:  pip install multicolorfits[web]"
        ) from exc
    return _embed(host=host, port=port, session=session, state=state,
                  files=files, colors=colors, labels=labels, **kwargs)


def detect_notebook_environment():
    """
    Detect Jupyter / Colab / VS Code / Binder, or ``None`` outside IPython.

    Used by :func:`gui_embed` when ``display='auto'``.
    """
    try:
        from .gui_web.embed import detect_notebook_environment as _detect
    except ImportError as exc:
        raise ImportError(
            "Notebook helpers require the 'web' extra:  pip install multicolorfits[web]"
        ) from exc
    return _detect()


def mcf_gui(*args, **kwargs):
    """Alias for :func:`gui` (backward compatible with v2.x scripts)."""
    gui(*args, **kwargs)


def gui_qt(session=None, state=None, files=None, colors=None, labels=None):
    """
    Launch the PySide6 desktop GUI.

    Parameters
    ----------
    session : McfSession or None
        Optional pre-configured session.  When ``state`` or ``files`` is given,
        a fresh session is built unless ``session`` is also passed.
    state : str or dict or None
        JSON session file or dict to restore on startup.
    files : sequence of str or None
        FITS paths to load into successive panels on startup.
    colors, labels : sequences or None
        Optional per-file colors and labels.

    Requires the qt extra:  pip install multicolorfits[qt]
    """
    try:
        from .qt_gui.app import main as qt_main
    except ImportError as exc:
        raise ImportError(
            "The desktop GUI requires the 'qt' extra dependencies. "
            "Install them with:  pip install multicolorfits[qt]"
        ) from exc
    if session is None:
        session = McfSession()
        if state:
            session.load_state(state)
        elif files:
            session.load_files(files, colors=colors, labels=labels)
    return qt_main(session=session)
