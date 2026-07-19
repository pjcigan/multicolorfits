"""
Backward-compatibility aliases for the historical (v2.x) flat API.

All old public names live here, in one place, so they are easy to corral
(and eventually remove).  Each alias is a thin wrapper around the renamed
function -- fully functional, no warnings emitted (yet) -- whose docstring
points at the new name.  ``LEGACY_ALIASES`` maps every old name to its new
name for introspection, migration tooling, and the docs migration table.

Old scripts keep working unchanged::

    import multicolorfits as mcf
    grey = mcf.greyRGBize_image(data, rescalefn='asinh')   # -> mcf.to_grey_rgb

Planned lifecycle: silent aliases in v3.0; a ``DeprecationWarning`` may be
added in a later 3.x release; removal no earlier than v4.0.
"""

import functools

from . import core as _core

#: Every renamed public symbol: old (v2.x) name -> new (v3) name.
LEGACY_ALIASES = {
    # colorize / color pipeline
    'greyRGBize_image': 'to_grey_rgb',
    'hexinv': 'hex_complement',
    'rgb_to_hsv_vectorized': 'rgb_to_hsv',
    'hsv_to_rgb_vectorized': 'hsv_to_rgb',
    # perceptual / subtractive colorspaces
    'rgb_to_hsl_vectorized': 'rgb_to_hsl',
    'hsl_to_rgb_vectorized': 'hsl_to_rgb',
    'rgb_to_ryb_vectorized': 'rgb_to_ryb',
    'ryb_to_rgb_vectorized': 'ryb_to_rgb',
    'rgb_to_cmyk_vectorized': 'rgb_to_cmyk',
    'cmyk_to_rgb_vectorized': 'cmyk_to_rgb',
    'image_value': 'layer_coverage',
    # output helpers
    'plotsinglemulticolorRGB': 'plot_combined_rgb',
    'comparemulticolorRGB_pureRGB': 'compare_multicolor_vs_rgb',
    'saveRGBfits': 'save_rgb_fits',
    # scaling
    'scaling_fns': 'stretch_functions',
    'drawProgressBar': 'draw_progress_bar',
    'nanpercofscore': 'nan_percentile_of_score',
    # header / WCS
    'makesimpleheader': 'make_simple_header',
    'getcdelts': 'get_cdelts',
    'getcdmatrix': 'get_cd_matrix',
    'getdegperpix': 'deg_per_pixel',
    'getasecperpix': 'arcsec_per_pixel',
    'getsteradperpix': 'sterad_per_pixel',
    'convsky2pix': 'sky_to_pixel',
    'convpix2sky': 'pixel_to_sky',
    'angulardistance': 'angular_distance',
    'beampars_asec_fromhdr': 'beam_params_arcsec',
    'pixperbeam_from_hdr': 'pixels_per_beam',
    'force_hdr_to_2D': 'force_header_2d',
    'force_hdr_to_3D': 'force_header_3d',
    'force_hdr_floats': 'force_header_floats',
    # cropping / reprojection
    'cropfits2D': 'crop_image',
    'cropfits3D': 'crop_cube',
    'cropfits2D_coords': 'crop_image_sky',
    'cropfits3D_coords': 'crop_cube_sky',
    'reproject2D': 'reproject_image',
    'reproject3D': 'reproject_cube',
}


def _make_alias(old_name, new_name):
    """Build a functional alias for ``new_name`` under the legacy ``old_name``."""
    target = getattr(_core, new_name)
    if not callable(target):
        # Plain objects (e.g. the scaling_fns dict) are re-exported directly.
        return target

    @functools.wraps(target)
    def alias(*args, **kwargs):
        return target(*args, **kwargs)

    alias.__name__ = old_name
    alias.__qualname__ = old_name
    alias.__doc__ = (
        "Deprecated alias for :func:`multicolorfits.%s` (renamed in v3.0).\n\n"
        "This alias remains fully functional for backward compatibility with\n"
        "v2.x scripts, but new code should call ``%s`` directly.\n\n"
        "%s" % (new_name, new_name, target.__doc__ or '')
    )
    return alias


for _old, _new in LEGACY_ALIASES.items():
    globals()[_old] = _make_alias(_old, _new)
del _old, _new

__all__ = ['LEGACY_ALIASES'] + list(LEGACY_ALIASES)
