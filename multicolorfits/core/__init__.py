"""
multicolorfits.core -- GUI-free computational core.

All functions here work headless with only numpy/scipy/astropy/scikit-image
(matplotlib is used lazily by the plotting helpers in io_output).
"""

from .scaling import (
    stretch_functions,
    SIGNED_STRETCHES,
    adjust_gamma,
    draw_progress_bar,
    nan_percentile_of_score,
    make_norm,
    rescale_image,
    zscale_limits,
)
from .wcs_tools import (
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
)
from .mcfheader import (
    McfHeader,
    as_mcfheader,
    save_header,
    load_header,
)
from .reproject_tools import (
    reproject_image,
    reproject_cube,
)
from .colorize import (
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
)
from .io_output import (
    plot_combined_rgb,
    compare_multicolor_vs_rgb,
    save_rgb_fits,
    annotate_provenance_header,
)
# Experimental additions (see module docstrings)
from .colormix import (
    rgb_to_hsl,
    hsl_to_rgb,
    rgb_to_lab,
    lab_to_rgb,
    combine_multicolor_colorspace,
    mix_colors_hex,
    greyscale_image,
    simulate_colorblindness,
)
from .skyframes import (
    header_frame_name,
    convert_header_frame,
    reproject_to_frame,
    reproject_to_galactic,
    optimal_common_header,
)
from .stack_tools import (
    reproject_stack_to_header,
    reproject_stack_to_reference,
    align_stack,
    downsample_for_preview,
)
from .paintmix import (
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
)
from .swatch import (
    combine_colorized_layers,
    swatch_layout,
    combo_swatch,
)
from .cutout import (
    luminance,
    alpha_from_intensity,
    make_transparent_cutout,
    save_transparent_cutout,
    make_checkerboard,
    preview_cutout_on_backgrounds,
    batch_transparent_cutouts,
    flatten_rgba,
    deblend_background,
)
