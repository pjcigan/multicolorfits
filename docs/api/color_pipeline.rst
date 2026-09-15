Color pipeline
==============

Stretch greyscale FITS data, tint each layer, and sum channels the classic
way.  Narrative: :doc:`/guide/intensity_scaling`,
:doc:`/guide/color_compositing`.

.. currentmodule:: multicolorfits

Stretch & rescale
-----------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   rescale_image
   zscale_limits
   suggest_levels
   describe_image
   describe_images
   make_norm
   adjust_gamma
   stretch_functions
   SIGNED_STRETCHES
   STRETCHES
   nan_percentile_of_score
   draw_progress_bar

The built-in :func:`~multicolorfits.draw_progress_bar` is a zero-dependency
ANSI stdout bar (used by :func:`~multicolorfits.reproject_cube` when
``print_progress=True``).  For richer progress in your own loops, use
`tqdm <https://tqdm.github.io/>`_ directly — see the function docstring
for a side-by-side example.  Multicolorfits does not depend on or wrap
tqdm.

:func:`~multicolorfits.suggest_levels` recommends a starting stretch and
absolute vmin/vmax from the pixel histogram. It is not a finished display:
a single stretch cannot show faint structure and bright peaks across many
decades. The returned ``reason`` states the rule that fired (negative
pixels, background on the floor, decade span, or a thin bright tail).
:func:`~multicolorfits.describe_image` returns that suggestion after an
optional header summary and prints both by default (``verbose=False`` to
keep the dict quiet). :func:`~multicolorfits.describe_images` does the
same for a set of layers and adds suggested colors plus parallel lists
for ``load_files`` / ``set_display``. The full rules are in the function
docstring and in :doc:`/guide/intensity_scaling`. The GUI **Auto levels**
button applies the same suggestion and leaves the fields editable.
**Zscale** still sets limits only.

Colorize
--------

.. autosummary::
   :toctree: generated
   :nosignatures:

   to_grey_rgb
   colorize_image
   colorize_image_direct_rgb
   smooth_image
   hex_to_rgb
   rgb_to_hex
   hex_to_hsv
   hex_complement
   to_hex

Combine (classic RGB)
---------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   combine_multicolor
   plot_combined_rgb
   compare_multicolor_vs_rgb
   save_rgb_fits
