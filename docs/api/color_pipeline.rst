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
