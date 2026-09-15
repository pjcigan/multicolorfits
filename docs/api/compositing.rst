Compositing
===========

Perceptual and subtractive mixers, backgrounds, and the shared dispatcher.
Concepts: :doc:`/guide/color_compositing`. Job snippets:
:doc:`/capabilities/index`. Recipes: ``mcf.recipes('lab')``,
``mcf.recipes('cutout')``.

Prefer :func:`~multicolorfits.combine_colorized_layers` when you want the
same ``mode`` / ``blend`` / ``background`` switch the GUI and session use.
Transparent cutouts and deblend helpers are listed at the bottom of this
page; narrative: :doc:`/guide/saving_outputs`,
:doc:`/examples/m74_transparent_cutout`.

.. currentmodule:: multicolorfits

Modes & dispatcher
------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   combine_colorized_layers
   combine_multicolor_colorspace
   combine_multicolor_alpha
   COMBINE_MODES
   COLORSPACES
   BLENDS

Backgrounds & paint mix
-----------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   composite_over_background
   layer_coverage
   rgb_to_ryb
   ryb_to_rgb
   rgb_to_cmyk
   cmyk_to_rgb

Lab / HSV / HSL helpers
-----------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   rgb_to_lab
   lab_to_rgb
   rgb_to_hsl
   hsl_to_rgb
   mix_colors_hex
   greyscale_image
   simulate_colorblindness

Swatches
--------

.. autosummary::
   :toctree: generated
   :nosignatures:

   combo_swatch
   swatch_layout

Transparent cutouts
-------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   make_transparent_cutout
   save_transparent_cutout
   preview_cutout_on_backgrounds
   batch_transparent_cutouts
   deblend_background
   flatten_rgba
   luminance
   alpha_from_intensity
   make_checkerboard
