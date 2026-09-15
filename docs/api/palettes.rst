Palettes
========

Curated color sets, generators, and color-vision helpers for layer tinting.

* Guide (how to generate colors, GUI map, ``plt.show`` tips):
  :doc:`/guide/palettes`
* Job snippets: :doc:`/capabilities/index`
* Recipe: ``mcf.recipes('palette')``

Layer hex colors are separate from compositing ``mode`` (``lab`` / ``hsv`` /
…). Prefer :func:`~multicolorfits.suggest_colors` (CIE LCh) for auto-N;
:func:`~multicolorfits.colors_from_hsv` is the classical HSV wheel when you
want it explicitly.

.. currentmodule:: multicolorfits

.. autosummary::
   :toctree: generated
   :nosignatures:

   get_palette
   list_palettes
   suggest_colors
   colors_from_hsv
   preview_palette

Related package utilities (also under ``mcf.palettes``):

.. currentmodule:: multicolorfits.palettes

.. autosummary::
   :toctree: generated
   :nosignatures:

   list_hue_patterns
   list_palette_menu
   colors_from_hue_angles
   colors_for_hue_pattern
   rotate_colors_from_base
   resolve_palette_colors
   check_palette_colorblind
   palette_colorblind_report
   CVD_KINDS

:func:`~multicolorfits.preview_palette` builds a pyplot figure with tiles, a
mode-aware combo swatch, and CVD tile rows / warnings. After the call,
``plt.show()`` or ``res.show()`` displays it.
