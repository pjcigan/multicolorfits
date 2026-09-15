Figures
=======

Combined WCS figures, publication styling, and component mosaics.
Guide: :doc:`/guide/figures_and_mosaics`, :doc:`/guide/saving_outputs`,
:doc:`/guide/overlays`. Job snippets: :doc:`/capabilities/index`.
Recipes: ``mcf.recipes('mosaic')``, ``mcf.recipes('figure')``.

``make_combined_figure`` returns a figure only (``fig.axes[0]`` is the hero).
``make_component_mosaic`` returns ``(fig, axes)`` with
``axes['combined']`` / ``axes['components']``.

.. currentmodule:: multicolorfits

Builders
--------

.. autosummary::
   :toctree: generated
   :nosignatures:

   make_combined_figure
   make_component_mosaic
   setup_combined_axes
   apply_bare_plot_style
   apply_bare_axes

Annotations (in core; no ``[overlays]`` extra)
----------------------------------------------

.. currentmodule:: multicolorfits.figures

.. autosummary::
   :toctree: generated
   :nosignatures:

   add_swatch_inset
   add_band_labels
