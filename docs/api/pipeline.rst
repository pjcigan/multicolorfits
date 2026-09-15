High-level pipeline
===================

One-shot helpers for scripts that skip the GUI / session model.

* Concepts: :doc:`/guide/color_compositing`, :doc:`/guide/preparing_images`,
  :doc:`/guide/wcs_and_alignment`
* Job snippets: :doc:`/capabilities/index`
* Recipes: ``mcf.recipes('combine_layers')``, ``mcf.recipes('align')``,
  ``mcf.recipes('prep_layers')``

``combine_layers`` / ``combine_from_files`` run stretch → colorize → combine.
``align_stack`` / ``prep_layers`` live here in the API index (stack tools);
north-up / beam / crop prose is in the preparing-images guide.
``save_combined`` writes an RGB FITS cube with optional provenance HISTORY.
``downsample_for_preview`` is for interactive previews only.

.. currentmodule:: multicolorfits

.. autosummary::
   :toctree: generated
   :nosignatures:

   combine_layers
   combine_from_files
   save_combined
   align_stack
   prep_layers
   reproject_stack_to_header
   reproject_stack_to_reference
   downsample_for_preview
