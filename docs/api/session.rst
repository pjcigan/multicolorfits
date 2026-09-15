Session
=======

GUI-agnostic editing state shared by the browser and Qt apps.

* Guide (method → task map): :doc:`/guide/session_model`
* Job snippets: :doc:`/capabilities/index`
* Recipes: ``mcf.recipes('session')``, ``mcf.recipes('export')``,
  ``mcf.recipes('align_panels')``

High-traffic ``McfSession`` methods (see class docstring / guide for full
lists): ``load_files``, ``set_display``, ``apply_suggested_levels``,
``align_panels``, ``render_combined``, ``export_script``, ``save_state`` /
``load_state``, ``export_transparent_cutout``, ``preview_palette``.

.. currentmodule:: multicolorfits

Core classes
------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   McfSession
   PanelState
   ComposeState

Panel count
-----------

.. autosummary::
   :toctree: generated
   :nosignatures:

   DEFAULT_N_PANELS
   MIN_PANELS
   MAX_PANELS
   ALIGN_FRAMES
   reproject_available
