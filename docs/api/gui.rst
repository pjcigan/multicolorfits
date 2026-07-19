GUIs
====

Browser and desktop front ends over :class:`~multicolorfits.McfSession`.
Tutorial: :doc:`/tutorials/gui_web_walkthrough`,
:doc:`/tutorials/gui_notebook_embed`.

The browser GUI opens your system default browser unless you pass
``browser='firefox'`` (etc.) to :func:`~multicolorfits.gui` or
``mcf-web --browser firefox``.  See the walkthrough for details.

.. currentmodule:: multicolorfits

.. autosummary::
   :toctree: generated
   :nosignatures:

   gui
   gui_qt
   mcf_gui
   gui_embed
   start_gui_server
   detect_notebook_environment
