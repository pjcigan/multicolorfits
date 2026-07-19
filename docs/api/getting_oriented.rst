Getting oriented
================

New here — human or AI agent? Start with the runtime helpers, then dig into
the topic pages.

.. code-block:: python

   import multicolorfits as mcf
   mcf.overview()              # mental model, conventions, task index
   mcf.recipes('cutout')       # copy-paste recipe matching a keyword
   cat = mcf.overview(as_dict=True)   # structured catalog for tools

The same catalog is published as machine-readable files at the docs site
root:

* https://multicolorfits.readthedocs.io/en/latest/llms.txt
* https://multicolorfits.readthedocs.io/en/latest/llms-full.txt
  (adds function signatures)

Both files are **generated** from ``multicolorfits/_overview.py`` by
``scripts/make_llms_txt.py``.  ``tests/test_overview.py`` fails the CI if
either file drifts from the catalog — regenerate, never hand-edit.

.. autofunction:: multicolorfits.overview
.. autofunction:: multicolorfits.recipes
