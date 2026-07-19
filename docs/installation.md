# Installation

```bash
pip install "multicolorfits[all]"
```

That pulls in both GUIs, alignment (`reproject`), overlays, and related
helpers.  For a smaller install:

| Extra | Provides |
|-------|----------|
| *(none)* | Core library (numpy, astropy, scipy, matplotlib, scikit-image) |
| `[web]` | Browser GUI |
| `[qt]` | PySide6 desktop GUI |
| `[reproject]` | Layer alignment / reprojection |
| `[overlays]` | Compass, beam, scale bar (skyplothelper) |
| `[pysymlog]` | Optional signed-stretch helper |
| `[all]` | Everything above |
| `[test]` | pytest stack for contributors |

```bash
pip install "multicolorfits[web,reproject]"   # example à la carte
```

## Conda environments and the `[qt]` / PySide6 extra

The browser GUI (`[web]`) needs no Qt and is safe in shared conda stacks.
The desktop GUI (`[qt]`) installs **PySide6** from PyPI — a full Qt6 runtime.
That is fine in a **dedicated** environment, but it is a common way to break
interactive matplotlib in a conda env that already has conda-forge's matched
**Qt5** stack (`pyqt` / `qt-main`).

**What goes wrong.**  A pip Qt6 wheel (`PySide6`, or leftover `PyQt6` /
`PyQt6-Qt6` / `PyQt6-sip` from older testing) ships its own `libQt6*.so` and
related libraries into `site-packages`.  Even when matplotlib's `QtAgg`
backend loads conda **PyQt5**, having a second Qt major with colliding
SONAMEs in-process can corrupt the heap once a real GUI event loop starts
(`plt.show()` → `malloc(): invalid size` / abort).  Headless **Agg**
rendering never trips it, so test suites often stay green.

**What to do.**

* Prefer `[web]` (or core-only mcf) inside shared conda Qt5 environments.
* Install `[qt]` / PySide6 only in a **throwaway or dedicated** conda/venv,
  not in a shared env others use for interactive plotting.
* Do not pip-install `PyQt6` into a conda Qt5 env (same class of poison).

**If interactive matplotlib dies** (repro: bare `plt.plot([0, 1]); plt.show()`
with no mcf involved):

```bash
conda list | grep -Ei 'pyqt|qt-main|pyside'
pip list | grep -iE 'pyqt|pyside|shiboken'
```

Look for a **pip (PyPI)** Qt6 row beside conda Qt5.  Remove the pip side,
then verify the modules are gone (a first `pip uninstall` sometimes leaves
files behind — run again if needed):

```bash
pip uninstall -y PyQt6 PyQt6-Qt6 PyQt6-sip
# if PySide6 was the pip install that mixed in:
# pip uninstall -y PySide6 shiboken6
python -c "import importlib.util as u; print(u.find_spec('PyQt6'), u.find_spec('PySide6'))"
```

Nuclear repair for the conda Qt5 stack if needed:

```bash
conda install -c conda-forge --force-reinstall pyqt qt-main qt-webengine
```

Docs build extras for this Sphinx tree:

```bash
pip install -r docs/requirements.txt
```

See the [README](https://github.com/pjcigan/multicolorfits) for citations and
gallery images, and [MIGRATION.md](https://github.com/pjcigan/multicolorfits/blob/master/MIGRATION.md)
for v2 → v3 notes.
