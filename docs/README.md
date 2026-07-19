# Documentation

Sphinx builds from this directory (`conf.py`). Prefer
`.readthedocs.yaml` at the repository root (legacy `readthedocs.yml`
kept for older RTD project settings).

```bash
pip install -r requirements.txt   # or: pip install "multicolorfits[docs]"
make html
# open _build/html/index.html
```

Theme: **pydata-sphinx-theme** with mcf accents (`_static/custom.css`), denim
Pygments styles, and `plot-theme.js` (navbar toggle for future paired
`.plot-light` / `.plot-dark` figures). On Read the Docs, a version
switcher reads `_static/switcher.json`.

| Tree | Contents |
|------|----------|
| `guide/` | Conceptual MyST pages (compositing, session, save, mosaics) |
| `tutorials/` | Jupytext percent `.py` ↔ `.ipynb` (nbsphinx) plus GUI markdown |
| `examples/` | Worked narratives (NGC 602, WLM, color-space suite) + gallery |
| `api/` | Autosummary reference grouped by topic (`api/generated/` is build output) |
| `_static/compositing/`, `mosaics/`, `gui/`, `examples/` | Light/dark pairs + classic result images |

Author notebooks as jupytext sources, then sync:

```bash
make tutorials-sync
```

Refresh committed plot assets in `_static/` (after regenerating PNGs into
local-only `examples/output/`):

```bash
make figures   # runs make_docs_figures.py
```

Browser GUI walkthrough screenshots (needs Chrome + local NGC 602 FITS):

```bash
python docs/capture_web_gui_screenshots.py
```

Compositing-mode montage on that page: regenerate with
`python hidden/testing/gui_compositing_demo.py`, then `make figures`.

Branding marks (navbar, favicon, Crab hero / social card) live under
`_static/logo/` and are regenerated with:

```bash
python docs/make_logo.py
```

Active navbar mark is the filled color swatch
(`logo_navbar-color.png`); mono outline variants
(`logo_navbar-mono-{light,dark}.png`) are generated too — swap them in
`_templates/navbar-logo.html` if you want the quieter look later.

`nbsphinx_execute = "never"` — tutorials ship stub/commented cells so RTD
does not need OpenFITS / LITTLE THINGS data. Plot light/dark follows the
navbar plot-color toggle (`plot-theme.js`).

