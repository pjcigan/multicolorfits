# Tutorials

Runnable tutorials are jupytext **percent** `.py` sources paired with `.ipynb`
(see the notebook download on Read the Docs / GitHub).  Screenshot-oriented
pages are MyST markdown. Fuller narrative walkthroughs with published figures
live under {doc}`../examples/index`.

```{toctree}
:maxdepth: 1

getting_started
ngc602_same_grid
wlm_align_reproject
colorspace_compositing
backgrounds_and_ryb
gui_web_walkthrough
gui_notebook_embed
save_and_export
component_mosaic
session_export_script
```

| Tutorial | Format | Topic |
|----------|--------|-------|
| Getting started | notebook | Minimal `McfSession` / Lab |
| NGC 602 | notebook | Same-grid OpenFITS → {doc}`../examples/ngc602` |
| WLM | notebook | Crop, align, combine → {doc}`../examples/wlm` |
| Color-space compositing | markdown | Pointer → {doc}`../examples/colorspace_comparison` |
| Backgrounds & RYB | notebook | White/transparent canvas, paint mix |
| GUI walkthrough | markdown | Browser GUI layout & workflow |
| GUI in notebooks | notebook | `gui_embed` / Colab / Binder |
| Save & export | notebook | PNG/FITS/cutout/script |
| Component mosaic | notebook | Hero + strip figure |
| Session & export script | notebook | JSON round-trip, `export_script` |
