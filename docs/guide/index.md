# User guide

Conceptual documentation for multicolorfits v3 — how the color pipeline fits
together before you dive into step-by-step tutorials.

```{toctree}
:maxdepth: 1

color_compositing
intensity_scaling
session_model
wcs_and_alignment
saving_outputs
figures_and_mosaics
```

| Page | Topics |
|------|--------|
| {doc}`color_compositing` | Compositing modes, `mode` vs `colorspace`, backgrounds, RYB/CMYK |
| {doc}`intensity_scaling` | Stretches, zscale, signed flux (`symlog` / `pysymlog`) |
| {doc}`session_model` | `McfSession`, panels, legend/swatch, preview, palettes |
| {doc}`wcs_and_alignment` | `McfHeader`, crop, `align_stack` / Align layers… |
| {doc}`saving_outputs` | Rasters, FITS, cutouts, export script |
| {doc}`figures_and_mosaics` | Combined WCS figures, mosaics, overlays |
