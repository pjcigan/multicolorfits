# User guide

Conceptual documentation for multicolorfits v3 — how the color pipeline fits
together before you dive into step-by-step tutorials.

For a job-oriented index with short snippets (mosaic, overlays, GUI, cutout,
…), start at {doc}`../capabilities/index`.

```{toctree}
:maxdepth: 1

color_compositing
intensity_scaling
session_model
palettes
wcs_and_alignment
preparing_images
saving_outputs
figures_and_mosaics
overlays
```

| Page | Topics |
|------|--------|
| {doc}`color_compositing` | Compositing modes, `mode` vs `colorspace`, backgrounds, RYB/CMYK |
| {doc}`intensity_scaling` | Stretches, `suggest_levels`, `describe_images`, Auto levels |
| {doc}`session_model` | `McfSession`, panels, palettes / `preview_palette`, legend/swatch |
| {doc}`palettes` | Generating layer colors (LCh / hue patterns / HSV), GUI map, `plt.show` |
| {doc}`wcs_and_alignment` | `McfHeader`, crop, `align_stack` / Align layers… |
| {doc}`preparing_images` | `tidy_header`, north-up, oversample, overlap crop, beam match |
| {doc}`saving_outputs` | Rasters, FITS, cutouts, export script |
| {doc}`figures_and_mosaics` | Combined WCS figures, mosaics, which axes is the hero |
| {doc}`overlays` | Scale bar, compass, beam, swatch, labels — before or after the figure |
