# API reference

Public names are re-exported from the top-level package
(`import multicolorfits as mcf`), matching the historical flat namespace.
Grouped by topic below; each page lists functions and classes with
signatures pulled from the source.

Every v2.x name (`greyRGBize_image`, `saveRGBfits`, `reproject2D`, ...)
still works as a silent alias for its renamed v3 counterpart; the full
old-to-new table is in {doc}`/migration`, and the mapping is importable
as `mcf.LEGACY_ALIASES`.

```{toctree}
:maxdepth: 2

getting_oriented
color_pipeline
compositing
wcs_fits
session
figures
palettes
overlays
pipeline
gui
```

::::{grid} 1 2 2 3
:gutter: 2

:::{grid-item-card} Getting oriented
:link: getting_oriented
:link-type: doc
`mcf.overview()` / `mcf.recipes()` and the generated `llms.txt` corpus.
Job map: {doc}`../capabilities/index`.
:::

:::{grid-item-card} Color pipeline
:link: color_pipeline
:link-type: doc
`to_grey_rgb`, `colorize_image`, `combine_multicolor`, `describe_images`.
:::

:::{grid-item-card} Compositing
:link: compositing
:link-type: doc
Lab / HSV / HSL, RYB / CMYK, backgrounds, `combine_colorized_layers`, swatches.
:::

:::{grid-item-card} WCS & FITS
:link: wcs_fits
:link-type: doc
`McfHeader`, squeeze/crop, coordinates, reprojection helpers.
:::

:::{grid-item-card} Session
:link: session
:link-type: doc
`McfSession`, `PanelState`, `ComposeState`, save/load, export script.
:::

:::{grid-item-card} Figures
:link: figures
:link-type: doc
Combined WCS figures, component mosaics, bare/publication styling.
:::

:::{grid-item-card} Palettes
:link: palettes
:link-type: doc
Curated palettes, LCh / HSV generators, `preview_palette`, color-vision checks.
:::

:::{grid-item-card} Overlays
:link: overlays
:link-type: doc
Optional compass / beam / scale bar. Guide: {doc}`../guide/overlays`.
:::

:::{grid-item-card} High-level pipeline
:link: pipeline
:link-type: doc
`combine_layers`, `combine_from_files`, `save_combined`, stack align.
:::

:::{grid-item-card} GUIs
:link: gui
:link-type: doc
`gui`, `gui_qt`, `gui_embed`, notebook helpers.
:::

::::
