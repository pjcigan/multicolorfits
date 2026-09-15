# Examples

Worked narratives on public datasets, plus a small image gallery. Prefer
{doc}`../tutorials/index` when you want a short notebook-oriented path;
these pages keep the fuller scripting story (session-first, then classic
API). Job-oriented snippets: {doc}`../capabilities/index`.

```{toctree}
:maxdepth: 1

ngc602
wlm
m74_transparent_cutout
crab_palette_gallery
colorspace_comparison
gallery
```

::::{grid} 1 2 2 2
:gutter: 2

:::{grid-item-card} NGC 602 — three bands, one decision
:link: ngc602
:link-type: doc

Custom palette, RGB vs Lab side by side, component mosaic, and a ship-it
footer. Chandra OpenFITS IR / R / B.
:::

:::{grid-item-card} WLM
:link: wlm
:link-type: doc

LITTLE THINGS HI / V / NUV — crop, reproject / `align_stack`, combine, save.
:::

:::{grid-item-card} M74 — transparent cutouts
:link: m74_transparent_cutout
:link-type: doc

RGB composite → transparent-sky RGBA stamp for slides and posters; fade
shaping and background previews.
:::

:::{grid-item-card} Crab — palette gallery
:link: crab_palette_gallery
:link-type: doc

Same four HST layers, eight looks (RGB / Lab / HSL, black vs white) —
how mode and background change the feel of an image.
:::

:::{grid-item-card} Color-space suite
:link: colorspace_comparison
:link-type: doc

Kepler SNR + NGC 602 grids for RGB vs Lab / HSV / HSL (light/dark figures).
:::

:::{grid-item-card} Gallery
:link: gallery
:link-type: doc

Classic result images (NGC 602, WLM, M51, Kepler, GUI).
:::

::::

Generator scripts and offline markdown pointers still live under the
repository [`examples/`](https://github.com/pjcigan/multicolorfits/tree/master/examples)
directory (`compare_colorspaces.py`, etc.).
