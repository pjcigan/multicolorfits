# Crab Nebula — palette & compositing gallery

Same four HST WFPC2 layers (F502N / F547M / F631N / F673N), many looks.
Changing colors, the compositing **mode** (RGB vs Lab / HSL), and the
background (black vs white) completely changes the feel of an image — useful
for artistic preference, or for figures that need to sit on a light
publication page or slide deck.

```{image} ../_static/showcase/crab_gallery.png
:class: mcf-plot
:alt: Eight Crab Nebula stamps with different palettes and compositing modes
```

| Panel | Mode | Notes |
|---|---|---|
| **red orange** | RGB / black | Classic `m1_dark1` / `big_M1_4color_dark1.jpg` |
| **purple orange** | RGB / black | Classic `m1_dark2` / `big_M1_4color_dark2.jpg` |
| **gui** | RGB / black | README / Qt splash palette |
| **ember** | Lab / max | Warm lava + violet/cyan accents (artistic) |
| **ice lava** | Lab / max | Hot filaments on an icy core (artistic) |
| **aurora** | RGB / black | Green / teal / violet / pink (artistic; classic additive is more vibrant than Lab) |
| **copper gold** | Lab / screen | Warm near-monochrome (artistic) |
| **sunset** | HSL / screen | Oranges / violet / blue (artistic) |
| **on white** | RGB / white | Same dark1 colors via the alpha **Background** path — for light slides / papers |

Each stamp is a transparent-sky RGBA cutout built on top of a
**black-background** composite (so classic RGB recipes keep the same colors as
``big_M1_4color_dark1.jpg`` / ``dark2.jpg``).  A soft sky matte
(``alpha_smooth=3``, ``alpha_lo≈45``, ``alpha_hi≈82``, ``alpha_gamma=0.7``)
clears mosaic-edge junk while keeping the inner filaments opaque.

```{note}
Earlier cutouts looked neon / blown-out because ``export_transparent_cutout``
forced a transparent compose background, which routes classic RGB through the
alpha compositor and changes the colors.  That helper now preserves
black-background math by default (``render_background='auto'``).
```


## Reproduce locally

The gallery mosaic above is shipped in ``docs/_static/showcase/``.  To explore
the same looks interactively, point a session at a local Crab HST WFPC2 crop
(same four layers as the README splash) and try the table’s colors / modes:

```python
import multicolorfits as mcf

s = mcf.McfSession(n_panels=4)
s.load_files(
    ['f502n.fits', 'f547m.fits', 'f631n.fits', 'f673n.fits'],
    colors=['#FF2200', '#FFAA00', '#FFE080', '#FFFFFF'],  # e.g. red_orange
    labels=['F502N', 'F547M', 'F631N', 'F673N'],
)
s.compose.combine_mode = 'rgb'          # or 'lab' / 'hsl' per the table
s.compose.combine_blend = 'screen'      # lab/hsl: screen / max / …
s.compose.combine_background = 'black'  # or 'white' for the slide-friendly look
mcf.gui(s)                              # or: fig = s.plot_combined()
```

For a transparent-sky stamp like the gallery cells, use
``s.export_transparent_cutout(...)`` (see {doc}`m74_transparent_cutout`).
