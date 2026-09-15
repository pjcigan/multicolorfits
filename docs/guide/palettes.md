# Choosing and generating layer colors

Layer **tints** (which hex each panel uses) are separate from **compositing
mode** (how those tinted layers are mixed). This page is about the first
problem — picking colors — and how the scripting API matches the GUI palette
menu.

Job snippets: {doc}`../capabilities/index`. API: {doc}`../api/palettes`.
Compose modes (RGB / Lab / HSV mixers): {doc}`color_compositing`.

```{contents}
:local:
```

## Two different “Lab / HSV” ideas

| Concept | Where | What it controls |
|---------|-------|------------------|
| **Layer colors** | `suggest_colors`, `get_palette`, … | The hex tint of each band before combine |
| **Compose mode** | `s.compose.combine_mode = 'lab'` / `'hsv'` | How overlapping tints are *mixed* |

CIE Lab / HSV as a **compose mode** does **not** generate palettes. Palette
generators below build hex lists; you then set `panels[i].color` or
`s.apply_palette(...)`, and optionally mix with Lab/HSV/RGB.

The recommended auto generator is **CIE LCh** (`suggest_colors`), not HSV.
HSV wheel steps look uneven (yellows dominate). Use `colors_from_hsv` only
when you explicitly want a classical HSV ramp.

## What the GUI palette menu maps to

| GUI entry | Scripting |
|-----------|-----------|
| Even (auto-N) / Perceptual | `mcf.suggest_colors(n)` |
| Complement / Triad / Split / Square / Analogous | `mcf.colors_for_hue_pattern('triad', n)` (etc.) |
| Pastel ring / Vivid ring | `mcf.get_palette('pastel')` / `'vivid'` (fixed LCh rings) |
| Wong / Tol / IBM / POB / … | `mcf.get_palette('cb_safe')`, `'tol'`, `'pob'`, … |
| Apply to layers | `s.apply_palette(name)` or assign hex list |
| Hue-wheel tuner | `mcf.rotate_colors_from_base(colors, delta_deg)` / `s.rotate_panel_hues(...)` |

```python
import multicolorfits as mcf

print(mcf.list_palettes())           # curated static names
print(mcf.list_hue_patterns())       # complement, triad, …
print(mcf.list_palette_menu())       # same groups as the GUI
```

## Curated palettes

```python
colors = mcf.get_palette('pob')              # purple / orange / blue
colors = mcf.get_palette('cb_safe', n=4)     # first 4 of Wong set
s.apply_palette('tol')                       # session: recolor loaded panels
```

## Perceptual auto (CIE LCh) — default for “N layers”

Equal hue steps on a constant lightness / chroma ring (GUI **Even (auto-N)**):

```python
colors = mcf.suggest_colors(4)                       # defaults L*=65, C=55
colors = mcf.suggest_colors(5, lightness=78, chroma=38, hue_start=20)  # pastel-ish
colors = mcf.suggest_colors(3, lightness=55, chroma=72, hue_start=10)  # vivid-ish
```

Custom angles on the same LCh ring:

```python
colors = mcf.colors_from_hue_angles([0, 120, 240], rotation=30,
                                    lightness=65, chroma=55)
```

## Hue-geometry patterns (Paletton-style)

```python
colors = mcf.colors_for_hue_pattern('triad', n=3, rotation=15)
colors = mcf.colors_for_hue_pattern('complement', n=2)
colors = mcf.colors_for_hue_pattern('analogous', n=4, lightness=60, chroma=50)
# Alias of suggest_colors:
colors = mcf.colors_for_hue_pattern('perceptual', n=6, rotation=25)
```

When `n` is larger than the pattern’s anchors, extra hues are filled in
around the wheel.

## Classical HSV wheel (optional)

```python
colors = mcf.colors_from_hsv(4, saturation=0.85, value=0.9, hue_start=0)
```

Prefer LCh for science plots; keep HSV for teaching demos or matching an
external HSV tool.

## Preview, CVD check, and `plt.show()`

`preview_palette` builds a **pyplot** figure (tiles + combo swatch + CVD
rows). The combo swatch’s `mode` / `blend` match **compositing**, not the
generator:

```python
import matplotlib.pyplot as plt
import multicolorfits as mcf

colors = mcf.suggest_colors(3)
res = mcf.preview_palette(
    colors, labels=['IR', 'R', 'B'],
    mode='lab', blend='screen',          # how overlaps mix in the swatch
    background='black')

plt.show()          # works — figure is pyplot-managed
# res.show()        # same
# In a notebook cell, `res.fig` as the last expression also embeds it.

res.fig.savefig('palette.png', dpi=120, facecolor=res.fig.get_facecolor())
print(res.cvd['ok'], res.colors)
```

Named palettes work the same way: `mcf.preview_palette('pob', n=3, mode='lab')`.
From a session: `s.preview_palette()` uses current panel colors and compose
settings.

```{note}
Earlier builds used a bare ``Figure()`` that was **not** registered with
pyplot, so ``plt.show()`` showed nothing. That is fixed — use a current
``multicolorfits`` with ``preview_palette`` from 3.1.0 onward.
```

## Apply colors to a session or pipeline

```python
s = mcf.McfSession(n_panels=3)
s.load_files(paths, colors=mcf.suggest_colors(3), labels=['A', 'B', 'C'])
# or later:
s.apply_palette('triad')                 # hue pattern or curated name
s.rotate_panel_hues(40)                  # spin current colors
print(s.colorblind_report())
```

Classic API: pass the hex list into `combine_layers(..., colors=...)` or
`colorize_image(..., color)`.

## See also

- Side-by-side looks on one target: {doc}`../examples/crab_palette_gallery`
- Session overview: {doc}`session_model`
- Recipe: `mcf.recipes('palette')`
