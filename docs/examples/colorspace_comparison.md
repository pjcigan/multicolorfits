# Color-space compositing suite

Side-by-side comparisons of classic RGB vs Lab / HSV / HSL on Kepler SNR and
NGC 602. Prefer **`combine_mode='lab'`** with **`blend='screen'`** as the
general-purpose alternative to classic RGB; treat HSV / HSL as experimental.

Guide: {doc}`../guide/color_compositing`. Short tutorial twin:
{doc}`../tutorials/colorspace_compositing`. Generator:
[`examples/compare_colorspaces.py`](https://github.com/pjcigan/multicolorfits/blob/master/examples/compare_colorspaces.py).

## Which mode?

| Goal | Recommendation |
|------|----------------|
| Match classic publications | **RGB** |
| Rich overlaps, avoid white clipping | **Lab / screen** (GUI default) |
| Strong per-layer hue separation | **Lab / max** |
| Paint-like mixing | **RYB** + white/custom background |
| Hue-vector experiments | HSV / HSL (watch palette balance) |

```python
import multicolorfits as mcf

combined = mcf.combine_colorized_layers(
    colorized_layers,
    mode='lab',
    blend='screen',
    gamma=2.2,
    background='black',
)
```

Use the navbar **plot-color** control to flip light / dark figure chrome.

---

## Kepler SNR — color spaces

Infrared / soft X-ray / hard X-ray — purple / orange / teal (classic linear
stretches; optical WFPC2 omitted — it only covers a tiny corner of the field).

```{image} ../_static/compositing/kepler_colorspaces_light.png
:class: mcf-plot plot-light
:alt: Kepler SNR — RGB vs Lab / HSV / HSL
```

```{image} ../_static/compositing/kepler_colorspaces_dark.png
:class: mcf-plot plot-dark
:alt: Kepler SNR — RGB vs Lab / HSV / HSL (dark chrome)
```

### Lab lightness blends

```{image} ../_static/compositing/kepler_lab_blends_light.png
:class: mcf-plot plot-light
:alt: Kepler SNR — Lab blends screen/sum/max/mean
```

```{image} ../_static/compositing/kepler_lab_blends_dark.png
:class: mcf-plot plot-dark
:alt: Kepler SNR — Lab blends (dark chrome)
```

### HSV lightness blends

```{image} ../_static/compositing/kepler_hsv_blends_light.png
:class: mcf-plot plot-light
:alt: Kepler SNR — HSV blends
```

```{image} ../_static/compositing/kepler_hsv_blends_dark.png
:class: mcf-plot plot-dark
:alt: Kepler SNR — HSV blends (dark chrome)
```

### Complementary stress test

Equal-brightness complementary hues cancel toward grey in additive models —
physics, not a bug. Palettes help when both bands must stay vivid. Soft + hard
X-ray (same stretch as above) so overlaps fill the remnant, not the tiny
optical corner.

```{image} ../_static/compositing/kepler_complementary_light.png
:class: mcf-plot plot-light
:alt: Complementary orange/blue X-ray stress test
```

```{image} ../_static/compositing/kepler_complementary_dark.png
:class: mcf-plot plot-dark
:alt: Complementary orange/blue X-ray stress test (dark chrome)
```

---

## NGC 602 — three bands

Same POB palette as {doc}`ngc602`.

```{image} ../_static/compositing/ngc602_colorspaces_light.png
:class: mcf-plot plot-light
:alt: NGC 602 — RGB vs Lab / HSV / HSL
```

```{image} ../_static/compositing/ngc602_colorspaces_dark.png
:class: mcf-plot plot-dark
:alt: NGC 602 — RGB vs Lab / HSV / HSL (dark chrome)
```

```{image} ../_static/compositing/ngc602_lab_blends_light.png
:class: mcf-plot plot-light
:alt: NGC 602 — Lab blend variants
```

```{image} ../_static/compositing/ngc602_lab_blends_dark.png
:class: mcf-plot plot-dark
:alt: NGC 602 — Lab blend variants (dark chrome)
```

---

## Subtractive paint (RYB)

Same NGC 602 linear stretch as the grids above; only the layer colors
(R / yellow / B) and compositing mode change — so the RYB panel is comparable
to classic additive, not a zscale retune.

```{image} ../_static/compositing/paintmix_demo_light.png
:class: mcf-plot plot-light
:alt: RGB vs RYB paint mixing
```

```{image} ../_static/compositing/paintmix_demo_dark.png
:class: mcf-plot plot-dark
:alt: RGB vs RYB paint mixing (dark chrome)
```

---

## Interpreting the results

- **`lab / screen`** — best general-purpose alternative to classic RGB;
  overlaps keep hue longer instead of clipping to white.
- **`lab / max`** — preserves per-layer hues most faithfully.
- **`lab / sum`** — closest to classic RGB inside Lab; rarely best.
- **`mean`** — dims the composite; better for `mix_colors_hex` than science images.
- **`hsv` / `hsl`** — can shift global cast when the palette is unbalanced
  around the hue wheel; keep experimental.

## Reproduce

```bash
python examples/compare_colorspaces.py   # light+dark pairs → examples/output/
python docs/make_docs_figures.py         # copy into docs/_static/compositing/
```

Helper pattern used by the generator (abridged):

```python
def composite(colorized, mode, blend='screen', gamma=2.2):
    if mode == 'rgb':
        return mcf.combine_multicolor(colorized, gamma=gamma)
    return mcf.combine_multicolor_colorspace(
        colorized, colorspace=mode, blend=blend, gamma=gamma)
```
