# Color-space compositing

Full figure suite (Kepler SNR, NGC 602, Lab blends, complementary stress
test, RYB demo) now lives under **Examples**:

→ {doc}`../examples/colorspace_comparison`

Guide: {doc}`../guide/color_compositing`.

```bash
python examples/compare_colorspaces.py   # light+dark theme pairs → examples/output/
python docs/make_docs_figures.py         # copies into docs/_static/compositing/
```

## Cheat sheet

| Goal | Recommendation |
|------|----------------|
| Match classic publications | RGB |
| Rich overlaps, avoid white clipping | Lab / screen |
| Strong per-layer hue separation | Lab / max |
| Paint-like mixing | RYB + white/custom background |
| Slides with transparency | `combine_background='transparent'` + PNG |

Quick visual (RYB vs additive): see the paint-mix figure on
{doc}`../guide/color_compositing`.
