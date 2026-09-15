# Preparing images

Before colorizing, layers often need a usable header, a common orientation,
and a shared pixel grid. This page is the processing path for that: tidy the
header, put north up in the **image's own frame**, optionally oversample, align,
crop to the overlap, and (when beams differ) match the resolution.

Job map: {doc}`../capabilities/index`. Alignment concepts:
{doc}`wcs_and_alignment`. Recipes: `mcf.recipes('prep_layers')`,
`mcf.recipes('match_beam')`, `mcf.recipes('north')`.

```{contents}
:local:
```

## Tidy the header

``force_header_2d`` only drops extra axes. ``make_simple_header`` copies WCS
cards into a minimal header. For everyday scripts, start with
{func}`multicolorfits.tidy_header` — a copy that also drops conflicting
``PC`` + ``CD`` cards and orphan ``CROTA`` when a matrix is present, and warns
if the WCS is mirrored (a parity flip, not a pure rotation).

{func}`multicolorfits.describe_header` returns the frame, center, pixel scale,
flip flag, and beam (and prints them by default; pass ``verbose=False`` to
keep the dict quiet). {func}`multicolorfits.describe_image` adds a suggested
stretch and vmin/vmax from the pixel distribution. For a whole set of layers
(plus suggested colors), use {func}`multicolorfits.describe_images` — see
{doc}`intensity_scaling`. {func}`multicolorfits.read_fits` applies
``tidy_header`` by default.

```python
import multicolorfits as mcf

data, hdr = mcf.read_fits('image.fits')
info = mcf.describe_header(hdr, name='image.fits')          # printed
# info = mcf.describe_header(hdr, name='image.fits', verbose=False)
if mcf.wcs_is_flipped(hdr):
    print('mirrored WCS — north-up will force east-left unless allow_flip=True')
```

A positive ``CDELT1`` does not by itself mean "flipped". A normal east-left
image has a negative determinant of the pixel-scale matrix; a mirror has
the opposite sign (a negative PC determinant can cancel a positive
``CDELT1``). {func}`multicolorfits.east_increases_right` probes the
world-to-pixel mapping the same way skyplothelper does.

## North-up, in the right frame

{func}`multicolorfits.make_rotated_header` builds a new grid. The default
``rotation_deg=0`` is **north-up in the header's celestial frame**: a Galactic
image stays Galactic (``+b`` up), it is not forced onto ICRS. Pass
``frame='galactic'`` only when you also want a frame change.

{func}`multicolorfits.reproject_north_up` is the one-shot resample (needs the
``reproject`` extra). Arbitrary angles use
{func}`multicolorfits.reproject_to_rotation`.

```python
north, hdr_n = mcf.reproject_north_up(data, hdr, oversample=2, order=1)
# or a custom angle from north-up:
turned, hdr_t = mcf.reproject_to_rotation(data, hdr, rotation_deg=15)
```

``oversample=2`` halves ``|CDELT|`` so a derotated, chunky image is less
jagged. It does not recover structure smaller than the original pixel.
``order`` is the ``reproject_interp`` spline order (1 bilinear, 3 cubic).
``method='exact'`` or ``order=0`` uses ``reproject_exact``.

The target header is rebuilt (east-left, clean ``CDELT``). Do not just zero
``CROTA`` or ``PC`` in place — a ``PC`` matrix can include a reflection, and
``CDELT1`` is not always the sky scale.

## Align, then crop to the overlap

{func}`multicolorfits.prep_layers` tidies, optionally north-ups and oversamples,
reprojects every layer onto that grid, and crops to pixels that are finite in
every layer (the NaN margin left by a rotated footprint).

```python
prepared = mcf.prep_layers(
    [(data_a, hdr_a), (data_b, hdr_b)],
    north_up=True,
    oversample=2,
    crop='overlap',
)
```

{func}`multicolorfits.crop_to_overlap` is the crop step on its own.
{func}`multicolorfits.blank_missing` turns non-finite values into NaN; pass
``zeros=True`` only when exact zero really means "no data" (chip gaps), not
when zero is a real measurement. ``McfSession.blank_panels`` and the Align
layers dialog (north-up / crop-to-overlap) call the same helpers.

In the browser GUI, scroll and drag zoom the combined preview only — that
does not resample. **Crop to view** is available once layers share a grid,
and it calls {func}`multicolorfits.crop_image` on the visible pixels of the
fast preview (not the WCS axes figure).

## Matching beams

{func}`multicolorfits.match_beam` degrades an image to a coarser elliptical
Gaussian beam. The kernel is the one in J.P. Wild, *Australian Journal of
Physics* **23**, 113–115 (1970): the desired beam is the convolution of a
kernel with the current beam. Convolution cannot sharpen. For Jy/beam maps
set ``per_beam=True`` so the array is scaled by the ratio of beam areas
(this function does not reproject, so both areas use the same pixel scale).

The same math is available under the older names ``convolve2Dgaus`` and
``convolve2Dgaus_matchhdr``.

```python
# Header already has BMAJ/BMIN/BPA (degrees):
matched, hdr_m = mcf.match_beam(data, hdr, bmaj_to_asec=8, bmin_to_asec=6, bpa_to_deg=20)
# Explicit source beam when those cards are missing:
matched, hdr_m = mcf.match_beam(
    data, hdr, bmaj_to_asec=12, bmin_to_asec=12, bpa_to_deg=0,
    bmaj_from_asec=5, bmin_from_asec=5, bpa_from_deg=0)
# Or match another header's BMAJ/BMIN/BPA (Jy/beam → per_beam=True):
matched, hdr_m = mcf.match_beam_to_header(data, hdr, hdr_coarse, per_beam=True)
```

Session users usually match beams **before** `load_files` / `set_data`, then
align with `s.align_panels(...)` if the grids still differ.

## Writing the result

{func}`multicolorfits.write_fits` writes the array and a copied header
(beam and ``BUNIT`` included) and a short provenance note.
{func}`multicolorfits.save_rgb_fits` is still the helper for an RGB cube.

See the {doc}`../tutorials/prepare_and_align` tutorial for a walk-through on
a local FITS crop.
