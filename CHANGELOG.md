# Changelog

Notable changes to multicolorfits. Version **3.0.0** is the traits-free package
rebuild; the v2.x scripting API (`greyRGBize_image`, `colorize_image`,
`combine_multicolor`, `saveRGBfits`, reprojection helpers, and related names)
remains available via `import multicolorfits as mcf`. See
[MIGRATION.md](https://github.com/pjcigan/multicolorfits/blob/master/MIGRATION.md)
for GUI and packaging changes.

---

## [3.0.0] — 2026-07-18

Traits-free package rebuild: browser + Qt GUIs over a shared session model,
modern compositing (Lab / RYB / CMYK, backgrounds, cutouts), and a Sphinx docs
site — while keeping the classic scripting API via silent aliases.

### Removed
- **traits / traitsui / pyface / PyQt** dependency stack. The computational core
  imports with the standard scientific stack (numpy, astropy, scipy, matplotlib,
  scikit-image); GUIs are optional extras.
- Monolithic root `multicolorfits.py` TraitsUI entry point (replaced by the
  `multicolorfits/` package).

### Changed — packaging & license
- Setuptools via **`pyproject.toml`** with extras: `[web]`, `[qt]`, `[reproject]`,
  `[overlays]`, `[pysymlog]`, `[all]`, `[test]`.
- Relicensed from MIT to **BSD 3-Clause** (aligned with skyplothelper). Past
  releases remain MIT; 3.0.0 and later are BSD-3-Clause.
- Version **3.0.0**.

### Changed — public API (renames with silent aliases)
- Most public functions gained clearer PEP 8 names; **every v2.x name still
  works** as a fully functional alias (no warnings yet). Highlights:
  `to_grey_rgb` (`greyRGBize_image`), `hex_complement` (`hexinv`),
  `save_rgb_fits` (`saveRGBfits`), `plot_combined_rgb`
  (`plotsinglemulticolorRGB`), `compare_multicolor_vs_rgb`
  (`comparemulticolorRGB_pureRGB`), `stretch_functions` (`scaling_fns`),
  `make_simple_header` (`makesimpleheader`), `sky_to_pixel` / `pixel_to_sky`,
  `crop_image[_sky]` / `crop_cube[_sky]`, `reproject_image` / `reproject_cube`,
  `layer_coverage` (`image_value`), and colorspace converters without the
  `_vectorized` suffix. Full table in
  [MIGRATION.md](https://github.com/pjcigan/multicolorfits/blob/master/MIGRATION.md);
  aliases live in `multicolorfits/compat.py` (`mcf.LEGACY_ALIASES`).
  `McfSession.export_script()` emits the new names.

### Changed — architecture
- Package layout: `core/` (scaling, colorize, colormix, paintmix, swatch, WCS,
  reproject, I/O), `session.py`, `figures.py`, `palettes.py`, `pipeline.py`,
  `overlays.py`, `gui_web/`, `gui_qt/`.
- Flat `__init__.py` re-exports the full public / compat API.
- Optional deps via lazy imports (`reproject`, `kapteyn`, `skyplothelper`).

### Added — browser & Qt GUIs
- **Browser GUI** (primary): FastAPI + vanilla JS (`mcf-web` / `mcf.gui()`).
- **PySide6 desktop GUI** (secondary): `mcf-qt` / `mcf.gui_qt()`.
- Shared **`McfSession`** controller; both GUIs stay in feature parity for the
  core workflow.
- **Dynamic panels:** start with 4 slots; **+ Add panel** / **Remove** (max 16).
  API: `add_panel()`, `remove_panel()`, `set_n_panels()`.
- **Notebook embed:** `mcf.gui_embed()`, `mcf.start_gui_server()`,
  `mcf.detect_notebook_environment()`.
- **Session JSON** Save/Load in both GUIs; web `POST /api/session/save|load`,
  `GET /api/session.json`.
- **Preload:** `mcf.gui(state=, files=, colors=, labels=)` and matching
  `gui_qt` / `mcf-web --state` / `--files`.
- Compose UI grouped into **Compositing / Colors / Axes / Canvas** tabs.
- Light / dark **theme toggle** (persisted; `?theme=light|dark` on web).
- **Fast preview** (default on): downsampled float32 combined updates while
  editing; **Plot Full Resolution** for the publication WCS figure (axes,
  overlays). Web also has a client-side blit path for classic RGB / Lab-on-black.
- Cursor **pixel readout** on the combined view (per-layer flux + sky coords).
- Per-panel labels, channel **legend**, **color-combo swatch**, **band labels**.
- **Palette** picker + live color-vision (CVD) status.
- **Align layers…** when grids differ (needs `[reproject]`).
- **Export Script** replays alignment and emits a faithful reproduce script.
- **Save FITS…** targets: annotated combined RGB cube and per-layer data/color
  components with provenance keywords (`MCFTYPE`, `MCFMODE`, `MCFBLEND`, …).
- **Export transparent…** cutout dialog (web + Qt).

### Added — compositing & color
- **CIE Lab** compositing (`combine_multicolor_colorspace`, blend `screen` /
  `sum` / `max` / `mean`) — stable, documented, in both GUIs. HSV/HSL remain
  available but marked experimental (hue drift under unbalanced palettes).
- **RYB (paint)** subtractive mixing and **CMYK** print-ink mode
  (`combine_multicolor_alpha`, `rgb_to_ryb` / `ryb_to_rgb`, …).
- **Background** control: black / white / transparent / custom color via
  `composite_over_background` (supersedes naive `inverse=True` for publication
  canvases; Inverse checkbox remains for the classic complement flip).
- Gamma applies to the combined image in RGB and in RYB/Lab/CMYK tone mixing;
  left-panel thumbnails intentionally ignore gamma.
- Accessibility helpers: `greyscale_image()`, `simulate_colorblindness()`.
- `mix_colors_hex()` for palette design; curated palettes + `suggest_colors(n)`.

### Added — cutouts, mosaics & figures
- Transparent cutouts: `make_transparent_cutout`, `save_transparent_cutout`,
  `preview_cutout_on_backgrounds`, `batch_transparent_cutouts`,
  `McfSession.export_transparent_cutout` — with `alpha_smooth`,
  `alpha_source='sat'|'dist'`, soft edges, crop/size.
- `deblend_background` / `flatten_rgba` to recover RGBA from a composite
  flattened onto a known solid color.
- `make_component_mosaic` / `session.plot_component_mosaic` — hero + colorized
  strip.
- Shared WCS figure helpers (`setup_combined_axes`, canvas facecolor,
  overlays via optional `[overlays]` / skyplothelper).
- High-level `pipeline.combine_layers` / `combine_from_files` / `save_combined`.

### Added — WCS, headers & alignment
- `McfHeader` — stateful header wrapper with cached WCS / frame / pixel scale
  and `.matches_grid()`.
- Celestial frame helpers and stack align: `align_stack`,
  `reproject_to_frame`, `reproject_to_galactic`, …
- Portable header sidecars: `save_header` / `load_header`,
  `McfHeader.minimal_wcs_header()`.
- Signed-data stretches: `symlog` and optional `symmetric_log` (`[pysymlog]`).
- WCS/coord fixes ported from skyplothelper (sexagesimal rollover, PC/CD
  cleanup, `getcdelts` / `getdegperpix`, Stokes-axis squeeze on load).

### Added — performance
- Faster `colorize_image` (direct-RGB gamma path) and `combine_multicolor`
  (pure-numpy, no large nansum stack) — roughly **~2×** end-to-end on large
  three-layer composites.
- Float32 + downsample **preview** path (`preview=True`, `max_size=`) for
  interactive work; full float64 remains the export default.

### Added — docs, orientation & CI
- Sphinx site (`docs/`) with guides, autosummary API, tutorials, Examples
  (NGC 602, WLM, M74 cutouts, color-space suite, gallery), light/dark figure
  pairs, and plot-theme toggle.
- `mcf.overview()` / `mcf.recipes()` capability map; committed `llms.txt` /
  `llms-full.txt` (CI fails if they drift).
- Teaching `TypeError`s on common mistakes (`colorize_image` on raw 2D;
  `combine_multicolor` on a single array).
- GitHub Actions test workflow on push/PR.
- README install matrix, migration notes, AI-agent section.

### Fixed
- `drawProgressBar` missing `sys` import.
- NumPy 2.0: `np.NZERO` → `-0.0`.
- `combine_multicolor` divide-by-zero when all channels saturate in inverse mode.
- scikit-image `multichannel` → `channel_axis` in `smooth_image`.
- Web combined preview now honors session `vmin`/`vmax` on first load; panel
  canvases match FITS origin (vertical flip).
- Compose gamma no longer cancelled in RYB/Lab/CMYK display-referred mixes.

---

# Earlier releases (v1–v2)

# 2.1

2025-06-02

* Removed orphaned matplotlib imports no longer used, some of which were preventing mcf to successfully load (AnchoredEllipse is now deprecated)
* Made some improvements to speed in colorize_image: treating hex color colorization separately, and also vectorizing hsv conversion
* Made improvements to behind-the-scenes testing of sphinx documentation building

2023-09-27

* Updated PyQt imports to use PyQt6 by default, to work with py 3.10 and M processor Macs.  PyQt5 is still there as a backup.  PyQt4 should still work for Python 3.4 and lower at the moment, but this may be phased out in future releases.

2022-09-04

* Added buttons in the left side panels to edit image headers
    - This may be useful for tweaking things (e.g. bump reference pixel values to shift an HST image alignment), to declutter HISTORY and COMMENT cards, etc.  You can also try deleting offending header cards if something is preventing proper image combination.
    - Currently, after making edits in the popup window the user must also click the Apply Header Changes button to apply the changes.  Clicking Edit Header again before applying will discard the changes. 
    - Headers can also be saved to / loaded from text files
* Fixed the Save Image button functionality 
* Added options for setting xaxislabel and yaxislabel in the convenience functions plotsinglemulticolorRGB() and comparemulticolorRGB_pureRGB()
* Added fields for x/y axis ticklabel formats.  (See [the astropy documentation](https://docs.astropy.org/en/stable/visualization/wcsaxes/ticks_labels_grid.html#tick-label-format for options) for available formats and syntax)
* Combined image x/y axis labels now automatically generate from header CTYPE1/2 values (when present)
* Added a button for turning minor ticks on/off
* Added simple gaussian smoothing function.  This can be useful for reducing noise, matching resolutions, etc.


# 2.0

2019-09-12

* Updated to use in python 3 or python 2. (Tested in 3.6.8 and 2.7.6)
* Updated to PyQt5 
* Preferred usage is now either 
    1. run multicolorfits.py from terminal, or 
    2. import and run GUI with multicolorfits.mcf_gui()  
* Added several functions to aid scripting (outside of the GUI)
    - reprojection to specified image header (requires reproject package or kapteyn package)
    - cropping by specifying coordinates and desired size
    - color conversions
    - some header wcs handling
    - Quick-look plotting functions to save out final combined images, and to compare to the simple RGB case
* Added an option to disable auto-refresh/update of the preview plots in the left panel after adjusting min/max values.  This will increase speed for large files, but with this option users will need to click the 'Plot' button (or the new 'Update plot' button) after each change to see the results.
* improved handling of image NaNs


# 1.2 

2018-12-12

* Roughly fixed incompatibility with matplotlib 1.5. 
    - defined a replacement to_hex() function
    - separated out rcParams keywords that don't exist in mpl 1.5


# 1.1 

2018-12-11

* Implemented the 'inverse plot' buttons


# 1.0

2018-12-10

* Initial version on github
