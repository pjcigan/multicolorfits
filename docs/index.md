---
sd_hide_title: true
---

# multicolorfits

```{toctree}
:hidden:

installation
guide/index
tutorials/index
examples/index
api/index
changelog
migration
```

# multicolorfits

**Colorize and combine FITS images for visually aesthetic scientific plots.**

```{image} _static/logo/logo_hero_crab.png
:alt: Crab Nebula multicolor composite (red–orange palette)
:class: mcf-hero-crab
:width: 420px
:align: center
```

Stretch and tint multiwavelength layers, mix them with classic RGB or Lab
compositing, and finish publication figures — from scripts or the browser /
desktop GUI.

```{tip}
**Human or AI agent?** Start with ``import multicolorfits as mcf;
mcf.overview()`` for the mental model and conventions, then
``mcf.recipes('<keyword>')`` for copy-paste code.  The same catalog is
served at the site root as <a href="llms.txt"><code>llms.txt</code></a> /
<a href="llms-full.txt"><code>llms-full.txt</code></a> (generated from
``multicolorfits/_overview.py`` — see {doc}`api/getting_oriented`).
```

```{code-block} python
import multicolorfits as mcf

s = mcf.McfSession(n_panels=3)
s.load_files(
    ['ir.fits', 'r.fits', 'b.fits'],
    colors=['#BE599E', '#DEA215', '#77C0F9'],
    labels=['IR', 'R', 'B'],
)
s.compose.combine_mode = 'lab'
s.compose.combine_blend = 'screen'
rgb = s.render_combined()
```

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {octicon}`download` Installation
:link: installation
:link-type: doc

`pip install "multicolorfits[all]"` for both GUIs and helpers, or pick extras
à la carte.
:::

:::{grid-item-card} {octicon}`book` User guide
:link: guide/index
:link-type: doc

Compositing modes, intensity scaling, the session model, saving outputs, and
component mosaics.
:::

:::{grid-item-card} {octicon}`mortar-board` Tutorials
:link: tutorials/index
:link-type: doc

Getting started, GUI walkthroughs, save/export, notebook embed, and short
stubs that pair with the Examples narratives.
:::

:::{grid-item-card} {octicon}`image` Examples
:link: examples/index
:link-type: doc

NGC 602, WLM, M74 transparent cutouts, color-space suite, and a gallery.
:::

:::{grid-item-card} {octicon}`browser` Browser GUI
:link: tutorials/gui_web_walkthrough
:link-type: doc

`mcf.gui()` / `mcf-web` — panels, Lab compositing, live preview, session JSON,
export script.
:::

:::{grid-item-card} {octicon}`git-compare` Migrating from v2
:link: migration
:link-type: doc

Traits-free GUIs, optional extras, and the shared `McfSession` controller.
:::

:::{grid-item-card} {octicon}`code-square` API reference
:link: api/index
:link-type: doc

Public functions and classes, grouped by topic, with autosummary signatures.
:::

:::{grid-item-card} {octicon}`mark-github` Package overview
:link: https://github.com/pjcigan/multicolorfits
:link-type: url

README, gallery, citation, and source on GitHub.
:::

::::

## Quick install

```bash
pip install "multicolorfits[all]"
mcf-web    # or: python -c "import multicolorfits as mcf; mcf.gui()"
```
