/* Independent light/dark control for gallery plot images.
 *
 * Plot figures are shipped in both modes, tagged with .plot-light /
 * .plot-dark classes, and shown/hidden by CSS keyed on the
 * html[data-plot-theme] attribute (see custom.css). By default plots
 * follow the site theme; a navbar button (inserted next to the site
 * light/dark switch) cycles auto -> light -> dark so users can, e.g.,
 * preview publication-style light figures while reading in dark mode.
 * The preference persists in localStorage.
 */
(function () {
  "use strict";

  var KEY = "sphPlotTheme";          // 'auto' | 'light' | 'dark'
  var ICONS = { auto: "A", light: "L", dark: "D" };
  var TITLES = {
    auto: "Plot colors: follow site theme (click to force light)",
    light: "Plot colors: light (click to force dark)",
    dark: "Plot colors: dark (click to follow site theme)",
  };

  function pref() {
    var v = null;
    try { v = localStorage.getItem(KEY); } catch (e) { /* private mode */ }
    return v === "light" || v === "dark" ? v : "auto";
  }

  function siteTheme() {
    var t = document.documentElement.dataset.theme;
    return t === "dark" ? "dark" : "light";
  }

  function apply() {
    var p = pref();
    var resolved = p === "auto" ? siteTheme() : p;
    document.documentElement.setAttribute("data-plot-theme", resolved);
    var btn = document.querySelector(".sph-plot-theme-toggle");
    if (btn) {
      btn.title = TITLES[p];
      btn.setAttribute("data-mode", p);
      var badge = btn.querySelector(".sph-plot-theme-badge");
      if (badge) { badge.textContent = ICONS[p]; }
    }
  }

  function cycle() {
    var next = { auto: "light", light: "dark", dark: "auto" }[pref()];
    try { localStorage.setItem(KEY, next); } catch (e) { /* ignore */ }
    apply();
  }

  function insertButton() {
    var end = document.querySelector(".bd-header .navbar-header-items__end");
    if (!end || document.querySelector(".sph-plot-theme-toggle")) { return; }
    var holder = document.createElement("div");
    holder.className = "navbar-item";
    var btn = document.createElement("button");
    btn.type = "button";
    btn.className = "btn btn-sm nav-link pst-navbar-icon sph-plot-theme-toggle";
    btn.setAttribute("aria-label", "Plot color mode");
    // Small custom glyph: an L-shaped plot frame with a 3-segment "data"
    // zigzag in the lower-left, leaving the upper-right open for the mode
    // letter so it stays clearly legible (not fused into the chart).
    btn.innerHTML =
      '<span class="sph-plot-theme-glyph" aria-hidden="true">' +
      '<svg viewBox="0 0 24 24">' +
      '<path class="sph-plot-axes" d="M5 3 V19 H21" fill="none" ' +
      'stroke="currentColor" stroke-width="1.5" stroke-linecap="round" ' +
      'stroke-linejoin="round"/>' +
      '<path class="sph-plot-trace" d="M5 8 L10 17 L13.5 11 L19 18" ' +
      'fill="none" stroke="currentColor" stroke-width="2" ' +
      'stroke-linecap="round" stroke-linejoin="round"/>' +
      '</svg>' +
      '<span class="sph-plot-theme-badge"></span></span>';
    btn.addEventListener("click", cycle);
    holder.appendChild(btn);
    // Place just before the site light/dark switcher when present.
    var themeSwitch = end.querySelector(".theme-switch-container");
    var anchor = themeSwitch ? themeSwitch.closest(".navbar-item") : null;
    end.insertBefore(holder, anchor || end.firstChild);
  }

  // Resolve the URL of the _static/ directory from any asset this page already
  // loaded (works at any page depth, no Sphinx globals required).
  function staticBase() {
    var el = document.querySelector(
      'link[href*="_static/"], script[src*="_static/"]');
    if (!el) { return null; }
    var url = el.href || el.src;
    var i = url.indexOf("_static/");
    return i < 0 ? null : url.slice(0, i + "_static/".length);
  }

  // Tutorial notebooks (nbsphinx) emit a single light figure per output. Pair
  // each with its dark counterpart under _static/nb_dark/ (produced by
  // docs/make_tutorial_dark_figs.py) so the same plot-color toggle that drives
  // the gallery also drives notebook figures. Each dark image is preloaded
  // first: only if it actually exists do we insert it and tag the light one
  // with .plot-light (the shared CSS then shows exactly one per
  // data-plot-theme). Notebooks without shipped dark figures keep their light
  // figure in every mode — no broken images.
  //
  // Naming: a per-notebook manifest maps figure order → a descriptive slug, so
  // dark figures are named nb_dark/<stem>__<slug>.png (self-documenting, no
  // content churn on reorder — only the manifest's order changes). The manifest
  // ships as a JS file nb_dark/<stem>.dark.js (NOT .json) that sets
  //   (window.__SPH_DARK = window.__SPH_DARK || {})["<stem>"] = ["slug", ...];
  // loaded via a <script> tag rather than fetch()/XHR, because those are blocked
  // for file:// resources (local builds) — the same reason the figures
  // themselves are loaded as <img> probes. Slugs are opaque to this code; all
  // slugify logic lives in the generator. A notebook with no manifest (or a
  // figure with no manifested slug) simply keeps its light figure in every mode.
  function pairNotebookFigures() {
    var imgs = document.querySelectorAll(".nboutput img");
    if (!imgs.length) { return; }
    var base = staticBase();
    if (!base) { return; }
    var stem = (location.pathname.split("/").pop() || "").replace(/\.html$/, "");
    if (!stem) { return; }
    var darkDir = base + "nb_dark/";

    function pairOne(light, i, slugs) {
      if (light.dataset.sphPaired) { return; }
      light.dataset.sphPaired = "1";
      var slug = slugs && slugs[i];
      if (!slug) { return; }   // no manifested dark counterpart — keep the light figure
      var probe = new Image();
      probe.onload = function () {
        probe.alt = light.alt;
        probe.className = "plot-dark dark-light";
        // Match the light image's layout (nbsphinx may set an inline width).
        if (light.getAttribute("style")) {
          probe.setAttribute("style", light.getAttribute("style"));
        }
        light.classList.add("plot-light", "dark-light");
        light.parentNode.insertBefore(probe, light.nextSibling);
      };
      // onerror: dark file missing despite the manifest — leave the light figure.
      probe.src = darkDir + stem + "__" + slug + ".png";
    }

    function pairAll(slugs) {
      imgs.forEach(function (light, i) { pairOne(light, i, slugs); });
    }

    // Load the slug manifest via a <script> tag (file://-safe), then pair.
    // Absent/erroring manifest → leave the light figures untouched.
    var reg = (window.__SPH_DARK = window.__SPH_DARK || {});
    if (reg[stem]) { pairAll(reg[stem]); return; }
    var s = document.createElement("script");
    s.onload = function () { pairAll(reg[stem] || null); };
    s.onerror = function () { pairAll(null); };
    s.src = darkDir + stem + ".dark.js";
    document.head.appendChild(s);
  }

  function ready(fn) {
    if (document.readyState !== "loading") { fn(); }
    else { document.addEventListener("DOMContentLoaded", fn); }
  }

  // Apply ASAP (before DOM ready) to minimize light/dark image flicker.
  apply();

  ready(function () {
    insertButton();
    pairNotebookFigures();
    apply();
    // In auto mode, follow site-theme changes live.
    new MutationObserver(function () {
      if (pref() === "auto") { apply(); }
      else { apply(); /* keep button state in sync regardless */ }
    }).observe(document.documentElement, {
      attributes: true, attributeFilter: ["data-theme"],
    });
  });
})();
