/* multicolorfits browser GUI -- vanilla JS, no build step */

'use strict';

const DEFAULT_N_PANELS = 4;
const MAX_PANELS = 16;
const MIN_PANELS = 1;
const $ = (sel, root) => (root || document).querySelector(sel);
const $$ = (sel, root) => Array.from((root || document).querySelectorAll(sel));

const state = {
  panels: [],          // server panel dicts
  compose: {},
  autorefresh: true,
  lastBrowseDir: null, // last directory used in the server file browser
  maxPanels: MAX_PANELS,
  minPanels: MIN_PANELS,
};

function nPanels() {
  return Math.max(state.panels.length, $$('#panels .panel').length);
}

function forEachPanelIndex(fn) {
  const n = nPanels();
  for (let i = 0; i < n; i++) fn(i);
}

try {
  state.lastBrowseDir = localStorage.getItem('mcf-last-browse-dir') || null;
} catch (e) { /* ignore */ }

function rememberBrowseDir(path) {
  if (!path) return;
  let dir = String(path);
  if (/\.(fits|fit|fts|json)(\.gz)?$/i.test(dir)) {
    const cut = Math.max(dir.lastIndexOf('/'), dir.lastIndexOf('\\'));
    if (cut > 0) dir = dir.slice(0, cut);
  }
  if (!dir) return;
  state.lastBrowseDir = dir;
  try { localStorage.setItem('mcf-last-browse-dir', dir); } catch (e) { /* ignore */ }
}

/** Prefer panel file dir → last browse dir → server CWD. */
function preferredBrowseStart(panelIdx) {
  if (panelIdx != null && state.panels[panelIdx] && state.panels[panelIdx].filepath) {
    const fp = state.panels[panelIdx].filepath;
    const cut = Math.max(fp.lastIndexOf('/'), fp.lastIndexOf('\\'));
    if (cut > 0) return fp.slice(0, cut);
  }
  for (const p of state.panels) {
    if (p && p.filepath) {
      const cut = Math.max(p.filepath.lastIndexOf('/'), p.filepath.lastIndexOf('\\'));
      if (cut > 0) return p.filepath.slice(0, cut);
    }
  }
  if (state.lastBrowseDir) return state.lastBrowseDir;
  return '.';  // server process CWD
}

// ---------------------------------------------------------------- utilities

function setStatus(msg, isError) {
  const el = $('#status');
  el.textContent = msg || '';
  el.classList.toggle('error', !!isError);
}

async function api(path, opts) {
  const resp = await fetch(path, opts);
  if (!resp.ok) {
    let detail = resp.statusText;
    try { detail = (await resp.json()).detail || detail; } catch (e) { /* not json */ }
    throw new Error(detail);
  }
  const ctype = resp.headers.get('content-type') || '';
  return ctype.includes('application/json') ? resp.json() : resp.text();
}

async function apiPost(path, payload) {
  return api(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload || {}),
  });
}

function debounce(fn, ms) {
  let t = null;
  return (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
}

function syncSwatchOffset(v) {
  const n = Math.min(1, Math.max(0, parseFloat(v) || 0));
  const slider = $('#swatch-label-offset');
  const num = $('#swatch-label-offset-num');
  if (slider) slider.value = n;
  if (num) num.value = n;
  return n;
}

function syncSwatchInsetScale(v) {
  const n = Math.min(0.45, Math.max(0.08, parseFloat(v) || 0.24));
  const slider = $('#swatch-inset-scale');
  const num = $('#swatch-inset-scale-num');
  if (slider) slider.value = n;
  if (num) num.value = n;
  return n;
}

function normalizePanelColor(c) {
  return PreviewClient.normalizeHexColor(c);
}

function syncPanelColorFields(panel, color) {
  const norm = normalizePanelColor(color);
  $('.color-hex', panel).value = norm.toUpperCase();
  $('.color-picker', panel).value = norm;
  return norm;
}

function syncPanelColorPicker(i) {
  const panel = panelEl(i);
  const p = state.panels[i];
  if (!panel) return;
  syncPanelColorFields(panel, $('.color-hex', panel).value || (p && p.color) || '#ffffff');
}

function syncSwatchSize(v) {
  const n = Math.min(640, Math.max(128, parseInt(v, 10) || 320));
  const el = $('#swatch-size');
  if (el) el.value = n;
  return n;
}

// mpl-style tick colors can be '0.9' greys or hex; convert for the color input
function colorToHexInput(c) {
  if (/^#[0-9a-fA-F]{6}$/.test(c)) return c.toLowerCase();
  const f = parseFloat(c);
  if (!isNaN(f) && f >= 0 && f <= 1) {
    const v = Math.round(f * 255).toString(16).padStart(2, '0');
    return '#' + v + v + v;
  }
  return '#e5e5e5';
}

// Resolve any CSS/matplotlib color string to '#rrggbb', or '' for transparent /
// unresolvable (e.g. 'none').  Handles named colors ('white') via the browser.
function cssColorToHex(c) {
  const s = String(c == null ? '' : c).trim().toLowerCase();
  if (!s || s === 'none' || s === 'transparent') return '';
  const f = parseFloat(s);
  if (!isNaN(f) && String(f) === s && f >= 0 && f <= 1) {
    const v = Math.round(f * 255).toString(16).padStart(2, '0');
    return '#' + v + v + v;
  }
  if (/^#[0-9a-fA-F]{6}$/.test(s)) return s;
  const probe = document.createElement('div');
  probe.style.color = s;
  document.body.appendChild(probe);
  const rgb = getComputedStyle(probe).color;
  probe.remove();
  const m = rgb && rgb.match(/\d+/g);
  if (!m) return '';
  return '#' + m.slice(0, 3).map(x => (+x).toString(16).padStart(2, '0')).join('');
}

// ---------------------------------------------------------------- modal

function showModal(title, bodyEl, buttons) {
  $('#modal-title').textContent = title;
  const body = $('#modal-body');
  body.innerHTML = '';
  body.appendChild(bodyEl);
  const btns = $('#modal-buttons');
  btns.innerHTML = '';
  (buttons || [{ label: 'Close' }]).forEach(b => {
    const btn = document.createElement('button');
    btn.textContent = b.label;
    if (b.primary) btn.classList.add('primary');
    btn.addEventListener('click', async () => {
      if (b.onclick) {
        const keepOpen = await b.onclick();
        if (keepOpen === true) return;
      }
      hideModal();
    });
    btns.appendChild(btn);
  });
  $('#modal-backdrop').hidden = false;
}

function hideModal() { $('#modal-backdrop').hidden = true; }
$('#modal-backdrop').addEventListener('click', e => {
  if (e.target === $('#modal-backdrop')) hideModal();
});

// ---------------------------------------------------------------- session save / load helpers

async function fetchSessionJson() {
  const resp = await fetch('/api/session.json');
  if (!resp.ok) {
    let detail = resp.statusText;
    try { detail = (await resp.json()).detail || detail; } catch (e) { /* not json */ }
    throw new Error(detail);
  }
  const body = await resp.text();
  if (!body || !body.trim()) throw new Error('Session state is empty');
  return body;
}

function triggerDownload(filename, body, mime) {
  const blob = new Blob([body], { type: mime || 'application/octet-stream' });
  const a = document.createElement('a');
  const url = URL.createObjectURL(blob);
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

async function saveSessionToClient(body, defaultName) {
  const handle = await window.showSaveFilePicker({
    suggestedName: defaultName,
    types: [{
      description: 'MultiColorFits session',
      accept: { 'application/json': ['.json'] },
    }],
  });
  const writable = await handle.createWritable();
  await writable.write(body);
  await writable.close();
  setStatus('Session saved (' + handle.name + ')');
}

async function saveSession() {
  const defaultName = 'multicolorfits_session.json';
  try {
    const body = await fetchSessionJson();
    if (typeof window.showSaveFilePicker === 'function') {
      try {
        await saveSessionToClient(body, defaultName);
        return;
      } catch (err) {
        if (err && err.name === 'AbortError') return;
        console.warn('Native save picker failed, opening folder browser', err);
      }
    }
    openSessionSaveBrowser(body, defaultName);
  } catch (err) {
    setStatus('Save session failed: ' + err.message, true);
  }
}

function openSessionSaveBrowser(jsonBody, defaultName) {
  const wrap = document.createElement('div');
  const pathEl = document.createElement('div');
  pathEl.className = 'browser-path';
  const list = document.createElement('ul');
  list.className = 'browser-list';
  const nameRow = document.createElement('div');
  nameRow.className = 'row';
  const nameLabel = document.createElement('label');
  nameLabel.textContent = 'Filename:';
  const nameInput = document.createElement('input');
  nameInput.type = 'text';
  nameInput.value = defaultName || 'multicolorfits_session.json';
  nameInput.size = 36;
  nameRow.append(nameLabel, nameInput);
  wrap.append(pathEl, list, nameRow);

  let currentDir = '';

  async function go(path) {
    let data;
    try { data = await api('/api/browse?path=' + encodeURIComponent(path)); }
    catch (err) { setStatus(err.message, true); return; }
    currentDir = data.path;
    pathEl.textContent = data.path;
    rememberBrowseDir(data.path);
    list.innerHTML = '';
    if (data.parent && data.parent !== data.path) {
      const up = document.createElement('li');
      up.className = 'dir';
      up.textContent = '\u2b06 ..';
      up.addEventListener('click', () => go(data.parent));
      list.appendChild(up);
    }
    data.dirs.forEach(d => {
      const li = document.createElement('li');
      li.className = 'dir';
      li.textContent = '\ud83d\udcc1 ' + d;
      li.addEventListener('click', () => go(data.path + '/' + d));
      list.appendChild(li);
    });
  }

  showModal('Save session JSON', wrap, [
    {
      label: 'Save here', primary: true,
      onclick: async () => {
        if (!currentDir) return true;
        let name = nameInput.value.trim() || defaultName || 'multicolorfits_session.json';
        if (!name.toLowerCase().endsWith('.json')) name += '.json';
        const full = currentDir.replace(/\/$/, '') + '/' + name;
        try {
          await apiPost('/api/session/save', { path: full });
          rememberBrowseDir(full);
          setStatus('Session saved to ' + full);
        } catch (err) { setStatus(err.message, true); return true; }
      },
    },
    {
      label: 'Save to this computer…',
      onclick: async () => {
        try {
          const body = jsonBody || await fetchSessionJson();
          triggerDownload(defaultName || 'multicolorfits_session.json', body, 'application/json');
          setStatus('Session downloaded (' + (defaultName || 'multicolorfits_session.json') + ')');
        } catch (err) {
          setStatus(err.message, true);
          return true;
        }
      },
    },
    { label: 'Cancel' },
  ]);
  go(preferredBrowseStart());
}

// ---------------------------------------------------------------- panels

function updatePanelActions() {
  const addBtn = $('#btn-add-panel');
  if (!addBtn) return;
  const atMax = state.panels.length >= (state.maxPanels || MAX_PANELS);
  addBtn.disabled = atMax;
  addBtn.title = atMax
    ? `At most ${state.maxPanels || MAX_PANELS} panels`
    : 'Add another image panel';
  $$('#panels .btn-remove').forEach((btn) => {
    btn.disabled = state.panels.length <= (state.minPanels || MIN_PANELS);
  });
}

function buildPanels(count) {
  const n = count == null
    ? Math.max(state.panels.length || 0, DEFAULT_N_PANELS)
    : count;
  const container = $('#panels');
  container.innerHTML = '';
  for (let i = 0; i < n; i++) {
    const node = $('#panel-template').content.cloneNode(true);
    const panel = node.querySelector('.panel');
    panel.dataset.idx = i;
    node.querySelector('.pnum').textContent = i + 1;
    if (i === 0) panel.classList.remove('collapsed');
    container.appendChild(node);
    wirePanel(panel, i);
  }
  updatePanelActions();
}

/** Rebuild the left panel strip when the server panel count changes. */
function syncPanelDomToState(preferredOpenIdx) {
  const want = state.panels.length || DEFAULT_N_PANELS;
  const have = $$('#panels .panel').length;
  if (want !== have) {
    let openIdx = preferredOpenIdx;
    if (openIdx == null) {
      const open = $('#panels .panel:not(.collapsed)');
      openIdx = open ? parseInt(open.dataset.idx, 10) : 0;
    }
    buildPanels(want);
    const panel = panelEl(Math.min(Math.max(0, openIdx), want - 1));
    if (panel) panel.classList.remove('collapsed');
  } else {
    updatePanelActions();
  }
}

function panelEl(i) { return $(`.panel[data-idx="${i}"]`); }

const panelSliderPreviewDebouncers = {};

/** Sync percentile limits to the server without touching the panel thumbnail. */
async function syncPanelLevelsToServer(i) {
  const panel = panelEl(i);
  const p = state.panels[i];
  if (!panel || !p || !p.loaded) return;
  try {
    await sendPanelParams(i, {
      percent_min: parseFloat($('.pmin', panel).value),
      percent_max: parseFloat($('.pmax', panel).value),
    }, { skipPreview: true });
  } catch (e) { /* ignore transient sync errors while dragging */ }
}

/**
 * Fetch the server PNG thumbnail. Keeps the client canvas visible until the
 * new image has loaded so we never flash a stale / blank frame.
 */
async function refreshServerPanelPreview(i) {
  const panel = panelEl(i);
  const p = state.panels[i];
  if (!panel || !p || !p.loaded) return;
  await syncPanelLevelsToServer(i);
  const img = $('.preview', panel);
  const empty = $('.preview-empty', panel);
  if (empty) empty.setAttribute('hidden', '');
  const inv = state.compose.inverse ? 'true' : 'false';
  const url = `/api/panel/${i}/preview.png?inverse=${inv}&_=${Date.now()}`;
  const token = (panel._previewToken = (panel._previewToken || 0) + 1);
  await new Promise((resolve) => {
    const done = () => {
      if (panel._previewToken !== token) { resolve(); return; }
      PreviewClient.showServerPanelPreview(panel);
      resolve();
    };
    img.onload = done;
    img.onerror = done;
    img.src = url;
  });
}

function panelLevelSyncDebounced(i) {
  if (!panelSliderPreviewDebouncers[i]) {
    panelSliderPreviewDebouncers[i] = debounce(() => syncPanelLevelsToServer(i), 120);
  }
  panelSliderPreviewDebouncers[i]();
}

function updatePanelLevelPreview(i, opts) {
  const opts_ = opts || {};
  // When Min%/Max% move, push percentile-derived absolute limits into the
  // inputs first so the client renderer (which prefers vmin/vmax) matches.
  if (opts_.fromPercentiles) updateVminVmaxReadout(i);
  const panelOpts = interactiveOptsForPanel(i);
  // Instant client-side thumbnail (and optional live combined). Never swap to
  // the server PNG here — that caused a saturated/flipped flash then blank.
  if (panelOpts.liveCombined) {
    PreviewClient.updatePanelLevels(i, panelOpts.panelEls, panelOpts.activeIndices,
                                    state.compose, true);
  } else if (PreviewClient.hasBuffer(i)) {
    PreviewClient.renderPanelThumbnail(panelEl(i), state.compose);
  }
  // While dragging, only the canvas updates; server state catches up on pointerup.
  // Keyboard / click adjustments (not dragging) debounce a params-only sync.
  if (!opts_.skipServerSync && !PreviewClient.isLevelDragging()) {
    panelLevelSyncDebounced(i);
  }
}

function updateVminVmaxReadout(i) {
  const panel = panelEl(i);
  // Always derive from percentile sliders on the preview buffer (not the
  // absolute inputs), so dragging Min%/Max% updates the readout correctly.
  const limits = PreviewClient.percentileLimitsForPanel(panel);
  if (!limits) return;
  $('.vmin', panel).value = Number(limits.vmin.toPrecision(6));
  $('.vmax', panel).value = Number(limits.vmax.toPrecision(6));
}

function interactiveOptsForPanel(panelIdx) {
  const active = state.panels.map((p, j) => (p && p.loaded) ? j : -1).filter(j => j >= 0);
  return {
    panelEls: state.panels.map((_, j) => panelEl(j)),
    activeIndices: active,
    liveCombined: livePreviewOn() && combinedIsVisible(),
    panelIdx,
  };
}

function serverPanelPreviewDebounced(i) {
  panelLevelSyncDebounced(i);
}

function wireComposeTabs() {
  const tabs = $$('.compose-tab');
  const panels = $$('.compose-tabpanel');
  tabs.forEach(tab => {
    tab.addEventListener('click', () => {
      const name = tab.dataset.tab;
      tabs.forEach(t => {
        const on = t === tab;
        t.classList.toggle('active', on);
        t.setAttribute('aria-selected', on ? 'true' : 'false');
      });
      panels.forEach(p => {
        const on = p.id === `compose-tab-${name}`;
        p.classList.toggle('active', on);
        if (on) p.removeAttribute('hidden');
        else p.setAttribute('hidden', '');
      });
    });
  });
}

function wirePanel(panel, i) {
  const q = sel => $(sel, panel);

  q('.panel-head').addEventListener('click', () => {
    const wasCollapsed = panel.classList.contains('collapsed');
    panel.classList.toggle('collapsed');
    if (wasCollapsed) syncPanelColorPicker(i);
  });

  q('.btn-load').addEventListener('click', async () => {
    const path = q('.filepath').value.trim();
    if (!path) { setStatus('Enter a FITS file path first', true); return; }
    await loadPanelFromPath(i, path);
  });
  q('.filepath').addEventListener('keydown', e => {
    if (e.key === 'Enter') q('.btn-load').click();
  });

  q('.btn-browse').addEventListener('click', () => openFileBrowser(i));

  q('.fileupload').addEventListener('change', async e => {
    const file = e.target.files[0];
    if (!file) return;
    const form = new FormData();
    form.append('file', file);
    try {
      setStatus(`Uploading ${file.name}…`);
      const p = await api(`/api/panel/${i}/upload`, { method: 'POST', body: form });
      onPanelLoaded(i, p);
    } catch (err) { setStatus('Upload failed: ' + err.message, true); }
    e.target.value = '';
  });

  q('.stretch').addEventListener('change', async () => {
    // Stretch is applied client-side to the existing float buffer — do not
    // invalidate/reload. (Invalidating forced a server-PNG fallback and hid
    // first-load client preview bugs behind a one-time "reselect stretch" fix.)
    await sendPanelParams(i, { stretch: q('.stretch').value });
  });

  // color picker <-> hex text, kept in sync
  function applyPanelColor(v) {
    const color = normalizePanelColor(v);
    syncPanelColorFields(panel, color);
    updatePanelLevelPreview(i, { skipServerSync: true });
    sendPanelParams(i, { color }, { skipPreview: true });
    if (combinedIsVisible()) refreshCombined();
    refreshCvd();
  }
  q('.color-picker').addEventListener('input', debounce(() => {
    applyPanelColor(q('.color-picker').value);
  }, 120));
  q('.color-hex').addEventListener('change', () => {
    let v = q('.color-hex').value.trim();
    if (!v.startsWith('#')) v = '#' + v;
    if (/^#[0-9a-fA-F]{6}$/.test(v)) applyPanelColor(v);
    else setStatus('Color must be a #RRGGBB hex string', true);
  });
  q('.color-hex').addEventListener('keydown', e => {
    if (e.key !== 'Enter') return;
    e.preventDefault();
    e.stopPropagation();
    let v = q('.color-hex').value.trim();
    if (!v.startsWith('#')) v = '#' + v;
    if (/^#[0-9a-fA-F]{6}$/.test(v)) {
      applyPanelColor(v);
      q('.color-hex').blur();
    } else {
      setStatus('Color must be a #RRGGBB hex string', true);
    }
  });

  q('.panel-label').addEventListener('change', () => {
    sendPanelParams(i, { label: q('.panel-label').value });
    if (state.compose.show_legend) refreshCombined();
  });

  q('.vmin').addEventListener('change', () => sendPanelParams(i, { vmin: parseFloat(q('.vmin').value) }));
  q('.vmax').addEventListener('change', () => sendPanelParams(i, { vmax: parseFloat(q('.vmax').value) }));

  function flushLevelSync() {
    PreviewClient.setLevelDragging(false);
    // Sync limits to the server; keep the client canvas (no PNG swap flash).
    syncPanelLevelsToServer(i);
    if (!PreviewClient.hasBuffer(i)) refreshServerPanelPreview(i);
    else if (livePreviewOn() && combinedIsVisible()) refreshCombined();
  }

  q('.pmin').addEventListener('pointerdown', () => PreviewClient.setLevelDragging(true));
  q('.pmax').addEventListener('pointerdown', () => PreviewClient.setLevelDragging(true));
  q('.pmin').addEventListener('pointerup', flushLevelSync);
  q('.pmax').addEventListener('pointerup', flushLevelSync);
  q('.pmin').addEventListener('pointercancel', flushLevelSync);
  q('.pmax').addEventListener('pointercancel', flushLevelSync);

  q('.pmin').addEventListener('input', () => {
    q('.pmin-val').textContent = q('.pmin').value;
    updatePanelLevelPreview(i, { fromPercentiles: true });
  });
  q('.pmax').addEventListener('input', () => {
    q('.pmax-val').textContent = q('.pmax').value;
    updatePanelLevelPreview(i, { fromPercentiles: true });
  });

  q('.btn-minmax').addEventListener('click', async () => {
    try { updatePanelUI(i, await apiPost(`/api/panel/${i}/minmax`)); refreshPanelPreview(i); }
    catch (err) { setStatus(err.message, true); }
  });
  q('.btn-zscale').addEventListener('click', async () => {
    try { updatePanelUI(i, await apiPost(`/api/panel/${i}/zscale`)); refreshPanelPreview(i); }
    catch (err) { setStatus(err.message, true); }
  });
  q('.btn-header').addEventListener('click', () => openHeaderEditor(i));
  q('.btn-clear').addEventListener('click', async () => {
    await apiPost(`/api/panel/${i}/clear`);
    PreviewClient.invalidatePanel(i);
    state.panels[i] = { loaded: false };
    panel.classList.remove('loaded');
    $('.panel-file', panel).textContent = '';
    $('.preview', panel).hidden = true;
    const canvas = $('.preview-canvas', panel);
    if (canvas) canvas.hidden = true;
    $('.preview-empty', panel).hidden = false;
    setStatus(`Panel ${i + 1} cleared`);
    checkGridStatus();
    refreshCvd();
  });
  q('.btn-remove').addEventListener('click', (ev) => {
    ev.stopPropagation();
    confirmRemovePanel(i);
  });

  q('.smooth').addEventListener('change', () => sendPanelParams(i, { smooth: q('.smooth').checked }));
  q('.smooth-sigma').addEventListener('change', () =>
    sendPanelParams(i, { smooth_sigma: parseFloat(q('.smooth-sigma').value), smooth: q('.smooth').checked }));
}

async function addPanelSlot() {
  try {
    const res = await apiPost('/api/panels/add', {});
    const idx = res.index;
    await applyServerState(res.state, { openPanelIdx: idx });
    setStatus(`Added image panel ${idx + 1}`);
  } catch (err) {
    setStatus(err.message, true);
  }
}

function confirmRemovePanel(i) {
  const p = state.panels[i];
  const loaded = p && p.loaded;
  const msg = document.createElement('p');
  msg.textContent = loaded
    ? `Remove Image ${i + 1} and unload its FITS data? Later panels will renumber.`
    : `Remove Image ${i + 1}? Later panels will renumber.`;
  showModal('Remove panel', msg, [
    {
      label: 'Remove', primary: true,
      onclick: async () => {
        try {
          const res = await apiPost(`/api/panels/${i}/remove`);
          forEachPanelIndex((j) => PreviewClient.invalidatePanel(j));
          await applyServerState(res.state, {
            openPanelIdx: Math.max(0, Math.min(i, (res.state.panels || []).length - 1)),
          });
          setStatus(`Removed panel ${i + 1}`);
        } catch (err) {
          setStatus('Remove failed: ' + err.message, true);
          return true;
        }
      },
    },
    { label: 'Cancel' },
  ]);
}

async function loadPanelFromPath(i, path) {
  try {
    setStatus(`Loading ${path}…`);
    const p = await apiPost(`/api/panel/${i}/load`, { path });
    onPanelLoaded(i, p);
  } catch (err) { setStatus('Load failed: ' + err.message, true); }
}

async function onPanelLoaded(i, p) {
  state.panels[i] = p;
  updatePanelUI(i, p);
  if (p.filepath) rememberBrowseDir(p.filepath);
  panelEl(i).classList.remove('collapsed');
  // Always drop the previous float buffer — ensurePanelBuffer() would keep it
  // and the thumbnail would still show the old image after a replace-load.
  PreviewClient.invalidatePanel(i);
  let bufOk = false;
  try {
    await PreviewClient.loadPanelBuffer(i);
    bufOk = true;
  } catch (e) { /* fall back to server PNG */ }
  if (bufOk && PreviewClient.renderPanel(panelEl(i), state.compose)) {
    setStatus(`Panel ${i + 1}: loaded ${p.filepath} (${p.shape[1]}\u00d7${p.shape[0]})`);
  } else {
    refreshPanelPreview(i);
    setStatus(`Panel ${i + 1}: loaded ${p.filepath} (${p.shape[1]}\u00d7${p.shape[0]})`);
  }
  checkGridStatus();
  refreshCvd();
  if (combinedIsVisible()) refreshCombined();
}

// ---------------------------------------------------------------- alignment

const ALIGN_LABELS = {
  reference: 'Reference panel grid',
  icrs: 'ICRS (equatorial)',
  galactic: 'Galactic',
  fk5: 'FK5 (J2000)',
  fk4: 'FK4 (B1950)',
  ecliptic: 'Ecliptic',
};

async function checkGridStatus() {
  let status;
  try { status = await api('/api/grid_status'); }
  catch (err) { return; }
  state.gridStatus = status;
  const banner = $('#grid-warning');
  $('#btn-align').disabled = !status.reproject_available;
  if (status.aligned || (status.mismatched || []).length === 0) {
    banner.hidden = true;
    return;
  }
  const which = status.mismatched.map(n => n + 1).join(', ');
  const refShape = status.shapes[status.reference];
  const plural = status.mismatched.length > 1;
  let msg = `Panel${plural ? 's' : ''} ${which} ${plural ? "don't" : "doesn't"} `
    + `share the reference grid (panel ${status.reference + 1}`
    + (refShape ? `, ${refShape[1]}\u00d7${refShape[0]}` : '') + `). `
    + 'Use Align layers… in the toolbar to reproject onto a common grid before combining.';
  if (!status.reproject_available) {
    msg += ' (Install the "reproject" package to enable alignment.)';
  }
  $('#grid-warning-text').textContent = msg;
  banner.hidden = false;
}

function showAlignDialog() {
  const status = state.gridStatus || {};
  const targets = status.align_targets || ['reference', 'icrs', 'galactic'];
  const active = status.active || state.panels.map((p, i) => (p && p.loaded) ? i : -1).filter(i => i >= 0);

  const form = document.createElement('div');
  form.className = 'save-form';

  const tLabel = document.createElement('label');
  tLabel.textContent = 'Reproject all loaded layers onto:';
  const tSel = document.createElement('select');
  targets.forEach(t => {
    const o = document.createElement('option');
    o.value = t; o.textContent = ALIGN_LABELS[t] || t;
    tSel.appendChild(o);
  });

  const rLabel = document.createElement('label');
  rLabel.textContent = 'Reference panel (defines center / pixel scale):';
  const rSel = document.createElement('select');
  active.forEach(i => {
    const o = document.createElement('option');
    o.value = i; o.textContent = `Panel ${i + 1}`;
    rSel.appendChild(o);
  });
  if (status.reference != null) rSel.value = String(status.reference);

  form.append(tLabel, tSel, rLabel, rSel);

  showModal('Align layers (reproject)', form, [
    {
      label: 'Align', primary: true,
      onclick: async () => {
        await doAlign(tSel.value, parseInt(rSel.value, 10));
      },
    },
    { label: 'Cancel' },
  ]);
}

async function doAlign(target, reference) {
  try {
    setStatus('Reprojecting layers…');
    const res = await apiPost('/api/align', { target, reference });
    (res.state.panels || []).forEach((p, i) => {
      state.panels[i] = p;
      if (p.loaded) {
        updatePanelUI(i, p);
        PreviewClient.invalidatePanel(i);
        PreviewClient.loadPanelBuffer(i).then(() => {
          if (!PreviewClient.renderPanel(panelEl(i), state.compose)) refreshPanelPreview(i);
        }).catch(() => refreshPanelPreview(i));
      }
    });
    setStatus(res.result.message || 'Layers aligned');
    await checkGridStatus();
    refreshCombined();
  } catch (err) {
    setStatus('Alignment failed: ' + err.message, true);
  }
}

// ---------------------------------------------------------------- palettes

async function wirePalettes() {
  const sel = $('#palette-select');
  try {
    const data = await api('/api/palettes');
    state.palettes = {};
    (data.palettes || []).forEach(p => { state.palettes[p.name] = p.colors; });
    const menu = data.menu || [];
    if (menu.length) {
      menu.forEach(({ group, items }) => {
        const og = document.createElement('optgroup');
        og.label = group;
        items.forEach(({ label, key }) => {
          const o = document.createElement('option');
          o.value = key;
          o.textContent = label;
          og.appendChild(o);
        });
        sel.appendChild(og);
      });
    } else {
      const opt = document.createElement('option');
      opt.value = 'perceptual';
      opt.textContent = 'Perceptual (auto)';
      sel.appendChild(opt);
      (data.palettes || []).forEach(p => {
        const o = document.createElement('option');
        o.value = p.name;
        o.textContent = p.name + '  (' + p.colors.length + ')';
        sel.appendChild(o);
      });
    }
  } catch (err) { /* palettes optional */ }
  sel.addEventListener('change', updatePaletteSwatches);
  updatePaletteSwatches();
  $('#btn-apply-palette').addEventListener('click', applyPalette);
  $('#btn-tune-colors').addEventListener('click', showHueTuneDialog);
}

function activePanelColors() {
  const colors = [];
  forEachPanelIndex((i) => {
    const p = state.panels[i];
    if (p && p.loaded) colors.push(p.color);
  });
  return colors;
}

async function fetchSwatchPreview(deltaDeg, baseColors) {
  const r = await fetch('/api/preview_swatch', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ delta_deg: deltaDeg, base_colors: baseColors }),
  });
  if (!r.ok) {
    let msg = r.statusText;
    try { msg = (await r.json()).detail || msg; } catch (_) { /* ignore */ }
    throw new Error(msg);
  }
  return r.blob();
}

async function showHueTuneDialog() {
  const baseColors = activePanelColors();
  if (!baseColors.length) {
    setStatus('Load at least one image first', true);
    return;
  }

  const form = document.createElement('div');
  form.className = 'hue-tune-dialog';
  form.innerHTML = `
    <p>Drag the slider to rotate hues while keeping their relative spacing.
    The overlapping-circle preview uses the same composite-mode mixing as your combined image.
    Layer colors apply when you click OK.</p>
    <div class="hue-swatch-wrap"><img id="hue-swatch-preview" alt="Color combination preview"></div>
    <div class="hue-rotation-row">
      <label style="flex:1">Hue rotation
        <input type="range" id="hue-rotation" min="0" max="360" step="1" value="0">
      </label>
      <span id="hue-rotation-val">0°</span>
    </div>
    <label class="chk" title="Refresh panel thumbnails while dragging (slower on large images)">
      <input type="checkbox" id="hue-live-panels"> Update panel previews while dragging
    </label>`;

  const slider = form.querySelector('#hue-rotation');
  const valEl = form.querySelector('#hue-rotation-val');
  const img = form.querySelector('#hue-swatch-preview');
  const liveCheck = form.querySelector('#hue-live-panels');
  let previewUrl = null;
  let liveApplied = false;
  let previewTimer = null;

  const setPreviewUrl = (blob) => {
    if (previewUrl) URL.revokeObjectURL(previewUrl);
    previewUrl = URL.createObjectURL(blob);
    img.src = previewUrl;
  };

  const refreshSwatch = async (deg) => {
    try {
      setPreviewUrl(await fetchSwatchPreview(deg, baseColors));
    } catch (err) {
      setStatus('Swatch preview failed: ' + err.message, true);
    }
  };

  const applyToPanels = async (deg, refreshPanels) => {
    const res = await apiPost('/api/rotate_hues', { delta_deg: deg, base_colors: baseColors });
    (res.state.panels || []).forEach((p, i) => {
      state.panels[i] = p;
      if (p.loaded) updatePanelUI(i, p);
    });
    renderCvdStatus(res.colorblind);
    if (refreshPanels && state.autorefresh) {
      forEachPanelIndex((i) => {
        if (state.panels[i] && state.panels[i].loaded) refreshPanelPreview(i);
      });
    }
    if (deg !== 0) liveApplied = true;
    return res;
  };

  const schedulePreview = (deg) => {
    clearTimeout(previewTimer);
    previewTimer = setTimeout(async () => {
      await refreshSwatch(deg);
      if (liveCheck.checked) {
        try {
          await applyToPanels(deg, true);
          if (livePreviewOn()) refreshCombined();
        } catch (err) {
          setStatus('Hue update failed: ' + err.message, true);
        }
      }
    }, 50);
  };

  slider.addEventListener('input', () => {
    const deg = parseInt(slider.value, 10);
    valEl.textContent = deg + '°';
    schedulePreview(deg);
  });

  showModal('Tune layer colors', form, [
    {
      label: 'Cancel',
      onclick: async () => {
        if (liveApplied) {
          try {
            await applyToPanels(0, state.autorefresh);
            if (livePreviewOn()) refreshCombined();
          } catch (err) {
            setStatus('Could not restore colors: ' + err.message, true);
          }
        }
        if (previewUrl) URL.revokeObjectURL(previewUrl);
        clearTimeout(previewTimer);
      },
    },
    {
      label: 'OK', primary: true,
      onclick: async () => {
        const deg = parseInt(slider.value, 10);
        try {
          await applyToPanels(deg, true);
          refreshCombined();
          if (previewUrl) URL.revokeObjectURL(previewUrl);
          clearTimeout(previewTimer);
          setStatus('Layer colors updated');
        } catch (err) {
          setStatus('Hue tune failed: ' + err.message, true);
          return true;
        }
      },
    },
  ]);

  await refreshSwatch(0);
}

function updatePaletteSwatches() {
  const sel = $('#palette-select');
  const box = $('#palette-swatches');
  box.innerHTML = '';
  const key = sel.value;
  const colors = (state.palettes || {})[key];
  if (!colors) {  // perceptual / hue geometry: preview depends on #loaded layers
    box.textContent = 'auto';
    return;
  }
  colors.forEach(c => {
    const sw = document.createElement('i');
    sw.style.background = c;
    sw.title = c;
    box.appendChild(sw);
  });
}

async function applyPalette() {
  const name = $('#palette-select').value;
  try {
    const res = await apiPost('/api/apply_palette', { name });
    (res.state.panels || []).forEach((p, i) => {
      state.panels[i] = p;
      if (p.loaded) { updatePanelUI(i, p); refreshPanelPreview(i); }
    });
    renderCvdStatus(res.colorblind);
    refreshCombined();
    const applied = (res.result.colors || []).length;
    setStatus(applied ? `Applied ${name} palette to ${applied} layer(s)` : 'No layers loaded');
  } catch (err) {
    setStatus('Palette failed: ' + err.message, true);
  }
}

async function refreshCvd() {
  try { renderCvdStatus(await api('/api/colorblind')); }
  catch (err) { /* ignore */ }
}

function renderCvdStatus(report) {
  const el = $('#cvd-status');
  if (!report || (report.colors || []).length < 2) { el.textContent = ''; return; }
  if (report.ok) {
    el.textContent = '\u2713 colorblind-safe';
    el.className = 'cvd-status ok';
    el.title = 'All layer colors stay distinguishable under protan/deutan/tritan simulation.';
    return;
  }
  const bits = [];
  Object.keys(report.kinds).forEach(kind => {
    const f = report.kinds[kind].failures;
    if (f && f.length) {
      const pairs = f.map(([a, b]) => `${a}&${b}`).join(', ');
      bits.push(`${kind.slice(0, 6)}: ${pairs}`);
    }
  });
  el.textContent = '\u26a0 confusable (' + bits.join('; ') + ')';
  el.className = 'cvd-status warn';
  el.title = 'These layer pairs may look alike under color-vision deficiency. '
    + 'Try a different palette or the Perceptual (auto) option.';
}

function updatePanelUI(i, p) {
  state.panels[i] = p;
  const panel = panelEl(i);
  if (!panel) return;
  panel.classList.toggle('loaded', !!p.loaded);
  $('.panel-file', panel).textContent = p.filepath ? p.filepath.split('/').pop() : '';
  $('.filepath', panel).value = p.filepath || '';
  // Percentiles and stretch before absolute limits so a bad vmin cannot abort
  // the rest of the sync (which previously left sliders at 0–100 defaults).
  if (p.stretch != null) $('.stretch', panel).value = p.stretch;
  if (p.percent_min != null && Number.isFinite(+p.percent_min)) {
    $('.pmin', panel).value = p.percent_min;
    $('.pmin-val', panel).textContent = (+p.percent_min).toFixed(2);
  }
  if (p.percent_max != null && Number.isFinite(+p.percent_max)) {
    $('.pmax', panel).value = p.percent_max;
    $('.pmax-val', panel).textContent = (+p.percent_max).toFixed(2);
  }
  if (p.vmin != null && Number.isFinite(+p.vmin)) {
    $('.vmin', panel).value = Number((+p.vmin).toPrecision(6));
  } else {
    $('.vmin', panel).value = '';
  }
  if (p.vmax != null && Number.isFinite(+p.vmax)) {
    $('.vmax', panel).value = Number((+p.vmax).toPrecision(6));
  } else {
    $('.vmax', panel).value = '';
  }
  if (p.color != null) syncPanelColorFields(panel, p.color);
  if (p.label !== undefined && document.activeElement !== $('.panel-label', panel)) {
    $('.panel-label', panel).value = p.label;
  }
  $('.smooth', panel).checked = p.smooth;
  $('.smooth-sigma', panel).value = p.smooth_sigma;
}

async function sendPanelParams(i, payload, opts) {
  if (!state.panels[i] || !state.panels[i].loaded) {
    // Allow setting color/stretch before an image is loaded; server keeps state
  }
  try {
    const p = await apiPost(`/api/panel/${i}/params`, payload);
    updatePanelUI(i, p);
    const skipPreview = opts && opts.skipPreview;
    if (!skipPreview && state.autorefresh && p.loaded) {
      if (!PreviewClient.renderPanel(panelEl(i), state.compose)) refreshPanelPreview(i);
    }
    if (!skipPreview && p.loaded && livePreviewOn() && combinedIsVisible()) {
      const active = state.panels.map((pp, j) => (pp && pp.loaded) ? j : -1).filter(j => j >= 0);
      const els = state.panels.map((_, j) => panelEl(j));
      if (!PreviewClient.canRenderCombined(active, els, state.compose)
          || !PreviewClient.renderCombined(active, els, state.compose)) {
        refreshCombined();
      }
    } else if (!skipPreview && p.loaded && livePreviewOn()) {
      refreshCombined();
    }
  } catch (err) { setStatus(err.message, true); }
}

function showCombinedPlaceholder() {
  $('#combined-img').hidden = true;
  $('#combined-canvas').hidden = true;
  $('#combined-placeholder').hidden = false;
  $('#combined-spinner').hidden = true;
  const el = $('#cursor-readout');
  if (el) {
    el.textContent = '';
    delete el.dataset.baseReadout;
  }
}

function confirmResetSession() {
  const msg = document.createElement('p');
  msg.textContent = 'Unload all FITS images and restore default compose settings? This cannot be undone.';
  showModal('Reset session', msg, [
    {
      label: 'Reset', primary: true,
      onclick: async () => {
        try {
          const r = await apiPost('/api/session/reset');
          forEachPanelIndex((i) => PreviewClient.invalidatePanel(i));
          $('#live-preview').checked = true;
          await applyServerState(r.state);
          forEachPanelIndex((i) => {
            const panel = panelEl(i);
            if (!panel) return;
            $('.preview', panel).hidden = true;
            const canvas = $('.preview-canvas', panel);
            if (canvas) canvas.hidden = true;
            $('.preview-empty', panel).hidden = false;
          });
          showCombinedPlaceholder();
          setStatus('Session reset — all panels cleared');
        } catch (err) {
          setStatus('Reset failed: ' + err.message, true);
          return true;
        }
      },
    },
    { label: 'Cancel' },
  ]);
}

function combinedIsVisible() {
  return !$('#combined-img').hidden || !$('#combined-canvas').hidden;
}

function refreshPanelPreview(i) {
  const p = state.panels[i];
  if (!p || !p.loaded) return;
  try {
    if (typeof PreviewClient !== 'undefined'
        && PreviewClient.renderPanel(panelEl(i), state.compose)) return;
  } catch (err) {
    console.warn('client panel preview failed; using server PNG', err);
  }
  if (typeof PreviewClient !== 'undefined') PreviewClient.showServerPanelPreview(panelEl(i));
  const img = $('.preview', panelEl(i));
  const inv = state.compose.inverse ? 'true' : 'false';
  img.src = `/api/panel/${i}/preview.png?inverse=${inv}&_=${Date.now()}`;
  img.removeAttribute('hidden');
  $('.preview-empty', panelEl(i)).setAttribute('hidden', '');
}

function refreshAllPreviews() {
  forEachPanelIndex((i) => refreshPanelPreview(i));
}

// ---------------------------------------------------------------- header editor

async function openHeaderEditor(i) {
  let text;
  try { text = await api(`/api/panel/${i}/header`); }
  catch (err) { setStatus(err.message, true); return; }
  const ta = document.createElement('textarea');
  ta.value = text;
  showModal(`FITS Header \u2014 Image ${i + 1}`, ta, [
    {
      label: 'Apply', primary: true,
      onclick: async () => {
        try { await apiPost(`/api/panel/${i}/header`, { text: ta.value }); setStatus('Header updated'); }
        catch (err) { setStatus(err.message, true); return true; }
      },
    },
    { label: 'Cancel' },
  ]);
}

// ---------------------------------------------------------------- file browser

async function openFileBrowser(i, startPath) {
  const wrap = document.createElement('div');
  const pathEl = document.createElement('div');
  pathEl.className = 'browser-path';
  const list = document.createElement('ul');
  list.className = 'browser-list';
  wrap.appendChild(pathEl);
  wrap.appendChild(list);

  async function go(path) {
    let data;
    try { data = await api('/api/browse?path=' + encodeURIComponent(path)); }
    catch (err) { setStatus(err.message, true); return; }
    pathEl.textContent = data.path;
    rememberBrowseDir(data.path);
    list.innerHTML = '';
    if (data.parent && data.parent !== data.path) {
      const up = document.createElement('li');
      up.className = 'dir';
      up.textContent = '\u2b06 ..';
      up.addEventListener('click', () => go(data.parent));
      list.appendChild(up);
    }
    data.dirs.forEach(d => {
      const li = document.createElement('li');
      li.className = 'dir';
      li.textContent = '\ud83d\udcc1 ' + d;
      li.addEventListener('click', () => go(data.path + '/' + d));
      list.appendChild(li);
    });
    data.files.forEach(f => {
      const li = document.createElement('li');
      li.className = 'fits';
      li.textContent = f;
      li.addEventListener('click', async () => {
        hideModal();
        rememberBrowseDir(data.path);
        await loadPanelFromPath(i, data.path + '/' + f);
      });
      list.appendChild(li);
    });
  }

  showModal(`Select FITS file \u2014 Image ${i + 1}`, wrap, [{ label: 'Cancel' }]);
  go(startPath || preferredBrowseStart(i));
}

// ---------------------------------------------------------------- combined view

function livePreviewOn() { return $('#live-preview').checked; }

// Blend only applies to the perceptual color-space modes (not RGB / subtractive).
function updateBlendEnabled(mode) {
  $('#combine-blend').disabled = (mode === 'rgb' || mode === 'ryb' || mode === 'cmyk');
}

function subtractiveModesDefaultWhite(mode) {
  return mode === 'ryb' || mode === 'cmyk';
}

// Current compositing background from the Background select (+ custom color).
function currentBackground() {
  const sel = $('#combine-bg').value;
  if (sel === 'custom') return $('#combine-bg-color').value.trim() || '#12243a';
  return sel;   // 'black' | 'white' | 'transparent'
}
function updateBgUI() {
  $('#combine-bg-pair').style.display = ($('#combine-bg').value === 'custom') ? '' : 'none';
}

// Fast, downsampled float32 pixel preview (no WCS axes) -- for live updates.
// Prefers the in-browser GPU/CPU path using the same buffers as left panels;
// falls back to a server-side downsampled PNG when the mode isn't client-capable.
function livePreview() {
  const anyLoaded = state.panels.some(p => p && p.loaded);
  if (!anyLoaded) return Promise.resolve(false);
  const active = state.panels.map((p, j) => (p && p.loaded) ? j : -1).filter(j => j >= 0);
  const els = state.panels.map((_, j) => panelEl(j));
  if (PreviewClient.canRenderCombined(active, els, state.compose)
      && PreviewClient.renderCombined(active, els, state.compose)) {
    const mode = PreviewClient.usesGpu() ? 'GPU' : 'CPU';
    setStatus(`Fast preview (${mode}, γ=${Number(state.compose.gamma).toFixed(2)})`);
    return Promise.resolve(true);
  }
  PreviewClient.showServerCombined();
  const img = $('#combined-img');
  const inv = state.compose.inverse ? 'true' : 'false';
  const g = state.compose.gamma;
  const mode = state.compose.combine_mode || 'rgb';
  // Include mode+gamma in the URL so caches cannot reuse a stale PNG when
  // only compose params change (Date.now alone is not always enough).
  const url = `/api/combined.png?preview=true&max_size=512&inverse=${inv}`
    + `&mode=${encodeURIComponent(mode)}&gamma=${encodeURIComponent(g)}&_=${Date.now()}`;
  return new Promise((resolve) => {
    img.onload = () => {
      img.hidden = false;
      $('#combined-placeholder').hidden = true;
      setStatus(`Fast preview (server, ${mode}, γ=${Number(g).toFixed(2)})`);
      resolve(true);
    };
    img.onerror = () => {
      setStatus('Fast preview failed to load', true);
      resolve(false);
    };
    img.src = url;
  });
}
const livePreviewDebounced = debounce(() => { livePreview(); }, 120);

// Refresh the combined view after a settings change.
// With Fast preview on: cheap update. Off: only refresh if a figure is already showing
// (avoid surprising full WCS renders on every tweak).
function refreshCombined(opts) {
  const force = !!(opts && opts.force);
  if (!force && $('#combined-img').hidden && $('#combined-canvas').hidden) return;
  if (livePreviewOn()) livePreviewDebounced();
  else if (force || combinedIsVisible()) plotCombined();
}

async function plotCombined() {
  const anyLoaded = state.panels.some(p => p && p.loaded);
  if (!anyLoaded) { setStatus('Load at least one image first', true); return; }
  // Always full-resolution WCS figure — Fast preview handles interactive updates.
  PreviewClient.showServerCombined();
  $('#combined-spinner').hidden = false;
  setStatus('Rendering full WCS plot…');
  const img = $('#combined-img');
  let plotUrl = null;
  try {
    const resp = await fetch(`/api/combined_plot.png?_=${Date.now()}`);
    if (!resp.ok) {
      let detail = resp.statusText;
      try { detail = (await resp.json()).detail || detail; } catch (e) { /* not json */ }
      throw new Error(detail);
    }
    const blob = await resp.blob();
    plotUrl = URL.createObjectURL(blob);
    await new Promise((resolve, reject) => {
      img.onload = () => resolve();
      img.onerror = () => reject(new Error('Could not display combined plot image'));
      img.src = plotUrl;
    });
    img.hidden = false;
    $('#combined-canvas').hidden = true;
    $('#combined-placeholder').hidden = true;
    setStatus('Full WCS plot updated');
  } catch (err) {
    setStatus('Plot failed: ' + err.message, true);
  } finally {
    $('#combined-spinner').hidden = true;
    if (plotUrl) URL.revokeObjectURL(plotUrl);
  }
}

const sendCompose = async (payload, replot) => {
  try {
    state.compose = await apiPost('/api/compose', payload);
    if (replot !== false) refreshCombined();
  } catch (err) { setStatus(err.message, true); }
};

function wireCombined() {
  $('#btn-plot-combined').addEventListener('click', plotCombined);
  $('#btn-align').addEventListener('click', showAlignDialog);

  function bindGammaRefresh() {
    const el = $('#gamma');
    if (!el) return;
    const apply = debounce(async () => {
      const g = parseFloat(el.value);
      if (!Number.isFinite(g) || g <= 0) return;
      try {
        state.compose = await apiPost('/api/compose', { gamma: g });
      } catch (err) {
        setStatus(err.message, true);
        return;
      }
      const applied = state.compose.gamma;
      // Left-panel thumbs intentionally ignore compose gamma (it cancels in
      // the display path). Always refresh the *combined* view here.
      if (livePreviewOn()) {
        setStatus(`Gamma ${Number(applied).toFixed(2)} — updating combined…`);
        await livePreview();
      } else if (combinedIsVisible()) {
        setStatus(`Gamma ${Number(applied).toFixed(2)} — rendering full plot…`);
        await plotCombined();
      } else {
        setStatus(`Gamma set to ${Number(applied).toFixed(2)} — enable Fast preview or Plot Full Resolution to see it`);
      }
    }, 150);
    el.addEventListener('input', apply);
    el.addEventListener('change', apply);
  }

  bindGammaRefresh();
  $('#inverse').addEventListener('change', () => sendCompose({ inverse: $('#inverse').checked }));
  $('#autorefresh').addEventListener('change', () => { state.autorefresh = $('#autorefresh').checked; });
  $('#panel-preview-hd').addEventListener('change', async () => {
    await sendCompose({ panel_preview_hd: $('#panel-preview-hd').checked }, false);
    if (state.autorefresh) refreshAllPreviews();
  });
  $('#live-preview').addEventListener('change', () => {
    try { localStorage.setItem('mcf-fast-preview', $('#live-preview').checked ? '1' : '0'); } catch (e) { /* ignore */ }
    const anyLoaded = state.panels.some(p => p && p.loaded);
    if (!anyLoaded) return;
    if (livePreviewOn()) livePreview();
    else plotCombined();
  });

  $('#tickcolor').addEventListener('change', () => {
    $('#tickcolor-picker').value = colorToHexInput($('#tickcolor').value.trim());
    sendCompose({ tickcolor: $('#tickcolor').value.trim() });
  });
  $('#tickcolor-picker').addEventListener('input', debounce(() => {
    $('#tickcolor').value = $('#tickcolor-picker').value;
    sendCompose({ tickcolor: $('#tickcolor-picker').value });
  }, 200));

  $('#facecolor').addEventListener('change', () => {
    const val = $('#facecolor').value.trim();
    const hex = cssColorToHex(val);
    if (hex) $('#facecolor-picker').value = hex;
    sendCompose({ facecolor: val || 'none' });
  });
  $('#facecolor-picker').addEventListener('input', debounce(() => {
    $('#facecolor').value = $('#facecolor-picker').value;
    sendCompose({ facecolor: $('#facecolor-picker').value });
  }, 200));

  $('#coord-style').addEventListener('change', () => {
    const style = $('#coord-style').value;
    sendCompose({
      coord_style: style,
      x_format: style === 'sexagesimal' ? 'hh:mm:ss.ss' : 'd.dddddd',
      y_format: style === 'sexagesimal' ? 'dd:mm:ss.ss' : 'd.dddddd',
    });
  });
  $('#minorticks').addEventListener('change', () => sendCompose({ minorticks: $('#minorticks').checked }));
  const tickNum = (id) => parseFloat($('#' + id).value);
  for (const [id, key] of [
    ['tick-major-size', 'tick_major_size'],
    ['tick-minor-size', 'tick_minor_size'],
    ['tick-major-width', 'tick_major_width'],
    ['tick-minor-width', 'tick_minor_width'],
  ]) {
    $('#' + id).addEventListener('change', () => sendCompose({ [key]: tickNum(id) }));
  }
  $('#tick-direction').addEventListener('change', () => sendCompose({ tick_direction: $('#tick-direction').value }));
  $('#plot-title').addEventListener('change', () => sendCompose({ title: $('#plot-title').value }));
  $('#plot-xlabel').addEventListener('change', () => sendCompose({ xlabel: $('#plot-xlabel').value }));
  $('#plot-ylabel').addEventListener('change', () => sendCompose({ ylabel: $('#plot-ylabel').value }));
  $('#bare-plot').addEventListener('change', () => sendCompose({ bare_plot: $('#bare-plot').checked }));
  $('#show-legend').addEventListener('change', async () => {
    await sendCompose({ show_legend: $('#show-legend').checked }, false);
    refreshCombined();
  });
  $('#legend-loc').addEventListener('change', async () => {
    await sendCompose({ legend_loc: $('#legend-loc').value }, false);
    if ($('#show-legend').checked) refreshCombined();
  });
  $('#show-swatch').addEventListener('change', async () => {
    await sendCompose({ show_combo_swatch: $('#show-swatch').checked }, false);
    refreshCombined();
  });
  $('#swatch-loc').addEventListener('change', async () => {
    await sendCompose({ combo_swatch_loc: $('#swatch-loc').value }, false);
    if ($('#show-swatch').checked) refreshCombined();
  });
  $('#swatch-labels').addEventListener('change', async () => {
    await sendCompose({ combo_swatch_labels: $('#swatch-labels').checked }, false);
    if ($('#show-swatch').checked) refreshCombined();
  });
  $('#swatch-label-offset').addEventListener('input', debounce(async () => {
    const n = syncSwatchOffset($('#swatch-label-offset').value);
    await sendCompose({ combo_swatch_label_offset: n }, false);
    if ($('#show-swatch').checked && $('#swatch-labels').checked) refreshCombined();
  }, 150));
  $('#swatch-label-offset-num').addEventListener('change', async () => {
    const n = syncSwatchOffset($('#swatch-label-offset-num').value);
    await sendCompose({ combo_swatch_label_offset: n }, false);
    if ($('#show-swatch').checked && $('#swatch-labels').checked) refreshCombined();
  });
  $('#swatch-inset-scale').addEventListener('input', debounce(async () => {
    const n = syncSwatchInsetScale($('#swatch-inset-scale').value);
    await sendCompose({ combo_swatch_inset_scale: n }, false);
    if ($('#show-swatch').checked) refreshCombined();
  }, 150));
  $('#swatch-inset-scale-num').addEventListener('change', async () => {
    const n = syncSwatchInsetScale($('#swatch-inset-scale-num').value);
    await sendCompose({ combo_swatch_inset_scale: n }, false);
    if ($('#show-swatch').checked) refreshCombined();
  });
  $('#swatch-size').addEventListener('change', async () => {
    const n = syncSwatchSize($('#swatch-size').value);
    await sendCompose({ combo_swatch_size: n }, false);
    if ($('#show-swatch').checked) refreshCombined();
  });
  $('#show-band-labels').addEventListener('change', async () => {
    await sendCompose({ show_band_labels: $('#show-band-labels').checked }, false);
    refreshCombined();
  });
  $('#band-labels-loc').addEventListener('change', async () => {
    await sendCompose({ band_labels_loc: $('#band-labels-loc').value }, false);
    if ($('#show-band-labels').checked) refreshCombined();
  });
  $('#show-compass').addEventListener('change', async () => {
    await sendCompose({ show_compass: $('#show-compass').checked }, false);
    refreshCombined();
  });
  $('#compass-loc').addEventListener('change', async () => {
    await sendCompose({ compass_loc: $('#compass-loc').value }, false);
    if ($('#show-compass').checked) refreshCombined();
  });
  $('#show-beam').addEventListener('change', async () => {
    await sendCompose({ show_beam: $('#show-beam').checked }, false);
    refreshCombined();
  });
  $('#beam-loc').addEventListener('change', async () => {
    await sendCompose({ beam_loc: $('#beam-loc').value }, false);
    if ($('#show-beam').checked) refreshCombined();
  });
  $('#beam-style').addEventListener('change', async () => {
    await sendCompose({ beam_style: $('#beam-style').value }, false);
    if ($('#show-beam').checked) refreshCombined();
  });
  $('#show-scale-bar').addEventListener('change', async () => {
    await sendCompose({ show_scale_bar: $('#show-scale-bar').checked }, false);
    refreshCombined();
  });
  $('#scale-bar-asec').addEventListener('change', async () => {
    await sendCompose({ scale_bar_asec: parseFloat($('#scale-bar-asec').value) || 0 }, false);
    if ($('#show-scale-bar').checked) refreshCombined();
  });
  $('#scale-bar-loc').addEventListener('change', async () => {
    await sendCompose({ scale_bar_loc: parseInt($('#scale-bar-loc').value) || 4 }, false);
    if ($('#show-scale-bar').checked) refreshCombined();
  });
  $('#scale-bar-color').addEventListener('change', async () => {
    await sendCompose({ scale_bar_color: $('#scale-bar-color').value.trim() }, false);
    if ($('#show-scale-bar').checked) refreshCombined();
  });
  $('#scale-bar-stroke').addEventListener('change', async () => {
    await sendCompose({ scale_bar_stroke_color: $('#scale-bar-stroke').value.trim() }, false);
    if ($('#show-scale-bar').checked) refreshCombined();
  });
  $('#scale-bar-stroke-lw').addEventListener('change', async () => {
    await sendCompose({ scale_bar_stroke_lw: parseFloat($('#scale-bar-stroke-lw').value) || 0 }, false);
    if ($('#show-scale-bar').checked) refreshCombined();
  });

  $('#combine-mode').addEventListener('change', () => {
    const mode = $('#combine-mode').value;
    updateBlendEnabled(mode);
    const payload = { combine_mode: mode };
    if (subtractiveModesDefaultWhite(mode) && $('#combine-bg').value === 'black') {
      $('#combine-bg').value = 'white';
      updateBgUI();
      payload.combine_background = 'white';
    }
    sendCompose(payload);
  });
  $('#combine-blend').addEventListener('change', () => sendCompose({ combine_blend: $('#combine-blend').value }));

  $('#combine-bg').addEventListener('change', async () => {
    updateBgUI();
    await sendCompose({ combine_background: currentBackground() }, false);
    refreshCombined();
  });
  $('#combine-bg-color').addEventListener('change', () => {
    const hex = cssColorToHex($('#combine-bg-color').value.trim());
    if (hex) $('#combine-bg-picker').value = hex;
    sendCompose({ combine_background: currentBackground() });
  });
  $('#combine-bg-picker').addEventListener('input', debounce(() => {
    $('#combine-bg-color').value = $('#combine-bg-picker').value;
    sendCompose({ combine_background: currentBackground() });
  }, 200));

  $('#btn-save-image').addEventListener('click', () => openSaveDialog('image'));
  $('#btn-save-fits').addEventListener('click', () => openSaveDialog('fits'));
  $('#btn-export-cutout').addEventListener('click', () => openCutoutDialog());

  $('#btn-save-session').addEventListener('click', () => { saveSession(); });

  $('#btn-load-session').addEventListener('click', () => {
    $('#session-file-input').click();
  });

  $('#btn-reset-session').addEventListener('click', confirmResetSession);

  $('#session-file-input').addEventListener('change', async (ev) => {
    const file = ev.target.files && ev.target.files[0];
    ev.target.value = '';
    if (!file) return;
    try {
      const text = await file.text();
      const st = JSON.parse(text);
      const r = await apiPost('/api/session/load', { state: st });
      await applyServerState(r.state);
      if (r.warnings && r.warnings.length) {
        let msg = 'Session loaded with warnings: ' + r.warnings.join('; ');
        if (r.warnings.some(w => w.includes('file not found'))) {
          msg += ' Missing FITS files must be in the same folder as the first panel '
            + 'or the session JSON, or re-loaded with Browse.';
        }
        setStatus(msg, true);
      } else {
        setStatus('Session loaded from ' + file.name);
      }
    } catch (err) { setStatus(err.message, true); }
  });

  wireCursorReadout();

  $('#btn-params').addEventListener('click', async () => {
    const pre = document.createElement('pre');
    pre.textContent = await api('/api/params');
    showModal('Current plot parameters', pre, [
      { label: 'Copy', onclick: () => { navigator.clipboard.writeText(pre.textContent); return true; } },
      { label: 'Close', primary: true },
    ]);
  });

  $('#btn-export').addEventListener('click', async () => {
    const script = await api('/api/export_script');
    const pre = document.createElement('pre');
    pre.textContent = script;
    showModal('Standalone script for the current state', pre, [
      { label: 'Copy', onclick: () => { navigator.clipboard.writeText(script); return true; } },
      {
        label: 'Download .py',
        onclick: () => {
          const a = document.createElement('a');
          a.href = URL.createObjectURL(new Blob([script], { type: 'text/x-python' }));
          a.download = 'multicolorfits_recreate.py';
          a.click();
          return true;
        },
      },
      { label: 'Close', primary: true },
    ]);
  });
}

async function openSaveDialog(kind) {
  const isImage = kind === 'image';
  const form = document.createElement('div');
  form.className = 'save-form';

  // FITS: let the user pick which product to save (combined / per-layer).
  let targetSelect = null;
  if (!isImage) {
    let targets = [];
    try { targets = (await api('/api/fits_targets')).targets || []; } catch (e) { targets = []; }
    const tLabel = document.createElement('label');
    tLabel.textContent = 'What to save:';
    targetSelect = document.createElement('select');
    targets.forEach((t) => {
      const opt = document.createElement('option');
      opt.value = `${t.kind}:${t.index}`;
      opt.textContent = t.label;
      targetSelect.appendChild(opt);
    });
    form.appendChild(tLabel);
    form.appendChild(targetSelect);
  }

  const pathInput = document.createElement('input');
  pathInput.type = 'text';
  pathInput.placeholder = isImage ? '~/my_multicolor_image.png  (.png .jpg .pdf .eps)' : '~/my_multicolor.fits';
  const pathLabel = document.createElement('label');
  pathLabel.textContent = 'Save to path (on the machine running the server):';
  form.appendChild(pathLabel);
  form.appendChild(pathInput);
  let dpiInput = null;
  if (isImage) {
    const dpiLabel = document.createElement('label');
    dpiLabel.textContent = 'DPI:';
    dpiInput = document.createElement('input');
    dpiInput.type = 'number';
    dpiInput.value = 300;
    form.appendChild(dpiLabel);
    form.appendChild(dpiInput);
  }
  showModal(isImage ? 'Save combined image' : 'Save FITS', form, [
    {
      label: 'Save', primary: true,
      onclick: async () => {
        const path = pathInput.value.trim();
        if (!path) return true;
        try {
          let payload;
          if (isImage) {
            payload = { path, dpi: parseInt(dpiInput.value) || 300 };
          } else {
            const [tkind, tindex] = (targetSelect ? targetSelect.value : 'combined:-1').split(':');
            payload = { path, kind: tkind, index: parseInt(tindex) };
          }
          const res = await apiPost(`/api/save/${isImage ? 'image' : 'fits'}`, payload);
          setStatus('Saved: ' + res.path);
        } catch (err) { setStatus('Save failed: ' + err.message, true); return true; }
      },
    },
    { label: 'Cancel' },
  ]);
  pathInput.focus();
}

function openCutoutDialog() {
  const form = document.createElement('div');
  form.className = 'save-form cutout-form';

  function field(labelText, el) {
    const lab = document.createElement('label');
    lab.textContent = labelText;
    form.appendChild(lab);
    form.appendChild(el);
    return el;
  }
  function num(val, step, min, max) {
    const inp = document.createElement('input');
    inp.type = 'number';
    inp.value = val;
    if (step != null) inp.step = step;
    if (min != null) inp.min = min;
    if (max != null) inp.max = max;
    return inp;
  }
  function sel(opts, value) {
    const s = document.createElement('select');
    opts.forEach(([v, t]) => {
      const o = document.createElement('option');
      o.value = v;
      o.textContent = t;
      if (v === value) o.selected = true;
      s.appendChild(o);
    });
    return s;
  }

  const sizeInput = field('Max size (px; blank = native)', num('', 1, 16, 8192));
  sizeInput.placeholder = 'e.g. 400';
  const loInput = field('Alpha low %', num(55, 0.1, 0, 100));
  const hiInput = field('Alpha high %', num(99.3, 0.1, 0, 100));
  const gammaInput = field('Alpha gamma (<1 = bolder)', num(0.5, 0.05, 0.05, 4));
  const srcInput = field('Alpha from', sel([
    ['luma', 'Luminance'],
    ['max', 'Max RGB'],
    ['mean', 'Mean RGB'],
  ], 'luma'));
  const cropInput = field('Crop', sel([
    ['auto', 'Auto (tight + pad)'],
    ['none', 'Full frame'],
  ], 'auto'));
  const padInput = field('Crop pad (fraction)', num(0.05, 0.01, 0, 0.5));
  const matteInput = field('Matte', sel([
    ['none', 'None (luminance only)'],
    ['circle', 'Circle'],
    ['ellipse', 'Ellipse'],
  ], 'none'));
  const softInput = field('Soft edge σ (px)', num(0, 0.5, 0, 50));
  const invertLab = document.createElement('label');
  invertLab.className = 'chk';
  const invertInput = document.createElement('input');
  invertInput.type = 'checkbox';
  invertLab.appendChild(invertInput);
  invertLab.appendChild(document.createTextNode(' Invert (dark = signal)'));
  form.appendChild(invertLab);

  const pathInput = field('Save to server path (optional)', document.createElement('input'));
  pathInput.type = 'text';
  pathInput.placeholder = '~/galaxy_stamp.png';

  const preview = document.createElement('img');
  preview.alt = '';
  preview.style.maxWidth = '100%';
  preview.style.maxHeight = '220px';
  preview.style.display = 'none';
  preview.style.marginTop = '8px';
  preview.style.background =
    'repeating-conic-gradient(#bbb 0% 25%, #888 0% 50%) 50% / 16px 16px';
  form.appendChild(preview);

  function payloadFromForm(extra) {
    const sizeVal = sizeInput.value.trim();
    return Object.assign({
      size: sizeVal === '' ? null : parseInt(sizeVal, 10),
      alpha_lo: parseFloat(loInput.value),
      alpha_hi: parseFloat(hiInput.value),
      alpha_gamma: parseFloat(gammaInput.value),
      alpha_source: srcInput.value,
      crop: cropInput.value,
      pad: parseFloat(padInput.value),
      matte: matteInput.value,
      soft_edge: parseFloat(softInput.value) || 0,
      invert: !!invertInput.checked,
      path: pathInput.value.trim(),
    }, extra || {});
  }

  showModal('Export transparent cutout', form, [
    {
      label: 'Preview',
      onclick: async () => {
        try {
          const resp = await fetch('/api/export/cutout', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payloadFromForm({ download: true })),
          });
          if (!resp.ok) {
            let detail = resp.statusText;
            try { detail = (await resp.json()).detail || detail; } catch (e) { /* */ }
            throw new Error(detail);
          }
          const blob = await resp.blob();
          if (preview.dataset.url) URL.revokeObjectURL(preview.dataset.url);
          const url = URL.createObjectURL(blob);
          preview.dataset.url = url;
          preview.src = url;
          preview.style.display = 'block';
        } catch (err) {
          setStatus('Cutout preview failed: ' + err.message, true);
        }
        return true; // keep modal open
      },
    },
    {
      label: 'Download PNG',
      primary: true,
      onclick: async () => {
        try {
          const resp = await fetch('/api/export/cutout', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payloadFromForm({ download: true })),
          });
          if (!resp.ok) {
            let detail = resp.statusText;
            try { detail = (await resp.json()).detail || detail; } catch (e) { /* */ }
            throw new Error(detail);
          }
          const blob = await resp.blob();
          triggerDownload('multicolorfits_cutout.png', blob, 'image/png');
          setStatus('Downloaded transparent cutout PNG');
        } catch (err) {
          setStatus('Cutout download failed: ' + err.message, true);
          return true;
        }
      },
    },
    {
      label: 'Save to path',
      onclick: async () => {
        const path = pathInput.value.trim();
        if (!path) {
          setStatus('Enter a server path, or use Download PNG', true);
          return true;
        }
        try {
          const res = await apiPost('/api/export/cutout', payloadFromForm({ download: false }));
          setStatus('Saved cutout: ' + res.path);
        } catch (err) {
          setStatus('Cutout save failed: ' + err.message, true);
          return true;
        }
      },
    },
    { label: 'Cancel' },
  ]);
}

// Map mouse position on the displayed combined image to data pixel coords.
function imgPixelFromEvent(img, evt) {
  const rect = img.getBoundingClientRect();
  const nw = img.naturalWidth || img.width;
  const nh = img.naturalHeight || img.height;
  if (!nw || !nh) return null;
  const scale = Math.min(rect.width / nw, rect.height / nh);
  const dispW = nw * scale, dispH = nh * scale;
  const offX = rect.left + (rect.width - dispW) / 2;
  const offY = rect.top + (rect.height - dispH) / 2;
  const fx = (evt.clientX - offX) / dispW;
  const fy = (evt.clientY - offY) / dispH;
  if (fx < 0 || fx > 1 || fy < 0 || fy > 1) return null;
  return {
    ix: Math.round(fx * (nw - 1)),
    iy: Math.round(fy * (nh - 1)),
    img_w: nw,
    img_h: nh,
  };
}

const cursorReadoutDebounced = debounce(async (ix, iy, img_w, img_h) => {
  const token = ++cursorReadoutToken;
  try {
    const info = await api(`/api/cursor?x=${ix}&y=${iy}&img_w=${img_w}&img_h=${img_h}`);
    if (token !== cursorReadoutToken) return;
    const el = $('#cursor-readout');
    if (!info.layers || !info.layers.length) {
      el.textContent = '';
      delete el.dataset.baseReadout;
      return;
    }
    let text = formatCursorReadout(info);
    if (info.sky) text += `  RA=${info.sky.ra.toFixed(5)}, Dec=${info.sky.dec.toFixed(5)}`;
    el.textContent = text;
    el.dataset.baseReadout = formatCursorReadout(info);
  } catch (e) { /* ignore transient errors while moving */ }
}, 120);

let cursorReadoutToken = 0;
let cursorRaf = null;
let cursorAbort = null;
let cursorPending = null;

function formatCursorReadout(info) {
  if (!info.layers || !info.layers.length) return '';
  const vals = info.layers.map(L => `${L.label}: ${L.value.toPrecision(5)}`).join('  ');
  return `x,y=${info.x},${info.y}   ${vals}`;
}

async function fetchCursorSky(ix, iy, img_w, img_h, token) {
  if (cursorAbort) cursorAbort.abort();
  cursorAbort = new AbortController();
  const signal = cursorAbort.signal;
  try {
    const resp = await fetch(`/api/cursor?x=${ix}&y=${iy}&img_w=${img_w}&img_h=${img_h}`, { signal });
    if (!resp.ok) return;
    const info = await resp.json();
    if (token !== cursorReadoutToken) return;
    const el = $('#cursor-readout');
    const base = el.dataset.baseReadout || '';
    if (!info.sky) {
      el.textContent = base;
      return;
    }
    el.textContent = base + `  RA=${info.sky.ra.toFixed(5)}, Dec=${info.sky.dec.toFixed(5)}`;
  } catch (e) {
    if (e.name === 'AbortError') return;
  }
}

function updateCursorReadout(p) {
  const el = $('#cursor-readout');
  state.panels.forEach((panel, i) => {
    if (panel && panel.loaded) PreviewClient.ensurePanelBuffer(i);
  });
  const shape = PreviewClient.dataShape(state.panels);
  if (!shape) {
    cursorReadoutDebounced(p.ix, p.iy, p.img_w, p.img_h);
    return;
  }
  const { x, y } = PreviewClient.displayToDataPixel(p.ix, p.iy, p.img_w, p.img_h, shape);
  const info = PreviewClient.sampleAtDataPixel(x, y, state.panels);
  if (!info.layers.length) {
    cursorReadoutDebounced(p.ix, p.iy, p.img_w, p.img_h);
    return;
  }
  const token = ++cursorReadoutToken;
  const text = formatCursorReadout(info);
  el.dataset.baseReadout = text;
  el.textContent = text;
  fetchCursorSky(p.ix, p.iy, p.img_w, p.img_h, token);
}

function wireCursorReadout() {
  const el = $('#cursor-readout');
  function onMove(evt) {
    const target = evt.currentTarget;
    if (target.hidden) return;
    const p = imgPixelFromEvent(target, evt);
    if (!p) { el.textContent = ''; delete el.dataset.baseReadout; return; }
    cursorPending = p;
    if (cursorRaf != null) return;
    cursorRaf = requestAnimationFrame(() => {
      cursorRaf = null;
      const next = cursorPending;
      cursorPending = null;
      if (!next) return;
      updateCursorReadout(next);
    });
  }
  function onLeave() {
    cursorPending = null;
    if (cursorAbort) cursorAbort.abort();
    el.textContent = '';
    delete el.dataset.baseReadout;
  }
  $('#combined-img').addEventListener('mousemove', onMove);
  $('#combined-img').addEventListener('mouseleave', onLeave);
  $('#combined-canvas').addEventListener('mousemove', onMove);
  $('#combined-canvas').addEventListener('mouseleave', onLeave);
}

async function applyServerState(s, opts) {
  opts = opts || {};
  state.compose = s.compose;
  if (s.max_panels != null) state.maxPanels = s.max_panels;
  if (s.min_panels != null) state.minPanels = s.min_panels;
  state.panels = Array.isArray(s.panels) ? s.panels.slice() : [];
  syncPanelDomToState(opts.openPanelIdx);
  $('#gamma').value = s.compose.gamma;
  $('#inverse').checked = s.compose.inverse;
  $('#tickcolor').value = s.compose.tickcolor;
  $('#tickcolor-picker').value = colorToHexInput(s.compose.tickcolor);
  if (s.compose.facecolor != null) {
    $('#facecolor').value = s.compose.facecolor;
    const fhex = cssColorToHex(s.compose.facecolor);
    if (fhex) $('#facecolor-picker').value = fhex;
  }
  $('#minorticks').checked = s.compose.minorticks;
  if (s.compose.tick_major_size != null) $('#tick-major-size').value = s.compose.tick_major_size;
  if (s.compose.tick_minor_size != null) $('#tick-minor-size').value = s.compose.tick_minor_size;
  if (s.compose.tick_major_width != null) $('#tick-major-width').value = s.compose.tick_major_width;
  if (s.compose.tick_minor_width != null) $('#tick-minor-width').value = s.compose.tick_minor_width;
  if (s.compose.tick_direction) $('#tick-direction').value = s.compose.tick_direction;
  if (s.compose.show_legend !== undefined) $('#show-legend').checked = s.compose.show_legend;
  if (s.compose.legend_loc) $('#legend-loc').value = s.compose.legend_loc;
  if (s.compose.show_combo_swatch !== undefined) $('#show-swatch').checked = s.compose.show_combo_swatch;
  if (s.compose.combo_swatch_loc) $('#swatch-loc').value = s.compose.combo_swatch_loc;
  if (s.compose.combo_swatch_labels !== undefined) $('#swatch-labels').checked = s.compose.combo_swatch_labels;
  if (s.compose.combo_swatch_label_offset != null) syncSwatchOffset(s.compose.combo_swatch_label_offset);
  if (s.compose.combo_swatch_inset_scale != null) syncSwatchInsetScale(s.compose.combo_swatch_inset_scale);
  if (s.compose.combo_swatch_size != null) syncSwatchSize(s.compose.combo_swatch_size);
  if (s.compose.show_band_labels !== undefined) $('#show-band-labels').checked = s.compose.show_band_labels;
  if (s.compose.band_labels_loc) $('#band-labels-loc').value = s.compose.band_labels_loc;
  if (s.compose.show_compass !== undefined) $('#show-compass').checked = s.compose.show_compass;
  if (s.compose.compass_loc) $('#compass-loc').value = s.compose.compass_loc;
  if (s.compose.show_beam !== undefined) $('#show-beam').checked = s.compose.show_beam;
  if (s.compose.beam_loc) $('#beam-loc').value = s.compose.beam_loc;
  if (s.compose.beam_style) $('#beam-style').value = s.compose.beam_style;
  if (s.compose.show_scale_bar !== undefined) $('#show-scale-bar').checked = s.compose.show_scale_bar;
  if (s.compose.scale_bar_asec != null) $('#scale-bar-asec').value = s.compose.scale_bar_asec;
  if (s.compose.scale_bar_loc != null) $('#scale-bar-loc').value = String(s.compose.scale_bar_loc);
  if (s.compose.scale_bar_color != null) $('#scale-bar-color').value = s.compose.scale_bar_color;
  if (s.compose.scale_bar_stroke_color != null) $('#scale-bar-stroke').value = s.compose.scale_bar_stroke_color;
  if (s.compose.scale_bar_stroke_lw != null) $('#scale-bar-stroke-lw').value = s.compose.scale_bar_stroke_lw;
  $('#coord-style').value = s.compose.coord_style;
  if (s.compose.combine_mode) $('#combine-mode').value = s.compose.combine_mode;
  if (s.compose.combine_blend) $('#combine-blend').value = s.compose.combine_blend;
  updateBlendEnabled(s.compose.combine_mode || 'rgb');
  if (s.compose.combine_background != null) {
    const bg = String(s.compose.combine_background).toLowerCase();
    if (['black', 'white', 'transparent', 'none'].includes(bg)) {
      $('#combine-bg').value = (bg === 'none') ? 'transparent' : bg;
    } else {
      $('#combine-bg').value = 'custom';
      $('#combine-bg-color').value = s.compose.combine_background;
      const bhex = cssColorToHex(s.compose.combine_background);
      if (bhex) $('#combine-bg-picker').value = bhex;
    }
    updateBgUI();
  }
  if (s.compose.title != null) $('#plot-title').value = s.compose.title;
  if (s.compose.xlabel != null) $('#plot-xlabel').value = s.compose.xlabel;
  if (s.compose.ylabel != null) $('#plot-ylabel').value = s.compose.ylabel;
  if (s.compose.bare_plot !== undefined) $('#bare-plot').checked = s.compose.bare_plot;
  if (s.compose.panel_preview_hd !== undefined) $('#panel-preview-hd').checked = s.compose.panel_preview_hd;
  forEachPanelIndex((i) => {
    const p = state.panels[i];
    if (!p) return;
    updatePanelUI(i, p);
    if (p.loaded) PreviewClient.invalidatePanel(i);
  });
  for (let i = 0; i < state.panels.length; i++) {
    const p = state.panels[i];
    if (!p || !p.loaded) continue;
    await PreviewClient.ensurePanelBuffer(i);
    // Re-sync UI immediately before paint so absolute vmin/vmax are in the DOM
    // (client preview prefers those over recomputed buffer percentiles).
    updatePanelUI(i, p);
    refreshPanelPreview(i);
  }
  checkGridStatus();
  refreshCvd();
  // Show a fast combined preview as soon as a session has loaded panels.
  if (state.panels.some(p => p && p.loaded)) {
    refreshCombined({ force: true });
  }
}

function wireTheme() {
  const btn = $('#theme-toggle');
  if (!btn) return;
  btn.addEventListener('click', () => {
    const cur = document.documentElement.getAttribute('data-theme') === 'light' ? 'light' : 'dark';
    const next = cur === 'light' ? 'dark' : 'light';
    document.documentElement.setAttribute('data-theme', next);
    try { localStorage.setItem('mcf-theme', next); } catch (e) { /* ignore */ }
  });
}

async function init() {
  wireTheme();
  wireComposeTabs();
  buildPanels(DEFAULT_N_PANELS);
  try {
    const pref = localStorage.getItem('mcf-fast-preview');
    if (pref === '0' || pref === '1') $('#live-preview').checked = pref === '1';
  } catch (e) { /* ignore */ }
  const addBtn = $('#btn-add-panel');
  if (addBtn) addBtn.addEventListener('click', addPanelSlot);
  wireCombined();
  wirePalettes();
  try {
    const s = await api('/api/state');
    await applyServerState(s);
  } catch (err) {
    setStatus('Could not reach backend: ' + err.message, true);
  }
}

init();
