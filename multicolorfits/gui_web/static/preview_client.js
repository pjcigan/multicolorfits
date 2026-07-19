/* Fast in-browser previews: WebGL2 GPU path + CPU fallback.
 *
 * Loads a downsampled float32 buffer once (binary API), uploads to GPU,
 * then stretch / colorize / combine run on the GPU while dragging levels.
 * Plot Full Resolution still uses the full server matplotlib render.
 */
'use strict';

const PreviewClient = (function () {
  const buffers = {};   // idx -> { data, ny, nx, sorted, dataTex }
  const STRETCHES = new Set(['linear', 'sqrt', 'squared', 'log', 'power', 'sinh', 'asinh']);
  const CLIENT_MODES = new Set(['rgb', 'lab']);
  // Match /api/combined.png?preview=true&max_size=… so RGB+black (client) and
  // Lab/white/etc (server) previews land at the same displayed size.
  const COMBINED_PREVIEW_MAX = 512;
  let useGpu = null;
  let rafPending = null;
  let pendingInteractive = null;
  let levelDragging = false;

  function gpuReady() {
    if (useGpu === null) useGpu = typeof PreviewGL !== 'undefined' && PreviewGL.available();
    return useGpu;
  }

  function parseHex(hex) {
    const h = String(hex || '#ffffff').replace('#', '');
    return [parseInt(h.slice(0, 2), 16) / 255,
            parseInt(h.slice(2, 4), 16) / 255,
            parseInt(h.slice(4, 6), 16) / 255];
  }

  function normalizeHexColor(c) {
    if (c == null) return '#ffffff';
    let s = String(c).trim();
    if (!s) return '#ffffff';
    if (!s.startsWith('#')) s = '#' + s;
    if (/^#[0-9a-fA-F]{6}$/.test(s)) return s.toLowerCase();
    const m3 = s.match(/^#([0-9a-fA-F])([0-9a-fA-F])([0-9a-fA-F])$/);
    if (m3) {
      return ('#' + m3[1] + m3[1] + m3[2] + m3[2] + m3[3] + m3[3]).toLowerCase();
    }
    return '#ffffff';
  }

  function panelColorFromDOM(panelEl) {
    const idx = parseInt(panelEl.dataset.idx, 10);
    const serverColor = (typeof state !== 'undefined' && state.panels && state.panels[idx])
      ? normalizeHexColor(state.panels[idx].color) : '';
    const q = sel => panelEl.querySelector(sel);
    const hexEl = q('.color-hex');
    const hexVal = hexEl && hexEl.value.trim() ? normalizeHexColor(hexEl.value) : '';
    const pickerEl = q('.color-picker');
    const pickerVal = pickerEl && pickerEl.value ? normalizeHexColor(pickerEl.value) : '';
    if (serverColor && hexVal === '#ffffff' && pickerVal === '#ffffff' && serverColor !== '#ffffff') {
      return serverColor;
    }
    if (hexVal) return hexVal;
    if (pickerVal) return pickerVal;
    return serverColor || '#ffffff';
  }

  function percentileSorted(sorted, p) {
    if (!sorted.length) return 0;
    const idx = (p / 100) * (sorted.length - 1);
    const lo = Math.floor(idx);
    const hi = Math.ceil(idx);
    if (lo === hi) return sorted[lo];
    const f = idx - lo;
    return sorted[lo] * (1 - f) + sorted[hi] * f;
  }

  function normalizeLimits(vmin, vmax) {
    if (vmin > vmax) { const t = vmin; vmin = vmax; vmax = t; }
    if (!(vmax > vmin)) {
      const eps = Math.max(Math.abs(vmin), 1e-12) * 1e-6 + 1e-12;
      vmax = vmin + eps;
    }
    return { vmin, vmax };
  }

  /** Percentile limits on the preview buffer (for slider readout only). */
  function percentileLimits(buf, percentMin, percentMax) {
    let vmin = percentileSorted(buf.sorted, percentMin);
    let vmax = percentileSorted(buf.sorted, percentMax);
    return normalizeLimits(vmin, vmax);
  }

  /**
   * Display limits for rendering. Prefer absolute vmin/vmax (same as the
   * server PanelState) so session load matches Plot Full Resolution. Fall back to
   * percentiles on the downsampled buffer only when absolute limits are missing.
   */
  function limitsFromParams(buf, params) {
    const avmin = +params.vmin;
    const avmax = +params.vmax;
    if (Number.isFinite(avmin) && Number.isFinite(avmax)) {
      return normalizeLimits(avmin, avmax);
    }
    return percentileLimits(buf, params.percent_min, params.percent_max);
  }

  function composeGamma(compose) {
    const g = compose && compose.gamma;
    return (g != null && Number.isFinite(+g) && +g > 0) ? +g : 2.2;
  }

  function limitsForPanel(panelEl) {
    const idx = parseInt(panelEl.dataset.idx, 10);
    const buf = buffers[idx];
    if (!buf) return null;
    return limitsFromParams(buf, panelParamsFromDOM(panelEl));
  }

  /** Percentile-derived limits for updating the vmin/vmax readout while dragging. */
  function percentileLimitsForPanel(panelEl) {
    const idx = parseInt(panelEl.dataset.idx, 10);
    const buf = buffers[idx];
    if (!buf) return null;
    const params = panelParamsFromDOM(panelEl);
    return percentileLimits(buf, params.percent_min, params.percent_max);
  }

  function panelParamsFromDOM(panelEl) {
    const q = sel => panelEl.querySelector(sel);
    const idx = parseInt(panelEl.dataset.idx, 10);
    const server = (typeof state !== 'undefined' && state.panels && state.panels[idx]) || {};
    const vminEl = q('.vmin');
    const vmaxEl = q('.vmax');
    let vmin = vminEl ? parseFloat(vminEl.value) : NaN;
    let vmax = vmaxEl ? parseFloat(vmaxEl.value) : NaN;
    // Session / server state is authoritative when the number inputs are empty
    // (e.g. mid-sync on first load before updatePanelUI finishes).
    if (!Number.isFinite(vmin) && Number.isFinite(+server.vmin)) vmin = +server.vmin;
    if (!Number.isFinite(vmax) && Number.isFinite(+server.vmax)) vmax = +server.vmax;
    const stretch = (q('.stretch') && q('.stretch').value) || server.stretch || 'linear';
    let percent_min = q('.pmin') ? parseFloat(q('.pmin').value) : NaN;
    let percent_max = q('.pmax') ? parseFloat(q('.pmax').value) : NaN;
    if (!Number.isFinite(percent_min) && Number.isFinite(+server.percent_min)) {
      percent_min = +server.percent_min;
    }
    if (!Number.isFinite(percent_max) && Number.isFinite(+server.percent_max)) {
      percent_max = +server.percent_max;
    }
    if (!Number.isFinite(percent_min)) percent_min = 0;
    if (!Number.isFinite(percent_max)) percent_max = 100;
    return {
      stretch,
      percent_min,
      percent_max,
      vmin,
      vmax,
      color: panelColorFromDOM(panelEl),
      smooth: q('.smooth') ? q('.smooth').checked : !!server.smooth,
    };
  }

  function ensureGpuTexture(buf) {
    if (!gpuReady()) return false;
    if (buf.dataTex) return true;
    const gl = PreviewGL.offscreen.getContext('webgl2');
    if (!gl) return false;
    buf.dataTex = gl.createTexture();
    PreviewGL.uploadDataTex(buf.dataTex, buf.data, buf.nx, buf.ny);
    return true;
  }

  async function loadPanelBuffer(idx, maxSize) {
    const previewMax = maxSize || 768;
    const url = `/api/panel/${idx}/preview_buffer.bin?max_size=${previewMax}`;
    const resp = await fetch(url);
    if (!resp.ok) throw new Error('preview_buffer failed');
    const raw = new Uint8Array(await resp.arrayBuffer());
    const view = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const ny = view.getUint32(0, true);
    const nx = view.getUint32(4, true);
    const data = new Float32Array(raw.buffer, raw.byteOffset + 8, (raw.byteLength - 8) / 4);
    const sorted = [];
    for (let i = 0; i < data.length; i++) {
      const v = data[i];
      if (Number.isFinite(v)) sorted.push(v);
    }
    sorted.sort((a, b) => a - b);
    const old = buffers[idx];
    if (old && old.dataTex && gpuReady()) {
      const gl = PreviewGL.offscreen.getContext('webgl2');
      gl.deleteTexture(old.dataTex);
    }
    const buf = { data, ny, nx, sorted, dataTex: null };
    buffers[idx] = buf;
    ensureGpuTexture(buf);
    return buf;
  }

  function invalidatePanel(idx) {
    const buf = buffers[idx];
    if (buf && buf.dataTex && gpuReady()) {
      const gl = PreviewGL.offscreen.getContext('webgl2');
      gl.deleteTexture(buf.dataTex);
    }
    delete buffers[idx];
  }

  function hasBuffer(idx) {
    return !!buffers[idx];
  }

  function canRenderPanel(panelEl) {
    const idx = parseInt(panelEl.dataset.idx, 10);
    if (!buffers[idx]) return false;
    const params = panelParamsFromDOM(panelEl);
    if (params.smooth) return false;
    if (!STRETCHES.has(params.stretch)) return false;
    return true;
  }

  function canRenderCombined(activeIndices, panelEls, compose) {
    if (!compose || !CLIENT_MODES.has(compose.combine_mode || 'rgb')) return false;
    const bg = String(compose.combine_background || 'black').toLowerCase();
    if (bg !== 'black' && bg !== '') return false;
    for (const i of activeIndices) {
      const el = panelEls[i];
      if (!el || !canRenderPanel(el)) return false;
    }
    return gpuReady() || true;
  }

  function showPanelCanvas(panelEl) {
    const canvas = panelEl.querySelector('.preview-canvas');
    const img = panelEl.querySelector('.preview');
    if (canvas) canvas.removeAttribute('hidden');
    if (img) img.setAttribute('hidden', '');
    const empty = panelEl.querySelector('.preview-empty');
    if (empty) empty.setAttribute('hidden', '');
  }

  /** CPU-only panel thumbnail (reliable for small preview boxes). */
  function renderPanelThumbnail(panelEl, compose) {
    const idx = parseInt(panelEl.dataset.idx, 10);
    if (!buffers[idx]) return false;
    return renderPanelCpu(panelEl, compose);
  }

  function renderLayersGpu(activeIndices, panelEls, compose, onlyPanelIdx, opts) {
    if (!gpuReady()) return null;
    const blitPanels = !(opts && opts.blitPanels === false);
    const inverse = !!(compose && compose.inverse);
    let nx = 0;
    let ny = 0;
    let slot = 0;
    const slotByPanel = {};
    for (const i of activeIndices) {
      if (onlyPanelIdx !== undefined && onlyPanelIdx !== null && i !== onlyPanelIdx) continue;
      const el = panelEls[i];
      if (!el || !canRenderPanel(el)) continue;
      const buf = buffers[i];
      if (!buf || !ensureGpuTexture(buf)) continue;
      const params = panelParamsFromDOM(el);
      const { vmin, vmax } = limitsFromParams(buf, params);
      const color = parseHex(params.color);
      const gamma = composeGamma(compose);
      nx = buf.nx;
      ny = buf.ny;
      const tex = PreviewGL.renderLayer(buf.dataTex, slot, nx, ny, vmin, vmax,
                                          params.stretch, color, false, gamma);
      if (!tex) continue;
      slotByPanel[i] = slot;
      if (blitPanels) {
        PreviewGL.blitLayerToCanvas(slot, el.querySelector('.preview-canvas'), inverse, gamma);
        showPanelCanvas(el);
      }
      slot += 1;
    }
    if (!slot) return null;
    return { count: slot, nx, ny, slotByPanel };
  }

  function renderPanelGpu(panelEl, compose) {
    try {
      const idx = parseInt(panelEl.dataset.idx, 10);
      const panelElsArr = [];
      panelElsArr[idx] = panelEl;
      return !!renderLayersGpu([idx], panelElsArr, compose, null);
    } catch (err) {
      console.warn('GPU panel preview failed:', err);
      return false;
    }
  }

  function renderCombinedGpu(activeIndices, panelEls, compose) {
    try {
      // Do not overwrite left-panel CPU thumbnails while composing.
      const layers = renderLayersGpu(activeIndices, panelEls, compose, null,
                                     { blitPanels: false });
      if (!layers) return false;
      const mode = compose.combine_mode || 'rgb';
      PreviewGL.combineLayers(layers.count, mode, compose.combine_blend || 'screen',
                              compose.gamma, compose.inverse);
      const canvas = document.getElementById('combined-canvas');
      PreviewGL.blitToCanvas(canvas, COMBINED_PREVIEW_MAX);
      canvas.hidden = false;
      const img = document.getElementById('combined-img');
      if (img) img.hidden = true;
      document.getElementById('combined-placeholder').hidden = true;
      return true;
    } catch (err) {
      console.warn('GPU combined preview failed:', err);
      return false;
    }
  }

  // ---- CPU fallback (no WebGL2) ----

  function stretch01(x, name) {
    switch (name) {
      case 'linear': return x;
      case 'sqrt': return Math.sqrt(x);
      case 'squared': return x * x;
      case 'log': {
        // Match astropy LogStretch(a=1000): log(a*x+1)/log(a+1)
        const a = 1000;
        return Math.log(a * x + 1) / Math.log(a + 1);
      }
      case 'power': {
        // Match astropy PowerDistStretch(a=1000): (a^x - 1)/(a - 1)
        const a = 1000;
        return (Math.pow(a, x) - 1) / (a - 1);
      }
      case 'sinh': {
        // Match astropy SinhStretch(a=1/3): sinh(x/a)/sinh(1/a)
        const a = 1 / 3;
        return Math.sinh(x / a) / Math.sinh(1 / a);
      }
      case 'asinh': {
        // Match astropy AsinhStretch(a=0.1): asinh(x/a)/asinh(1/a) on [0,1]
        const a = 0.1;
        return Math.asinh(x / a) / Math.asinh(1 / a);
      }
      default: return null;
    }
  }

  function scaledBuffer(buf, params) {
    const { vmin, vmax } = limitsFromParams(buf, params);
    const out = new Float32Array(buf.data.length);
    for (let i = 0; i < buf.data.length; i++) {
      const v = buf.data[i];
      if (!Number.isFinite(v)) { out[i] = 0; continue; }
      const span = vmax - vmin;
      let n = span > 0 ? (v - vmin) / span : 0;
      if (n < 0) n = 0; else if (n > 1) n = 1;
      const s = stretch01(n, params.stretch);
      out[i] = s == null ? 0 : s;
    }
    return { scaled: out, ny: buf.ny, nx: buf.nx };
  }

  function drawPanelCanvasCpu(canvas, scaled, color, inverse, gamma) {
    // Match PanelState.render_display: stretch → **gamma → color**gamma → **(1/gamma).
    // Vertically flip so FITS origin (lower-left) matches the server PNG / combined plot.
    const [cr, cg, cb] = parseHex(color);
    const g = composeGamma({ gamma });
    const invG = 1 / g;
    const crG = Math.pow(cr, g);
    const cgG = Math.pow(cg, g);
    const cbG = Math.pow(cb, g);
    const { nx, ny } = scaled;
    const img = canvas.getContext('2d').createImageData(nx, ny);
    const px = img.data;
    const arr = scaled.scaled;
    for (let y = 0; y < ny; y++) {
      const srcRow = y * nx;
      const dstRow = (ny - 1 - y) * nx;
      for (let x = 0; x < nx; x++) {
        let t = Math.max(0, Math.min(1, arr[srcRow + x]));
        t = Math.pow(t, g);
        let r = Math.pow(Math.min(1, Math.max(0, t * crG)), invG);
        let gv = Math.pow(Math.min(1, Math.max(0, t * cgG)), invG);
        let b = Math.pow(Math.min(1, Math.max(0, t * cbG)), invG);
        if (inverse) { r = 1 - r; gv = 1 - gv; b = 1 - b; }
        const j = (dstRow + x) * 4;
        px[j] = (r * 255) | 0;
        px[j + 1] = (gv * 255) | 0;
        px[j + 2] = (b * 255) | 0;
        px[j + 3] = 255;
      }
    }
    if (canvas.width !== nx || canvas.height !== ny) {
      canvas.width = nx;
      canvas.height = ny;
    }
    canvas.getContext('2d').putImageData(img, 0, 0);
  }

  function renderPanelCpu(panelEl, compose) {
    const idx = parseInt(panelEl.dataset.idx, 10);
    const buf = buffers[idx];
    if (!buf) return false;
    const params = panelParamsFromDOM(panelEl);
    const scaled = scaledBuffer(buf, params);
    const canvas = panelEl.querySelector('.preview-canvas');
    if (!canvas) return false;
    drawPanelCanvasCpu(canvas, scaled, params.color, !!(compose && compose.inverse),
                       composeGamma(compose));
    showPanelCanvas(panelEl);
    return true;
  }

  function renderPanel(panelEl, compose) {
    const idx = parseInt(panelEl.dataset.idx, 10);
    if (!buffers[idx]) return false;
    // Prefer CPU for panel thumbnails so gamma matches server render_display
    // (GPU panel path is a shortcut; combined live preview still uses GPU).
    return renderPanelCpu(panelEl, compose);
  }

  function renderCombined(activeIndices, panelEls, compose) {
    if (!canRenderCombined(activeIndices, panelEls, compose)) return false;
    if (gpuReady() && renderCombinedGpu(activeIndices, panelEls, compose)) return true;
    if (compose.combine_mode !== 'rgb') return false;
    // CPU RGB only
    const layers = [];
    for (const i of activeIndices) {
      const el = panelEls[i];
      const params = panelParamsFromDOM(el);
      const scaled = scaledBuffer(buffers[i], params);
      layers.push({ ...scaled, color: params.color });
    }
    const combined = combineRgbLayersCpu(layers, compose.gamma, compose.inverse);
    const canvas = document.getElementById('combined-canvas');
    drawRgbCanvasCpu(canvas, combined, COMBINED_PREVIEW_MAX);
    canvas.hidden = false;
    document.getElementById('combined-img').hidden = true;
    document.getElementById('combined-placeholder').hidden = true;
    return true;
  }

  function combineRgbLayersCpu(layers, gamma, inverse) {
    // Match combine_multicolor inputs: intensity**gamma * color**gamma, then mix + **(1/gamma).
    const n = layers[0].scaled.length;
    const g = gamma || 2.2;
    const sum = new Float32Array(n * 3);
    for (const layer of layers) {
      const [cr, cg, cb] = parseHex(layer.color);
      const crG = Math.pow(cr, g);
      const cgG = Math.pow(cg, g);
      const cbG = Math.pow(cb, g);
      const arr = layer.scaled;
      for (let i = 0; i < n; i++) {
        const t = Math.pow(Math.max(0, Math.min(1, arr[i])), g);
        sum[i * 3] += t * crG;
        sum[i * 3 + 1] += t * cgG;
        sum[i * 3 + 2] += t * cbG;
      }
    }
    let maxR = 0, maxG = 0, maxB = 0, minR = 1, minG = 1, minB = 1;
    for (let i = 0; i < n; i++) {
      const r = sum[i * 3], gv = sum[i * 3 + 1], b = sum[i * 3 + 2];
      if (r > maxR) maxR = r; if (gv > maxG) maxG = gv; if (b > maxB) maxB = b;
      if (r < minR) minR = r; if (gv < minG) minG = gv; if (b < minB) minB = b;
    }
    const maxOfMax = inverse ? Math.max(1 - maxR, 1 - maxG, 1 - maxB) : Math.max(maxR, maxG, maxB);
    const out = new Float32Array(n * 3);
    const invG = 1 / g;
    for (let i = 0; i < n; i++) {
      let r = sum[i * 3], gv = sum[i * 3 + 1], b = sum[i * 3 + 2];
      if (maxOfMax > 0) {
        const dr = maxR - minR, dg = maxG - minG, db = maxB - minB;
        if (dr > 0) r = (r - minR) * (maxR / maxOfMax) / dr;
        if (dg > 0) gv = (gv - minG) * (maxG / maxOfMax) / dg;
        if (db > 0) b = (b - minB) * (maxB / maxOfMax) / db;
      }
      r = Math.pow(Math.min(1, Math.max(0, r)), invG);
      gv = Math.pow(Math.min(1, Math.max(0, gv)), invG);
      b = Math.pow(Math.min(1, Math.max(0, b)), invG);
      if (inverse) { r = 1 - r; gv = 1 - gv; b = 1 - b; }
      out[i * 3] = r; out[i * 3 + 1] = gv; out[i * 3 + 2] = b;
    }
    return { rgb: out, ny: layers[0].ny, nx: layers[0].nx };
  }

  function drawRgbCanvasCpu(canvas, combined, maxSide) {
    // Flip vertically to match server PNG / matplotlib (FITS origin lower-left).
    // Optionally downsample so client RGB previews match server max_size=512.
    const { ny, nx, rgb } = combined;
    let dw = nx;
    let dh = ny;
    if (maxSide && maxSide > 0) {
      const longest = Math.max(dw, dh);
      if (longest > maxSide) {
        const s = maxSide / longest;
        dw = Math.max(1, Math.round(dw * s));
        dh = Math.max(1, Math.round(dh * s));
      }
    }
    const img = canvas.getContext('2d').createImageData(dw, dh);
    const px = img.data;
    for (let y = 0; y < dh; y++) {
      // Display row 0 = top of sky = last FITS row
      const srcY = Math.min(ny - 1, Math.floor((y + 0.5) * ny / dh));
      const srcRow = (ny - 1 - srcY) * nx;
      for (let x = 0; x < dw; x++) {
        const srcX = Math.min(nx - 1, Math.floor((x + 0.5) * nx / dw));
        const i = srcRow + srcX;
        const j = (y * dw + x) * 4;
        px[j] = (rgb[i * 3] * 255) | 0;
        px[j + 1] = (rgb[i * 3 + 1] * 255) | 0;
        px[j + 2] = (rgb[i * 3 + 2] * 255) | 0;
        px[j + 3] = 255;
      }
    }
    if (canvas.width !== dw || canvas.height !== dh) {
      canvas.width = dw;
      canvas.height = dh;
    }
    canvas.getContext('2d').putImageData(img, 0, 0);
  }

  function dataShape(panels) {
    for (const p of panels) {
      if (p && p.loaded && p.shape && p.shape.length === 2) return p.shape;
    }
    return null;
  }

  function displayToDataPixel(ix, iy, imgW, imgH, shape) {
    const ny = shape[0];
    const nx = shape[1];
    const fx = ix / Math.max(imgW - 1, 1);
    const fy = iy / Math.max(imgH - 1, 1);
    return {
      x: Math.round(fx * (nx - 1)),
      y: Math.round((1.0 - fy) * (ny - 1)),
    };
  }

  function sampleAtDataPixel(x, y, panels) {
    const layers = [];
    for (let i = 0; i < panels.length; i++) {
      const p = panels[i];
      if (!p || !p.loaded || !p.shape) continue;
      const ny = p.shape[0];
      const nx = p.shape[1];
      if (x < 0 || y < 0 || x >= nx || y >= ny) continue;
      const buf = buffers[i];
      let value = null;
      if (buf && buf.nx > 0 && buf.ny > 0) {
        const bx = Math.min(buf.nx - 1, Math.round(x * (buf.nx - 1) / Math.max(nx - 1, 1)));
        const by = Math.min(buf.ny - 1, Math.round(y * (buf.ny - 1) / Math.max(ny - 1, 1)));
        value = buf.data[by * buf.nx + bx];
      }
      if (value == null || !Number.isFinite(value)) continue;
      layers.push({
        index: i,
        label: p.label || (`Image ${i + 1}`),
        color: p.color,
        value,
      });
    }
    return { x, y, layers };
  }

  function showServerPanelPreview(panelEl) {
    const canvas = panelEl.querySelector('.preview-canvas');
    const img = panelEl.querySelector('.preview');
    if (canvas) canvas.setAttribute('hidden', '');
    if (img) img.removeAttribute('hidden');
  }

  function showServerCombined() {
    const canvas = document.getElementById('combined-canvas');
    const img = document.getElementById('combined-img');
    if (canvas) canvas.hidden = true;
    if (img) img.hidden = false;
  }

  async function ensurePanelBuffer(idx) {
    if (buffers[idx]) return true;
    try {
      await loadPanelBuffer(idx);
      return true;
    } catch (e) {
      return false;
    }
  }

  function updatePanelLevels(panelIdx, panelEls, activeIndices, compose, liveCombined) {
    const panel = panelEls[panelIdx];
    if (!panel) return false;
    // Panel thumbnail: CPU is reliable in the small preview box; GPU is for combined.
    if (!renderPanelThumbnail(panel, compose)) return false;
    if (liveCombined && activeIndices.length) {
      if (gpuReady() && renderCombinedGpu(activeIndices, panelEls, compose)) return true;
      return renderCombined(activeIndices, panelEls, compose);
    }
    return true;
  }

  function renderInteractiveGpu(activeIndices, panelEls, compose, liveCombined, panelIdx) {
    if (!activeIndices.length) return false;
    if (panelIdx !== undefined && panelIdx !== null) {
      return updatePanelLevels(panelIdx, panelEls, activeIndices, compose, liveCombined);
    }
    if (liveCombined) {
      return renderCombinedGpu(activeIndices, panelEls, compose);
    }
    return !!renderLayersGpu(activeIndices, panelEls, compose, null);
  }

  function renderInteractiveCpu(panelIdx, activeIndices, panelEls, compose, liveCombined) {
    const panel = panelEls[panelIdx];
    if (!renderPanelCpu(panel, compose)) return false;
    for (const i of activeIndices) {
      if (i === panelIdx) continue;
      const el = panelEls[i];
      if (el && canRenderPanel(el)) renderPanelCpu(el, compose);
    }
    if (liveCombined && activeIndices.length) {
      return renderCombined(activeIndices, panelEls, compose);
    }
    return true;
  }

  function runInteractiveRender(spec) {
    const { panelIdx, compose, panelEls, activeIndices, liveCombined } = spec;
    if (!activeIndices.length) return false;
    if (gpuReady()) {
      return renderInteractiveGpu(activeIndices, panelEls, compose, liveCombined, panelIdx);
    }
    return renderInteractiveCpu(panelIdx, activeIndices, panelEls, compose, liveCombined);
  }

  function renderPanelImmediate(panelIdx, compose, opts) {
    const panelEls = opts.panelEls;
    const activeIndices = opts.activeIndices;
    const liveCombined = opts.liveCombined;
    return runInteractiveRender({ panelIdx, compose, panelEls, activeIndices, liveCombined });
  }

  function scheduleInteractive(panelIdx, compose, opts) {
    const job = {
      panelIdx,
      compose,
      panelEls: opts.panelEls,
      activeIndices: opts.activeIndices,
      liveCombined: opts.liveCombined,
    };
    // Draw immediately so panel thumbnails track the slider without waiting for rAF.
    renderPanelImmediate(panelIdx, compose, opts);
    pendingInteractive = job;
    if (rafPending) return;
    rafPending = requestAnimationFrame(() => {
      rafPending = null;
      const next = pendingInteractive;
      pendingInteractive = null;
      if (!next) return;
      if (runInteractiveRender(next)) {
        const mode = gpuReady() ? 'GPU' : 'CPU';
        if (typeof setStatus === 'function') setStatus(`Fast preview (${mode}, instant)`);
      }
    });
  }

  function setLevelDragging(on) {
    levelDragging = on;
  }

  function isLevelDragging() {
    return levelDragging;
  }

  function usesGpu() {
    return gpuReady();
  }

  return {
    loadPanelBuffer,
    ensurePanelBuffer,
    invalidatePanel,
    hasBuffer,
    limitsForPanel,
    percentileLimitsForPanel,
    canRenderPanel,
    renderPanel,
    renderPanelThumbnail,
    sampleAtDataPixel,
    displayToDataPixel,
    dataShape,
    updatePanelLevels,
    canRenderCombined,
    renderCombined,
    showServerPanelPreview,
    showServerCombined,
    scheduleInteractive,
    renderPanelImmediate,
    setLevelDragging,
    isLevelDragging,
    usesGpu,
    STRETCHES,
    COMBINED_PREVIEW_MAX,
    normalizeHexColor,
    panelColorFromDOM,
  };
})();
