/* WebGL2 preview engine: GPU stretch, colorize, RGB/Lab combine. */
'use strict';

const PreviewGL = (function () {
  const STRETCH_ID = {
    linear: 0, sqrt: 1, squared: 2, log: 3, power: 4, sinh: 5, asinh: 6,
  };
  const MAX_LAYERS = 4;

  let gl = null;
  let offscreen = null;
  let programs = null;
  let quadVao = null;
  let layerFbos = [];
  let width = 0;
  let height = 0;

  const FS_BLIT = `#version 300 es
precision highp float;
precision highp sampler2D;
uniform sampler2D u_tex;
uniform float u_gamma;
uniform float u_inverse;
in vec2 v_uv;
out vec4 outColor;
void main() {
  // Match render_display on a colorized layer: **(1/gamma), then optional inverse.
  vec3 c = texture(u_tex, v_uv).rgb;
  float invG = 1.0 / max(u_gamma, 0.01);
  c = pow(clamp(c, 0.0, 1.0), vec3(invG));
  if (u_inverse > 0.5) c = 1.0 - c;
  outColor = vec4(c, 1.0);
}`;

  const VS = `#version 300 es
in vec2 a_pos;
out vec2 v_uv;
void main() {
  v_uv = a_pos * 0.5 + 0.5;
  v_uv.y = 1.0 - v_uv.y;
  gl_Position = vec4(a_pos, 0.0, 1.0);
}`;

  const FS_PANEL = `#version 300 es
precision highp float;
precision highp sampler2D;
uniform sampler2D u_data;
uniform float u_vmin, u_vmax;
uniform int u_stretch;
uniform vec3 u_color;
uniform float u_gamma;
uniform float u_inverse;
in vec2 v_uv;
out vec4 outColor;

float norm01(float v) {
  float span = u_vmax - u_vmin;
  if (span <= 0.0) return 0.0;
  return clamp((v - u_vmin) / span, 0.0, 1.0);
}

float stretch(float x) {
  if (u_stretch == 1) return sqrt(max(x, 0.0));
  if (u_stretch == 2) return x * x;
  if (u_stretch == 3) {
    // LogStretch(a=1000)
    float a = 1000.0;
    return log(a * x + 1.0) / log(a + 1.0);
  }
  if (u_stretch == 4) {
    // PowerDistStretch(a=1000)
    float a = 1000.0;
    return (pow(a, x) - 1.0) / (a - 1.0);
  }
  if (u_stretch == 5) {
    // SinhStretch(a=1/3)
    float a = 1.0 / 3.0;
    return sinh(x / a) / sinh(1.0 / a);
  }
  if (u_stretch == 6) {
    // AsinhStretch(a=0.1) — must match astropy (NOT a centered mapping)
    float a = 0.1;
    return asinh(x / a) / asinh(1.0 / a);
  }
  return x;
}

void main() {
  float v = texture(u_data, v_uv).r;
  if (v != v) v = 0.0;
  float t = stretch(norm01(v));
  // Match render_color_rgb: intensity**gamma * color**gamma (pre-display).
  float g = max(u_gamma, 0.01);
  t = pow(clamp(t, 0.0, 1.0), g);
  vec3 col = pow(clamp(u_color, 0.0, 1.0), vec3(g));
  outColor = vec4(t * col, 1.0);
}`;

  const FS_COMBINE_RGB = `#version 300 es
precision highp float;
precision highp sampler2D;
uniform sampler2D u_layers[4];
uniform int u_count;
uniform float u_gamma;
uniform float u_inverse;
in vec2 v_uv;
out vec4 outColor;

void main() {
  vec3 sum = vec3(0.0);
  vec3 mx = vec3(0.0);
  vec3 mn = vec3(1.0);
  for (int i = 0; i < 4; i++) {
    if (i >= u_count) break;
    vec3 c = texture(u_layers[i], v_uv).rgb;
    sum += c;
    if (i == 0) { mx = c; mn = c; }
    else { mx = max(mx, c); mn = min(mn, c); }
  }
  float maxOfMax = u_inverse > 0.5
    ? max(max(1.0 - mx.r, 1.0 - mx.g), 1.0 - mx.b)
    : max(max(mx.r, mx.g), mx.b);
  vec3 outv = sum;
  if (maxOfMax > 0.0) {
    vec3 denom = mx - mn;
    outv.r = denom.r > 0.0 ? (sum.r - mn.r) * (mx.r / maxOfMax) / denom.r : sum.r;
    outv.g = denom.g > 0.0 ? (sum.g - mn.g) * (mx.g / maxOfMax) / denom.g : sum.g;
    outv.b = denom.b > 0.0 ? (sum.b - mn.b) * (mx.b / maxOfMax) / denom.b : sum.b;
  }
  outv = clamp(outv, 0.0, 1.0);
  float invG = 1.0 / max(u_gamma, 0.01);
  outv = pow(outv, vec3(invG));
  if (u_inverse > 0.5) outv = 1.0 - outv;
  outColor = vec4(clamp(outv, 0.0, 1.0), 1.0);
}`;

  const FS_COMBINE_LAB = `#version 300 es
precision highp float;
precision highp sampler2D;
uniform sampler2D u_layers[4];
uniform int u_count;
uniform int u_blend;
uniform float u_gamma;
uniform float u_inverse;
in vec2 v_uv;
out vec4 outColor;

vec3 srgbToLinear(vec3 c) {
  return mix(c / 12.92, pow((c + 0.055) / 1.055, vec3(2.4)), step(0.04045, c));
}

vec3 linearToSrgb(vec3 c) {
  c = clamp(c, 0.0, 1.0);
  return mix(c * 12.92, 1.055 * pow(c, vec3(1.0 / 2.4)) - 0.055, step(0.0031308, c));
}

vec3 rgb2lab(vec3 rgb) {
  vec3 lin = srgbToLinear(rgb);
  float x = lin.r * 0.4124564 + lin.g * 0.3575761 + lin.b * 0.1804375;
  float y = lin.r * 0.2126729 + lin.g * 0.7151522 + lin.b * 0.0721750;
  float z = lin.r * 0.0193339 + lin.g * 0.1191920 + lin.b * 0.9503041;
  x /= 0.95047; z /= 1.08883;
  vec3 f = vec3(
    x > 0.008856 ? pow(x, 1.0/3.0) : (7.787 * x + 16.0/116.0),
    y > 0.008856 ? pow(y, 1.0/3.0) : (7.787 * y + 16.0/116.0),
    z > 0.008856 ? pow(z, 1.0/3.0) : (7.787 * z + 16.0/116.0)
  );
  return vec3(116.0 * f.y - 16.0, 500.0 * (f.x - f.y), 200.0 * (f.y - f.z));
}

vec3 lab2rgb(vec3 lab) {
  float y = (lab.x + 16.0) / 116.0;
  float x = lab.y / 500.0 + y;
  float z = y - lab.z / 200.0;
  vec3 f = vec3(x, y, z);
  vec3 f3 = f * f * f;
  vec3 xyz = vec3(
    0.95047 * mix((f.x - 16.0/116.0) / 7.787, f3.x, step(0.008856, f3.x)),
    mix((f.y - 16.0/116.0) / 7.787, f3.y, step(0.008856, f3.y)),
    1.08883 * mix((f.z - 16.0/116.0) / 7.787, f3.z, step(0.008856, f3.z))
  );
  vec3 lin = vec3(
    xyz.x *  3.2404542 + xyz.y * -1.5371385 + xyz.z * -0.4985314,
    xyz.x * -0.9692660 + xyz.y *  1.8760108 + xyz.z *  0.0415560,
    xyz.x *  0.0556434 + xyz.y * -0.2040259 + xyz.z *  1.0572252
  );
  return linearToSrgb(clamp(lin, 0.0, 1.0));
}

float blendL(float L) {
  if (u_blend == 1) return clamp(L, 0.0, 1.0);
  if (u_blend == 2) return L;
  if (u_blend == 3) return L;
  return L;
}

void main() {
  float Ls[4];
  float as[4];
  float bs[4];
  int n = u_count;
  for (int i = 0; i < 4; i++) {
    if (i >= n) break;
    vec3 rgb = texture(u_layers[i], v_uv).rgb;
    float invG = 1.0 / max(u_gamma, 0.01);
    rgb = pow(clamp(rgb, 0.0, 1.0), vec3(invG));
    vec3 lab = rgb2lab(rgb);
    Ls[i] = lab.x / 100.0;
    as[i] = lab.y;
    bs[i] = lab.z;
  }
  float Lout = 0.0;
  if (u_blend == 1) {
    for (int i = 0; i < 4; i++) { if (i >= n) break; Lout += Ls[i]; }
    Lout = clamp(Lout, 0.0, 1.0);
  } else if (u_blend == 2) {
  float prod = 1.0;
    for (int i = 0; i < 4; i++) { if (i >= n) break; prod *= (1.0 - Ls[i]); }
    Lout = 1.0 - prod;
  } else if (u_blend == 3) {
    Lout = Ls[0];
    for (int i = 1; i < 4; i++) { if (i >= n) break; Lout = max(Lout, Ls[i]); }
  } else {
    Lout = 0.0;
    for (int i = 0; i < 4; i++) { if (i >= n) break; Lout += Ls[i]; }
    Lout /= float(n);
  }
  float wsum = 0.0;
  float aout = 0.0;
  float bout = 0.0;
  for (int i = 0; i < 4; i++) {
    if (i >= n) break;
    float w = Ls[i];
    wsum += w;
    aout += w * as[i];
    bout += w * bs[i];
  }
  if (wsum > 1e-9) { aout /= wsum; bout /= wsum; }
  if (u_inverse > 0.5) Lout = 1.0 - Lout;
  vec3 rgb = lab2rgb(vec3(Lout * 100.0, aout, bout));
  outColor = vec4(clamp(rgb, 0.0, 1.0), 1.0);
}`;

  function compile(glCtx, type, src) {
    const sh = glCtx.createShader(type);
    glCtx.shaderSource(sh, src);
    glCtx.compileShader(sh);
    if (!glCtx.getShaderParameter(sh, glCtx.COMPILE_STATUS)) {
      console.warn('shader compile:', glCtx.getShaderInfoLog(sh));
      return null;
    }
    return sh;
  }

  function link(glCtx, vs, fs) {
    const prog = glCtx.createProgram();
    glCtx.attachShader(prog, vs);
    glCtx.attachShader(prog, fs);
    glCtx.linkProgram(prog);
    if (!glCtx.getProgramParameter(prog, glCtx.LINK_STATUS)) {
      console.warn('program link:', glCtx.getProgramInfoLog(prog));
      return null;
    }
    return prog;
  }

  function init() {
    if (gl && programs && layerFbos.length) return true;
    gl = null;
    programs = null;
    layerFbos = [];
    width = 0;
    height = 0;

    offscreen = document.createElement('canvas');
    gl = offscreen.getContext('webgl2', {
      alpha: false, antialias: false, preserveDrawingBuffer: true,
    });
    if (!gl) return false;

    const vs = compile(gl, gl.VERTEX_SHADER, VS);
    if (!vs) { gl = null; return false; }

    const fsPanel = compile(gl, gl.FRAGMENT_SHADER, FS_PANEL);
    const fsRgb = compile(gl, gl.FRAGMENT_SHADER, FS_COMBINE_RGB);
    const fsLab = compile(gl, gl.FRAGMENT_SHADER, FS_COMBINE_LAB);
    const fsBlit = compile(gl, gl.FRAGMENT_SHADER, FS_BLIT);
    if (!fsPanel || !fsRgb || !fsLab || !fsBlit) { gl = null; return false; }

    programs = {
      panel: link(gl, vs, fsPanel),
      combineRgb: link(gl, vs, fsRgb),
      combineLab: link(gl, vs, fsLab),
      blit: link(gl, vs, fsBlit),
    };
    if (!programs.panel || !programs.combineRgb || !programs.combineLab || !programs.blit) {
      gl = null;
      return false;
    }

    quadVao = gl.createVertexArray();
    gl.bindVertexArray(quadVao);
    const buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([
      -1, -1, 1, -1, -1, 1, 1, 1,
    ]), gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(programs.panel, 'a_pos');
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
    gl.bindVertexArray(null);

    for (let i = 0; i < MAX_LAYERS; i++) {
      const tex = gl.createTexture();
      const fbo = gl.createFramebuffer();
      layerFbos.push({ tex, fbo });
    }
    return true;
  }

  function ensureSize(nx, ny) {
    if (width === nx && height === ny) return;
    width = nx;
    height = ny;
    offscreen.width = nx;
    offscreen.height = ny;
    gl.viewport(0, 0, nx, ny);
    for (const lf of layerFbos) {
      gl.bindTexture(gl.TEXTURE_2D, lf.tex);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, nx, ny, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
      gl.bindFramebuffer(gl.FRAMEBUFFER, lf.fbo);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, lf.tex, 0);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  }

  function uploadDataTex(tex, data, nx, ny) {
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, nx, ny, 0, gl.RED, gl.FLOAT, data);
  }

  function renderPanelScreen(dataTex, nx, ny, vmin, vmax, stretch, color, inverse, gamma) {
    ensureSize(nx, ny);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.useProgram(programs.panel);
    gl.bindVertexArray(quadVao);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, dataTex);
    gl.uniform1i(gl.getUniformLocation(programs.panel, 'u_data'), 0);
    gl.uniform1f(gl.getUniformLocation(programs.panel, 'u_vmin'), vmin);
    gl.uniform1f(gl.getUniformLocation(programs.panel, 'u_vmax'), vmax);
    gl.uniform1i(gl.getUniformLocation(programs.panel, 'u_stretch'), STRETCH_ID[stretch] || 0);
    gl.uniform3f(gl.getUniformLocation(programs.panel, 'u_color'), color[0], color[1], color[2]);
    gl.uniform1f(gl.getUniformLocation(programs.panel, 'u_gamma'), gamma || 2.2);
    gl.uniform1f(gl.getUniformLocation(programs.panel, 'u_inverse'), inverse ? 1 : 0);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }

  function drawPanelToLayer(dataTex, nx, ny, vmin, vmax, stretch, color, inverse, gamma) {
    return renderLayer(dataTex, 0, nx, ny, vmin, vmax, stretch, color, inverse, gamma);
  }

  function renderLayer(dataTex, layerSlot, nx, ny, vmin, vmax, stretch, color, inverse, gamma) {
    if (!init()) return null;
    if (layerSlot < 0 || layerSlot >= layerFbos.length) return null;
    ensureSize(nx, ny);
    const lf = layerFbos[layerSlot];
    if (!lf) return null;
    gl.bindFramebuffer(gl.FRAMEBUFFER, lf.fbo);
    gl.useProgram(programs.panel);
    gl.bindVertexArray(quadVao);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, dataTex);
    gl.uniform1i(gl.getUniformLocation(programs.panel, 'u_data'), 0);
    gl.uniform1f(gl.getUniformLocation(programs.panel, 'u_vmin'), vmin);
    gl.uniform1f(gl.getUniformLocation(programs.panel, 'u_vmax'), vmax);
    gl.uniform1i(gl.getUniformLocation(programs.panel, 'u_stretch'), STRETCH_ID[stretch] || 0);
    gl.uniform3f(gl.getUniformLocation(programs.panel, 'u_color'), color[0], color[1], color[2]);
    gl.uniform1f(gl.getUniformLocation(programs.panel, 'u_gamma'), gamma || 2.2);
    gl.uniform1f(gl.getUniformLocation(programs.panel, 'u_inverse'), inverse ? 1 : 0);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    return lf.tex;
  }

  function combineLayers(count, mode, blend, gamma, inverse) {
    if (!init() || count < 1 || count > layerFbos.length) return;
    const prog = mode === 'lab' ? programs.combineLab : programs.combineRgb;
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.useProgram(prog);
    gl.bindVertexArray(quadVao);
    for (let i = 0; i < count; i++) {
      const lf = layerFbos[i];
      if (!lf) return;
      gl.activeTexture(gl.TEXTURE0 + i);
      gl.bindTexture(gl.TEXTURE_2D, lf.tex);
      gl.uniform1i(gl.getUniformLocation(prog, `u_layers[${i}]`), i);
    }
    gl.uniform1i(gl.getUniformLocation(prog, 'u_count'), count);
    gl.uniform1f(gl.getUniformLocation(prog, 'u_gamma'), gamma || 2.2);
    gl.uniform1f(gl.getUniformLocation(prog, 'u_inverse'), inverse ? 1 : 0);
    if (mode === 'lab') {
      const blendId = { sum: 1, screen: 2, max: 3, mean: 4 }[blend] || 2;
      gl.uniform1i(gl.getUniformLocation(prog, 'u_blend'), blendId);
    }
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }

  function blitLayerToScreen(layerSlot, inverse, gamma) {
    if (!init() || layerSlot < 0 || layerSlot >= layerFbos.length) return;
    const lf = layerFbos[layerSlot];
    if (!lf) return;
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.useProgram(programs.blit);
    gl.bindVertexArray(quadVao);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, lf.tex);
    gl.uniform1i(gl.getUniformLocation(programs.blit, 'u_tex'), 0);
    gl.uniform1f(gl.getUniformLocation(programs.blit, 'u_gamma'), gamma || 2.2);
    gl.uniform1f(gl.getUniformLocation(programs.blit, 'u_inverse'), inverse ? 1 : 0);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }

  function blitLayerToCanvas(layerSlot, targetCanvas, inverse, gamma) {
    if (!targetCanvas || layerSlot < 0) return;
    blitLayerToScreen(layerSlot, inverse, gamma);
    blitToCanvas(targetCanvas);
  }

  function blitToCanvas(targetCanvas, maxSide) {
    if (!targetCanvas || !width || !height) return;
    let dw = width;
    let dh = height;
    if (maxSide && maxSide > 0) {
      const longest = Math.max(dw, dh);
      if (longest > maxSide) {
        const s = maxSide / longest;
        dw = Math.max(1, Math.round(dw * s));
        dh = Math.max(1, Math.round(dh * s));
      }
    }
    targetCanvas.width = dw;
    targetCanvas.height = dh;
    const ctx = targetCanvas.getContext('2d');
    if (ctx) ctx.drawImage(offscreen, 0, 0, dw, dh);
  }

  function available() {
    return init();
  }

  return {
    available,
    STRETCH_ID,
    uploadDataTex,
    renderPanelScreen,
    renderLayer,
    combineLayers,
    blitLayerToCanvas,
    blitToCanvas,
    get offscreen() { return offscreen; },
  };
})();
