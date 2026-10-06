// hex_globe.js
// An interactive globe for hexify grids, drawn with WebGPU.
//
// The scene arrives from R as binary buffers (base64): triangle meshes for
// the faces and land, and polylines for coastlines and face edges. Every
// point carries its place on the solid and on the sphere, and the
// shader folds one into the other. An ISEA grid comes as its frame and the
// IDs and values of its cells, and the shader finds each pixel's cell on the
// faces mesh; an H3 grid comes as cell meshes and boundary polylines. The
// camera follows hexify's plot() method: surface_view(), view_frame() and
// project() are ported below.

(function () {
  "use strict";

  const D2R = Math.PI / 180;

  // ---------------------------------------------------------------------------
  // Vectors
  // ---------------------------------------------------------------------------

  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1],
                           a[2] * b[0] - a[0] * b[2],
                           a[0] * b[1] - a[1] * b[0]];
  const norm = (a) => Math.sqrt(dot(a, a));
  const unit = (a) => mul(a, 1 / norm(a));

  function unitVec(lonDeg, latDeg) {
    const lon = lonDeg * D2R, lat = latDeg * D2R;
    return [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)];
  }

  // ---------------------------------------------------------------------------
  // Camera (R/plot_surface.R)
  // ---------------------------------------------------------------------------

  // surface_view(): the camera looking at `center`, turned by `rotation`
  // about the line of sight, swung by `tilt` about the point it looks at.
  function surfaceView(lonDeg, latDeg, distance, tilt, rotation) {
    const lon = lonDeg * D2R;
    const e3 = unitVec(lonDeg, latDeg);
    const e1 = [-Math.sin(lon), Math.cos(lon), 0];
    const e2 = cross(e3, e1);
    const r = rotation * D2R;
    const u = add(mul(e1, Math.cos(r)), mul(e2, Math.sin(r)));
    const v = add(mul(e1, -Math.sin(r)), mul(e2, Math.cos(r)));
    const t = tilt * D2R;
    const back = sub(mul(e3, Math.cos(t)), mul(v, Math.sin(t)));
    const up = add(mul(v, Math.cos(t)), mul(e3, Math.sin(t)));
    const l = unit([-0.45, 0.55, 0.70]);
    const light = add(add(mul(u, l[0]), mul(up, l[1])), mul(back, l[2]));
    const cam = { right: u, up: up, back: back };
    if (!isFinite(distance)) {
      return { cam: cam, eye: null, scale: 1, near: NaN, dir: e3, horizon: 0,
               u: u, v: v, light: light };
    }
    const eye = add(e3, mul(back, distance - 1));
    const reach = norm(eye);
    const dir = mul(eye, 1 / reach);
    return { cam: cam, eye: eye, scale: Math.sqrt(distance * distance - 1),
             near: 0.01 * (reach - 1), far: reach + 1.5, dir: dir,
             horizon: 1 / reach, u: u, v: cross(dir, u), light: light };
  }

  // project(): screen coordinates of a scene point, null behind the near plane.
  function project(p, view) {
    const c = view.cam;
    if (!view.eye) return [dot(p, c.right), dot(p, c.up)];
    const rel = sub(p, view.eye);
    const z = -dot(rel, c.back);
    if (z < view.near) return null;
    return [view.scale * dot(rel, c.right) / z, view.scale * dot(rel, c.up) / z];
  }

  function horizonRing(view, n) {
    const h = view.horizon, s = Math.sqrt(1 - h * h);
    const out = [];
    for (let k = 0; k < n; k++) {
      const a = 2 * Math.PI * k / (n - 1);
      out.push(add(mul(view.dir, h),
                   add(mul(view.u, s * Math.cos(a)), mul(view.v, s * Math.sin(a)))));
    }
    return out;
  }

  // view_frame(): centre and half width of the square of screen shown. The
  // whole sphere is framed, which holds the solid inscribed in it at
  // every stage of the fold.
  function viewFrame(view, fov) {
    if (!view.eye) return [0, 0, 1.02];
    if (fov != null) return [0, 0, 1.02 * view.scale * Math.tan(fov / 2 * D2R)];
    const widest = view.scale * Math.tan(60 * D2R);
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const p of horizonRing(view, 721)) {
      const s = project(p, view);
      if (!s) continue;
      const x = Math.min(Math.max(s[0], -widest), widest);
      const y = Math.min(Math.max(s[1], -widest), widest);
      x0 = Math.min(x0, x); x1 = Math.max(x1, x);
      y0 = Math.min(y0, y); y1 = Math.max(y1, y);
    }
    return [(x0 + x1) / 2, (y0 + y1) / 2, 1.02 * Math.max(x1 - x0, y1 - y0) / 2];
  }

  // ---------------------------------------------------------------------------
  // Data from R
  // ---------------------------------------------------------------------------

  function decode(b64, Type) {
    const bin = atob(b64);
    const bytes = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    return new Type(bytes.buffer);
  }

  function gpuBuffer(device, data, usage) {
    const size = Math.max(16, Math.ceil(data.byteLength / 4) * 4);
    const buf = device.createBuffer({ size: size, usage: usage, mappedAtCreation: true });
    new Uint8Array(buf.getMappedRange()).set(new Uint8Array(data.buffer, data.byteOffset, data.byteLength));
    buf.unmap();
    return buf;
  }

  // ---------------------------------------------------------------------------
  // The globe
  // ---------------------------------------------------------------------------

  const SAMPLES = 4;

  // Fractions of the radius each layer is lifted off the surface, so that a
  // later layer lies on top: more than a surface triangle's sag below the
  // sphere, which is under 3e-4 radii at the mesh spacing R uses.
  const LIFT = { ocean: 0, land: 6e-4, cells: 1.2e-3, grid: 2e-3, coast: 2.2e-3, edges: 2.4e-3 };

  class HexGlobe {
    constructor(el, x) {
      this.el = el;
      this.x = x;
      this.state = this.initialState();
      this.layers = [];
      this.meshGpu = new Map();
      this.frameRequested = false;
      el.innerHTML = "";
      delete el.dataset.state;
      el.classList.add("hexify-globe");
      el.hexGlobe = this;
      this.canvas = document.createElement("canvas");
      this.canvas.style.cssText = "width:100%;height:100%;display:block;touch-action:none;cursor:grab;";
      el.appendChild(this.canvas);
      this.ready = this.start().catch((err) => this.fail(String(err && err.message || err)));
    }

    initialState() {
      const c = this.x.camera;
      return {
        lon: c.center[0], lat: c.center[1],
        projection: c.projection, distance: c.distance == null ? 3 : c.distance,
        tilt: c.tilt, rotation: c.rotation, fov: c.fov, zoom: 1,
        fold: this.x.fold
      };
    }

    // `data-state` on the widget reads "drawn" once the first frame is on
    // screen and "failed" when the globe cannot be drawn, for tools that
    // capture the page, such as hex_globe_png().
    fail(message) {
      this.el.dataset.state = "failed";
      this.el.dataset.message = message;
      this.el.innerHTML = "";
      const box = document.createElement("div");
      box.style.cssText = "display:flex;align-items:center;justify-content:center;" +
        "height:100%;padding:16px;box-sizing:border-box;font:14px sans-serif;" +
        "color:#555;text-align:center;background:#F4F6F8;";
      box.textContent = message;
      this.el.appendChild(box);
    }

    async start() {
      if (!navigator.gpu) {
        this.fail("This viewer has no WebGPU, which hex_globe() draws with. " +
                  "Open the page in a browser with WebGPU enabled.");
        return;
      }
      const adapter = await navigator.gpu.requestAdapter();
      if (!adapter) {
        this.fail("WebGPU found no graphics adapter in this viewer.");
        return;
      }
      this.device = await adapter.requestDevice({
        requiredLimits: {
          maxStorageBufferBindingSize: adapter.limits.maxStorageBufferBindingSize,
          maxBufferSize: adapter.limits.maxBufferSize
        }
      });
      this.context = this.canvas.getContext("webgpu");
      this.format = navigator.gpu.getPreferredCanvasFormat();
      this.context.configure({ device: this.device, format: this.format, alphaMode: "premultiplied" });
      this.buildPipelines();
      this.buildLayers();
      this.addControls();
      this.addInteraction();
      this.resizeObserver = new ResizeObserver(() => this.resize());
      this.resizeObserver.observe(this.el);
      this.resize();
    }

    buildPipelines() {
      const device = this.device;
      const module = device.createShaderModule({ code: this.x.shader });
      this.module = module;
      const blend = {
        color: { srcFactor: "one", dstFactor: "one-minus-src-alpha", operation: "add" },
        alpha: { srcFactor: "one", dstFactor: "one-minus-src-alpha", operation: "add" }
      };
      const meshBuffers = [
        { arrayStride: 24, attributes: [
          { shaderLocation: 0, offset: 0, format: "float32x3" },
          { shaderLocation: 1, offset: 12, format: "float32x3" }] },
        { arrayStride: 4, attributes: [{ shaderLocation: 2, offset: 0, format: "uint32" }] }
      ];
      const make = (vs, fs, buffers, blended, depthWrite) => device.createRenderPipeline({
        layout: "auto",
        vertex: { module: module, entryPoint: vs, buffers: buffers },
        fragment: { module: module, entryPoint: fs,
                    targets: [{ format: this.format, blend: blended ? blend : undefined }] },
        primitive: { topology: "triangle-list", cullMode: "none" },
        depthStencil: { format: "depth24plus", depthWriteEnabled: depthWrite,
                        depthCompare: depthWrite ? "less" : "less-equal" },
        multisample: { count: SAMPLES }
      });
      const gridBuffers = [
        meshBuffers[0], meshBuffers[1],
        { arrayStride: 8, attributes: [{ shaderLocation: 3, offset: 0, format: "float32x2" }] }
      ];
      this.pipelines = {
        surface: make("vs_mesh", "fs_mesh", meshBuffers, false, true),
        overlay: make("vs_mesh", "fs_mesh", meshBuffers, true, false),
        grid: make("vs_grid", "fs_grid", gridBuffers, true, false),
        pick: device.createRenderPipeline({
          layout: "auto",
          vertex: { module: module, entryPoint: "vs_grid", buffers: gridBuffers },
          fragment: { module: module, entryPoint: "fs_pick", targets: [{ format: "rgba32uint" }] },
          primitive: { topology: "triangle-list", cullMode: "none" },
          depthStencil: { format: "depth24plus", depthWriteEnabled: true, depthCompare: "less" }
        }),
        line: make("vs_line", "fs_line",
                   [{ arrayStride: 4, stepMode: "instance",
                      attributes: [{ shaderLocation: 0, offset: 0, format: "uint32" }] }],
                   true, false)
      };
      this.cameraBuffer = device.createBuffer({
        size: 112, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
      this.cameraGroups = {};
      for (const name in this.pipelines) {
        this.cameraGroups[name] = device.createBindGroup({
          layout: this.pipelines[name].getBindGroupLayout(0),
          entries: [{ binding: 0, resource: { buffer: this.cameraBuffer } }]
        });
      }
    }

    // A layer's colour and settings, kept in its own uniform buffer.
    layerUniform(color, lift, width, shaded, ramped, fade) {
      const buf = this.device.createBuffer({
        size: 48, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
      const dpr = window.devicePixelRatio || 1;
      this.device.queue.writeBuffer(buf, 0, new Float32Array([
        color[0], color[1], color[2], color[3], lift, width * dpr, shaded ? 1 : 0, ramped ? 1 : 0,
        (fade || 0) * dpr, 0, 0, 0]));
      return buf;
    }

    // The vertex and index buffers of a mesh, made once and shared by the
    // layers drawn on it.
    meshBuffers(mesh) {
      if (!this.meshGpu.has(mesh)) {
        const device = this.device;
        this.meshGpu.set(mesh, {
          pos: gpuBuffer(device, decode(mesh.position, Float32Array), GPUBufferUsage.VERTEX),
          item: gpuBuffer(device, decode(mesh.item, Uint32Array), GPUBufferUsage.VERTEX),
          index: gpuBuffer(device, decode(mesh.index, Uint32Array), GPUBufferUsage.INDEX),
          tri: mesh.tri ? gpuBuffer(device, decode(mesh.tri, Float32Array), GPUBufferUsage.VERTEX) : null
        });
      }
      return this.meshGpu.get(mesh);
    }

    meshLayer(name, mesh, color, shaded, pipeline, values, ramp) {
      if (!mesh || !color) return;
      const device = this.device;
      const S = GPUBufferUsage.STORAGE;
      const { pos, item, index } = this.meshBuffers(mesh);
      const ramped = Boolean(values);
      const valueBuf = gpuBuffer(device, values || new Float32Array(4), S);
      const rampBuf = gpuBuffer(device, ramp || new Float32Array(4), S);
      const uniform = this.layerUniform(color, LIFT[name], 0, shaded, ramped);
      const group = device.createBindGroup({
        layout: this.pipelines[pipeline].getBindGroupLayout(1),
        entries: [
          { binding: 0, resource: { buffer: uniform } },
          { binding: 1, resource: { buffer: valueBuf } },
          { binding: 2, resource: { buffer: rampBuf } }
        ]
      });
      this.layers.push({ kind: "mesh", pipeline: pipeline, pos: pos, item: item, index: index,
                         count: mesh.n_index, group: group });
    }

    // The Grid uniform of globe.wgsl: the projection, the grid's frame, and
    // the fill of NA cells.
    gridUniform(g, naFill) {
      const buf = new ArrayBuffer(3152);
      const f = new Float32Array(buf), u = new Uint32Array(buf), i = new Int32Array(buf);
      const c = this.x.projection.constants;
      f.set(c.slice(0, 12), 0);
      u.set([g.dim, g.index, g.c, g.n_keys], 12);
      i.set(g.generator, 16);
      u.set([g.all ? 1 : 0, g.values ? 1 : 0, g.dense ? 1 : 0,
             this.x.projection.n_faces], 20);
      u.set(g.per_quad, 24);
      f.set(naFill || [0, 0, 0, 0], 28);
      f.set(g.ramp_map || [0, 0, 0, 0], 32);
      f.set(this.x.projection.faces, 36);
      i.set(this.x.projection.edges, 356);
      const out = this.device.createBuffer({
        size: buf.byteLength, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
      this.device.queue.writeBuffer(out, 0, buf);
      return out;
    }

    // ISEA cells found per pixel on the faces mesh: fill and borders at once.
    gridLayer(g, mesh, color, width, naFill) {
      const device = this.device;
      const S = GPUBufferUsage.STORAGE;
      const { pos, item, index, tri } = this.meshBuffers(mesh);
      this.gridBuffer = this.gridUniform(g, naFill);
      const keys = gpuBuffer(device, g.n_keys > 0 ? decode(g.keys, Uint32Array) : new Uint32Array(4), S);
      this.valueArray = g.values ? decode(g.values, Float32Array) : null;
      const values = gpuBuffer(device, this.valueArray || new Float32Array(4), S);
      const ramp = gpuBuffer(device, new Float32Array(this.x.palette), S);
      const uniform = this.layerUniform(color || [0, 0, 0, 0], LIFT.cells, color ? width : 0,
                                        true, Boolean(g.values));
      const group = device.createBindGroup({
        layout: this.pipelines.grid.getBindGroupLayout(1),
        entries: [
          { binding: 0, resource: { buffer: uniform } },
          { binding: 1, resource: { buffer: values } },
          { binding: 2, resource: { buffer: ramp } },
          { binding: 4, resource: { buffer: this.gridBuffer } },
          { binding: 5, resource: { buffer: keys } }
        ]
      });
      const layer = { kind: "mesh", pipeline: "grid", pos: pos, item: item, tri: tri,
                      index: index, count: mesh.n_index, group: group };
      this.layers.push(layer);
      this.pickLayer = Object.assign({}, layer, {
        pipeline: "pick",
        group: device.createBindGroup({
          layout: this.pipelines.pick.getBindGroupLayout(1),
          entries: [
            { binding: 0, resource: { buffer: uniform } },
            { binding: 4, resource: { buffer: this.gridBuffer } },
            { binding: 5, resource: { buffer: keys } }
          ]
        })
      });
    }

    // The cell under device pixel (px, py), drawn alone into a one-pixel
    // scissor of an integer target: { id, value } or null off the grid.
    async pick(px, py) {
      const device = this.device, w = this.canvas.width, h = this.canvas.height;
      if (!this.pickLayer || px < 0 || py < 0 || px >= w || py >= h) return null;
      if (!this.pickTarget || this.pickTarget.width !== w || this.pickTarget.height !== h) {
        if (this.pickTarget) {
          this.pickTarget.color.destroy();
          this.pickTarget.depth.destroy();
        }
        this.pickTarget = {
          width: w, height: h,
          color: device.createTexture({ size: [w, h], format: "rgba32uint",
            usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC }),
          depth: device.createTexture({ size: [w, h], format: "depth24plus",
            usage: GPUTextureUsage.RENDER_ATTACHMENT }),
          read: device.createBuffer({ size: 256, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST })
        };
      }
      const t = this.pickTarget, layer = this.pickLayer;
      this.writeCamera();
      const encoder = device.createCommandEncoder();
      const pass = encoder.beginRenderPass({
        colorAttachments: [{ view: t.color.createView(), clearValue: { r: 0, g: 0, b: 0, a: 0 },
                             loadOp: "clear", storeOp: "store" }],
        depthStencilAttachment: { view: t.depth.createView(), depthClearValue: 1,
                                  depthLoadOp: "clear", depthStoreOp: "discard" }
      });
      pass.setScissorRect(px, py, 1, 1);
      pass.setPipeline(this.pipelines.pick);
      pass.setBindGroup(0, this.cameraGroups.pick);
      pass.setBindGroup(1, layer.group);
      pass.setVertexBuffer(0, layer.pos);
      pass.setVertexBuffer(1, layer.item);
      pass.setVertexBuffer(2, layer.tri);
      pass.setIndexBuffer(layer.index, "uint32");
      pass.drawIndexed(layer.count);
      pass.end();
      encoder.copyTextureToBuffer({ texture: t.color, origin: { x: px, y: py } },
                                  { buffer: t.read, bytesPerRow: 256 }, [1, 1]);
      device.queue.submit([encoder.finish()]);
      await t.read.mapAsync(GPUMapMode.READ);
      const r = new Uint32Array(t.read.getMappedRange().slice(0, 16));
      t.read.unmap();
      if (r[3] === 0 || (!this.x.grid.all && r[2] === 0)) return null;
      const out = { id: r[0] * 4294967296 + r[1] };
      if (this.valueArray && r[2] > 0) out.value = this.valueArray[r[2] - 1];
      return out;
    }

    // The cell IDs of directions given as base64 xyz float triples, found by
    // the shader's own lookup.
    async locate(b64) {
      const device = this.device;
      const xyz = decode(b64, Float32Array);
      const n = xyz.length / 3;
      if (n > 64 * 65535) throw new Error("at most " + 64 * 65535 + " points at once");
      const probe = new Float32Array(4 * n);
      for (let k = 0; k < n; k++) probe.set(xyz.subarray(3 * k, 3 * k + 3), 4 * k);
      const pipeline = device.createComputePipeline({
        layout: "auto", compute: { module: this.module, entryPoint: "cs_locate" } });
      const size = Math.max(16, 8 * n);
      const out = device.createBuffer({ size: size, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC });
      const read = device.createBuffer({ size: size, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
      const groups = [
        device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [] }),
        device.createBindGroup({ layout: pipeline.getBindGroupLayout(1),
                                 entries: [{ binding: 4, resource: { buffer: this.gridBuffer } }] }),
        device.createBindGroup({ layout: pipeline.getBindGroupLayout(2), entries: [
          { binding: 0, resource: { buffer: gpuBuffer(device, probe, GPUBufferUsage.STORAGE) } },
          { binding: 1, resource: { buffer: out } }] })
      ];
      const encoder = device.createCommandEncoder();
      const pass = encoder.beginComputePass();
      pass.setPipeline(pipeline);
      groups.forEach((g, k) => pass.setBindGroup(k, g));
      pass.dispatchWorkgroups(Math.ceil(n / 64));
      pass.end();
      encoder.copyBufferToBuffer(out, 0, read, 0, size);
      device.queue.submit([encoder.finish()]);
      await read.mapAsync(GPUMapMode.READ);
      const words = new Uint32Array(read.getMappedRange().slice(0, 8 * n));
      read.unmap();
      const ids = new Float64Array(n);
      for (let k = 0; k < n; k++) ids[k] = words[2 * k] * 4294967296 + words[2 * k + 1];
      return ids;
    }

    lineLayer(name, lines, color, width, fade) {
      if (!lines || !color || lines.n_segment === 0 || !(width > 0)) return;
      const device = this.device;
      const points = gpuBuffer(device, decode(lines.position, Float32Array), GPUBufferUsage.STORAGE);
      const segment = gpuBuffer(device, decode(lines.segment, Uint32Array), GPUBufferUsage.VERTEX);
      const uniform = this.layerUniform(color, LIFT[name], width, false, false, fade);
      const group = device.createBindGroup({
        layout: this.pipelines.line.getBindGroupLayout(1),
        entries: [
          { binding: 0, resource: { buffer: uniform } },
          { binding: 3, resource: { buffer: points } }
        ]
      });
      this.layers.push({ kind: "line", pipeline: "line", segment: segment,
                         count: lines.n_segment, group: group });
    }

    buildLayers() {
      const x = this.x, s = x.style;
      this.meshLayer("ocean", x.surface, s.ocean_fill, true, "surface");
      this.meshLayer("land", x.land, s.land_fill, true, "surface");
      if (x.grid) {
        this.gridLayer(x.grid, x.surface, s.grid_border, s.grid_lwd, s.na_fill);
      }
      if (x.cells) {
        const ramp = new Float32Array(x.palette);
        this.meshLayer("cells", x.cells, s.na_fill || [0, 0, 0, 0], true, "overlay",
                       decode(x.cells.values, Float32Array), ramp);
      }
      this.lineLayer("grid", x.grid_lines, s.grid_border, s.grid_lwd, x.cell_width);
      this.lineLayer("coast", x.land_lines, s.land_border, s.land_lwd);
      this.lineLayer("edges", x.edge_lines, s.edge_col, s.edge_lwd);
    }

    resize() {
      if (!this.device) return;
      const dpr = window.devicePixelRatio || 1;
      const w = Math.max(1, Math.round(this.el.clientWidth * dpr));
      const h = Math.max(1, Math.round(this.el.clientHeight * dpr));
      if (this.canvas.width === w && this.canvas.height === h && this.colorTexture) {
        this.requestFrame();
        return;
      }
      this.canvas.width = w;
      this.canvas.height = h;
      if (this.colorTexture) this.colorTexture.destroy();
      if (this.depthTexture) this.depthTexture.destroy();
      this.colorTexture = this.device.createTexture({
        size: [w, h], sampleCount: SAMPLES, format: this.format,
        usage: GPUTextureUsage.RENDER_ATTACHMENT });
      this.depthTexture = this.device.createTexture({
        size: [w, h], sampleCount: SAMPLES, format: "depth24plus",
        usage: GPUTextureUsage.RENDER_ATTACHMENT });
      this.requestFrame();
    }

    currentView() {
      const s = this.state;
      const distance = s.projection === "perspective" ? s.distance : Infinity;
      const tilt = s.projection === "perspective" ? s.tilt : 0;
      const view = surfaceView(s.lon, s.lat, distance, tilt, s.rotation);
      const frame = viewFrame(view, s.projection === "perspective" ? s.fov : null);
      if (s.projection !== "perspective") frame[2] /= s.zoom;
      return { view: view, frame: frame };
    }

    requestFrame() {
      if (this.frameRequested) return;
      this.frameRequested = true;
      requestAnimationFrame(() => {
        this.frameRequested = false;
        this.draw();
      });
    }

    writeCamera() {
      const { view, frame } = this.currentView();
      const c = view.cam, e = view.eye || [0, 0, 0];
      this.device.queue.writeBuffer(this.cameraBuffer, 0, new Float32Array([
        c.right[0], c.right[1], c.right[2], 0,
        c.up[0], c.up[1], c.up[2], 0,
        c.back[0], c.back[1], c.back[2], 0,
        e[0], e[1], e[2], view.eye ? 1 : 0,
        frame[0], frame[1], frame[2], view.scale,
        this.canvas.width, this.canvas.height, view.eye ? view.near : 0, view.eye ? view.far : 1,
        view.light[0], view.light[1], view.light[2], this.state.fold
      ]));
    }

    draw() {
      if (!this.device || !this.colorTexture) return;
      this.writeCamera();

      const encoder = this.device.createCommandEncoder();
      const pass = encoder.beginRenderPass({
        colorAttachments: [{
          view: this.colorTexture.createView(),
          resolveTarget: this.context.getCurrentTexture().createView(),
          clearValue: { r: 0, g: 0, b: 0, a: 0 }, loadOp: "clear", storeOp: "discard"
        }],
        depthStencilAttachment: {
          view: this.depthTexture.createView(),
          depthClearValue: 1, depthLoadOp: "clear", depthStoreOp: "discard"
        }
      });
      for (const layer of this.layers) {
        pass.setPipeline(this.pipelines[layer.pipeline]);
        pass.setBindGroup(0, this.cameraGroups[layer.pipeline]);
        pass.setBindGroup(1, layer.group);
        if (layer.kind === "mesh") {
          pass.setVertexBuffer(0, layer.pos);
          pass.setVertexBuffer(1, layer.item);
          if (layer.tri) pass.setVertexBuffer(2, layer.tri);
          pass.setIndexBuffer(layer.index, "uint32");
          pass.drawIndexed(layer.count);
        } else {
          pass.setVertexBuffer(0, layer.segment);
          pass.draw(6, layer.count);
        }
      }
      pass.end();
      this.device.queue.submit([encoder.finish()]);
      if (this.el.dataset.state !== "drawn") {
        this.device.queue.onSubmittedWorkDone().then(() => { this.el.dataset.state = "drawn"; });
      }
    }

    // -------------------------------------------------------------------------
    // Controls
    // -------------------------------------------------------------------------

    addControls() {
      const panel = document.createElement("div");
      panel.className = "hexify-globe-controls";
      panel.style.cssText = "position:absolute;left:8px;top:8px;display:flex;gap:10px;" +
        "align-items:center;font:12px sans-serif;color:#333;background:rgba(255,255,255,0.85);" +
        "padding:4px 8px;border-radius:4px;user-select:none;";
      if (this.x.foldable) {
        const label = document.createElement("label");
        label.textContent = "Icosahedron ";
        const slider = document.createElement("input");
        slider.type = "range";
        slider.min = "0";
        slider.max = "1";
        slider.step = "0.001";
        slider.value = String(this.state.fold);
        slider.style.width = "110px";
        slider.addEventListener("input", () => {
          this.state.fold = Number(slider.value);
          this.requestFrame();
        });
        label.appendChild(slider);
        label.appendChild(document.createTextNode(" Sphere"));
        panel.appendChild(label);
        this.foldSlider = slider;
      }
      const select = document.createElement("select");
      for (const p of ["orthographic", "perspective"]) {
        const opt = document.createElement("option");
        opt.value = p;
        opt.textContent = p;
        select.appendChild(opt);
      }
      select.value = this.state.projection;
      select.addEventListener("change", () => {
        this.state.projection = select.value;
        this.requestFrame();
      });
      panel.appendChild(select);
      const reset = document.createElement("button");
      reset.textContent = "Reset";
      reset.addEventListener("click", () => {
        this.state = this.initialState();
        select.value = this.state.projection;
        if (this.foldSlider) this.foldSlider.value = String(this.state.fold);
        this.requestFrame();
      });
      panel.appendChild(reset);
      if (getComputedStyle(this.el).position === "static") this.el.style.position = "relative";
      this.el.appendChild(panel);
      this.readout = document.createElement("div");
      this.readout.className = "hexify-globe-controls";
      this.readout.style.cssText = "position:absolute;display:none;pointer-events:none;" +
        "font:12px sans-serif;color:#222;background:rgba(255,255,255,0.9);" +
        "padding:3px 6px;border-radius:3px;white-space:pre;";
      this.el.appendChild(this.readout);
      this.canvas.addEventListener("pointerleave", () => { this.readout.style.display = "none"; });
    }

    // Drag turns the globe under the pointer; shift-drag tilts (up and down)
    // and turns the view (left and right); the wheel moves the camera closer
    // or, in the orthographic view, zooms.
    addInteraction() {
      const canvas = this.canvas;
      let last = null;
      canvas.addEventListener("pointerdown", (ev) => {
        last = [ev.clientX, ev.clientY];
        canvas.setPointerCapture(ev.pointerId);
        canvas.style.cursor = "grabbing";
      });
      canvas.addEventListener("pointerup", (ev) => {
        last = null;
        canvas.releasePointerCapture(ev.pointerId);
        canvas.style.cursor = "grab";
      });
      canvas.addEventListener("pointermove", (ev) => {
        if (!last) {
          this.hover(ev);
          return;
        }
        this.readout.style.display = "none";
        const dx = ev.clientX - last[0], dy = ev.clientY - last[1];
        last = [ev.clientX, ev.clientY];
        if (ev.shiftKey) this.turn(dx, dy);
        else this.drag(dx, dy);
        this.requestFrame();
      });
      canvas.addEventListener("wheel", (ev) => {
        ev.preventDefault();
        const s = this.state, k = Math.exp(ev.deltaY * 0.001);
        if (s.projection === "perspective") {
          s.distance = Math.min(Math.max(1 + (s.distance - 1) * k, 1.02), 60);
        } else {
          s.zoom = Math.min(Math.max(s.zoom / k, 0.5), 500);
        }
        this.requestFrame();
      }, { passive: false });
    }

    // The cell under the pointer, and its value, shown beside it. One pick
    // runs at a time; a move during it is picked when it returns.
    hover(ev) {
      if (!this.pickLayer) return;
      const rect = this.canvas.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      this.hoverAt = [ev.clientX - rect.left, ev.clientY - rect.top];
      if (this.picking) return;
      this.picking = true;
      const [x, y] = this.hoverAt;
      this.pick(Math.floor(x * dpr), Math.floor(y * dpr)).then((cell) => {
        this.picking = false;
        const box = this.readout;
        if (!cell) {
          box.style.display = "none";
        } else {
          let text = "cell " + cell.id;
          if ("value" in cell) {
            text += "\n" + (Number.isNaN(cell.value) ? "NA" : String(Number(cell.value.toPrecision(6))));
          }
          box.textContent = text;
          box.style.left = (x + 14) + "px";
          box.style.top = (y + 14) + "px";
          box.style.display = "block";
        }
        if (this.hoverAt[0] !== x || this.hoverAt[1] !== y) this.hover({
          clientX: this.hoverAt[0] + rect.left, clientY: this.hoverAt[1] + rect.top });
      }).catch(() => { this.picking = false; });
    }

    drag(dx, dy) {
      const { view, frame } = this.currentView();
      const side = Math.min(this.el.clientWidth, this.el.clientHeight);
      let perScreen = 2 * frame[2] / side;
      if (view.eye) perScreen *= (this.state.distance - 1) / view.scale;
      const d = Math.hypot(dx, dy);
      if (d === 0) return;
      const a = d * perScreen;
      const e3 = unitVec(this.state.lon, this.state.lat);
      const dir = unit(add(mul(view.u, -dx), mul(view.v, dy)));
      const c = unit(add(mul(e3, Math.cos(a)), mul(dir, Math.sin(a))));
      this.state.lon = Math.atan2(c[1], c[0]) / D2R;
      this.state.lat = Math.asin(Math.max(-1, Math.min(1, c[2]))) / D2R;
    }

    turn(dx, dy) {
      const s = this.state;
      s.rotation += dx * 0.3;
      if (s.projection === "perspective") {
        s.tilt = Math.min(Math.max(s.tilt + dy * 0.3, -85), 85);
      }
    }

    destroy() {
      if (this.resizeObserver) this.resizeObserver.disconnect();
      if (this.device) this.device.destroy();
    }
  }

  HTMLWidgets.widget({
    name: "hex_globe",
    type: "output",
    factory: function (el) {
      let globe = null;
      return {
        renderValue: function (x) {
          if (globe) globe.destroy();
          globe = new HexGlobe(el, x);
        },
        resize: function () {
          if (globe) globe.resize();
        }
      };
    }
  });
})();
