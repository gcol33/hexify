// globe.wgsl
// A grid on the icosahedron and the sphere, folded between the two.
//
// Every vertex carries its place on the flat face of the icosahedron and on
// the unit sphere; `fold` blends them, 0 showing the icosahedron and 1 the
// sphere. The camera is the one of hexify's plot() method: an orthographic
// view, or a pinhole at `eye`, looking along -back, with the view framed by
// a square of screen of centre `frame.xy` and half width `frame.z`.

struct Camera {
  right: vec4f,      // camera axes in the scene
  up: vec4f,
  back: vec4f,
  eye: vec4f,        // xyz: camera position; w: 1 perspective, 0 orthographic
  frame: vec4f,      // centre x, y and half width of the framed square; scale
  viewport: vec4f,   // width, height in device pixels; near, far
  light: vec4f,      // xyz: direction light comes from; w: fold
};

struct Layer {
  color: vec4f,      // fill or line colour, straight alpha
  params: vec4f,     // x: lift, y: line width in device pixels, z: shaded, w: coloured by value
  fade: vec4f,       // x: cell width in sphere radii, 0 for lines that never fade
};

@group(0) @binding(0) var<uniform> camera: Camera;
@group(1) @binding(0) var<uniform> layer: Layer;

// The point between the icosahedron and the sphere, lifted off the surface
// by a fraction of its radius so that layers drawn later lie on top.
fn fold_point(solid: vec3f, sphere: vec3f) -> vec3f {
  return mix(solid, sphere, camera.light.w) * (1.0 + layer.params.x);
}

// Clip-space position of a scene point.
fn to_clip(p: vec3f) -> vec4f {
  let side = min(camera.viewport.x, camera.viewport.y);
  let kx = side / camera.viewport.x / camera.frame.z;
  let ky = side / camera.viewport.y / camera.frame.z;
  if (camera.eye.w > 0.5) {
    let rel = p - camera.eye.xyz;
    let z = -dot(rel, camera.back.xyz);
    let x = (camera.frame.w * dot(rel, camera.right.xyz) - camera.frame.x * z) * kx;
    let y = (camera.frame.w * dot(rel, camera.up.xyz) - camera.frame.y * z) * ky;
    let near = camera.viewport.z;
    let far = camera.viewport.w;
    return vec4f(x, y, far * (z - near) / (far - near), z);
  }
  let x = (dot(p, camera.right.xyz) - camera.frame.x) * kx;
  let y = (dot(p, camera.up.xyz) - camera.frame.y) * ky;
  return vec4f(x, y, 0.5 - 0.25 * dot(p, camera.back.xyz), 1.0);
}

// ---------------------------------------------------------------------------
// Triangle meshes: faces, land and cells
// ---------------------------------------------------------------------------

@group(1) @binding(1) var<storage, read> values: array<f32>;
@group(1) @binding(2) var<storage, read> ramp: array<vec4f>;

struct MeshOut {
  @builtin(position) position: vec4f,
  @location(0) world: vec3f,
  @location(1) @interpolate(flat) item: u32,
};

@vertex
fn vs_mesh(@location(0) solid: vec3f, @location(1) sphere: vec3f,
           @location(2) item: u32) -> MeshOut {
  var out: MeshOut;
  let p = fold_point(solid, sphere);
  out.position = to_clip(p);
  out.world = p;
  out.item = item;
  return out;
}

@fragment
fn fs_mesh(in: MeshOut) -> @location(0) vec4f {
  var color = layer.color;
  if (layer.params.w > 0.5) {
    let v = values[in.item];
    if (v >= 0.0) {
      let n = arrayLength(&ramp);
      color = ramp[min(u32(round(v * f32(n - 1u))), n - 1u)];
    }
  }
  // The flat icosahedron is shaded by the turn of each face to the light,
  // as plot() shades it; the shading fades as the faces fold into the sphere.
  if (layer.params.z > 0.5) {
    var n = normalize(cross(dpdx(in.world), dpdy(in.world)));
    if (dot(n, in.world) < 0.0) {
      n = -n;
    }
    let shade = 0.80 + 0.20 * max(0.0, dot(n, camera.light.xyz));
    color = vec4f(color.rgb * mix(shade, 1.0, camera.light.w), color.a);
  }
  return vec4f(color.rgb * color.a, color.a);
}

// ---------------------------------------------------------------------------
// Lines: cell boundaries, coastlines and face edges
// ---------------------------------------------------------------------------
//
// Each segment between two consecutive points is one instance, drawn as a
// quad of the line's width in screen space with ends extended by half that
// width, so that the segments of a path join. The fragment fades the last
// pixel across the line. Cell boundaries fade out where cells shrink to a
// few pixels on screen, which they would otherwise cover; the cell's width
// on screen is read from the segment's own length there.

@group(1) @binding(3) var<storage, read> points: array<f32>;

struct LineOut {
  @builtin(position) position: vec4f,
  @location(0) across: f32,
  @location(1) alpha: f32,
};

fn line_point(k: u32) -> vec3f {
  let o = 6u * k;
  return fold_point(vec3f(points[o], points[o + 1u], points[o + 2u]),
                    vec3f(points[o + 3u], points[o + 4u], points[o + 5u]));
}

@vertex
fn vs_line(@builtin(vertex_index) vi: u32, @location(0) start: u32) -> LineOut {
  var out: LineOut;
  let pa = line_point(start);
  let pb = line_point(start + 1u);
  let ca = to_clip(pa);
  let cb = to_clip(pb);
  if (ca.w <= 0.0 || cb.w <= 0.0) {
    out.position = vec4f(2.0, 2.0, 2.0, 1.0);
    out.across = 0.0;
    out.alpha = 0.0;
    return out;
  }
  let half_px = 0.5 * camera.viewport.xy;
  let sa = ca.xy / ca.w * half_px;
  let sb = cb.xy / cb.w * half_px;
  var dir = sb - sa;
  let len = length(dir);
  dir = select(vec2f(1.0, 0.0), dir / max(len, 1e-6), len > 1e-6);
  let nrm = vec2f(-dir.y, dir.x);
  let half_width = 0.5 * layer.params.y + 1.0;

  // Corners of the quad: which end, and which side of the line.
  let at_b = array<f32, 6>(0.0, 1.0, 1.0, 0.0, 1.0, 0.0)[vi];
  let side = array<f32, 6>(-1.0, -1.0, 1.0, -1.0, 1.0, 1.0)[vi];
  let end_dir = select(-1.0, 1.0, at_b > 0.5);
  let s = mix(sa, sb, at_b) + nrm * side * half_width +
          dir * end_dir * 0.5 * layer.params.y;
  let c = select(ca, cb, at_b > 0.5);
  out.position = vec4f(s / half_px * c.w, c.z, c.w);
  out.across = side * half_width;
  out.alpha = 1.0;
  let world = length(pb - pa);
  if (layer.fade.x > 0.0 && world > 0.0) {
    let cell_px = len / world * layer.fade.x;
    out.alpha = smoothstep(3.0, 10.0, cell_px);
  }
  return out;
}

@fragment
fn fs_line(in: LineOut) -> @location(0) vec4f {
  let coverage = clamp(0.5 * layer.params.y + 0.5 - abs(in.across), 0.0, 1.0);
  let a = layer.color.a * coverage * in.alpha;
  return vec4f(layer.color.rgb * a, a);
}
