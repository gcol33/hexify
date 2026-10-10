// globe.wgsl
// A grid on the icosahedron and the sphere, folded between the two.
//
// Every vertex carries its place on the flat face of the icosahedron and on
// the unit sphere, and its face; `fold` blends the two places, 0 showing the
// icosahedron and 1 the sphere. Instead of the fold, `stage` can move every
// point through Lambert's construction of Snyder's projection
// (construction_point()). The camera is the one of hexify's plot() method:
// an orthographic view, or a pinhole at `eye`, looking along -back, with the
// view framed by a square of screen of centre `frame.xy` and half width
// `frame.z`.

struct Camera {
  right: vec4f,      // camera axes in the scene
  up: vec4f,
  back: vec4f,
  eye: vec4f,        // xyz: camera position; w: 1 perspective, 0 orthographic
  frame: vec4f,      // centre x, y and half width of the framed square; scale
  viewport: vec4f,   // width, height in device pixels; near, far
  light: vec4f,      // xyz: direction light comes from; w: fold
  stage: vec4f,      // x: step of Lambert's construction, 0 to 3; y: 0 the fold
                     // instead, 1 the construction step by step, 2 as one path;
                     // z: R'; w: R' over the cosine of the arc from a face
                     // centre to its vertices
  centres: array<vec4f, 20>,  // xyz: each face's centre on the unit sphere
};

struct Layer {
  color: vec4f,      // fill or line colour, straight alpha
  params: vec4f,     // x: lift, y: line width in device pixels, z: shaded, w: coloured by value
  fade: vec4f,       // x: cell width in sphere radii, 0 for lines that never fade
};

@group(0) @binding(0) var<uniform> camera: Camera;
@group(1) @binding(0) var<uniform> layer: Layer;

// The point of a vertex on face f between the icosahedron and the sphere,
// by the fold or along Lambert's construction, lifted off the surface by a
// fraction of its radius so that layers drawn later lie on top.
fn fold_point(solid: vec3f, sphere: vec3f, f: u32) -> vec3f {
  var p = mix(solid, sphere, camera.light.w);
  if (camera.stage.y > 0.5) {
    p = construction_point(solid, sphere, f);
  }
  return p * (1.0 + layer.params.x);
}

// How round the surface is: 1 on the sphere, 0 on flat faces. It fades the
// shading of the faces.
fn roundness() -> f32 {
  if (camera.stage.y < 0.5) {
    return camera.light.w;
  }
  if (camera.stage.y > 1.5) {
    return 1.0 - clamp(camera.stage.x / 3.0, 0.0, 1.0);
  }
  return 1.0 - clamp(camera.stage.x, 0.0, 1.0);
}

// How far a fragment's cell is read from its direction on the sphere rather
// than from its triangle coordinates on the face, which are linear on the
// flat face and, along Lambert's construction, on the planes from the nudged
// point on.
fn sphere_weight() -> f32 {
  if (camera.stage.y < 0.5) {
    return camera.light.w;
  }
  if (camera.stage.y > 1.5) {
    return 1.0 - clamp(camera.stage.x - 2.0, 0.0, 1.0);
  }
  return 1.0 - clamp(camera.stage.x - 1.0, 0.0, 1.0);
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
           @location(2) item: u32, @location(4) face: f32) -> MeshOut {
  var out: MeshOut;
  let p = fold_point(solid, sphere, u32(face));
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
    color = vec4f(color.rgb * mix(shade, 1.0, roundness()), color.a);
  }
  return vec4f(color.rgb * color.a, color.a);
}

// ---------------------------------------------------------------------------
// ISEA cells, found per pixel
// ---------------------------------------------------------------------------
//
// The faces mesh carries each vertex's face and triangle coordinates. A
// fragment's triangle coordinates are the interpolated ones on the flat
// solid, where they are linear, and the face projection (Snyder's, Fuller's
// or IVEA) of its direction on the sphere, blended by sphere_weight(). The
// point
// goes into its face's quad, is scaled to the substrate, and its cell is the
// nearest multiple of the grid's generator there.
//
// A Hex9 grid's cells are a coset of their sublattice, shifted along j in
// each quad by the coset its Level carries, and no cell sits at a vertex of
// the solid; its keys count the cells from the first diamond quad's.
//
// The cells' values sit in a table laid out like the quads (see
// src/globe_table.cpp): per resolution, a block per diamond quad whose row u
// holds substrate column u, padded past the quad's box with the cells that
// own those points across its edges. The value of the cell nearest a point,
// and of the corners of the lattice triangle around it, is read where the
// shader finds them in the quad of the point's face, with no move into the
// quad that owns them. Only a cell's ID, for the readout under the pointer,
// takes that move: the solid's edge maps carry the centre into its own quad,
// and the ID follows hexify's numbering.

struct Face {
  centre: vec4f,     // xyz: face centre on the unit sphere
  az_a: vec4f,       // the face's azimuth is atan2(p . az_b, p . az_a)
  az_b: vec4f,
  quad: vec4f,       // quad, 60-degree turns into it, offset x, offset y
};

// One resolution of the grid and its part of the table
struct Level {
  frame: vec4u,      // quad side in substrate steps, sublattice index, c, rows of a block
  shape: vec4i,      // generator a + b omega, padding rows and columns of a block
  slots: vec4u,      // slots per block row, first slot of the part (high, low 32 bits),
                     // each quad's coset in two bits
  size: vec4f,       // x: mean cell width in radians
};

struct Grid {
  snyder0: vec4f,    // tan, cos of the edge angle, cot 30 degrees, sin G
  snyder1: vec4f,    // cos G, G, R', R'^2
  snyder2: vec4f,    // face-plane origin x, y, face edge, projection (0 ISEA, 1 Fuller,
                     // 2 IVEA, 3 AK, 4 AK with Hex9's warp)
  flags: vec4u,      // every cell drawn, values given, smooth fill, faces
  ids: vec4u,        // cells per quad (high, low 32 bits), levels in the table, keys stored
  table: vec4u,      // slots M, hash seed, offsets, cells before the first diamond
                     // quad (1, the vertex quad's, or 0 for Hex9)
  na_fill: vec4f,    // fill of a cell whose value is NA
  ramp_map: vec4f,   // a value v sits at clamp((v - x) * y + z, 0, 1) along the ramp
  faces: array<Face, 20>,
  edges: array<vec4i, 108>,  // per quad nine rows: its far corner's quad, whether
                             // it is a vertex quad, its fold axis and whether
                             // any quad folds, then each edge's map in two rows
  levels: array<Level, 32>,  // the grid's resolution, then each coarser one
};

@group(1) @binding(4) var<uniform> grid: Grid;
@group(1) @binding(6) var cell_words: texture_2d_array<u32>;
@group(1) @binding(7) var cell_keys: texture_2d_array<u32>;
@group(1) @binding(8) var cell_offsets: texture_2d<u32>;

const PI = 3.14159265358979;
const SIN60 = 0.866025403784439;

// Trigonometry from +, *, / and sqrt. WGSL bounds its built-in sin and cos
// only to 2^-11 absolute and atan2 to 4096 ulp, and some GPUs (Intel's) use
// that room, which moves points across cells. These polynomials (Cephes'
// single-precision ones) hold every result to about one ulp on any GPU.
const PIO2_HI = 1.5703125;
const PIO2_MID = 4.837512969970703125e-4;
const PIO2_LO = 7.54978995489188216e-8;
const TAN_PI8 = 0.414213562373095;

fn sin_poly(x: f32) -> f32 {
  let z = x * x;
  return ((-1.9515295891e-4 * z + 8.3321608736e-3) * z - 1.6666654611e-1) * z * x + x;
}

fn cos_poly(x: f32) -> f32 {
  let z = x * x;
  return ((2.443315711809948e-5 * z - 1.388731625493765e-3) * z
          + 4.166664568298827e-2) * z * z - 0.5 * z + 1.0;
}

// sin and cos of x, for |x| up to a few turns: x is reduced by the nearest
// multiple of pi / 2, split in three parts so the reduction is exact.
fn sincos(x: f32) -> vec2f {
  let j = round(x * (2.0 / PI));
  let y = ((x - j * PIO2_HI) - j * PIO2_MID) - j * PIO2_LO;
  let s = sin_poly(y);
  let c = cos_poly(y);
  let q = i32(j) & 3;
  if (q == 0) {
    return vec2f(s, c);
  } else if (q == 1) {
    return vec2f(c, -s);
  } else if (q == 2) {
    return vec2f(-s, -c);
  }
  return vec2f(-c, s);
}

// atan of a in [0, 1].
fn atan_unit(a: f32) -> f32 {
  var x = a;
  var base = 0.0;
  if (x > TAN_PI8) {
    x = (x - 1.0) / (x + 1.0);
    base = 0.25 * PI;
  }
  let z = x * x;
  return base + ((((8.05374449538e-2 * z - 1.38776856032e-1) * z + 1.99777106478e-1) * z
                  - 3.33329491539e-1) * z * x + x);
}

fn atan2_hx(y: f32, x: f32) -> f32 {
  let ax = abs(x);
  let ay = abs(y);
  let hi = max(ax, ay);
  if (hi == 0.0) {
    return 0.0;
  }
  var r = atan_unit(min(ax, ay) / hi);
  if (ay > ax) {
    r = 0.5 * PI - r;
  }
  if (x < 0.0) {
    r = PI - r;
  }
  return select(r, -r, y < 0.0);
}

fn atan_hx(t: f32) -> f32 {
  return atan2_hx(t, 1.0);
}

fn acos_hx(c: f32) -> f32 {
  return atan2_hx(sqrt(max((1.0 - c) * (1.0 + c), 0.0)), c);
}

// Lambert's construction of Snyder's projection (hexify's
// projection_stages()) for the point s of the unit sphere on face f, whose
// place on the inscribed solid is `solid`, at step t = stage.x. With T the
// face centre, z the arc TS and d the direction of s on the tangent plane at
// T, s = T + |TS| (cos(z / 2) d - sin(z / 2) T), |TS| = 2 sin(z / 2). From
// t = 0 to 1 the point swings about T down onto the tangent plane, on the
// circle of radius |TS|, to the Lambert point T + |TS| d. From 1 to 2 it
// turns from the azimuth of d to Snyder's azimuth and moves out to the
// nudged point N = P / R', P being Snyder's point, R' / cos(g) times the
// inscribed solid's point. From 2 to 3 the tangent plane scales about the
// sphere's centre by R' onto the face plane, taking N to P. As one path
// (stage.y = 2) the angle below the plane, the azimuth, the radius and the
// plane's height move together, from s at t = 0 to P at t = 3.
fn construction_point(solid: vec3f, s: vec3f, f: u32) -> vec3f {
  let t = camera.stage.x;
  let r1 = camera.stage.z;
  let c = camera.centres[f].xyz;
  let n = camera.stage.w / r1 * solid;
  let r = length(s - c);
  let off = s - dot(s, c) * c;
  let len = length(off);
  var d = vec3f(0.0);
  if (len > 0.0) {
    d = off / len;
  }
  let half_z = atan2_hx(0.5 * r, sqrt(max(1.0 - 0.25 * r * r, 0.0)));
  let b = n - c;
  let turn = atan2_hx(dot(cross(d, b), c), dot(d, b));
  let across = cross(c, d);
  if (camera.stage.y > 1.5) {
    let u = clamp(t / 3.0, 0.0, 1.0);
    let sc_turn = sincos(u * turn);
    let sc_down = sincos((1.0 - u) * half_z);
    let dir = sc_turn.y * d + sc_turn.x * across;
    let rho = mix(r, r1 * length(b), u);
    return mix(1.0, r1, u) * c + rho * (sc_down.y * dir - sc_down.x * c);
  }
  if (t <= 1.0) {
    let sc = sincos((1.0 - max(t, 0.0)) * half_z);
    return c + r * (sc.y * d - sc.x * c);
  }
  if (t <= 2.0) {
    let u = t - 1.0;
    let sc_turn = sincos(u * turn);
    return c + mix(r, length(b), u) * (sc_turn.y * d + sc_turn.x * across);
  }
  return mix(1.0, r1, min(t - 2.0, 1.0)) * n;
}

// Snyder's forward projection of the unit vector p onto face f, as triangle
// coordinates (a face edge is 1). sin(z / 2) is half the chord from the face
// centre, which keeps precision near the centre.
fn snyder_face(p: vec3f, f: u32) -> vec2f {
  let face = grid.faces[f];
  let tan_el = grid.snyder0.x;
  let cos_el = grid.snyder0.y;
  let cot30 = grid.snyder0.z;
  let sin_g = grid.snyder0.w;
  let cos_g = grid.snyder1.x;
  let g = grid.snyder1.y;
  let r1 = grid.snyder1.z;
  let r1sq = grid.snyder1.w;
  let half_chord = 0.5 * length(p - face.centre.xyz);

  var az = atan2_hx(dot(p, face.az_b.xyz), dot(p, face.az_a.xyz));
  if (az < 0.0) {
    az += 2.0 * PI;
  }
  var sector = 0.0;
  if (az >= 4.0 * PI / 3.0) {
    sector = 2.0;
  } else if (az >= 2.0 * PI / 3.0) {
    sector = 1.0;
  }
  az -= sector * 2.0 * PI / 3.0;

  let sc_az = sincos(az);
  let dz = atan2_hx(tan_el, sc_az.y + cot30 * sc_az.x);
  let h = acos_hx(clamp(sc_az.x * sin_g * cos_el - sc_az.y * cos_g, -1.0, 1.0));
  let ag = az + g + h - PI;
  let azt = atan2_hx(2.0 * ag, r1sq * tan_el * tan_el - 2.0 * ag * cot30);
  let sc_azt = sincos(azt);
  let denom = 2.0 * (sc_azt.y + cot30 * sc_azt.x) * sincos(0.5 * dz).x;
  let rho = 2.0 * r1 * tan_el / denom * half_chord;
  let sc_out = sincos(azt + sector * 2.0 * PI / 3.0);
  return (vec2f(rho * sc_out.x, rho * sc_out.y) + grid.snyder2.xy) / grid.snyder2.z;
}

// Fuller's projection (Gray 1995) of the unit vector p onto face f, as
// triangle coordinates. The point on the face's plane triangle is the
// gnomonic one; its distances along the triangle's edges become arc lengths
// a1, a2, a3 along the spherical face, and those lengths place the point on
// the plane triangle of edge ARC.
const FULLER_ARC = 1.10714871779409;      // atan(2)
const FULLER_EL = 1.05146222423827;       // sqrt(8) / sqrt(5 + sqrt(5))
const FULLER_DVE = 0.85065080835204;      // sqrt(3 + sqrt(5)) / sqrt(5 + sqrt(5))
const FULLER_Z0 = 0.794654472291766;      // sqrt(5 + 2 sqrt(5)) / sqrt(15)
const SQRT3 = 1.73205080756888;

fn fuller_face(p: vec3f, f: u32) -> vec2f {
  let face = grid.faces[f];
  let s = FULLER_Z0 / dot(p, face.centre.xyz);
  let xs = s * dot(p, face.az_b.xyz);
  let ys = s * dot(p, face.az_a.xyz);
  let b1 = atan_hx((2.0 * ys / SQRT3 - FULLER_EL / 6.0) / FULLER_DVE);
  let b2 = atan_hx((xs - ys / SQRT3 - FULLER_EL / 6.0) / FULLER_DVE);
  let b3 = atan_hx((-xs - ys / SQRT3 - FULLER_EL / 6.0) / FULLER_DVE);
  let xy = vec2f(0.5 * (b2 - b3), (2.0 * b1 - b2 - b3) / (2.0 * SQRT3));
  return xy / FULLER_ARC + vec2f(0.5, 0.5 / SQRT3);
}

// van Leeuwen and Strebe's (2006) vertex-oriented equal-area projection
// (IVEA) of the unit vector p onto face f, as triangle coordinates; see
// src/projection_ivea.cpp. The face splits into six right triangles (edge
// midpoint A, vertex B, centre C); the great circle from B through p maps to
// the line from B' to D' on C'A', placed by area, and p sits along it by
// (B'P' / B'D') = sin(x / 2) / sin(BD / 2), x = BP.
fn ivea_face(p: vec3f, f: u32) -> vec2f {
  let face = grid.faces[f];
  let tan_bc = grid.snyder0.x;
  let cos_bc = grid.snyder0.y;
  let sin_bc = tan_bc * cos_bc;
  let beta = grid.snyder1.y;
  let tan_ab = grid.snyder1.x * tan_bc;
  let cos_ab = 1.0 / sqrt(1.0 + tan_ab * tan_ab);
  let excess = beta - PI / 6.0;

  // Arc z from the centre, from the half chord
  let h = 0.5 * length(p - face.centre.xyz);
  let sz = 2.0 * h * sqrt(max(1.0 - h * h, 0.0));
  let cz = 1.0 - 2.0 * h * h;

  var az = atan2_hx(dot(p, face.az_b.xyz), dot(p, face.az_a.xyz));
  if (az < 0.0) {
    az += 2.0 * PI;
  }
  var sector = 0.0;
  if (az >= 4.0 * PI / 3.0) {
    sector = 2.0;
  } else if (az >= 2.0 * PI / 3.0) {
    sector = 1.0;
  }
  var a = az - sector * 2.0 * PI / 3.0;
  let mirror = a > PI / 3.0;
  if (mirror) {
    a = 2.0 * PI / 3.0 - a;
  }

  // The triangle's vertex B, and sin(x / 2) as half the chord to it
  let az_v = select(sector, sector + 1.0, mirror) * 2.0 * PI / 3.0;
  let sc_v = sincos(az_v);
  let vb = cos_bc * face.centre.xyz + sin_bc * (sc_v.y * face.az_a.xyz + sc_v.x * face.az_b.xyz);
  let sin_half_x = 0.5 * length(p - vb);

  let sc_a = sincos(a);
  let theta = atan2_hx(sc_a.x * sz, sin_bc * cz - cos_bc * sz * sc_a.y);
  let rho = beta - theta;
  let sc_rho = sincos(rho);
  let delta = acos_hx(sc_rho.x * cos_ab);
  let frac = (beta + PI / 3.0 - rho - delta) / excess;
  let bd = atan_hx(tan_ab / sc_rho.y);
  let s = sin_half_x / sincos(0.5 * bd).x;

  let rv = 1.0 / SQRT3;
  let rin = 0.5 / SQRT3;
  var q = vec2f(s * frac * rin * SIN60, rv + s * (frac * rin * 0.5 - rv));
  if (mirror) {
    let u = vec2f(SIN60, 0.5);
    q = 2.0 * dot(q, u) * u - q;
  }
  let sc_k = sincos(sector * 2.0 * PI / 3.0);
  return vec2f(q.x * sc_k.y + q.y * sc_k.x + 0.5, q.y * sc_k.y - q.x * sc_k.x + rin);
}

// Hex9's warp (projection "akw"; src/hex9_warp.cpp), in 32-bit floats.
// The field texture holds libhex9's level-6 field as hex9_warp.cpp reads it:
// word 0 is n, the lattice steps along a face edge; words 1 .. n + 1 the
// first point of each wedge column i, the column's points running j = i,
// i + 2, ...; then six floats per wedge point, (dx, dy) and the gradients
// d dx/dx, d dx/dy, d dy/dx, d dy/dy. All of it lives in libhex9's chart of
// a face, a triangle of edge sqrt(2) pointing down, its lower corner C at
// (0, -2H/3), H = sqrt(6) / 2. The displacement elsewhere on the face follows
// from the face's six symmetries and the reflection across its right edge;
// between lattice points it is the Clough-Tocher cubic.

@group(1) @binding(9) var warp_field: texture_2d<u32>;

const WARP_W = 1.41421356237310;     // the chart's edge
const WARP_H = 1.22474487139159;     // its height
const WARP_VF = -0.816496580927726;  // its lower corner C, (0, -2H/3)
const WARP_INR = 0.288675134594813;  // the unit triangle's inradius
const WARP_TOL = 1e-7;

fn warp_word(k: u32) -> u32 {
  return textureLoad(warp_field, vec2u(k & 8191u, k >> 13u), 0).x;
}

// The face's symmetries in hex9_warp.cpp's fold order: the identity, the two
// turns, then each after the mirror x -> -x
fn warp_sym(k: u32) -> mat2x2f {
  switch k {
    case 1u: { return mat2x2f(vec2f(-0.5, SIN60), vec2f(-SIN60, -0.5)); }
    case 2u: { return mat2x2f(vec2f(-0.5, -SIN60), vec2f(SIN60, -0.5)); }
    case 3u: { return mat2x2f(vec2f(-1.0, 0.0), vec2f(0.0, 1.0)); }
    case 4u: { return mat2x2f(vec2f(0.5, SIN60), vec2f(SIN60, -0.5)); }
    case 5u: { return mat2x2f(vec2f(0.5, -SIN60), vec2f(-SIN60, -0.5)); }
    default: { return mat2x2f(vec2f(1.0, 0.0), vec2f(0.0, 1.0)); }
  }
}

// The wedge's side of the median, and the face's side of the right edge
const WARP_N2 = vec2f(-0.5, -0.866025403784439);
const WARP_N3 = vec2f(-0.866025403784439, 0.5);

// The symmetry carrying p into the wedge's cone
fn warp_fold(p: vec2f) -> u32 {
  for (var k = 0u; k < 6u; k++) {
    let c = warp_sym(k) * p;
    if (c.x >= -WARP_TOL && dot(c, WARP_N2) >= -WARP_TOL) {
      return k;
    }
  }
  return 0u;
}

// A field value and its gradients at a lattice point: d = (dx, dy), and J
// with J * v the change of d along v
struct WarpDatum {
  d: vec2f,
  J: mat2x2f,
};

// The datum of the wedge point nearest q, zero off the wedge
fn wedge_datum(q: vec2f, n: i32) -> WarpDatum {
  var out: WarpDatum;
  out.d = vec2f(0.0);
  out.J = mat2x2f(vec2f(0.0), vec2f(0.0));
  let ux = WARP_W / (2.0 * f32(n));
  let uy = WARP_H / f32(n);
  let i = i32(round(q.x / ux));
  let j = i32(round((q.y - WARP_VF) / uy));
  if (i < 0 || i > n || j < i || ((i + j) & 1) != 0 || i + 3 * j > 2 * n) {
    return out;
  }
  let start = warp_word(1u + u32(i));
  let at = 2u + u32(n) + 6u * (start + u32((j - i) / 2));
  var r: array<f32, 6>;
  for (var k = 0u; k < 6u; k++) {
    r[k] = bitcast<f32>(warp_word(at + k));
  }
  out.d = vec2f(r[0], r[1]);
  out.J = mat2x2f(vec2f(r[2], r[4]), vec2f(r[3], r[5]));
  return out;
}

// The datum at a lattice point of the face or its edge band: carried into
// the wedge by the face's symmetries and, past the right edge, the
// reflection there; with p = M q, d(p) = M d(q) and J(p) = M J(q) M^T.
fn lattice_datum(p_in: vec2f, n: i32) -> WarpDatum {
  let corner = vec2f(0.0, WARP_VF);
  let s3 = mat2x2f(vec2f(-0.5, SIN60), vec2f(SIN60, 0.5));
  var M = mat2x2f(vec2f(1.0, 0.0), vec2f(0.0, 1.0));
  var p = p_in;
  for (var turn = 0; turn < 3; turn++) {
    let T = warp_sym(warp_fold(p));
    let q = T * p;
    M = M * transpose(T);
    if (dot(q - corner, WARP_N3) >= -WARP_TOL) {
      let w = wedge_datum(q, n);
      var out: WarpDatum;
      out.d = M * w.d;
      out.J = M * w.J * transpose(M);
      return out;
    }
    p = corner + s3 * (q - corner);
    M = M * s3;
  }
  var zero: WarpDatum;
  zero.d = vec2f(0.0);
  zero.J = mat2x2f(vec2f(0.0), vec2f(0.0));
  return zero;
}

// One component's Clough-Tocher cubic over the lattice triangle P, from its
// values v and gradients gr at the corners, at barycentric b (Alfeld's
// split at the centroid, the cross-boundary weights G)
fn clough_tocher(P: array<vec2f, 3>, v: vec3f, gr: array<vec2f, 3>, G: vec3f, b: vec3f) -> f32 {
  let e01 = P[1] - P[0];
  let e02 = P[2] - P[0];
  let e12 = P[2] - P[1];
  let d01 = dot(gr[0], e01);
  let d02 = dot(gr[0], e02);
  let d10 = -dot(gr[1], e01);
  let d12 = dot(gr[1], e12);
  let d20 = -dot(gr[2], e02);
  let d21 = -dot(gr[2], e12);
  let c3000 = v.x;
  let c0300 = v.y;
  let c0030 = v.z;
  let c2100 = (d01 + 3.0 * c3000) / 3.0;
  let c1200 = (d10 + 3.0 * c0300) / 3.0;
  let c2010 = (d02 + 3.0 * c3000) / 3.0;
  let c0210 = (d12 + 3.0 * c0300) / 3.0;
  let c1020 = (d20 + 3.0 * c0030) / 3.0;
  let c0120 = (d21 + 3.0 * c0030) / 3.0;
  let c2001 = (c2100 + c2010 + c3000) / 3.0;
  let c0201 = (c1200 + c0300 + c0210) / 3.0;
  let c0021 = (c1020 + c0120 + c0030) / 3.0;
  let c0111 = (G.x * (-c0300 + 3.0 * c0210 - 3.0 * c0120 + c0030) +
               (-c0300 + 2.0 * c0210 - c0120 + c0021 + c0201)) / 2.0;
  let c1011 = (G.y * (-c0030 + 3.0 * c1020 - 3.0 * c2010 + c3000) +
               (-c0030 + 2.0 * c1020 - c2010 + c2001 + c0021)) / 2.0;
  let c1101 = (G.z * (-c3000 + 3.0 * c2100 - 3.0 * c1200 + c0300) +
               (-c3000 + 2.0 * c2100 - c1200 + c2001 + c0201)) / 2.0;
  let c1002 = (c1101 + c1011 + c2001) / 3.0;
  let c0102 = (c1101 + c0111 + c0201) / 3.0;
  let c0012 = (c1011 + c0111 + c0021) / 3.0;
  let c0003 = (c1002 + c0102 + c0012) / 3.0;
  let mn = min(b.x, min(b.y, b.z));
  let s1 = b.x - mn;
  let s2 = b.y - mn;
  let s3 = b.z - mn;
  let s4 = 3.0 * mn;
  return s1 * s1 * s1 * c3000 + 3.0 * s1 * s1 * s2 * c2100 + 3.0 * s1 * s1 * s3 * c2010 +
         3.0 * s1 * s1 * s4 * c2001 + 3.0 * s1 * s2 * s2 * c1200 +
         6.0 * s1 * s2 * s4 * c1101 + 3.0 * s1 * s3 * s3 * c1020 + 6.0 * s1 * s3 * s4 * c1011 +
         3.0 * s1 * s4 * s4 * c1002 + s2 * s2 * s2 * c0300 + 3.0 * s2 * s2 * s3 * c0210 +
         3.0 * s2 * s2 * s4 * c0201 + 3.0 * s2 * s3 * s3 * c0120 + 6.0 * s2 * s3 * s4 * c0111 +
         3.0 * s2 * s4 * s4 * c0102 + s3 * s3 * s3 * c0030 + 3.0 * s3 * s3 * s4 * c0021 +
         3.0 * s3 * s4 * s4 * c0012 + s4 * s4 * s4 * c0003;
}

// The displacement at a chart point of the wedge's cone, from the lattice
// triangle holding it, and the corners' gradients blended by its
// barycentric weights, which the solve below steps with
fn wedge_delta(x: vec2f) -> WarpDatum {
  let n = i32(warp_word(0u));
  let ux = WARP_W / (2.0 * f32(n));
  let uy = WARP_H / f32(n);
  // Lattice coordinates along the 60- and 120-degree steps; the unit rhombus
  // at (U, V) splits into the triangles on its lower and upper sides
  let s = x.x / ux;
  let t = (x.y - WARP_VF) / uy;
  let u = 0.5 * (s + t);
  let v = 0.5 * (t - s);
  let U = floor(u);
  let V = floor(v);
  let lower = (u - U) + (v - V) < 1.0;
  let cu = select(vec3f(U + 1.0, U, U + 1.0), vec3f(U, U + 1.0, U), lower);
  let cv = select(vec3f(V, V + 1.0, V + 1.0), vec3f(V, V, V + 1.0), lower);
  var P: array<vec2f, 3>;
  var D: array<WarpDatum, 3>;
  for (var k = 0; k < 3; k++) {
    P[k] = vec2f((cu[k] - cv[k]) * ux, WARP_VF + (cu[k] + cv[k]) * uy);
    D[k] = lattice_datum(P[k], n);
  }
  let T = mat2x2f(P[0] - P[2], P[1] - P[2]);
  let det = T[0].x * T[1].y - T[1].x * T[0].y;
  let e = x - P[2];
  let b0 = (T[1].y * e.x - T[1].x * e.y) / det;
  let b1 = (-T[0].y * e.x + T[0].x * e.y) / det;
  let b = vec3f(b0, b1, 1.0 - b0 - b1);
  // Cross-boundary weights from the centroid of the triangle across each
  // edge, the reflection of this one's third corner through the edge
  let V4 = (P[0] + P[1] + P[2]) / 3.0;
  var G: vec3f;
  for (var k = 0; k < 3; k++) {
    let A = P[(k + 1) % 3];
    let B = P[(k + 2) % 3];
    let nb = (2.0 * A + 2.0 * B - P[k]) / 3.0;
    let dd = nb - V4;
    let a = V4 - A;
    let ab = B - A;
    G[k] = (dd.y * a.x - dd.x * a.y) / (dd.x * ab.y - dd.y * ab.x);
  }
  var out: WarpDatum;
  for (var r = 0; r < 2; r++) {
    let vals = vec3f(D[0].d[r], D[1].d[r], D[2].d[r]);
    let grs = array<vec2f, 3>(vec2f(D[0].J[0][r], D[0].J[1][r]),
                              vec2f(D[1].J[0][r], D[1].J[1][r]),
                              vec2f(D[2].J[0][r], D[2].J[1][r]));
    out.d[r] = clough_tocher(P, vals, grs, G, b);
  }
  out.J = b.x * D[0].J + b.y * D[1].J + b.z * D[2].J;
  return out;
}

// The displacement at chart point p: folded into the wedge's cone, read
// there and unfolded, d(p) = T^T d(T p), J(p) = T^T J(T p) T
fn chart_delta(p: vec2f) -> WarpDatum {
  let T = warp_sym(warp_fold(p));
  let w = wedge_delta(T * p);
  var out: WarpDatum;
  out.d = transpose(T) * w.d;
  out.J = transpose(T) * w.J * T;
  return out;
}

// The lattice point L of face triangle coordinates whose L + d(L) is the
// face point t, which AK^-1 gives: by Newton's method in the chart from t,
// stepping with the blended gradients
fn warp_solve(t: vec2f) -> vec2f {
  let goal = vec2f(WARP_W * (t.x - 0.5), -WARP_W * (t.y - WARP_INR));
  var c = goal;
  for (var it = 0; it < 5; it++) {
    let w = chart_delta(c);
    let e = c + w.d - goal;
    let A = mat2x2f(vec2f(1.0, 0.0), vec2f(0.0, 1.0)) + w.J;
    let det = A[0].x * A[1].y - A[1].x * A[0].y;
    c -= vec2f(A[1].y * e.x - A[1].x * e.y, A[0].x * e.y - A[0].y * e.x) / det;
  }
  return vec2f(c.x / WARP_W + 0.5, -c.y / WARP_W + WARP_INR);
}

// Kaseorg's octahedral projection (AK; src/projection_ak.cpp) of the unit
// vector p onto face f, as triangle coordinates: the weights b on the face's
// vertices whose direction sum_k t_k (t_i^2 + t_j^2 + alpha t_i^2 t_j^2)^(1/4)
// V_k, t_k = tan(pi b_k / 2), is p's, by Gauss-Newton on b_0 and b_2 from
// the gnomonic start.
const AK_ALPHA = 3.22780623714388;

fn ak_face(p: vec3f, f: u32) -> vec2f {
  let face = grid.faces[f];
  let cos_el = grid.snyder0.y;
  let sin_el = grid.snyder0.x * cos_el;
  let ctr = face.centre.xyz;
  let v0 = cos_el * ctr + sin_el * face.az_a.xyz;
  let v1 = cos_el * ctr + sin_el * (-0.5 * face.az_a.xyz - SIN60 * face.az_b.xyz);
  let v2 = cos_el * ctr + sin_el * (-0.5 * face.az_a.xyz + SIN60 * face.az_b.xyz);
  let c = vec3f(dot(p, v0), dot(p, v1), dot(p, v2));
  var b = c / (c.x + c.y + c.z);
  for (var it = 0; it < 8; it++) {
    var t: vec3f;
    var dt: vec3f;
    for (var k = 0; k < 3; k++) {
      let sc = sincos(0.5 * PI * b[k]);
      t[k] = sc.x / sc.y;
      dt[k] = 0.5 * PI * (1.0 + t[k] * t[k]);
    }
    // X_k and its derivatives along b_0 and b_2 (b_1 = 1 - b_0 - b_2)
    var x: vec3f;
    var du: vec3f;
    var dw: vec3f;
    for (var k = 0; k < 3; k++) {
      let i = (k + 1) % 3;
      let j = (k + 2) % 3;
      let ti = t[i] * t[i];
      let tj = t[j] * t[j];
      let s = ti + tj + AK_ALPHA * ti * tj;
      let q = sqrt(sqrt(s));
      x[k] = t[k] * q;
      var g: vec3f;
      g[k] = dt[k] * q;
      g[i] = t[k] * 0.25 * q / s * 2.0 * t[i] * dt[i] * (1.0 + AK_ALPHA * tj);
      g[j] = t[k] * 0.25 * q / s * 2.0 * t[j] * dt[j] * (1.0 + AK_ALPHA * ti);
      du[k] = g[0] - g[1];
      dw[k] = g[2] - g[1];
    }
    // The direction of X is c's where X x c vanishes.
    let r = cross(x, c);
    let ju = cross(du, c);
    let jw = cross(dw, c);
    let a = dot(ju, ju);
    let bb = dot(ju, jw);
    let d = dot(jw, jw);
    let det = a * d - bb * bb;
    if (det == 0.0) {
      break;
    }
    let gu = dot(ju, r);
    let gw = dot(jw, r);
    let su = (d * gu - bb * gw) / det;
    let sw = (a * gw - bb * gu) / det;
    b = vec3f(b.x - su, b.y + su + sw, b.z - sw);
  }
  return vec2f(0.5 * b.x + b.z, SIN60 * b.x);
}

// The grid's face projection of the unit vector p onto face f: snyder2.w is
// 0 for ISEA, 1 for Fuller, 2 for IVEA, 3 for AK and 4 for AK with Hex9's
// warp.
fn face_xy(p: vec3f, f: u32) -> vec2f {
  if (grid.snyder2.w > 3.5) {
    return warp_solve(ak_face(p, f));
  }
  if (grid.snyder2.w > 2.5) {
    return ak_face(p, f);
  }
  if (grid.snyder2.w > 1.5) {
    return ivea_face(p, f);
  }
  if (grid.snyder2.w > 0.5) {
    return fuller_face(p, f);
  }
  return snyder_face(p, f);
}

// The face whose centre is nearest: the face the point lies on.
fn nearest_face(p: vec3f) -> u32 {
  var best = 0u;
  var best_dot = -2.0;
  for (var f = 0u; f < grid.flags.w; f++) {
    let d = dot(p, grid.faces[f].centre.xyz);
    if (d > best_dot) {
      best_dot = d;
      best = f;
    }
  }
  return best;
}

// 64-bit unsigned arithmetic on (high, low) pairs.
fn mul_wide(a: u32, b: u32) -> vec2u {
  let a0 = a & 0xffffu;
  let a1 = a >> 16u;
  let b0 = b & 0xffffu;
  let b1 = b >> 16u;
  let p00 = a0 * b0;
  let p01 = a0 * b1;
  let p10 = a1 * b0;
  let mid = (p00 >> 16u) + (p01 & 0xffffu) + (p10 & 0xffffu);
  return vec2u(a1 * b1 + (p01 >> 16u) + (p10 >> 16u) + (mid >> 16u),
               (p00 & 0xffffu) | (mid << 16u));
}

fn add64(a: vec2u, b: vec2u) -> vec2u {
  let lo = a.y + b.y;
  return vec2u(a.x + b.x + select(0u, 1u, lo < a.y), lo);
}

// Whether quad `quad` holds the point (i, j) of its frame: a vertex quad its
// origin; a diamond quad its half-open box, or along its fold axis
// (src/polyhedron.h, quad_holds()) the far edge's interior in place of the
// near edge's.
fn quad_holds(top: i32, quad: i32, i: i32, j: i32) -> bool {
  let head = grid.edges[9 * quad];
  if (head.y == 1) {
    return i == 0 && j == 0;
  }
  if (head.z < 0) {
    return i >= 0 && j >= 0 && i < top && j < top;
  }
  let x = select(j, i, head.z == 0);
  let y = select(i, j, head.z == 0);
  if (y < 0 || y >= top) {
    return false;
  }
  if (y == 0) {
    return x >= 0 && x < top;
  }
  return x > 0 && x <= top;
}

// One crossing of canonicalize(): DgQ2DDtoIConverter's reassignment, through
// the solid's edge maps. Across an edge, with d the point's distance past it
// along the crossed axis and 'along' the other coordinate, each new
// coordinate is k0 * top + k_along * along + k_d * d. A point beyond the far
// corner goes to that corner's vertex, a far edge's starting corner to its
// vertex quad when that corner is one, and a vertex quad has no box to
// leave. Along a quad's fold axis the near edge's interior lies across it,
// d = 0, and a point on the line of the other far edge past the held far
// edge crosses the held one.
fn canonicalize_step(top: i32, quad: i32, i: i32, j: i32) -> vec3i {
  let head = grid.edges[9 * quad];   // far-corner quad, vertex quad, fold axis, folds
  if (head.y == 1 || quad_holds(top, quad, i, j)) {
    return vec3i(quad, i, j);
  }
  var under_i = i < 0;
  var under_j = j < 0;
  var over_i = i >= top;
  var over_j = j >= top;
  if (head.z >= 0) {
    let x = select(j, i, head.z == 0);
    let y = select(i, j, head.z == 0);
    if (x == 0 && y > 0 && y < top) {
      under_i = under_i || head.z == 0;
      under_j = under_j || head.z == 1;
    }
    if (over_i && over_j) {
      if (i == top && j > top) {
        over_i = false;
      } else if (j == top && i > top) {
        over_j = false;
      }
    }
  }
  let n_over = u32(under_i) + u32(under_j) + u32(over_i) + u32(over_j);
  if (n_over == 0u) {
    return vec3i(quad, i, j);
  }
  if (over_i && over_j) {
    return vec3i(head.x, 0, 0);
  }
  if (n_over > 1u) {
    return vec3i(quad, i, j);
  }
  var e = 3;
  var along = i;
  var d = j - top;
  if (under_i) {
    e = 0;
    along = j;
    d = i;
  } else if (under_j) {
    e = 1;
    d = j;
  } else if (over_i) {
    e = 2;
    along = j;
    d = i - top;
  }
  let m0 = grid.edges[9 * quad + 1 + 2 * e];   // quad across, vertex quad, k[0][0], k[0][1]
  let m1 = grid.edges[9 * quad + 2 + 2 * e];   // k[0][2], k[1][0], k[1][1], k[1][2]
  if (m0.y >= 0 && along == 0) {
    return vec3i(m0.y, 0, 0);
  }
  return vec3i(m0.x, m0.z * top + m0.w * along + m1.x * d,
               m1.y * top + m1.z * along + m1.w * d);
}

// The point (i, j) of a quad that has stepped outside it, in the quad that
// owns it (src/coordinate_transforms.cpp, canonicalize_q2d()): one crossing,
// or on a solid whose quads fold up to three, until a quad holds it.
fn canonicalize(top: i32, quad_in: u32, i_in: i32, j_in: i32) -> vec3i {
  var p = canonicalize_step(top, i32(quad_in), i_in, j_in);
  if (grid.edges[0].w == 1) {
    for (var k = 0; k < 2 && !quad_holds(top, p.x, p.y, p.z); k++) {
      p = canonicalize_step(top, p.x, p.y, p.z);
    }
  }
  return p;
}

// A point's place among the cells of one level, in the quad of its face
struct Spot {
  quad: u32,
  centre: vec2i,     // its cell's centre in the quad's substrate, possibly past the box
  lattice: vec2f,    // the point as za + zb omega, in units of the generator
  spot: vec2f,       // the point in the plane of the cell lattice, centres 1 apart
  offset: vec2f,     // the point from its cell's centre, in the same plane
};

// The coset of a level's cells in a quad: 0 but on Hex9
fn level_coset(level: Level, quad: u32) -> i32 {
  return i32((level.slots.w >> (2u * quad)) & 3u);
}

// The cell centre qa + qb omega times the level's generator, shifted by the
// quad's coset
fn lattice_point(level: Level, quad: u32, qa: i32, qb: i32) -> vec2i {
  let ga = level.shape.x;
  let gb = level.shape.y;
  return vec2i(ga * qa - gb * qb, ga * qb + gb * qa - gb * qb + level_coset(level, quad));
}

// Hex9 has no cell at a vertex of the solid. Inside a face the plane's
// nearest centre is always one of the solid's, but a point that rounding puts
// past a vertex can be nearest the centre that would lie in the plane beyond
// it, off the solid, and a point at the vertex is as near that one as the
// two real ones. So on Hex9 the point is held inside its face, as
// lonlat_to_cell() holds it, and 4e-6 face edges (about 6e-6 radians) from
// its edges: some thirty float steps, which the rounding cannot undo, and
// within the shader's error.
fn held_in_face(t: vec2f) -> vec2f {
  if (grid.table.w != 0u) {
    return t;
  }
  let least = 4e-6;
  let top = max(t.y / SIN60, least);
  let right = max(t.x - 0.5 * t.y / SIN60, least);
  let left = max(1.0 - t.y / SIN60 - (t.x - 0.5 * t.y / SIN60), least);
  let sum = top + right + left;
  return vec2f((0.5 * top + right) / sum, SIN60 * top / sum);
}

// The point at triangle coordinates t of face f among the cells of level lv.
fn locate(lv: u32, f: u32, t_in: vec2f) -> Spot {
  let t = held_in_face(t_in);
  let level = grid.levels[lv];
  let place = grid.faces[f].quad;
  let turn = u32(place.y) % 6u;
  let cs = array<vec2f, 6>(vec2f(1.0, 0.0), vec2f(0.5, SIN60), vec2f(-0.5, SIN60),
                           vec2f(-1.0, 0.0), vec2f(-0.5, -SIN60), vec2f(0.5, -SIN60))[turn];
  let q = vec2f(cs.x * t.x - cs.y * t.y, cs.y * t.x + cs.x * t.y) - place.zw;

  // The substrate point as a + b omega, divided by the generator g: times
  // conj(g) = (ga - gb) - gb omega, over the norm N, with omega^2 = -1 - omega.
  let n = level.frame.y;
  let p = q * f32(level.frame.x);
  let b = p.y / SIN60 - f32(level_coset(level, u32(place.x)));
  let a = p.x + 0.5 * p.y / SIN60;
  let ca = f32(level.shape.x - level.shape.y);
  let cb = f32(-level.shape.y);
  let za = (a * ca - b * cb) / f32(n);
  let zb = (a * cb + b * ca - b * cb) / f32(n);

  // Nearest Eisenstein integer, by cube rounding the axial coordinates
  // (x, y, -x - y) of the 60-degree basis 1, 1 + omega: za + zb omega =
  // (za - zb) + zb (1 + omega).
  let x = za - zb;
  let y = zb;
  let z = -x - y;
  var rx = round(x);
  var ry = round(y);
  let rz = round(z);
  let dx = abs(rx - x);
  let dy = abs(ry - y);
  let dz = abs(rz - z);
  if (dx > dy && dx > dz) {
    rx = -ry - rz;
  } else if (dy > dz) {
    ry = -rx - rz;
  }
  let qa = i32(rx + ry);
  let qb = i32(ry);

  var s: Spot;
  s.quad = u32(place.x);
  s.centre = lattice_point(level, s.quad, qa, qb);
  s.lattice = vec2f(za, zb);
  s.spot = vec2f(za - 0.5 * zb, SIN60 * zb);
  let la = za - f32(qa);
  let lb = zb - f32(qb);
  s.offset = vec2f(la - 0.5 * lb, SIN60 * lb);
  return s;
}

// The ID of a spot's cell at the grid's resolution, high and low 32 bits:
// its centre moved into the quad that owns it, numbered as hexify numbers
// the cells.
fn spot_id(s: Spot) -> vec2u {
  let level = grid.levels[0];
  let dim = level.frame.x;
  let n = level.frame.y;
  let own = canonicalize(i32(dim), s.quad, s.centre.x, s.centre.y);
  if (own.x == 0) {
    return vec2u(0u, 1u);
  }
  // A cell on the far edge of a fold axis is numbered at the near edge's place
  var place = vec2i(own.y, own.z);
  let axis = grid.edges[9 * own.x].z;
  if (axis == 0 && place.x == i32(dim) && place.y > 0) {
    place.x = 0;
  }
  if (axis == 1 && place.y == i32(dim) && place.x > 0) {
    place.y = 0;
  }
  let u = u32(place.x);
  let v = u32(place.y);
  let residue = (level.frame.z * u + u32(level_coset(level, u32(own.x)))) % n;
  let within = add64(mul_wide(u, dim / n), vec2u(0u, (v - residue) / n));
  let k = u32(own.x) - 1u;
  let before = vec2u(grid.ids.x * k + mul_wide(grid.ids.y, k).x,
                     mul_wide(grid.ids.y, k).y);
  return add64(add64(before, within), vec2u(0u, 1u + grid.table.w));
}

// The word of a slot no cell is given at
const ABSENT = 0xffffffffu;

fn mix32(a: u32) -> u32 {
  var x = a;
  x ^= x >> 16u;
  x *= 0x7feb352du;
  x ^= x >> 15u;
  x *= 0x846ca68bu;
  x ^= x >> 16u;
  return x;
}

// The word of table slot s: texel (s mod 2^13, s / 2^13 mod 2^13) of layer
// s / 2^26.
fn slot_word(s: u32) -> u32 {
  return textureLoad(cell_words, vec2u(s & 8191u, (s >> 13u) & 8191u), s >> 26u, 0).x;
}

// The word of slot `key` (high, low 32 bits) of a table packed by the perfect
// hash: at slot (key mod M + offset) mod M, with the offset of the key's
// bucket, and the key stored beside it tells whether it is there.
fn table_word(key: vec2u) -> u32 {
  let bucket = mix32(key.y ^ mix32(key.x ^ grid.table.y)) % grid.table.z;
  let offset = textureLoad(cell_offsets, vec2u(bucket & 8191u, bucket >> 13u), 0).x;
  let m = grid.table.x;
  let slot = (key.y % m + offset) % m;
  if (any(textureLoad(cell_keys, vec2u(slot & 8191u, (slot >> 13u) & 8191u), slot >> 26u, 0).xy != key)) {
    return ABSENT;
  }
  return slot_word(slot);
}

// The word of the cell centred at substrate point c of diamond quad `quad`
// at level lv, c in or past the quad's box.
// A lattice point's place along its row is (v - residue(u)) / index, which is
// floor(v / index). A table laid out slot by slot has fewer than 2^31 slots.
fn cell_word(lv: u32, quad: u32, c: vec2i) -> u32 {
  let level = grid.levels[lv];
  let n = i32(level.frame.y);
  var place = c.y;
  if (n != 1) {
    place = c.y / n;
    if (place * n > c.y) {
      place -= 1;
    }
  }
  let row = c.x + level.shape.z;
  let col = place + level.shape.w;
  if (row < 0 || row >= i32(level.frame.w) || col < 0 || col >= i32(level.slots.x)) {
    return ABSENT;
  }
  let at = (quad - 1u) * level.frame.w + u32(row);
  if (grid.ids.w == 0u) {
    return slot_word(level.slots.z + at * level.slots.x + u32(col));
  }
  return table_word(add64(add64(mul_wide(at, level.slots.x), vec2u(0u, u32(col))),
                          level.slots.yz));
}

fn is_na(word: u32) -> bool {
  return (word & 0x7fffffffu) > 0x7f800000u;
}

// A cell's place on the colour ramp, negative for NA, and the share of it
// the cells given cover. The grid's own resolution holds each cell's value,
// a coarser one both in 16 bits each (src/globe_table.cpp).
fn cell_sample(lv: u32, word: u32) -> vec2f {
  if (lv == 0u) {
    if (is_na(word)) {
      return vec2f(-1.0, 1.0);
    }
    let v = bitcast<f32>(word);
    return vec2f(clamp((v - grid.ramp_map.x) * grid.ramp_map.y + grid.ramp_map.z, 0.0, 1.0), 1.0);
  }
  let pos = word >> 16u;
  let cover = f32(word & 0xffffu) / 65534.0;
  return vec2f(select(f32(pos) / 65534.0, -1.0, pos == 0xffffu), cover);
}

// The colour at place `pos` along the ramp; with the smooth fill, blended
// between its two nearest colours.
fn ramp_colour(pos: f32) -> vec4f {
  let n = arrayLength(&ramp);
  if (grid.flags.z == 1u) {
    let x = pos * f32(n - 1u);
    let k = min(u32(floor(x)), n - 2u);
    return mix(ramp[k], ramp[k + 1u], x - f32(k));
  }
  return ramp[min(u32(round(pos * f32(n - 1u))), n - 1u)];
}

// The fill at a spot of level lv whose cell holds `word`: the cell's place
// on the ramp, or with the smooth fill, the places of the three cell centres
// around the spot blended linearly (barycentric) on the lattice triangle
// they make, over those that have a value; its opacity is the share of the
// cell covered. The lattice points 0, 1 and 1 + omega make one equilateral
// triangle, 0, omega and 1 + omega the other.
fn spot_fill(lv: u32, s: Spot, word: u32) -> vec4f {
  let own = cell_sample(lv, word);
  var pos = own.x;
  if (own.x >= 0.0 && grid.flags.z == 1u) {
    let x0 = floor(s.lattice.x);
    let y0 = floor(s.lattice.y);
    let fa = s.lattice.x - x0;
    let fb = s.lattice.y - y0;
    let upper = fb > fa;
    let corners = array<vec2i, 3>(vec2i(i32(x0), i32(y0)),
                                  select(vec2i(i32(x0) + 1, i32(y0)), vec2i(i32(x0), i32(y0) + 1), upper),
                                  vec2i(i32(x0) + 1, i32(y0) + 1));
    let weights = select(vec3f(1.0 - fa, fa - fb, fb), vec3f(1.0 - fb, fb - fa, fa), upper);
    var sum = 0.0;
    var weight = 0.0;
    for (var k = 0; k < 3; k++) {
      let c = lattice_point(grid.levels[lv], s.quad, corners[k].x, corners[k].y);
      let w = cell_word(lv, s.quad, c);
      if (w != ABSENT) {
        let corner = cell_sample(lv, w);
        if (corner.x >= 0.0) {
          sum += weights[k] * corner.x;
          weight += weights[k];
        }
      }
    }
    if (weight > 0.0) {
      pos = sum / weight;
    }
  }
  let colour = select(ramp_colour(pos), grid.na_fill, own.x < 0.0);
  return vec4f(colour.rgb, colour.a * own.y);
}

// The level of the table to fill from where a pixel spans `reach` radians:
// the finest whose cells are at least a pixel wide.
fn table_level(reach: f32) -> u32 {
  var lv = 0u;
  while (lv + 1u < grid.ids.z && grid.levels[lv].size.x < reach) {
    lv++;
  }
  return lv;
}

// A spot's cell at the grid's resolution, for the readout under the pointer
// and for testing: its ID, its word, and flags 1 (a surface is here), 2 (the
// cell is drawn) and 4 (it has a value).
fn spot_record(s: Spot) -> vec4u {
  var word = ABSENT;
  if (grid.ids.z > 0u) {
    word = cell_word(0u, s.quad, s.centre);
  }
  let drawn = grid.flags.x == 1u || word != ABSENT;
  let valued = grid.flags.y == 1u && word != ABSENT;
  return vec4u(spot_id(s), word, 1u | select(0u, 2u, drawn) | select(0u, 4u, valued));
}

struct GridOut {
  @builtin(position) position: vec4f,
  @location(0) world: vec3f,
  @location(1) sphere: vec3f,
  @location(2) tri: vec2f,
  @location(3) @interpolate(flat) face: u32,
};

@vertex
fn vs_grid(@location(0) solid: vec3f, @location(1) sphere: vec3f,
           @location(3) tri: vec2f, @location(4) face: f32) -> GridOut {
  var out: GridOut;
  let p = fold_point(solid, sphere, u32(face));
  out.position = to_clip(p);
  out.world = p;
  out.sphere = sphere;
  out.tri = tri;
  out.face = u32(face);
  return out;
}

@fragment
fn fs_grid(in: GridOut) -> @location(0) vec4f {
  let t = mix(in.tri, face_xy(normalize(in.sphere), in.face), sphere_weight());
  let reach = max(length(dpdx(in.sphere)), length(dpdy(in.sphere)));
  let lv = table_level(reach);
  let s = locate(lv, in.face, t);
  let gx = dpdx(s.spot);
  let gy = dpdy(s.spot);
  var normal = normalize(cross(dpdx(in.world), dpdy(in.world)));

  // A cell is drawn when the grid is drawn whole or the cell is given.
  var word = ABSENT;
  if (grid.ids.z > 0u) {
    word = cell_word(lv, s.quad, s.centre);
  }
  if (grid.flags.x == 0u && word == ABSENT) {
    discard;
  }
  var fill = vec4f(0.0);
  if (grid.flags.y == 1u && word != ABSENT) {
    fill = spot_fill(lv, s, word);
  }

  if (layer.params.z > 0.5) {
    if (dot(normal, in.world) < 0.0) {
      normal = -normal;
    }
    let shade = 0.80 + 0.20 * max(0.0, dot(normal, camera.light.xyz));
    fill = vec4f(fill.rgb * mix(shade, 1.0, roundness()), fill.a);
  }

  // Distance to the cell's edge in pixels: the edges lie halfway to the six
  // neighbours, along three directions of the lattice plane. Borders are
  // drawn at the grid's own resolution only; a coarser level is read where
  // its cells are below a pixel, and there they have faded out.
  var line = 0.0;
  if (lv == 0u) {
    let dirs = array<vec2f, 3>(vec2f(1.0, 0.0), vec2f(-0.5, SIN60), vec2f(-0.5, -SIN60));
    var edge_px = 1e9;
    for (var e = 0; e < 3; e++) {
      let d = dirs[e];
      let per_px = length(vec2f(dot(gx, d), dot(gy, d)));
      edge_px = min(edge_px, (0.5 - abs(dot(s.offset, d))) / max(per_px, 1e-12));
    }
    let cell_px = 1.0 / max(max(length(gx), length(gy)), 1e-12);
    let width = layer.params.y;
    line = layer.color.a * clamp(0.5 * width + 0.5 - edge_px, 0.0, 1.0) *
           smoothstep(3.0, 10.0, cell_px);
  }

  let a = line + fill.a * (1.0 - line);
  return vec4f(layer.color.rgb * line + fill.rgb * fill.a * (1.0 - line), a);
}

// The cell under one pixel, for the readout under the pointer (spot_record()).
@fragment
fn fs_pick(in: GridOut) -> @location(0) vec4u {
  let t = mix(in.tri, face_xy(normalize(in.sphere), in.face), sphere_weight());
  return spot_record(locate(0u, in.face, t));
}

// The cell of every probe point (a direction in xyz), for testing the lookup
// against the C++ one.
@group(2) @binding(0) var<storage, read> probe: array<vec4f>;
@group(2) @binding(1) var<storage, read_write> probed: array<vec4u>;

@compute @workgroup_size(64)
fn cs_locate(@builtin(global_invocation_id) gid: vec3u) {
  if (gid.x >= arrayLength(&probe)) {
    return;
  }
  let p = normalize(probe[gid.x].xyz);
  let f = nearest_face(p);
  probed[gid.x] = spot_record(locate(0u, f, face_xy(p, f)));
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

// Point k of a layer's lines: its solid and sphere positions and its face,
// seven numbers.
fn line_point(k: u32) -> vec3f {
  let o = 7u * k;
  return fold_point(vec3f(points[o], points[o + 1u], points[o + 2u]),
                    vec3f(points[o + 3u], points[o + 4u], points[o + 5u]),
                    u32(points[o + 6u]));
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
