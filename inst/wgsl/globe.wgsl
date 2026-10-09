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
// ISEA cells, found per pixel
// ---------------------------------------------------------------------------
//
// The faces mesh carries each vertex's face and triangle coordinates. A
// fragment's triangle coordinates are the interpolated ones on the flat
// solid, where they are linear, and the face projection (Snyder's, Fuller's
// or IVEA) of its direction on the sphere, blended by the fold. The point
// goes into its face's quad, is scaled to the substrate, and its cell is the
// nearest multiple of the grid's generator there.
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
  slots: vec4u,      // slots per block row, first slot of the part (high, low 32 bits)
  size: vec4f,       // x: mean cell width in radians
};

struct Grid {
  snyder0: vec4f,    // tan, cos of the edge angle, cot 30 degrees, sin G
  snyder1: vec4f,    // cos G, G, R', R'^2
  snyder2: vec4f,    // face-plane origin x, y, face edge, projection (0 ISEA, 1 Fuller, 2 IVEA)
  flags: vec4u,      // every cell drawn, values given, smooth fill, faces
  ids: vec4u,        // cells per quad (high, low 32 bits), levels in the table, keys stored
  table: vec4u,      // slots M, texture width and height, hash seed
  buckets: vec4u,    // offsets, offset texture width
  na_fill: vec4f,    // fill of a cell whose value is NA
  ramp_map: vec4f,   // a value v sits at clamp((v - x) * y + z, 0, 1) along the ramp
  faces: array<Face, 20>,
  edges: array<vec4i, 108>,  // per quad nine rows: its far corner's quad and whether
                             // it is a vertex quad, then each edge's map in two rows
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

// The grid's face projection of the unit vector p onto face f: snyder2.w is
// 0 for ISEA, 1 for Fuller and 2 for IVEA.
fn face_xy(p: vec3f, f: u32) -> vec2f {
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

// The point (i, j) of a quad that has stepped outside it, in the quad that
// owns it: DgQ2DDtoIConverter's reassignment, through the solid's edge maps.
// Across an edge, with d the point's distance past it along the crossed axis
// and 'along' the other coordinate, each new coordinate is
// k0 * top + k_along * along + k_d * d. A point beyond the far corner goes to
// that corner's vertex, a far edge's starting corner to its vertex quad when
// that corner is one, and a vertex quad has no box to leave.
fn canonicalize(top: i32, quad_in: u32, i_in: i32, j_in: i32) -> vec3i {
  let quad = i32(quad_in);
  let i = i_in;
  let j = j_in;
  let under_i = i < 0;
  let under_j = j < 0;
  let over_i = i >= top;
  let over_j = j >= top;
  let n_over = u32(under_i) + u32(under_j) + u32(over_i) + u32(over_j);
  let head = grid.edges[9 * quad];   // far-corner quad, vertex quad
  if (n_over == 0u || head.y == 1) {
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

// A point's place among the cells of one level, in the quad of its face
struct Spot {
  quad: u32,
  centre: vec2i,     // its cell's centre in the quad's substrate, possibly past the box
  lattice: vec2f,    // the point as za + zb omega, in units of the generator
  spot: vec2f,       // the point in the plane of the cell lattice, centres 1 apart
  offset: vec2f,     // the point from its cell's centre, in the same plane
};

// The substrate point qa + qb omega times the level's generator
fn lattice_point(level: Level, qa: i32, qb: i32) -> vec2i {
  let ga = level.shape.x;
  let gb = level.shape.y;
  return vec2i(ga * qa - gb * qb, ga * qb + gb * qa - gb * qb);
}

// The point at triangle coordinates t of face f among the cells of level lv.
fn locate(lv: u32, f: u32, t: vec2f) -> Spot {
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
  let b = p.y / SIN60;
  let a = p.x + 0.5 * b;
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
  s.centre = lattice_point(level, qa, qb);
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
  let u = u32(own.y);
  let v = u32(own.z);
  let residue = (level.frame.z * u) % n;
  let within = add64(mul_wide(u, dim / n), vec2u(0u, (v - residue) / n));
  let k = u32(own.x) - 1u;
  let before = vec2u(grid.ids.x * k + mul_wide(grid.ids.y, k).x,
                     mul_wide(grid.ids.y, k).y);
  return add64(add64(before, within), vec2u(0u, 2u));
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

// The word of table slot `key` (high, low 32 bits). Where the table is laid
// out slot by slot it sits at the key's own place; where it is packed by the
// perfect hash it sits at (key mod M + offset) mod M, with the offset of the
// key's bucket, and the key stored beside it tells whether it is there.
fn table_word(key: vec2u) -> u32 {
  let bucket = mix32(key.y ^ mix32(key.x ^ grid.table.w)) % grid.buckets.x;
  let ow = grid.buckets.y;
  let offset = textureLoad(cell_offsets, vec2u(bucket % ow, bucket / ow), 0).x;
  let m = grid.table.x;
  let slot = (key.y % m + offset) % m;
  let w = grid.table.y;
  let h = grid.table.z;
  let at = vec2u(slot % w, (slot / w) % h);
  let layer = slot / (w * h);
  let word = textureLoad(cell_words, at, layer, 0).x;
  if (grid.ids.w == 1u && any(textureLoad(cell_keys, at, layer, 0).xy != key)) {
    return ABSENT;
  }
  return word;
}

// The word of the cell centred at substrate point c of diamond quad `quad`
// at level lv, c in or past the quad's box.
fn cell_word(lv: u32, quad: u32, c: vec2i) -> u32 {
  let level = grid.levels[lv];
  let n = i32(level.frame.y);
  let row = c.x + level.shape.z;
  let col = (c.y - ((i32(level.frame.z) * c.x) % n + n) % n) / n + level.shape.w;
  if (row < 0 || row >= i32(level.frame.w) || col < 0 || col >= i32(level.slots.x)) {
    return ABSENT;
  }
  let at = (quad - 1u) * level.frame.w + u32(row);
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
      let c = lattice_point(grid.levels[lv], corners[k].x, corners[k].y);
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
           @location(2) face: u32, @location(3) tri: vec2f) -> GridOut {
  var out: GridOut;
  let p = fold_point(solid, sphere);
  out.position = to_clip(p);
  out.world = p;
  out.sphere = sphere;
  out.tri = tri;
  out.face = face;
  return out;
}

@fragment
fn fs_grid(in: GridOut) -> @location(0) vec4f {
  let fold = camera.light.w;
  let t = mix(in.tri, face_xy(normalize(in.sphere), in.face), fold);
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
    fill = vec4f(fill.rgb * mix(shade, 1.0, fold), fill.a);
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
  let t = mix(in.tri, face_xy(normalize(in.sphere), in.face), camera.light.w);
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
