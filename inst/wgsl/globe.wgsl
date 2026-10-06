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
// solid, where they are linear, and the face projection (Snyder's or
// Fuller's) of its direction on the sphere, blended by the fold. The point goes into its
// face's quad, is scaled to the substrate, and its cell is the nearest
// multiple of the grid's generator; the solid's edge maps move that centre
// into the quad that owns it, and the cell ID follows hexify's numbering.

struct Face {
  centre: vec4f,     // xyz: face centre on the unit sphere
  az_a: vec4f,       // the face's azimuth is atan2(p . az_b, p . az_a)
  az_b: vec4f,
  quad: vec4f,       // quad, 60-degree turns into it, offset x, offset y
};

struct Grid {
  snyder0: vec4f,    // tan, cos of the edge angle, cot 30 degrees, sin G
  snyder1: vec4f,    // cos G, G, R', R'^2
  snyder2: vec4f,    // face-plane origin x, y, face edge, projection (0 ISEA, 1 Fuller)
  frame: vec4u,      // quad side in substrate steps, sublattice index, c, number of keys
  generator: vec4i,  // generator a + b omega
  flags: vec4u,      // every cell drawn, values given, values indexed by cell ID - 1, faces
  per_quad: vec4u,   // cells per quad: high and low 32 bits
  na_fill: vec4f,    // fill of a cell whose value is NA
  ramp_map: vec4f,   // a value v sits at clamp((v - x) * y + z, 0, 1) along the ramp
  faces: array<Face, 20>,
  edges: array<vec4i, 108>,  // per quad nine rows: its far corner's quad and whether
                             // it is a vertex quad, then each edge's map in two rows
};

@group(1) @binding(4) var<uniform> grid: Grid;
@group(1) @binding(5) var<storage, read> keys: array<vec2u>;

const PI = 3.14159265358979;
const SIN60 = 0.866025403784439;

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

  var az = atan2(dot(p, face.az_b.xyz), dot(p, face.az_a.xyz));
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

  let dz = atan2(tan_el, cos(az) + cot30 * sin(az));
  let h = acos(clamp(sin(az) * sin_g * cos_el - cos(az) * cos_g, -1.0, 1.0));
  let ag = az + g + h - PI;
  var azt = atan2(2.0 * ag, r1sq * tan_el * tan_el - 2.0 * ag * cot30);
  let denom = 2.0 * (cos(azt) + cot30 * sin(azt)) * sin(0.5 * dz);
  let rho = 2.0 * r1 * tan_el / denom * half_chord;
  azt += sector * 2.0 * PI / 3.0;
  return (vec2f(rho * sin(azt), rho * cos(azt)) + grid.snyder2.xy) / grid.snyder2.z;
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
  let b1 = atan((2.0 * ys / SQRT3 - FULLER_EL / 6.0) / FULLER_DVE);
  let b2 = atan((xs - ys / SQRT3 - FULLER_EL / 6.0) / FULLER_DVE);
  let b3 = atan((-xs - ys / SQRT3 - FULLER_EL / 6.0) / FULLER_DVE);
  let xy = vec2f(0.5 * (b2 - b3), (2.0 * b1 - b2 - b3) / (2.0 * SQRT3));
  return xy / FULLER_ARC + vec2f(0.5, 0.5 / SQRT3);
}

// The grid's face projection of the unit vector p onto face f.
fn face_xy(p: vec3f, f: u32) -> vec2f {
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

struct Cell {
  id: vec2u,         // cell ID, high and low 32 bits
  spot: vec2f,       // the point in the plane of the cell lattice, centres 1 apart
  offset: vec2f,     // the point from its cell's centre, in the same plane
};

// The cell at triangle coordinates t of face f.
fn cell_at(f: u32, t: vec2f) -> Cell {
  let place = grid.faces[f].quad;
  let turn = u32(place.y) % 6u;
  let cs = array<vec2f, 6>(vec2f(1.0, 0.0), vec2f(0.5, SIN60), vec2f(-0.5, SIN60),
                           vec2f(-1.0, 0.0), vec2f(-0.5, -SIN60), vec2f(0.5, -SIN60))[turn];
  let q = vec2f(cs.x * t.x - cs.y * t.y, cs.y * t.x + cs.x * t.y) - place.zw;

  // The substrate point as a + b omega, divided by the generator g: times
  // conj(g) = (ga - gb) - gb omega, over the norm N, with omega^2 = -1 - omega.
  let dim = grid.frame.x;
  let n = grid.frame.y;
  let p = q * f32(dim);
  let b = p.y / SIN60;
  let a = p.x + 0.5 * b;
  let ga = grid.generator.x;
  let gb = grid.generator.y;
  let ca = f32(ga - gb);
  let cb = f32(-gb);
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

  var cell: Cell;
  cell.spot = vec2f(za - 0.5 * zb, SIN60 * zb);
  let la = za - f32(qa);
  let lb = zb - f32(qb);
  cell.offset = vec2f(la - 0.5 * lb, SIN60 * lb);

  // The centre g * q in the substrate, in the quad that owns it.
  let ci = ga * qa - gb * qb;
  let cj = ga * qb + gb * qa - gb * qb;
  let own = canonicalize(i32(dim), u32(place.x), ci, cj);
  if (own.x == 0) {
    cell.id = vec2u(0u, 1u);
    return cell;
  }
  let u = u32(own.y);
  let v = u32(own.z);
  let residue = (grid.frame.z * u) % n;
  let within = add64(mul_wide(u, dim / n), vec2u(0u, (v - residue) / n));
  let k = u32(own.x) - 1u;
  let before = vec2u(grid.per_quad.x * k + mul_wide(grid.per_quad.y, k).x,
                     mul_wide(grid.per_quad.y, k).y);
  cell.id = add64(add64(before, within), vec2u(0u, 2u));
  return cell;
}

// Index of the cell's value: its ID - 1 when every cell has one, else its
// place among the keys (sorted cell IDs); -1 for a cell not given.
fn value_index(id: vec2u) -> i32 {
  if (grid.flags.z == 1u) {
    return select(-1, i32(id.y) - 1, id.x == 0u);
  }
  if (grid.frame.w == 0u) {
    return -1;
  }
  return key_index(id);
}

fn key_index(id: vec2u) -> i32 {
  var lo = 0u;
  var hi = grid.frame.w;
  while (lo < hi) {
    let mid = (lo + hi) / 2u;
    let k = keys[mid];
    if (k.x < id.x || (k.x == id.x && k.y < id.y)) {
      lo = mid + 1u;
    } else {
      hi = mid;
    }
  }
  if (lo < grid.frame.w && all(keys[lo] == id)) {
    return i32(lo);
  }
  return -1;
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
  let cell = cell_at(in.face, t);
  let gx = dpdx(cell.spot);
  let gy = dpdy(cell.spot);
  var normal = normalize(cross(dpdx(in.world), dpdy(in.world)));

  // A cell is drawn when the grid is drawn whole or the cell is given. NaN
  // stands for NA.
  var fill = vec4f(0.0);
  let k = value_index(cell.id);
  let drawn = grid.flags.x == 1u || k >= 0;
  if (k >= 0 && grid.flags.y == 1u) {
    let v = values[k];
    if ((bitcast<u32>(v) & 0x7fffffffu) > 0x7f800000u) {
      fill = grid.na_fill;
    } else {
      let pos = clamp((v - grid.ramp_map.x) * grid.ramp_map.y + grid.ramp_map.z, 0.0, 1.0);
      let n_ramp = arrayLength(&ramp);
      fill = ramp[min(u32(round(pos * f32(n_ramp - 1u))), n_ramp - 1u)];
    }
  }
  if (!drawn) {
    discard;
  }

  if (layer.params.z > 0.5) {
    if (dot(normal, in.world) < 0.0) {
      normal = -normal;
    }
    let shade = 0.80 + 0.20 * max(0.0, dot(normal, camera.light.xyz));
    fill = vec4f(fill.rgb * mix(shade, 1.0, fold), fill.a);
  }

  // Distance to the cell's edge in pixels: the edges lie halfway to the six
  // neighbours, along three directions of the lattice plane.
  let dirs = array<vec2f, 3>(vec2f(1.0, 0.0), vec2f(-0.5, SIN60), vec2f(-0.5, -SIN60));
  var edge_px = 1e9;
  for (var e = 0; e < 3; e++) {
    let d = dirs[e];
    let per_px = length(vec2f(dot(gx, d), dot(gy, d)));
    edge_px = min(edge_px, (0.5 - abs(dot(cell.offset, d))) / max(per_px, 1e-12));
  }
  let cell_px = 1.0 / max(max(length(gx), length(gy)), 1e-12);
  let width = layer.params.y;
  let line = layer.color.a * clamp(0.5 * width + 0.5 - edge_px, 0.0, 1.0) *
             smoothstep(3.0, 10.0, cell_px);

  let a = line + fill.a * (1.0 - line);
  return vec4f(layer.color.rgb * line + fill.rgb * fill.a * (1.0 - line), a);
}

// The cell under one pixel, for the readout under the pointer: its ID, the
// index of its value plus one (0 for none), and whether a surface is there.
@fragment
fn fs_pick(in: GridOut) -> @location(0) vec4u {
  let t = mix(in.tri, face_xy(normalize(in.sphere), in.face), camera.light.w);
  let cell = cell_at(in.face, t);
  return vec4u(cell.id, u32(value_index(cell.id) + 1), 1u);
}

// The cell of every probe point (a direction in xyz), for testing the lookup
// against the C++ one.
@group(2) @binding(0) var<storage, read> probe: array<vec4f>;
@group(2) @binding(1) var<storage, read_write> probed: array<vec2u>;

@compute @workgroup_size(64)
fn cs_locate(@builtin(global_invocation_id) gid: vec3u) {
  if (gid.x >= arrayLength(&probe)) {
    return;
  }
  let p = normalize(probe[gid.x].xyz);
  let f = nearest_face(p);
  probed[gid.x] = cell_at(f, face_xy(p, f)).id;
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
