// hex9_libhex9_dump.cpp - read Hex9 cells and points out of libhex9
//
// Links libhex9 (https://github.com/MrBenGriffin/libhex9, Apache-2.0) and
// writes what hexify's Hex9 tests and benches compare against. Build against
// a libhex9 checkout (tested at 1eb0be7, version 2.3.1), e.g. with MinGW:
//
//   g++ -std=gnu++17 -O2 -D_USE_MATH_DEFINES -ffp-contract=off \
//       -DH9_WARP_ENABLE=1 -DH9_WARP_BLOB="\"<libhex9>/core/Sphere_l6_fund.f64g.h9warp\"" \
//       -I<libhex9> -I<libhex9>/core -c <libhex9>/{hex9_c.cpp,core/h9_warp_runtime.cpp,
//       core/h9_warp_embedded.cpp}; ar rcs libhex9.a *.o
//   g++ -std=gnu++17 -O2 -D_USE_MATH_DEFINES -I<libhex9> hex9_libhex9_dump.cpp libhex9.a
//
// Modes:
//   points N SEED LEVELS...  N points uniform on the sphere: their place on
//                            the octahedron (libhex9's face chart lifted to
//                            x, y, z with |x| + |y| + |z| = 1) and their
//                            labels at the given levels
//   cells L N SEED           every cell of level L reached by N sphere points:
//                            label, neighbours, parent and children labels
//   lonlat LEVELS...         lon lat pairs on stdin (WGS84 degrees): the
//                            labels libhex9 gives them, and the decoded
//                            centre of the deepest level
#include "hex9_c.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <set>
#include <string>
#include <vector>

static std::string label(const uint8_t u[16], int layer) {
  char buf[64];
  if (hex9_label_key(u, layer, buf, sizeof buf) < 0) return "NA";
  return buf;
}

static void sphere_point(std::mt19937_64& rng, double& lon, double& lat) {
  std::normal_distribution<double> N(0.0, 1.0);
  double x, y, z, r;
  do {
    x = N(rng); y = N(rng); z = N(rng);
    r = std::sqrt(x * x + y * y + z * z);
  } while (r < 1e-9);
  lon = std::atan2(y, x) * 180.0 / M_PI;
  lat = std::asin(z / r) * 180.0 / M_PI;
}

int main(int argc, char** argv) {
  char err[256];
  if (hex9_init(err, sizeof err)) { std::fprintf(stderr, "hex9_init: %s\n", err); return 1; }
  if (argc < 2) { std::fprintf(stderr, "mode?\n"); return 1; }
  const std::string mode = argv[1];
  if (mode == "points") {
    const long n = std::atol(argv[2]);
    std::mt19937_64 rng(std::strtoull(argv[3], nullptr, 10));
    std::vector<int> levels;
    for (int a = 4; a < argc; a++) levels.push_back(std::atoi(argv[a]));
    std::printf("x,y,z");
    for (int l : levels) std::printf(",L%d", l);
    std::printf("\n");
    for (long k = 0; k < n; k++) {
      double lon, lat, cx, cy, w[3];
      int oid;
      sphere_point(rng, lon, lat);
      hex9_project_sphere(lon, lat, &cx, &cy, &oid);
      hex9_boct_to_woct(cx, cy, oid, w);
      uint8_t full[16];
      hex9_encode_boct(cx, cy, oid, full);
      std::printf("%.17g,%.17g,%.17g", w[0], w[1], w[2]);
      for (int l : levels) std::printf(",%s", label(full, l).c_str());
      std::printf("\n");
    }
  } else if (mode == "cells") {
    const int layer = std::atoi(argv[2]);
    const long n = std::atol(argv[3]);
    std::mt19937_64 rng(std::strtoull(argv[4], nullptr, 10));
    std::set<std::string> seen;
    std::printf("label,neighbours,parent,children\n");
    for (long k = 0; k < n; k++) {
      double lon, lat;
      sphere_point(rng, lon, lat);
      uint8_t full[16], bin[16];
      hex9_encode_sphere(lon, lat, full);
      hex9_bin(full, layer, bin);
      const std::string lab = label(bin, layer);
      if (!seen.insert(lab).second) continue;
      uint8_t nb[6 * 16];
      const int nn = hex9_neighbors(full, layer, nb);
      std::string nbs;
      for (int i = 0; i < nn; i++) nbs += (i ? ";" : "") + label(nb + 16 * i, layer);
      std::string par = "NA";
      uint8_t p[16];
      if (layer > 0 && !hex9_cell_parent(bin, p)) par = label(p, layer - 1);
      uint8_t ch[9 * 16];
      std::string chs;
      if (!hex9_cell_children(bin, ch)) {
        for (int i = 0; i < 9; i++) chs += (i ? ";" : "") + label(ch + 16 * i, layer + 1);
      }
      std::printf("%s,%s,%s,%s\n", lab.c_str(), nbs.c_str(), par.c_str(), chs.c_str());
    }
  } else if (mode == "lonlat") {
    std::vector<int> levels;
    for (int a = 2; a < argc; a++) levels.push_back(std::atoi(argv[a]));
    std::printf("lon,lat");
    for (int l : levels) std::printf(",L%d", l);
    std::printf(",centre_lon,centre_lat\n");
    double lon, lat;
    while (std::scanf("%lf %lf", &lon, &lat) == 2) {
      uint8_t full[16], bin[16];
      hex9_encode(lon, lat, full);
      std::printf("%.17g,%.17g", lon, lat);
      for (int l : levels) std::printf(",%s", label(full, l).c_str());
      double clon, clat;
      hex9_bin(full, levels.back(), bin);
      hex9_decode(bin, &clon, &clat);
      std::printf(",%.17g,%.17g\n", clon, clat);
    }
  }
  return 0;
}
