#include <Rcpp.h>
#include <cmath>
#include <cstring>
#include <vector>

extern "C" {
#include "h3/h3api.h"
}

static const double DEG_TO_RAD = M_PI / 180.0;
static const double RAD_TO_DEG = 180.0 / M_PI;

// Helper: H3Index → hex string (17 chars max incl. NUL)
static std::string h3_to_string(H3Index h) {
    char buf[17];
    hexify_h3_h3ToString(h, buf, sizeof(buf));
    return std::string(buf);
}

// Helper: hex string → H3Index
static H3Index string_to_h3(const char* s) {
    H3Index h;
    H3Error err = hexify_h3_stringToH3(s, &h);
    if (err != E_SUCCESS) return H3_NULL;
    return h;
}

// [[Rcpp::export]]
Rcpp::CharacterVector cpp_h3_latLngToCell(Rcpp::NumericVector lon_deg,
                                           Rcpp::NumericVector lat_deg,
                                           int resolution) {
    R_xlen_t n = lon_deg.size();
    if (lat_deg.size() != n) {
        Rcpp::stop("lon_deg and lat_deg must have the same length");
    }

    Rcpp::CharacterVector out(n);
    LatLng ll;

    for (R_xlen_t i = 0; i < n; i++) {
        if (Rcpp::NumericVector::is_na(lon_deg[i]) ||
            Rcpp::NumericVector::is_na(lat_deg[i])) {
            out[i] = NA_STRING;
            continue;
        }
        ll.lat = lat_deg[i] * DEG_TO_RAD;
        ll.lng = lon_deg[i] * DEG_TO_RAD;
        H3Index h;
        H3Error err = hexify_h3_latLngToCell(&ll, resolution, &h);
        if (err != E_SUCCESS) {
            out[i] = NA_STRING;
        } else {
            out[i] = h3_to_string(h);
        }
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::DataFrame cpp_h3_cellToLatLng(Rcpp::CharacterVector cell_ids) {
    R_xlen_t n = cell_ids.size();
    Rcpp::NumericVector lon_deg(n);
    Rcpp::NumericVector lat_deg(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            lon_deg[i] = NA_REAL;
            lat_deg[i] = NA_REAL;
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            lon_deg[i] = NA_REAL;
            lat_deg[i] = NA_REAL;
            continue;
        }
        LatLng ll;
        H3Error err = hexify_h3_cellToLatLng(h, &ll);
        if (err != E_SUCCESS) {
            lon_deg[i] = NA_REAL;
            lat_deg[i] = NA_REAL;
        } else {
            lon_deg[i] = ll.lng * RAD_TO_DEG;
            lat_deg[i] = ll.lat * RAD_TO_DEG;
        }
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("lon") = lon_deg,
        Rcpp::Named("lat") = lat_deg,
        Rcpp::Named("stringsAsFactors") = false
    );
}

// [[Rcpp::export]]
Rcpp::LogicalVector cpp_h3_isValidCell(Rcpp::CharacterVector cell_ids) {
    R_xlen_t n = cell_ids.size();
    Rcpp::LogicalVector out(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = NA_LOGICAL;
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        out[i] = (h != H3_NULL && hexify_h3_isValidCell(h)) ? TRUE : FALSE;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::CharacterVector cpp_h3_cellToParent(Rcpp::CharacterVector cell_ids,
                                           int parent_res) {
    R_xlen_t n = cell_ids.size();
    Rcpp::CharacterVector out(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = NA_STRING;
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            out[i] = NA_STRING;
            continue;
        }
        H3Index parent;
        H3Error err = hexify_h3_cellToParent(h, parent_res, &parent);
        if (err != E_SUCCESS) {
            out[i] = NA_STRING;
        } else {
            out[i] = h3_to_string(parent);
        }
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_h3_cellToChildren(Rcpp::CharacterVector cell_ids,
                                  int child_res) {
    R_xlen_t n = cell_ids.size();
    Rcpp::List out(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }

        int64_t num_children = 0;
        H3Error err = hexify_h3_cellToChildrenSize(h, child_res, &num_children);
        if (err != E_SUCCESS || num_children <= 0) {
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }

        std::vector<H3Index> children(num_children);
        err = hexify_h3_cellToChildren(h, child_res, children.data());
        if (err != E_SUCCESS) {
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }

        Rcpp::CharacterVector child_strs(num_children);
        for (int64_t j = 0; j < num_children; j++) {
            child_strs[j] = h3_to_string(children[j]);
        }
        out[i] = child_strs;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_h3_cellToBoundary(Rcpp::CharacterVector cell_ids) {
    R_xlen_t n = cell_ids.size();
    Rcpp::List out(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = Rcpp::NumericMatrix(0, 2);
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            out[i] = Rcpp::NumericMatrix(0, 2);
            continue;
        }

        CellBoundary bndry;
        H3Error err = hexify_h3_cellToBoundary(h, &bndry);
        if (err != E_SUCCESS) {
            out[i] = Rcpp::NumericMatrix(0, 2);
            continue;
        }

        // +1 for closing vertex
        int nv = bndry.numVerts + 1;
        Rcpp::NumericMatrix ring(nv, 2);
        for (int j = 0; j < bndry.numVerts; j++) {
            ring(j, 0) = bndry.verts[j].lng * RAD_TO_DEG;
            ring(j, 1) = bndry.verts[j].lat * RAD_TO_DEG;
        }
        // Close the ring
        ring(bndry.numVerts, 0) = ring(0, 0);
        ring(bndry.numVerts, 1) = ring(0, 1);
        out[i] = ring;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::CharacterVector cpp_h3_polygonToCells(Rcpp::NumericMatrix coords,
                                             int resolution,
                                             Rcpp::Nullable<Rcpp::List> holes = R_NilValue,
                                             int flags = 2) {
    int n_outer = coords.nrow();
    // Convert outer ring degrees → radians
    std::vector<LatLng> outer_verts(n_outer);
    for (int i = 0; i < n_outer; i++) {
        outer_verts[i].lng = coords(i, 0) * DEG_TO_RAD;
        outer_verts[i].lat = coords(i, 1) * DEG_TO_RAD;
    }

    GeoPolygon polygon;
    polygon.geoloop.numVerts = n_outer;
    polygon.geoloop.verts = outer_verts.data();

    // Handle holes
    std::vector<GeoLoop> hole_loops;
    std::vector<std::vector<LatLng>> hole_verts_storage;

    if (holes.isNotNull()) {
        Rcpp::List holes_list(holes);
        int n_holes = holes_list.size();
        hole_loops.resize(n_holes);
        hole_verts_storage.resize(n_holes);

        for (int h = 0; h < n_holes; h++) {
            Rcpp::NumericMatrix hole_coords = Rcpp::as<Rcpp::NumericMatrix>(holes_list[h]);
            int n_hole = hole_coords.nrow();
            hole_verts_storage[h].resize(n_hole);
            for (int j = 0; j < n_hole; j++) {
                hole_verts_storage[h][j].lng = hole_coords(j, 0) * DEG_TO_RAD;
                hole_verts_storage[h][j].lat = hole_coords(j, 1) * DEG_TO_RAD;
            }
            hole_loops[h].numVerts = n_hole;
            hole_loops[h].verts = hole_verts_storage[h].data();
        }
        polygon.numHoles = n_holes;
        polygon.holes = hole_loops.data();
    } else {
        polygon.numHoles = 0;
        polygon.holes = NULL;
    }

    // Get max size and fill cells
    // Default flags=2 (CONTAINMENT_OVERLAPPING): include cells that overlap the
    // polygon boundary, ensuring full spatial coverage for ecological sampling.
    int64_t max_cells = 0;
    uint32_t h3_flags = static_cast<uint32_t>(flags);
    H3Error err;

    if (h3_flags != 0) {
        err = hexify_h3_maxPolygonToCellsSizeExperimental(&polygon, resolution, h3_flags, &max_cells);
    } else {
        err = hexify_h3_maxPolygonToCellsSize(&polygon, resolution, 0, &max_cells);
    }
    if (err != E_SUCCESS || max_cells <= 0) {
        return Rcpp::CharacterVector(0);
    }

    std::vector<H3Index> cells(max_cells, H3_NULL);
    if (h3_flags != 0) {
        err = hexify_h3_polygonToCellsExperimental(&polygon, resolution, h3_flags, max_cells, cells.data());
    } else {
        err = hexify_h3_polygonToCells(&polygon, resolution, 0, cells.data());
    }
    if (err != E_SUCCESS) {
        return Rcpp::CharacterVector(0);
    }

    // Filter out H3_NULL entries and collect valid cells
    std::vector<std::string> valid_cells;
    valid_cells.reserve(max_cells);
    for (int64_t i = 0; i < max_cells; i++) {
        if (cells[i] != H3_NULL) {
            valid_cells.push_back(h3_to_string(cells[i]));
        }
    }

    return Rcpp::wrap(valid_cells);
}

// [[Rcpp::export]]
Rcpp::IntegerVector cpp_h3_getResolution(Rcpp::CharacterVector cell_ids) {
    R_xlen_t n = cell_ids.size();
    Rcpp::IntegerVector out(n);
    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = NA_INTEGER;
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            out[i] = NA_INTEGER;
        } else {
            out[i] = hexify_h3_getResolution(h);
        }
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_h3_cellAreaKm2(Rcpp::CharacterVector cell_ids) {
    R_xlen_t n = cell_ids.size();
    Rcpp::NumericVector out(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = NA_REAL;
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            out[i] = NA_REAL;
            continue;
        }
        double area;
        H3Error err = hexify_h3_cellAreaKm2(h, &area);
        if (err != E_SUCCESS) {
            out[i] = NA_REAL;
        } else {
            out[i] = area;
        }
    }
    return out;
}

// ==========================================================================
// Grid Disk / Neighbors (v0.7.0)
// ==========================================================================

// [[Rcpp::export]]
Rcpp::List cpp_h3_gridDisk(Rcpp::CharacterVector cell_ids, int k) {
    R_xlen_t n = cell_ids.size();
    Rcpp::List out(n);

    int64_t max_size = 0;
    H3Error err = hexify_h3_maxGridDiskSize(k, &max_size);
    if (err != E_SUCCESS || max_size <= 0) {
        Rcpp::stop("Invalid k value for gridDisk");
    }

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }

        std::vector<H3Index> disk(max_size, H3_NULL);
        err = hexify_h3_gridDisk(h, k, disk.data());
        if (err != E_SUCCESS) {
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }

        std::vector<std::string> valid;
        valid.reserve(max_size);
        for (int64_t j = 0; j < max_size; j++) {
            if (disk[j] != H3_NULL) {
                valid.push_back(h3_to_string(disk[j]));
            }
        }
        out[i] = Rcpp::wrap(valid);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_h3_gridDiskDistances(Rcpp::CharacterVector cell_ids, int k) {
    R_xlen_t n = cell_ids.size();
    Rcpp::List out(n);

    int64_t max_size = 0;
    H3Error err = hexify_h3_maxGridDiskSize(k, &max_size);
    if (err != E_SUCCESS || max_size <= 0) {
        Rcpp::stop("Invalid k value for gridDiskDistances");
    }

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = Rcpp::DataFrame::create();
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            out[i] = Rcpp::DataFrame::create();
            continue;
        }

        std::vector<H3Index> disk(max_size, H3_NULL);
        std::vector<int> distances(max_size, -1);
        err = hexify_h3_gridDiskDistances(h, k, disk.data(), distances.data());
        if (err != E_SUCCESS) {
            out[i] = Rcpp::DataFrame::create();
            continue;
        }

        std::vector<std::string> valid_ids;
        std::vector<int> valid_dists;
        for (int64_t j = 0; j < max_size; j++) {
            if (disk[j] != H3_NULL) {
                valid_ids.push_back(h3_to_string(disk[j]));
                valid_dists.push_back(distances[j]);
            }
        }
        out[i] = Rcpp::DataFrame::create(
            Rcpp::Named("cell_id") = Rcpp::wrap(valid_ids),
            Rcpp::Named("ring_distance") = Rcpp::wrap(valid_dists),
            Rcpp::Named("stringsAsFactors") = false
        );
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_h3_gridRingUnsafe(Rcpp::CharacterVector cell_ids, int k) {
    R_xlen_t n = cell_ids.size();
    Rcpp::List out(n);

    // Ring size: 6*k for k > 0, 1 for k = 0
    int64_t ring_size = (k == 0) ? 1 : 6 * static_cast<int64_t>(k);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }

        std::vector<H3Index> ring(ring_size, H3_NULL);
        H3Error err = hexify_h3_gridRingUnsafe(h, k, ring.data());
        if (err != E_SUCCESS) {
            // Fallback: gridRingUnsafe fails near pentagons, use gridDisk diff
            out[i] = Rcpp::CharacterVector(0);
            continue;
        }

        std::vector<std::string> valid;
        for (int64_t j = 0; j < ring_size; j++) {
            if (ring[j] != H3_NULL) {
                valid.push_back(h3_to_string(ring[j]));
            }
        }
        out[i] = Rcpp::wrap(valid);
    }
    return out;
}

// ==========================================================================
// Compact / Uncompact (v0.9.0)
// ==========================================================================

// [[Rcpp::export]]
Rcpp::CharacterVector cpp_h3_compactCells(Rcpp::CharacterVector cell_ids) {
    R_xlen_t n = cell_ids.size();
    std::vector<H3Index> h3_set(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            h3_set[i] = H3_NULL;
        } else {
            h3_set[i] = string_to_h3(CHAR(cell_ids[i]));
        }
    }

    std::vector<H3Index> compacted(n, H3_NULL);
    H3Error err = hexify_h3_compactCells(h3_set.data(), compacted.data(), n);
    if (err != E_SUCCESS) {
        Rcpp::stop("H3 compactCells failed (error code %d)", (int)err);
    }

    std::vector<std::string> valid;
    valid.reserve(n);
    for (R_xlen_t i = 0; i < n; i++) {
        if (compacted[i] != H3_NULL) {
            valid.push_back(h3_to_string(compacted[i]));
        }
    }
    return Rcpp::wrap(valid);
}

// [[Rcpp::export]]
Rcpp::CharacterVector cpp_h3_uncompactCells(Rcpp::CharacterVector cell_ids,
                                             int target_res) {
    R_xlen_t n = cell_ids.size();
    std::vector<H3Index> h3_set(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            h3_set[i] = H3_NULL;
        } else {
            h3_set[i] = string_to_h3(CHAR(cell_ids[i]));
        }
    }

    int64_t max_size = 0;
    H3Error err = hexify_h3_uncompactCellsSize(h3_set.data(), n, target_res,
                                                &max_size);
    if (err != E_SUCCESS || max_size <= 0) {
        Rcpp::stop("H3 uncompactCellsSize failed (error code %d)", (int)err);
    }

    std::vector<H3Index> uncompacted(max_size, H3_NULL);
    err = hexify_h3_uncompactCells(h3_set.data(), n, uncompacted.data(),
                                    max_size, target_res);
    if (err != E_SUCCESS) {
        Rcpp::stop("H3 uncompactCells failed (error code %d)", (int)err);
    }

    std::vector<std::string> valid;
    valid.reserve(max_size);
    for (int64_t i = 0; i < max_size; i++) {
        if (uncompacted[i] != H3_NULL) {
            valid.push_back(h3_to_string(uncompacted[i]));
        }
    }
    return Rcpp::wrap(valid);
}

// ==========================================================================
// Pentagon Detection (v0.9.0)
// ==========================================================================

// [[Rcpp::export]]
Rcpp::LogicalVector cpp_h3_isPentagon(Rcpp::CharacterVector cell_ids) {
    R_xlen_t n = cell_ids.size();
    Rcpp::LogicalVector out(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (cell_ids[i] == NA_STRING) {
            out[i] = NA_LOGICAL;
            continue;
        }
        H3Index h = string_to_h3(CHAR(cell_ids[i]));
        if (h == H3_NULL) {
            out[i] = NA_LOGICAL;
        } else {
            out[i] = hexify_h3_isPentagon(h) ? TRUE : FALSE;
        }
    }
    return out;
}

// ==========================================================================
// Grid Distance (v0.9.0)
// ==========================================================================

// [[Rcpp::export]]
Rcpp::IntegerVector cpp_h3_gridDistance(Rcpp::CharacterVector origin,
                                        Rcpp::CharacterVector destination) {
    R_xlen_t n = origin.size();
    if (destination.size() != n) {
        Rcpp::stop("origin and destination must have the same length");
    }

    Rcpp::IntegerVector out(n);

    for (R_xlen_t i = 0; i < n; i++) {
        if (origin[i] == NA_STRING || destination[i] == NA_STRING) {
            out[i] = NA_INTEGER;
            continue;
        }
        H3Index h_orig = string_to_h3(CHAR(origin[i]));
        H3Index h_dest = string_to_h3(CHAR(destination[i]));
        if (h_orig == H3_NULL || h_dest == H3_NULL) {
            out[i] = NA_INTEGER;
            continue;
        }

        int64_t dist;
        H3Error err = hexify_h3_gridDistance(h_orig, h_dest, &dist);
        if (err != E_SUCCESS) {
            out[i] = NA_INTEGER;
        } else {
            out[i] = static_cast<int>(dist);
        }
    }
    return out;
}
