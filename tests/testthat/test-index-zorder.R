
# tests/testthat/test-index_zorder.R
# Comprehensive tests for Z-Order indexing (generic, all apertures)

# =============================================================================
# Helper Functions
# =============================================================================

get_max_coord <- function(aperture, resolution) {
  # Calculate maximum valid coordinate for given aperture/resolution
  if (aperture == 3) {
    # For aperture 3, effective resolution
    eff_res <- (resolution + 1) %/% 2
    return(3^eff_res - 1)
  } else if (aperture == 4) {
    # For aperture 4, uses radix-2 (binary)
    return(2^resolution - 1)
  } else if (aperture == 7) {
    # For aperture 7, uses radix-7
    return(7^resolution - 1)
  } else {
    return(aperture^resolution - 1)
  }
}

# =============================================================================
# Aperture 3 Z-Order Tests
# =============================================================================

test_that("Z-Order aperture 3: basic encode/decode", {
  test_cases <- list(
    list(face = 0, i = 0, j = 0, res = 1),
    list(face = 0, i = 1, j = 2, res = 1),  # Changed from (1,0) - i%3=1 requires j%3=2
    list(face = 0, i = 0, j = 0, res = 1),  # Changed from (0,1) - invalid for Class II
    list(face = 5, i = 1, j = 2, res = 1),  # Changed from (1,1) - i%3=1 requires j%3=2
    list(face = 10, i = 2, j = 2, res = 2),
    list(face = 19, i = 5, j = 7, res = 3)
  )
  
  for (tc in test_cases) {
    index <- cpp_cell_to_index(tc$face, tc$i, tc$j, tc$res, 3, "zorder")
    result <- cpp_index_to_cell(index, 3, "zorder")
    
    expect_equal(result$face, tc$face, 
                 label = sprintf("face[%d,%d,%d]", tc$i, tc$j, tc$res))
    expect_equal(result$i, tc$i,
                 label = sprintf("i[%d,%d,%d]", tc$i, tc$j, tc$res))
    expect_equal(result$j, tc$j,
                 label = sprintf("j[%d,%d,%d]", tc$i, tc$j, tc$res))
    expect_equal(result$resolution, tc$res,
                 label = sprintf("res[%d,%d,%d]", tc$i, tc$j, tc$res))
  }
})

test_that("Z-Order aperture 3: digit alternation pattern", {
  # Z-Order for aperture 3 uses simple alternation
  # For res=2 (Class I), eff_res=1: alternates i,j digits
  
  # i=1, j=2 → radix-3: "1", "2" → alternation: "12"
  index <- cpp_cell_to_index(0, 1, 2, 2, 3, "zorder")
  expect_match(index, "^0012$")
  
  # i=2, j=1 → radix-3: "2", "1" → alternation: "21"
  index <- cpp_cell_to_index(0, 2, 1, 2, 3, "zorder")
  expect_match(index, "^0021$")
})

test_that("Z-Order aperture 3: Class II j-digit inference", {
  # For Class II (odd resolution), last j digit is inferred from last i digit
  # Rule: i=0→j=0, i=1→j=2, i=2→j=1
  
  test_cases <- list(
    list(i = 0, j = 0, res = 1),  # Last i=0 → inferred j=0
    list(i = 1, j = 2, res = 1),  # Last i=1 → inferred j=2
    list(i = 2, j = 1, res = 1),  # Last i=2 → inferred j=1
    list(i = 3, j = 6, res = 3),  # i ends in 0 → j ends in 0
    list(i = 4, j = 8, res = 3)   # i ends in 1 → j ends in 2
  )
  
  for (tc in test_cases) {
    index <- cpp_cell_to_index(0, tc$i, tc$j, tc$res, 3, "zorder")
    result <- cpp_index_to_cell(index, 3, "zorder")
    
    expect_equal(result$i, tc$i,
                 label = sprintf("i inference[%d,%d]", tc$i, tc$j))
    expect_equal(result$j, tc$j,
                 label = sprintf("j inference[%d,%d]", tc$i, tc$j))
  }
})

test_that("Z-Order aperture 3: Class I vs Class II string length", {
  # Class I (even res): even number of digits after face
  # Class II (odd res): odd number of digits after face
  
  for (res in 1:6) {
    index <- cpp_cell_to_index(0, 0, 0, res, 3, "zorder")
    digits_after_face <- nchar(index) - 2
    
    if (res %% 2 == 0) {
      # Class I: even number of digits
      expect_true(digits_after_face %% 2 == 0,
                  label = sprintf("Class I res=%d", res))
    } else {
      # Class II: odd number of digits
      expect_true(digits_after_face %% 2 == 1,
                  label = sprintf("Class II res=%d", res))
    }
  }
})

test_that("Z-Order aperture 3: comprehensive roundtrip", {
  # Reduced test set for speed
  test_cases <- expand.grid(
    face = c(0, 10),
    i = c(0, 2, 5),
    j = c(0, 2, 5),
    res = 1:2
  )

  for (idx in seq_len(nrow(test_cases))) {
    tc <- test_cases[idx, ]

    # Skip invalid coordinates
    max_coord <- get_max_coord(3, tc$res)
    if (tc$i > max_coord || tc$j > max_coord) next

    # Skip invalid Class II coordinates
    if (tc$res %% 2 == 1) {  # Class II (odd resolution)
      j_required <- c(0, 2, 1)[tc$i %% 3 + 1]
      if (tc$j %% 3 != j_required) next
    }

    # Roundtrip test
    index <- cpp_cell_to_index(tc$face, tc$i, tc$j, tc$res, 3, "zorder")
    result <- cpp_index_to_cell(index, 3, "zorder")

    expect_equal(result$face, tc$face, label = sprintf("idx=%d", idx))
    expect_equal(result$i, tc$i, label = sprintf("idx=%d", idx))
    expect_equal(result$j, tc$j, label = sprintf("idx=%d", idx))
    expect_equal(result$resolution, tc$res, label = sprintf("idx=%d", idx))
  }
})

# =============================================================================
# Aperture 4 Z-Order Tests
# =============================================================================

test_that("Z-Order aperture 4: bit interleaving formula", {
  # DGGRID formula: digit = i_bit * 2 + j_bit
  # This is tested by checking known cases
  
  test_cases <- list(
    # i=0 (binary: 0), j=0 (binary: 0) → digit: 0
    list(i = 0, j = 0, res = 1, expected = "000"),
    
    # i=1 (binary: 1), j=0 (binary: 0) → digit: 2
    list(i = 1, j = 0, res = 1, expected = "002"),
    
    # i=0 (binary: 0), j=1 (binary: 1) → digit: 1
    list(i = 0, j = 1, res = 1, expected = "001"),
    
    # i=1 (binary: 1), j=1 (binary: 1) → digit: 3
    list(i = 1, j = 1, res = 1, expected = "003"),
    
    # i=5 (binary: 101), j=6 (binary: 110) → digits: 3,1,2
    list(i = 5, j = 6, res = 3, expected = "00312"),
    
    # i=3 (binary: 11), j=3 (binary: 11) → digits: 3,3
    list(i = 3, j = 3, res = 2, expected = "0033"),
    
    # i=15 (binary: 1111), j=0 (binary: 0000) → digits: 2,2,2,2
    list(i = 15, j = 0, res = 4, expected = "002222")
  )
  
  for (tc in test_cases) {
    index <- cpp_cell_to_index(0, tc$i, tc$j, tc$res, 4, "zorder")
    expect_equal(index, tc$expected,
                 label = sprintf("i=%d,j=%d,res=%d", tc$i, tc$j, tc$res))
  }
})

test_that("Z-Order aperture 4: basic encode/decode", {
  test_cases <- list(
    list(face = 0, i = 0, j = 0, res = 1),
    list(face = 5, i = 1, j = 2, res = 2),
    list(face = 10, i = 5, j = 6, res = 3),
    list(face = 11, i = 10, j = 15, res = 4)  # Changed from face 19 to 11
  )
  
  for (tc in test_cases) {
    index <- cpp_cell_to_index(tc$face, tc$i, tc$j, tc$res, 4, "zorder")
    result <- cpp_index_to_cell(index, 4, "zorder")
    
    expect_equal(result$face, tc$face)
    expect_equal(result$i, tc$i)
    expect_equal(result$j, tc$j)
    expect_equal(result$resolution, tc$res)
  }
})

test_that("Z-Order aperture 4: comprehensive roundtrip", {
  # Reduced test set for speed
  test_cases <- expand.grid(
    face = c(0, 10),
    i = c(0, 3, 7),
    j = c(0, 3, 7),
    res = 1:3
  )

  for (idx in seq_len(nrow(test_cases))) {
    tc <- test_cases[idx, ]

    # Skip invalid coordinates
    max_coord <- get_max_coord(4, tc$res)
    if (tc$i > max_coord || tc$j > max_coord) next

    # Roundtrip test
    index <- cpp_cell_to_index(tc$face, tc$i, tc$j, tc$res, 4, "zorder")
    result <- cpp_index_to_cell(index, 4, "zorder")

    expect_equal(result$face, tc$face, label = sprintf("idx=%d", idx))
    expect_equal(result$i, tc$i, label = sprintf("idx=%d", idx))
    expect_equal(result$j, tc$j, label = sprintf("idx=%d", idx))
    expect_equal(result$resolution, tc$res, label = sprintf("idx=%d", idx))
  }
})

# =============================================================================
# Aperture 7 Z-Order Tests
# =============================================================================

test_that("Z-Order aperture 7: basic encode/decode", {
  test_cases <- list(
    list(face = 0, i = 0, j = 0, res = 1),
    list(face = 5, i = 1, j = 2, res = 1),
    list(face = 10, i = 3, j = 5, res = 2),
    list(face = 11, i = 10, j = 20, res = 2)  # Changed from face 19 to 11
  )
  
  for (tc in test_cases) {
    index <- cpp_cell_to_index(tc$face, tc$i, tc$j, tc$res, 7, "zorder")
    result <- cpp_index_to_cell(index, 7, "zorder")
    
    expect_equal(result$face, tc$face)
    expect_equal(result$i, tc$i)
    expect_equal(result$j, tc$j)
    expect_equal(result$resolution, tc$res)
  }
})

test_that("Z-Order aperture 7: comprehensive roundtrip", {
  # Reduced test set for speed
  test_cases <- expand.grid(
    face = c(0, 10),
    i = c(0, 3, 6),
    j = c(0, 3, 6),
    res = 1:2
  )

  for (idx in seq_len(nrow(test_cases))) {
    tc <- test_cases[idx, ]

    # Skip invalid coordinates
    max_coord <- get_max_coord(7, tc$res)
    if (tc$i > max_coord || tc$j > max_coord) next

    # Roundtrip test
    index <- cpp_cell_to_index(tc$face, tc$i, tc$j, tc$res, 7, "zorder")
    result <- cpp_index_to_cell(index, 7, "zorder")

    expect_equal(result$face, tc$face, label = sprintf("idx=%d", idx))
    expect_equal(result$i, tc$i, label = sprintf("idx=%d", idx))
    expect_equal(result$j, tc$j, label = sprintf("idx=%d", idx))
    expect_equal(result$resolution, tc$res, label = sprintf("idx=%d", idx))
  }
})

# =============================================================================
# Cross-Aperture Tests
# =============================================================================

test_that("Z-Order: resolution 0 works for all apertures", {
  for (face in c(0, 5, 10, 11)) {  # Changed from c(0, 5, 10, 19) to only valid faces
    for (ap in c(3, 4, 7)) {
      index <- cpp_cell_to_index(face, 0, 0, 0, ap, "zorder")
      
      # Should be just face number
      expect_equal(nchar(index), 2)
      expect_equal(substr(index, 1, 2), sprintf("%02d", face))
      
      # Decode should work
      result <- cpp_index_to_cell(index, ap, "zorder")
      expect_equal(result$face, face)
      expect_equal(result$resolution, 0)
      expect_equal(result$i, 0)
      expect_equal(result$j, 0)
    }
  }
})

test_that("Z-Order: face encoding correct for all apertures", {
  for (ap in c(3, 4, 7)) {
    # Single digit faces
    expect_equal(substr(cpp_cell_to_index(0, 0, 0, 1, ap, "zorder"), 1, 2), "00")
    expect_equal(substr(cpp_cell_to_index(5, 0, 0, 1, ap, "zorder"), 1, 2), "05")
    expect_equal(substr(cpp_cell_to_index(9, 0, 0, 1, ap, "zorder"), 1, 2), "09")
    
    # Double digit faces
    expect_equal(substr(cpp_cell_to_index(10, 0, 0, 1, ap, "zorder"), 1, 2), "10")
    # Face 19 is only valid for aperture 3 at resolution > 0
    if (ap == 3) {
      expect_equal(substr(cpp_cell_to_index(19, 0, 0, 1, ap, "zorder"), 1, 2), "19")
    }
  }
})

test_that("Z-Order: index comparison works", {
  idx1 <- cpp_cell_to_index(0, 0, 0, 1, 4, "zorder")
  idx2 <- cpp_cell_to_index(0, 1, 0, 1, 4, "zorder")
  idx3 <- cpp_cell_to_index(1, 0, 0, 1, 4, "zorder")
  
  expect_equal(cpp_compare_indices(idx1, idx1), 0)
  expect_true(cpp_compare_indices(idx1, idx2) != 0)
  expect_true(cpp_compare_indices(idx1, idx3) != 0)
})

test_that("Z-Order: parent/child relationships", {
  # Test with valid coordinates for each aperture
  test_cases <- list(
    list(ap = 3, i = 5, j = 4, res = 3),  # Valid: i%3=2, j%3=1 ✓ (was 5,5)
    list(ap = 4, i = 5, j = 5, res = 3),  # Valid: ap4 has no Class II constraint
    list(ap = 7, i = 0, j = 0, res = 2)   # Changed to (0,0) which is always valid
  )
  
  for (tc in test_cases) {
    # Create child at higher resolution
    child_idx <- cpp_cell_to_index(0, tc$i, tc$j, tc$res, tc$ap, "zorder")
    
    # Get parent
    parent_idx <- cpp_get_parent_index(child_idx, tc$ap, "zorder")
    
    # Parent should be prefix of child
    expect_true(substr(child_idx, 1, nchar(parent_idx)) == parent_idx,
                label = sprintf("ap=%d", tc$ap))
    
    # Get children
    children <- cpp_get_children_indices(parent_idx, tc$ap, "zorder")
    
    # Child should be in children list
    expect_true(child_idx %in% children,
                label = sprintf("ap=%d child (%d,%d) in list", tc$ap, tc$i, tc$j))
    
    # Should have correct number of children
    expect_equal(length(children), tc$ap,
                 label = sprintf("ap=%d children count", tc$ap))
  }
})

test_that("Z-Order: resolution extraction works", {
  for (ap in c(3, 4, 7)) {
    for (res in 0:3) {
      index <- cpp_cell_to_index(0, 0, 0, res, ap, "zorder")
      extracted <- cpp_get_index_resolution(index, ap, "zorder")
      expect_equal(extracted, res,
                   label = sprintf("ap=%d, res=%d", ap, res))
    }
  }
})
