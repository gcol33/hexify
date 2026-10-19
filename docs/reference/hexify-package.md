# hexify: Equal-Area Hex Grids on the Snyder ISEA Icosahedron

Provides functions to build and use hexagonal discrete global grids
using the Icosahedral Snyder Equal Area ('ISEA') projection (Snyder 1992
[doi:10.3138/27H7-8K88-4882-1752](https://doi.org/10.3138/27H7-8K88-4882-1752)
) and the 'H3' hierarchical hexagonal system (Uber Technologies).
Implements the 'ISEA' discrete global grid system (Sahr, White and
Kimerling 2003
[doi:10.1559/152304003100011090](https://doi.org/10.1559/152304003100011090)
) and the 'Hex9' octahedral grid (Griffin 2026
[arXiv:2608.00022](https://arxiv.org/abs/2608.00022)). Includes a fast
C++ core for 'ISEA' projection and aperture quantization, an included
'H3' v4.4.1 C library for native 'H3' grid operations, and
'sf'/'terra'-compatible R wrappers for grid generation and coordinate
assignment. Output is compatible with 'dggridR' for interoperability.

## See also

Useful links:

- <https://gillescolling.com/hexify/>

- Report bugs at <https://github.com/gcol33/hexify/issues>

## Author

**Maintainer**: Gilles Colling <gilles.colling051@gmail.com>
([ORCID](https://orcid.org/0000-0003-3070-6066)) \[copyright holder\]

Authors:

- Gilles Colling <gilles.colling051@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3070-6066)) \[copyright holder\]

Other contributors:

- Uber Technologies, Inc. (H3 C library, src/h3 (Apache-2.0))
  \[contributor, copyright holder\]

- Ben Griffin (Hex9 address tables in src/hex9.cpp, from libhex9
  (Apache-2.0)) \[contributor, copyright holder\]
