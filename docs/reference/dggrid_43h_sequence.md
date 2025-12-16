# Create DGGRID 43H aperture sequence

Create an aperture sequence following DGGRID's 43H pattern: first
num_ap4 resolutions use aperture 4, then aperture 3.

## Usage

``` r
dggrid_43h_sequence(num_ap4, num_ap3)
```

## Arguments

- num_ap4:

  Number of aperture-4 resolutions

- num_ap3:

  Number of aperture-3 resolutions

## Value

Integer vector of aperture sequence

## Examples

``` r
# DGGRID 43H with 2 ap4 resolutions, then 3 ap3 resolutions
seq <- dggrid_43h_sequence(2, 3)  # c(4, 4, 3, 3, 3)
```
