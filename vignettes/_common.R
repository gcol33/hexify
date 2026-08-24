# Figure defaults shared by every vignette.
#
# A figure is shown at the full width of the reader's viewport and scales to it,
# so its text arrives at the authored point size times viewport width over
# figure width. A 390 px phone showing a 7 inch figure scales by about 0.77,
# which is what puts the device default of 12 pt under the size the text stays
# legible at. Sizing against that viewport keeps every figure readable, and the
# site's stylesheet scales the text back down as the column widens.
#
# Figures are one width throughout, so a single point size serves all of them.
# A chunk that needs a different shape sets fig.height alone.

FIG_WIDTH <- 7
FIG_HEIGHT <- 5
FIG_POINTSIZE <- 16
FIG_BASE_SIZE <- 18

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = FIG_WIDTH,
  fig.height = FIG_HEIGHT,
  dev.args = list(pointsize = FIG_POINTSIZE)
)

# ggplot text is absolute points carried by the theme, so the device point size
# above does not reach it. hexify's own plots read this option for their base.
old_base_size <- options(hexify.base_size = FIG_BASE_SIZE)

# Put the option back once the article has been knit.
knitr::knit_hooks$set(document = function(x) { options(old_base_size); x })

# Neutral grey for secondary strokes and annotations. This is the grey of the
# site figure palette, which the stylesheet swaps for a lighter tint in dark
# mode; the fixed R greys stay invisible against the dark background.
GREY <- "#5C6166"
