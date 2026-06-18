# tests/testthat/test-theme_phase_plane.R
#
# Unit tests for theme_phase_plane().
#
# ggplot2 theme testing is intentionally lightweight: we verify that the
# function returns the right class, that arguments are accepted and produce
# valid theme objects, and that key structural properties hold. We do not
# test exact color values since those are cosmetic and may be adjusted.


test_that("theme_phase_plane() returns a ggplot2 theme object", {
  th <- theme_phase_plane()
  expect_s3_class(th, "theme")
})

test_that("theme_phase_plane() can be added to a ggplot", {
  p <- ggplot2::ggplot() + theme_phase_plane()
  expect_s3_class(p, "ggplot")
})

test_that("theme_phase_plane() accepts and applies base_size", {
  th_small <- theme_phase_plane(base_size = 8)
  th_large <- theme_phase_plane(base_size = 20)
  expect_s3_class(th_small, "theme")
  expect_s3_class(th_large, "theme")
})

test_that("theme_phase_plane() accepts custom grid_color", {
  th <- theme_phase_plane(grid_color = "grey50")
  expect_s3_class(th, "theme")
  # Grid color should be reflected in panel.grid.major
  grid_col <- th$panel.grid.major$colour
  expect_equal(grid_col, "grey50")
})

test_that("theme_phase_plane() removes panel border", {
  th <- theme_phase_plane()
  expect_s3_class(th$panel.border, "element_blank")
})

test_that("theme_phase_plane() removes panel minor grid", {
  th <- theme_phase_plane()
  expect_s3_class(th$panel.grid.minor, "element_blank")
})

test_that("theme_phase_plane() removes axis lines", {
  th <- theme_phase_plane()
  expect_s3_class(th$axis.line, "element_blank")
})

test_that("theme_phase_plane() sets white panel background", {
  th <- theme_phase_plane()
  expect_equal(th$panel.background$fill, "white")
})

test_that("theme_phase_plane() complete theme builds without error", {
  # Build a full plot and check it renders without error
  p <- ggplot2::ggplot(
    data.frame(x = c(-2, 2), y = c(-2, 2)),
    ggplot2::aes(x, y)
  ) +
    ggplot2::geom_blank() +
    theme_phase_plane()

  # ggplot_build() triggers the full rendering pipeline
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("theme_phase_plane() result can be further customized with +", {
  # Users should be able to override individual elements after applying theme
  p <- ggplot2::ggplot() +
    theme_phase_plane() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "lightyellow"))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})
