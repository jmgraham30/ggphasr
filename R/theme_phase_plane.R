# theme_phase_plane.R
#
# A ggplot2 theme designed for phase plane and phase portrait plots.
#
# Design goals:
#   - White background for clean printing and projection
#   - Light grey grid lines to help students read coordinates
#   - No panel border (axes through the origin replace the border)
#   - No default axis lines on the panel edges (origin lines added as layers)
#   - Comfortable text sizing for classroom projection
#
# The axis lines through the origin (x = 0 and y = 0) are NOT part of the
# theme itself — they are added as geom_hline() / geom_vline() layers by
# gg_flow_field(), gg_nullclines(), and gg_phase_portrait(). This means
# users can remove them with standard ggplot2 layer removal if they prefer
# border-style axes.


#' Phase plane ggplot2 theme
#'
#' A clean [ggplot2::theme()] designed for phase plane and phase portrait
#' plots. Key features:
#'
#' - White panel background for clean printing and projection
#' - Light grey major grid lines to help read off coordinates
#' - No panel border (axis lines through the origin are added separately
#'   as [ggplot2::geom_hline()] / [ggplot2::geom_vline()] layers by the
#'   `gg_*` plotting functions)
#' - Slightly larger base font size than the ggplot2 default, suitable
#'   for classroom projection
#'
#' @param base_size Numeric. Base font size in points. Default: `13`.
#' @param base_family Character. Base font family. Default: `""` (device
#'   default).
#' @param grid_color Character. Color of the major grid lines.
#'   Default: `"grey88"`.
#' @param axis_text_color Character. Color of axis tick labels.
#'   Default: `"grey30"`.
#'
#' @return A [ggplot2::theme()] object that can be added to any ggplot with
#'   `+`.
#'
#' @examples
#' library(ggplot2)
#'
#' # Apply to any ggplot
#' ggplot(data.frame(x = c(-2, 2), y = c(-2, 2)), aes(x, y)) +
#'   geom_blank() +
#'   theme_phase_plane()
#'
#' # Override font size for a smaller plot
#' ggplot(data.frame(x = c(-2, 2), y = c(-2, 2)), aes(x, y)) +
#'   geom_blank() +
#'   theme_phase_plane(base_size = 10)
#'
#' # Customize grid color
#' ggplot(data.frame(x = c(-2, 2), y = c(-2, 2)), aes(x, y)) +
#'   geom_blank() +
#'   theme_phase_plane(grid_color = "grey75")
#'
#' @seealso [ggplot2::theme()], [ggplot2::theme_bw()]
#' @export
theme_phase_plane <- function(base_size      = 13,
                               base_family    = "",
                               grid_color     = "grey88",
                               axis_text_color = "grey30") {

  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(

      # ── Panel ────────────────────────────────────────────────────────────
      # White background; no border (axis lines through origin replace it)
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border     = ggplot2::element_blank(),

      # Light grey major grid lines; no minor grid lines
      panel.grid.major = ggplot2::element_line(color = grid_color,
                                               linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),

      # ── Axes ─────────────────────────────────────────────────────────────
      # No axis lines on the panel edges (origin lines added as layers)
      axis.line  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = axis_text_color,
                                         linewidth = 0.4),
      axis.text  = ggplot2::element_text(color = axis_text_color,
                                         size  = ggplot2::rel(0.85)),
      axis.title = ggplot2::element_text(color = "grey20",
                                         size  = ggplot2::rel(0.95)),

      # ── Plot labels ───────────────────────────────────────────────────────
      plot.title    = ggplot2::element_text(color    = "grey10",
                                            size     = ggplot2::rel(1.1),
                                            face     = "bold",
                                            margin   = ggplot2::margin(
                                              b = 6)),
      plot.subtitle = ggplot2::element_text(color  = "grey30",
                                            size   = ggplot2::rel(0.9),
                                            margin = ggplot2::margin(b = 4)),
      plot.caption  = ggplot2::element_text(color  = "grey50",
                                            size   = ggplot2::rel(0.75),
                                            hjust  = 1),

      # ── Legend ───────────────────────────────────────────────────────────
      legend.background = ggplot2::element_rect(fill  = "white",
                                                color = "grey80",
                                                linewidth = 0.3),
      legend.key        = ggplot2::element_rect(fill = "white"),
      legend.title      = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.text       = ggplot2::element_text(size = ggplot2::rel(0.85)),

      # ── Spacing ──────────────────────────────────────────────────────────
      plot.margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
    )
}
