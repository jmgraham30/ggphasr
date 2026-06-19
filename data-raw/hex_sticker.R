# data-raw/hex_sticker.R

library(ggplot2)
library(hexSticker)
library(showtext)
library(ggphasr)

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------

col_bg     <- "#1a5276"
col_border <- "#aed6f1"
col_cycle  <- "#ffffff"
col_inner  <- "#aed6f1"
col_outer  <- "#f9e79f"
col_text   <- "#ffffff"
col_url    <- "#aed6f1"


# ---------------------------------------------------------------------------
# Step 1: Integrate trajectories
# ---------------------------------------------------------------------------
#
# mu = 0.5 gives a nearly circular, smooth limit cycle that spans
# approximately x in [-2.05, 2.05], y in [-2.1, 2.1].
# This fits comfortably in a tight window with no clipping artifacts.

vdp_params <- c(mu = 0.5)
norm_vdp   <- ggphasr:::.normalize_ode(ode_van_der_pol, "two.dim")

# Limit cycle: integrate from a point near the attractor for a long time,
# then keep only the final two periods (~12 time units each for mu=0.5)
traj_cycle_full <- ggphasr:::.integrate_trajectory(
  norm_vdp, y0 = c(2.0, 0.0), t_start = 0, t_end = 60,
  t_steps = 6000L, parameters = vdp_params
)
# Keep only the final settled portion — guarantees a clean closed curve
n_cyc      <- nrow(traj_cycle_full)
traj_cycle <- traj_cycle_full[round(n_cyc * 0.85):n_cyc, ]

# Inner spiral: starts near origin, spirals out toward limit cycle
traj_inner <- ggphasr:::.integrate_trajectory(
  norm_vdp, y0 = c(0.4, 0.0), t_start = 0, t_end = 20,
  t_steps = 2000L, parameters = vdp_params
)

# Outer spiral: starts just outside the limit cycle, spirals inward
# For mu=0.5 the limit cycle reaches ~x=2.05, so start at x=2.6
traj_outer <- ggphasr:::.integrate_trajectory(
  norm_vdp, y0 = c(2.6, 0.0), t_start = 0, t_end = 10,
  t_steps = 1000L, parameters = vdp_params
)

# Axis limits that contain all curves with a small margin —
# no hard clipping needed since mu=0.5 keeps everything small
lim_x <- c(-2.8, 2.8)
lim_y <- c(-2.8, 2.8)


# ---------------------------------------------------------------------------
# Step 2: Build the subplot
# ---------------------------------------------------------------------------

sticker_plot <- ggplot() +
  # Outer spiral (gold)
  geom_path(data      = traj_outer,
            mapping   = aes(x = x, y = y),
            color     = col_outer,
            linewidth = 0.7,
            alpha     = 0.9,
            lineend   = "round",
            linejoin  = "round") +
  # Inner spiral (light blue)
  geom_path(data      = traj_inner,
            mapping   = aes(x = x, y = y),
            color     = col_inner,
            linewidth = 0.7,
            alpha     = 0.9,
            lineend   = "round",
            linejoin  = "round") +
  # Limit cycle (white, dominant)
  geom_path(data      = traj_cycle,
            mapping   = aes(x = x, y = y),
            color     = col_cycle,
            linewidth = 1.1,
            lineend   = "round",
            linejoin  = "round") +
  # Use coord_cartesian (not coord_fixed) so hexSticker can scale freely;
  # the ratio is enforced by equal lim_x and lim_y ranges
  coord_cartesian(xlim   = lim_x,
                  ylim   = lim_y,
                  expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA)
  )


# ---------------------------------------------------------------------------
# Step 3: Font
# ---------------------------------------------------------------------------

font_add_google("Nunito", "nunito")
showtext_auto()


# ---------------------------------------------------------------------------
# Step 4: Assemble sticker
# ---------------------------------------------------------------------------

sticker(
  subplot    = sticker_plot,
  package    = "ggphasr",

  # Subplot: centered, moderate size so curves stay well inside the hex
  s_x        = 1.00,
  s_y        = 1.15,
  s_width    = 1.35,
  s_height   = 1.35,

  # Package name: lower band, clear of curves
  p_x        = 1.00,
  p_y        = 0.33,
  p_size     = 17,
  p_color    = col_text,
  p_family   = "nunito",
  p_fontface = "bold",

  # Hex
  h_fill     = col_bg,
  h_color    = col_border,
  h_size     = 1.4,

  # URL
  url        = "github.com/jmgraham30/ggphasr",
  u_size     = 3.5,
  u_color    = col_url,
  u_family   = "nunito",

  filename   = "man/figures/logo.png",
  dpi        = 600
)

showtext_auto(FALSE)
message("Hex sticker saved to man/figures/logo.png")
