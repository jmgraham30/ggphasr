# Phase plane ggplot2 theme

A clean
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
designed for phase plane and phase portrait plots. Key features:

## Usage

``` r
theme_phase_plane(
  base_size = 13,
  base_family = "",
  grid_color = "grey88",
  axis_text_color = "grey30"
)
```

## Arguments

- base_size:

  Numeric. Base font size in points. Default: `13`.

- base_family:

  Character. Base font family. Default: `""` (device default).

- grid_color:

  Character. Color of the major grid lines. Default: `"grey88"`.

- axis_text_color:

  Character. Color of axis tick labels. Default: `"grey30"`.

## Value

A
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
object that can be added to any ggplot with `+`.

## Details

- White panel background for clean printing and projection

- Light grey major grid lines to help read off coordinates

- No panel border (axis lines through the origin are added separately as
  [`ggplot2::geom_hline()`](https://ggplot2.tidyverse.org/reference/geom_abline.html)
  /
  [`ggplot2::geom_vline()`](https://ggplot2.tidyverse.org/reference/geom_abline.html)
  layers by the `gg_*` plotting functions)

- Slightly larger base font size than the ggplot2 default, suitable for
  classroom projection

## See also

[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html),
[`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)

## Examples

``` r
library(ggplot2)

# Apply to any ggplot
ggplot(data.frame(x = c(-2, 2), y = c(-2, 2)), aes(x, y)) +
  geom_blank() +
  theme_phase_plane()


# Override font size for a smaller plot
ggplot(data.frame(x = c(-2, 2), y = c(-2, 2)), aes(x, y)) +
  geom_blank() +
  theme_phase_plane(base_size = 10)


# Customize grid color
ggplot(data.frame(x = c(-2, 2), y = c(-2, 2)), aes(x, y)) +
  geom_blank() +
  theme_phase_plane(grid_color = "grey75")

```
