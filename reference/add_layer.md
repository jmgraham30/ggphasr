# Add ggplot2 layers to a ggphasr_result object

Adds a ggplot2 layer, scale, or theme element to the `$plot` component
of a
[`gg_phase_plane()`](https://jmgraham30.github.io/ggphasr/reference/gg_phase_plane.md)
result, returning an updated `ggphasr_result` object. This is the
recommended way to further customize a phase plane plot while preserving
the `$equilibria` table:

## Usage

``` r
add_layer(result, layer)
```

## Arguments

- result:

  A `ggphasr_result` object from
  [`gg_phase_plane()`](https://jmgraham30.github.io/ggphasr/reference/gg_phase_plane.md).

- layer:

  A ggplot2 layer, scale, theme, or list thereof.

## Value

A `ggphasr_result` object with the updated `$plot`.

## Details

    result <- gg_phase_plane(ode_lotka_volterra, ...)

    # Use add_layer() to customize and keep the ggphasr_result structure:
    result2 <- add_layer(result, ggplot2::labs(title = "My plot"))
    result2$plot        # updated plot
    result2$equilibria  # equilibria preserved

    # Or extract $plot and add layers with + directly:
    result$plot + ggplot2::labs(title = "My plot")
