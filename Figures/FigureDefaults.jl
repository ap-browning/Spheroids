#=
    Figure Defaults
=#

gr()
default()
default(
    framestyle      =:axis,
    grid            =:y,
    axis            =:x,
    tick_direction  = :out,
    foreground_color_border = "#aaa",
    foreground_color_axis = "#aaa",
    msw=0.0,
    lw=2.0,
    fontfamily      ="Helvetica",
    guidefontsize   = 9,
    titlefontsize   = 10,
)

# 3D defaults (need to include manually)
defaults_3D = (grid=:all, axis=:all, axes=:all, box=:on, camera=(40,50))

# Colorscheme
colors = Dict(
        "983b" => Dict(2500 => "#fdbb84", 5000 => "#ef6548", 10000 => "#990000"),
        "793b" => Dict(2500 => "#99d8c9", 5000 => "#41ae76", 10000 => "#00441b")
    )