set_theme!(
    # rowgap=5,
    # colgap=5,
    fontsize=16,
    # font=:dejavu,
    Axis=(
        xgridcolor=:black,
        ygridcolor=:black,
        backgroundcolor=:white,
        xlabelsize=22, ylabelsize=22,
        xgridwidth=0.15, ygridwidth=0.15,
        xticksmirrored=true, yticksmirrored=true,
        xminorticksvisible=false, yminorticksvisible=false, xminortickalign=1, yminortickalign=1,
        xtickalign=1, ytickalign=1,
        xgridvisible=true, ygridvisible=true),
    Colorbar=(
        topspinevisible=false,
        rightspinevisible=false,
        bottomspinevisible=false,
        leftspinevisible=false,
        width=12,
        height=Relative(1),
        tickalign=1,
        labelpadding=-5),
    Legend=(
        tellheight=false,
        tellwidth=false,
        halign=:left,
        valign=:bottom,
        labelsize=14,
        linewidth=2,
        margin=(10, 10, 10, 10),
        framevisible=true)
)