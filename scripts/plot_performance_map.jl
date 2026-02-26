#!/usr/bin/env julia

using TurboMachineModel
using Plots

const TM = TurboMachineModel.Physics.Turbomachine.Compressor

"""
    plot_performance_map(map::TM.TabulatedCompressorPerformanceMap; title_prefix="Performance Map")

Create a 2-panel contour plot for `PR` and `eta` over corrected flow/speed.
"""
function plot_performance_map(map::TM.TabulatedCompressorPerformanceMap; title_prefix::String="Performance Map")
    p_pr = contour(
        map.mdot_corr_grid,
        map.omega_corr_grid,
        map.pr_table;
        xlabel="mdot_corr",
        ylabel="omega_corr",
        title="$(title_prefix): Pressure Ratio",
        colorbar_title="PR",
        linewidth=2,
    )
    scatter!(p_pr, vec(repeat(map.mdot_corr_grid', length(map.omega_corr_grid))), vec(repeat(map.omega_corr_grid, 1, length(map.mdot_corr_grid))); ms=2, label=false)

    p_eta = contour(
        map.mdot_corr_grid,
        map.omega_corr_grid,
        map.eta_table;
        xlabel="mdot_corr",
        ylabel="omega_corr",
        title="$(title_prefix): Isentropic Efficiency",
        colorbar_title="eta",
        linewidth=2,
    )
    scatter!(p_eta, vec(repeat(map.mdot_corr_grid', length(map.omega_corr_grid))), vec(repeat(map.omega_corr_grid, 1, length(map.mdot_corr_grid))); ms=2, label=false)

    return plot(p_pr, p_eta; layout=(1, 2), size=(1100, 450))
end

function _main()
    output_path = length(ARGS) >= 1 ? ARGS[1] : "compressor_performance_map.png"
    map = TM.demo_compressor_performance_map()
    title_prefix = "Demo Compressor Map"

    fig = plot_performance_map(map; title_prefix=title_prefix)
    savefig(fig, output_path)
    println("Saved map plot to: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main()
end
