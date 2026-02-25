#!/usr/bin/env julia

using TurboMachineModel
using Plots

const TM = TurboMachineModel.Physics.Turbomachine

"""
    plot_performance_map(map::TM.TabulatedPerformanceMap; title_prefix="Performance Map")

Create a 2-panel contour plot for `PR` and `eta` over corrected flow/speed.
"""
function plot_performance_map(map::TM.TabulatedPerformanceMap; title_prefix::String="Performance Map")
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
    map_name = length(ARGS) >= 1 ? lowercase(ARGS[1]) : "compressor"
    output_path = length(ARGS) >= 2 ? ARGS[2] : "performance_map.png"

    map, title_prefix =
        if map_name == "compressor"
            TM.demo_compressor_performance_map(), "Demo Compressor Map"
        elseif map_name == "turbine"
            TM.demo_turbine_performance_map(), "Demo Turbine Map"
        else
            error("Unknown map kind '$map_name'. Use 'compressor' or 'turbine'.")
        end

    fig = plot_performance_map(map; title_prefix=title_prefix)
    savefig(fig, output_path)
    println("Saved map plot to: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main()
end
