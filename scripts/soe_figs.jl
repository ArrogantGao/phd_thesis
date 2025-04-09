include("settings.jl")

using CairoMakie, LaTeXStrings
using JLD

# fig 3.3
begin
    @load "data/kLz_lim.jld" 

    fig = Figure(size = (600, 450), fontsize = 20)
    ax = Axis(fig[1, 1], yticks = ([-5:-3:-17...], [L"$10^{-5}$", L"$10^{-8}$", L"$10^{-11}$", L"$10^{-14}$", L"$10^{-17}$"]), xlabel = L"L_z", ylabel = L"$\mathcal{E}$")
	labels = [L"s = 3,~M = 4", L"s = 4,~M = 8", L"s = 5,~M = 16"]

	for i in 1:3
        l = labels[i]
        ms = markerstyle[i]
        c = colors[i]
        sw = strokewidth
        msw = markerstrokewidth
        m = markersize
		scatter!(ax, Lz_array, log10.(abs.(error_soe_array[i])), label = l, marker = ms, markersize = m, color = c, strokecolor = strokecolor, strokewidth = sw)
	end

	lines!(ax, Lz_array, log10.(abs.(error_dir_array)), label = L"\text{Ewald2D}", linewidth = linewidth, color = :black, linestyle = :solid)
    xlims!(ax, (50, 1050))
    ylims!(ax, (-19.5, -3.5))

    axislegend(ax, position = :rb)

    save("../figs/fig_error_Lz.pdf", fig)

    fig
end

begin
	data_error_fixn = load("data/error_fixn.jld")
	s_array = data_error_fixn["s_array"]
	error_fix_n = data_error_fixn["error_fix_n"]
	error_dir = data_error_fixn["error_dir"]

    data_error_N = load("data/error_N.jld")
	error_4 = data_error_N["error_4"]
	error_8 = data_error_N["error_8"]

    fig = Figure(size = (1000, 400), fontsize = 20)
    ax1 = Axis(fig[1, 1], yticks = ([0:-3:-15...], [L"$10^{0}$", L"$10^{-3}$", L"$10^{-6}$", L"$10^{-9}$", L"$10^{-12}$", L"$10^{-15}$"]), xlabel = L"s", ylabel = L"$\mathcal{E}$")
    ylims!(ax1, (-16.5, 1.5))

    M_label = [L"M = 4", L"M = 8", L"M = 16"]
	
	for i in 1:3
        l = M_label[i]
        ms = markerstyle[i]
        c = colors[i]
        sw = strokewidth
        msw = markerstrokewidth
        m = markersize
		scatter!(ax1, s_array, log10.(abs.(error_fix_n[i])), label = l, marker = ms, markersize = m, color = c, strokecolor = strokecolor, strokewidth = sw)
	end
    lines!(ax1, s_array, log10.(abs.(error_dir)), label = L"\text{Ewald2D}", linewidth = linewidth, color = :black, linestyle = :solid)
    axislegend(ax1, position = :lb)

    ax2 = Axis(fig[1, 2], xlabel = L"N", ylabel = L"$\mathcal{E}_r$", xticks = [100:200:900...], yticks = ([-4:-2:-12...], [L"$10^{-4}$", L"$10^{-6}$", L"$10^{-8}$", L"$10^{-10}$", L"$10^{-12}$"]))
    xlims!(ax2, (50, 600))
    ylims!(ax2, (-12.5, -3.5))

    n_atoms_array = [100:50:550...]
	
    error_s = [log10.(abs.(error_4)), log10.(abs.(error_8))]
    labels = [ L"s = 5,~M = 4", L"s = 5,~M = 8"]

    for i in 1:2
        l = labels[i]
        ms = markerstyle[i]
        c = colors[i]
        sw = strokewidth
        msw = markerstrokewidth
        m = markersize

        @show n_atoms_array
        @show error_s[i]

        scatter!(ax2, n_atoms_array, error_s[i], label = l, marker = ms, markersize = m, color = c, strokecolor = strokecolor, strokewidth = sw)
	end

    axislegend(ax2, position = :lb)

    text!(ax1, 1.5, 0, text = "(a)", fontsize = 30, align = (:right, :center))
    text!(ax2, 130, -4.2, text = "(b)", fontsize = 30, align = (:right, :center))

    save("../figs/fig_error_fixn.pdf", fig)

    fig
end