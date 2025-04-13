using CairoMakie, SpecialFunctions, LaTeXStrings, GaussQuadrature

function fi_singular(k, k0)
    return besselj0(k0) * exp(-k0) / (exp(-2(k - k0)) - 1)
end

function fi_smooth(k, k0)
    return (besselj0(k) * exp(-k) - besselj0(k0) * exp(-k0))  / (exp(-2(k - k0)) - 1) - besselj0(k0) * exp(-k0)
end

function fd(k, k0)
    return besselj0(k0) * exp(-k0) / (exp(-2(k - k0)) + 1)
end

function plot_fi()
    k = range(0.0, 4.0, length=1000)
    ks1 = range(0.0, 0.990, length=1000)
    ks2 = range(1.001, 4.0, length=1000)
    k0 = 1.0
    f2 = x -> fi_singular(x, k0)
    fi_m = [fi_smooth(k[i], k0) for i in 1:length(k)]
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1]; xlabel = L"h", ylabel = L"f(h)")
    xlims!(ax, 0, 4)
    ylims!(ax, -5, 5)
    lines!(ax, ks1, f2.(ks1); label = "singular", color = :red)
    lines!(ax, ks2, f2.(ks2); color = :red)
    lines!(ax, k, fi_m; label = "smooth", color = :blue)
    axislegend(ax)
    save("../figs/fi.pdf", fig)
    fig
end

fig = plot_fi()

function truncated_int(f::Function, a::T, b::T, n::Int) where{T}
    x, w = legendre(n)
    return sum(w .* f.((b+a)/2 .+ (b-a)/2 .* x)) * (b-a)/2
end

function fo(k, k0)
    if k ≤ 2 * k0
        return besselj0(k) * exp(-k) / (exp(-2(k - k0)) - 1) + besselj0(k0) * exp(-k0) / (2 * (k - k0))
    else
        return besselj0(k) * exp(-k) / (exp(-2(k - k0)) - 1)
    end
end

function plot_int()
    k_fs = [10.0, 20.0, 30.0]
    colors = [:red, :blue, :green]
    exact_int = truncated_int(x -> fo(x, 1.0), 0.0, 2.0, 2000) + truncated_int(x -> fo(x, 1.0), 2.0, 40.0, 1000)
    results = []
    for k_f in k_fs
        rt = []
        for n in 1:100
            t = truncated_int(x -> fi_smooth(x, 1.0), 0.0, k_f, n) + besselj0(1.0) * exp(-1.0) / 2 * log(exp(2) - 1)
            push!(rt, t)
        end
        push!(results, rt)
    end

    rel_error = [abs.(ri .- exact_int) ./ abs(exact_int) for ri in results]

    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1]; xlabel = L"n", ylabel = L"\mathcal{E}_r", yscale = log10)
    xlims!(ax, 0, 100)
    ylims!(ax, 1e-16, 1e1)
    for (i, k_f) in enumerate(k_fs)
        # scatter!(ax, [1:100...], rel_error[i], color = colors[i], label = "a")
        lines!(ax, [1:100...], rel_error[i], color = colors[i], label = L"k_f = %$(k_f)")
    end
    axislegend(ax)

    save("../figs/int_convergence.pdf", fig)

    fig
end

fig = plot_int()