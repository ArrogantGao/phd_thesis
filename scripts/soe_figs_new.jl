using Plots, ExTinyMD, StatsPlots, SoEwald2D, LaTeXStrings, SpecialFunctions, JLD, Statistics, ForwardDiff

f=font(14,"Computer Modern")
fl=font(10,"Computer Modern")

begin
	k_array = [-10.00:0.1:10.0...]
	soexp_4 = [soexp_mul_erfc(k, 1.0, 1.0, SoePara4()) for k in k_array]
	soexp_8 = [soexp_mul_erfc(k, 1.0, 1.0, SoePara8()) for k in k_array]
	soexp_16 = [soexp_mul_erfc(k, 1.0, 1.0, SoePara16()) for k in k_array]
	exact_exp = [exp(k) * erfc(0.5 + k) for k in k_array]
	fig_exp = plot(size = [500,400], xlabel = L"z", ylabel = L"$\mathcal{E}$", xlim = [-10, 10], ylim = [-18.5, 0.5], framestyle = :box, yticks = ([0:-3:-18...], [L"$10^{0}$", L"$10^{-3}$", L"$10^{-6}$", L"$10^{-9}$", L"$10^{-12}$", L"$10^{-15}$", L"$10^{-18}$"]), titlefont=f, tickfont=f, legendfont=fl, guidefont = f)
	for (soeresult, name, m) in zip([soexp_4, soexp_8, soexp_16], ["M = 4", "M = 8", "M = 16"], [:solid, :dash, :dashdot])
		plot!(k_array, log10.(abs.(soeresult .- exact_exp)), label = :none, linewidth = 2, linestyle = m)
	end
	annotate!([(-9., -0.6, text(L"(a) $\xi^{\pm}(h, z)$", 20, :left, :top, :black, "Computer Modern"))])
	# title = L"(a)  $\xi^{\pm}(k, z)$"
	fig_exp
end

begin
	soe_4 = [soerfc(k, SoePara4()) for k in k_array]
	soe_8 = [soerfc(k, SoePara8()) for k in k_array]
	soe_16 = [soerfc(k, SoePara16()) for k in k_array]
	exact_erfc = [erfc(k) for k in k_array]
	fig_erf = plot(size = [500,400], xlabel = L"z", ylabel = L"$\mathcal{E}$", xlim = [0, 10], ylim = [-18.5, 0.5], framestyle = :box, yticks = ([0:-3:-18...], [L"$10^{0}$", L"$10^{-3}$", L"$10^{-6}$", L"$10^{-9}$", L"$10^{-12}$", L"$10^{-15}$", L"$10^{-18}$"]), titlefont=f, tickfont=f, legendfont=fl, guidefont = f)
	for (soeresult, name, m) in zip([soe_4, soe_8, soe_16], [L"M = 4", L"M = 8", L"M = 16"], [:solid, :dash, :dashdot])
		plot!(k_array, log10.(abs.(soeresult .- exact_erfc)), label = name, linewidth = 2, linestyle = m)
	end
	annotate!([(0.5, -1, text(L"(b) erf($\alpha z$)", 20, :left, :top, :black, "Computer Modern"))])
	fig_erf
end

begin
	k_array_positive = [0.0:0.1:10.0...]
	soexp_erfc_4 = [soexp_mul_erfc(1.0, 1.0, k, SoePara4()) for k in k_array_positive]
	soexp_erfc_8 = [soexp_mul_erfc(1.0, 1.0, k, SoePara8()) for k in k_array_positive]
	soexp_erfc_16 = [soexp_mul_erfc(1.0, 1.0, k, SoePara16()) for k in k_array_positive]
	exact_exp_erfc = [exp(k) * erfc(k / 2 + 1.0) for k in k_array_positive]
	fig_exp_erfc  = plot(
		size = [500,400], xlim = [0, 100], xlabel = L"h^2", ylabel = L"$\mathcal{E}$", framestyle = :box,
		yticks = ([-0:-5:-25...], [L"$10^{0}$", L"$10^{-5}$", L"$10^{-10}$", L"$10^{-15}$", L"$10^{-20}$", L"$10^{-25}$"]), ylim = [-26, 1],
		titlefont=f, tickfont=f, legendfont=fl, guidefont = f)
	for (soeresult, name, m) in zip([soexp_erfc_4, soexp_erfc_8, soexp_erfc_16], ["M = 4", "M = 8", "M = 16"], [:solid, :dash, :dashdot])
		plot!((k_array_positive).^2, log10.(abs.(soeresult .- exact_exp_erfc)), margin=4Plots.mm, linewidth = 2, linestyle = m, label = :none)
	end
	annotate!([(4, 0, text(L"(c) $\xi^{+}(h, z)$", 20, :left, :top, :black, "Computer Modern"))])
	fig_exp_erfc
end

# ╔═╡ 1e4fe9b7-90ee-4a8d-b43c-6afd33acbd68
begin
	xim_4 = [soexp_mul_erfc(- 1.0, 1.0, k, SoePara4()) for k in k_array_positive]
	xim_8 = [soexp_mul_erfc(- 1.0, 1.0, k, SoePara8()) for k in k_array_positive]
	xim_16 = [soexp_mul_erfc(- 1.0, 1.0, k, SoePara16()) for k in k_array_positive]
	exact_xim = [exp(-k) * erfc(k / 2 - 1.0) for k in k_array_positive]
	fig_xim  = plot(size = [500,400], xlim = [0, 100], xlabel = L"h^2", ylabel = L"$\mathcal{E}$", framestyle = :box, yticks = ([-0:-5:-25...], [L"$10^{0}$", L"$10^{-5}$", L"$10^{-10}$", L"$10^{-15}$", L"$10^{-20}$", L"$10^{-25}$"]), ylim = [-26, 1], titlefont=f, tickfont=f, legendfont=fl, guidefont = f, legend = :none)
	for (soeresult, name, m) in zip([xim_4, xim_8, xim_16], ["M = 4", "M = 8", "M = 16"], [:solid, :dash, :dashdot])
		plot!((k_array_positive).^2, log10.(abs.(soeresult .- exact_xim)), label = name, margin=4Plots.mm, linewidth = 2, linestyle = m)
	end
	annotate!([(4, 0, text(L"(d) $\xi^{-}(h, z)$", 20, :left, :top, :black, "Computer Modern"))])
	fig_xim
end

# ╔═╡ 44ff8a81-a802-4315-b726-ef8194fcc13b
begin
	fig_error_Uxi = plot(fig_exp, fig_erf, fig_exp_erfc, fig_xim, layout = (2, 2), size = [1000,800], margin=3Plots.mm)
	savefig(fig_error_Uxi, "../figs/fig_error_Uxi.pdf")
end

# ╔═╡ 08263031-b7ad-49a1-a5fa-8b8c9647650c
function dz_ξ(z::T, k::T) where{T}
	f = z -> exp(k * z) * erfc(k / 2.0 + z)
	return ForwardDiff.derivative(f, z)
end

# ╔═╡ 18633c5a-7b4e-4bb9-a83d-685108299eb9
function dz_erf(z::T) where{T}
	f = z -> erf(z)
	return ForwardDiff.derivative(f, z)
end

# ╔═╡ 2e520d9f-b754-4441-89c2-aba48ba9d990
function dz_soexp_mul_erfc(z::T, α::T, k::T, soepara::SoePara{T2}) where{T<:Real, T2}
    k = T2(k)
    sum = zero(ComplexF64)
    if z ≥ 0
        for (s, w) in soepara.sw
            sum += -s * α * α * w * 2 / sqrt(π) / (s * α + k) * exp(-s * α * z)
        end
    else
        for (s, w) in soepara.sw
            sum += α * w * 2 / sqrt(π) * ( 
                - s * α * exp(s * α * z) / (s * α - k) 
                + 2.0 * s * α / ((s * α)^2 - k^2) * k * exp(k * z))
        end
    end
    return T(real(exp(- k^2 / (4 * α^2)) * sum))
end

# ╔═╡ c136f707-5489-4e52-a2fb-a05514ddf689
function dz_soerfc(x::T1, sw::SoePara{T2}) where{T1<:Real, T2}
    sum = zero(T2)
    if x > 0.0
        for (si, wi) in sw.sw
            sum += - 2 / sqrt(π) * wi * exp( - si * x)
        end
    else
        for (si, wi) in sw.sw
            sum += - 2 / sqrt(π) * wi * exp(si * x)
        end
    end
    return T1(real(sum))
end

# ╔═╡ 5fe21efe-f780-46c2-9b4c-c8a869709068
function dz_soerf(x::T1, sw::SoePara{T2}) where{T1<:Real, T2}
    return - dz_soerfc(x, sw)
end

# ╔═╡ ab4bb87a-d4ff-417f-9186-d4d180e26d06
begin
	z_array = [-10.0:0.01:10.0...]
	dξz = [dz_ξ(z, 1.0) for z in z_array]
	dξM4z = [dz_soexp_mul_erfc(z, 1.0, 1.0, SoePara4()) for z in z_array]
	dξM8z = [dz_soexp_mul_erfc(z, 1.0, 1.0, SoePara8()) for z in z_array]
	dξM16z = [dz_soexp_mul_erfc(z, 1.0, 1.0, SoePara16()) for z in z_array]
end

# ╔═╡ b46982b0-c125-438c-aa2f-227a318ab7dd
begin
	fig_error_dzξ = plot(size = [500,400], xlabel = L"z", ylabel = L"$\mathcal{E}$", xlim = [-10, 10], ylim = [-18.5, 0.5], framestyle = :box, yticks = ([0:-3:-18...], [L"$10^{0}$", L"$10^{-3}$", L"$10^{-6}$", L"$10^{-9}$", L"$10^{-12}$", L"$10^{-15}$", L"$10^{-18}$"]), titlefont=f, tickfont=f, legendfont=fl, guidefont = f)

	for (dξzM, name, m) in zip([dξM4z, dξM8z, dξM16z], ["M = 4", "M = 8", "M = 16"], [:solid, :dash, :dashdot])
		plot!(z_array, log10.(abs.(dξz - dξzM)), label = :none, linewidth = 2, linestyle = m)
	end

	annotate!([(-9., -0.6, text(L"(a) $\partial_z \xi^{\pm}(h, z)$", 20, :left, :top, :black, "Computer Modern"))])

	fig_error_dzξ
end

# ╔═╡ d0434048-37b0-4a0b-ab47-f951e2af070d
begin
	z_array_positive = [0.0:0.01:10.0...]
	dz_erf_r = [dz_erf(z) for z in z_array_positive]
	dz_soerf_r4 = [dz_soerf(z, SoePara4()) for z in z_array_positive]
	dz_soerf_r8 = [dz_soerf(z, SoePara8()) for z in z_array_positive]
	dz_soerf_r16 = [dz_soerf(z, SoePara16()) for z in z_array_positive]
end

# ╔═╡ 9291ef76-a92e-4d57-85fc-23f604509e6d
begin
	fig_error_dz_erf = plot(size = [500,400], xlabel = L"z", ylabel = L"$\mathcal{E}$", xlim = [0, 10], ylim = [-18.5, 0.5], framestyle = :box, yticks = ([0:-3:-18...], [L"$10^{0}$", L"$10^{-3}$", L"$10^{-6}$", L"$10^{-9}$", L"$10^{-12}$", L"$10^{-15}$", L"$10^{-18}$"]), titlefont=f, tickfont=f, legendfont=fl, guidefont = f)

	for (dz_soerf, name, m) in zip([dz_soerf_r4, dz_soerf_r8, dz_soerf_r16], [L"M = 4", L"M = 8", L"M = 16"], [:solid, :dash, :dashdot])
		plot!(z_array_positive, log10.(abs.(dz_erf_r - dz_soerf)), label = name, linewidth = 2, linestyle = m)
	end

	annotate!([(0.5, -0.8, text(L"(b) $\partial_z \mathrm{erf}(\alpha z)$", 20, :left, :top, :black, "Computer Modern"))])

	fig_error_dz_erf
end

# ╔═╡ e997f531-104e-4e11-9b9a-487ee0b37c6f
begin
	dξk = [dz_ξ(1.0, k) for k in k_array_positive]
	dξM4k = [dz_soexp_mul_erfc(1.0, 1.0, k, SoePara4()) for k in k_array_positive]
	dξM8k = [dz_soexp_mul_erfc(1.0, 1.0, k, SoePara8()) for k in k_array_positive]
	dξM16k = [dz_soexp_mul_erfc(1.0, 1.0, k, SoePara16()) for k in k_array_positive]
end

# ╔═╡ 9e9b9c6b-af24-4a3e-9d07-965e5d391fb8
begin
	fig_error_dz_kp = plot(size = [500,400], xlim = [0, 100], xlabel = L"h^2", ylabel = L"$\mathcal{E}$", framestyle = :box, yticks = ([-0:-5:-25...], [L"$10^{0}$", L"$10^{-5}$", L"$10^{-10}$", L"$10^{-15}$", L"$10^{-20}$", L"$10^{-25}$"]), ylim = [-26, 1], titlefont=f, tickfont=f, legendfont=fl, guidefont = f, legend = :none)

	for (dz_ξk, name, m) in zip([dξM4k, dξM8k, dξM16k], ["M = 4", "M = 8", "M = 16"], [:solid, :dash, :dashdot])
		plot!(k_array_positive.^2, log10.(abs.(dξk .- dz_ξk)), label = :none, linewidth = 2, linestyle = m)
	end

	annotate!([(4, 0, text(L"(c) $\partial_z \xi^{+}(h, z)$", 20, :left, :top, :black, "Computer Modern"))])

	fig_error_dz_kp
end

# ╔═╡ 5d721aa9-84c9-4839-b1cd-c6a029ee706a
begin
	dξkn = [dz_ξ(- 1.0, k) for k in k_array_positive]
	dξM4kn = [dz_soexp_mul_erfc(- 1.0, 1.0, k, SoePara4()) for k in k_array_positive]
	dξM8kn = [dz_soexp_mul_erfc(- 1.0, 1.0, k, SoePara8()) for k in k_array_positive]
	dξM16kn = [dz_soexp_mul_erfc(- 1.0, 1.0, k, SoePara16()) for k in k_array_positive]
end

# ╔═╡ 179758bb-ce16-4a96-a431-a480be317443
begin
	fig_error_dz_kn = plot(size = [500,400], xlim = [0, 100], xlabel = L"h^2", ylabel = L"$\mathcal{E}$", framestyle = :box, yticks = ([-0:-5:-25...], [L"$10^{0}$", L"$10^{-5}$", L"$10^{-10}$", L"$10^{-15}$", L"$10^{-20}$", L"$10^{-25}$"]), ylim = [-26, 1], titlefont=f, tickfont=f, legendfont=fl, guidefont = f, legend = :none)

	for (dz_ξk, name, m) in zip([dξM4kn, dξM8kn, dξM16kn], ["M = 4", "M = 8", "M = 16"], [:solid, :dash, :dashdot])
		plot!(k_array_positive.^2, log10.(abs.(dξkn .- dz_ξk)), label = :none, linewidth = 2, linestyle = m)
	end

	annotate!([(4, 0, text(L"(d) $\partial_z \xi^{-}(h, z)$", 20, :left, :top, :black, "Computer Modern"))])

	fig_error_dz_kn
end

# ╔═╡ 1c442239-7623-4056-8457-5002c0e2ab1f
begin
	fig_error_diff = plot(fig_error_dzξ, fig_error_dz_erf, fig_error_dz_kp, fig_error_dz_kn, layout = (2, 2), size = [1000,800], margin=3Plots.mm)
	savefig(fig_error_diff, "../figs/fig_error_diff.pdf")
end

begin
	f=font(14,"Computer Modern")
	fl=font(10,"Computer Modern")
end

begin

	data = load("./data/kLz_lim.jld")
	Lz_array = data["Lz_array"]
	error_soe_array = data["error_soe_array"]
	error_dir_array = data["error_dir_array"]

	fig_error_Lz = plot(size = [600, 400], ylim = [-18.5, -3.5], titlefont=f, tickfont=f, legendfont=fl, guidefont = f, framestyle = :box, yticks = ([-5:-3:-17...], [L"$10^{-5}$", L"$10^{-8}$", L"$10^{-11}$", L"$10^{-14}$", L"$10^{-17}$"]), xlabel = L"H", ylabel = L"$\mathcal{E}$", xlim = [50, 1050])
	ms = [:square, :utriangle, :circle]
	labels = [L"s = 3,~M = 4", L"s = 4,~M = 8", L"s = 5,~M = 16"]
	for i in 1:3
		scatter!(Lz_array, log10.(abs.(error_soe_array[i])), marker = ms[i], label = labels[i], linewidth = 1, markerstrokewidth = 0.5, markersize = 4)
	end
	plot!(Lz_array, log10.(abs.(error_dir_array)), label = "Ewald2D", linewidth = 2, markerstrokewidth = 0.5, color = :red, linestyle = :solid)
	fig_error_Lz
end

# ╔═╡ 3c892de9-1534-4244-80c0-5501648830c5
savefig(fig_error_Lz, "../figs/fig_error_Lz.pdf")