# Thesis review findings

Full line-by-line pass over `main.tex`, `00II/00III`, `chapters/*.tex`, `appendix/*.tex`,
`references.bib`. Equations were re-derived independently where possible; "verified" below means
I reproduced the algebra. Line numbers refer to the files **as of the original review**, before
any fixes were applied — they will have shifted slightly in `chp_introduction.tex` and `main.tex`.

## Status

Branch `review-fixes`, five commits. Everything in sections **A, B (except the items listed as
still open below), C, D-partial, E-partial and F has been applied**; the build now reports
0 errors, 0 undefined references, 0 undefined citations, 0 amsmath warnings, 216 pages
(214 before — the restored abstract adds two).

**Still open — these change results or narrative, so they are yours to decide:**

| Item | Why it was left |
|---|---|
| **B26** odd-`M` gap in `eq:beta_ij` | The bound needs restating or restricting to even `M`; either choice alters a lemma. |
| **B28** the `\lfloor (M+1)/2\rfloor^{-1}` in `eq:U_gamma_bound` | Restoring it makes the energy and force rates equal, which contradicts the sentence at `chp_icmewald2d.tex:204` and the surrounding discussion. |
| **B27** missing `L_xL_y/(2\pi)^2` in Ch. 3's lattice-sum→integral steps | Fixing it changes every explicit prefactor in `eq:U_gamma_bound` and the force bound. |
| **B34** sign of the `4\pi H/\max\{L_x,L_y\}` term in `eq::38` | The printed form is one image level more conservative; may well be intended. |
| **B36** `|g_{\rm u}g_{\rm d}|<1` claimed for `\gamma=0.6, H=1` (actually 1.265) | Correcting it flips which case of the parameter-selection strategy applies, and Table 3.1's `\gamma=0.6` column needs recomputing. |
| **B37** trapezoidal remainder `e^{-\alpha^2(L_z-H)^2}` with images present | Should probably be `(L_z-(M+1)H)`; affects Eq. `eq:error_icmelc` and Step 3. |
| **B48** the `\alpha^3` in Eq. (4.x) | I could not reproduce it; needs re-derivation, not a guess. |
| **B49** dimensionally inconsistent `3+\alpha/\sqrt\pi` and `(1+2\alpha)L_z` | Needs the correct grouping, which follows from B89. |
| **B50** Eq. `eq::estiZ` | The stated chain does not hold; the repair needs the intended argument. |
| **B51 / B73** factor-2 and factor-2 questions in the variance proofs | Loose, not wrong; fixing changes theorem constants. |
| **B87** self-energy treatment in Ch. 6 | Needs an explicit statement of what is kept, which only you can supply. |
| **B89** the `\lambda_D^2` in `eq::E.4` | Load-bearing: every DH bound in Chs. 4 and 6 carries `\lambda_D^2`/`\lambda_D^4`. |
| **B101** the `z_i` in the ion–wall force | Looks like the energy expression; please confirm before I touch it. |
| **B104** non-neutral rows (C), (D) of `Table::Dielec` | The tabulated numbers need checking against the actual runs. |
| **B66 / D** `\alpha` convention differs in Ch. 6; the rest of the symbol collisions | Unifying notation across chapters is an editorial decision with wide blast radius. |
| **E** `chp_icmewald2d.tex:123` (`\bm h=\bm 0` mode claim), `chp_applications.tex:111` (`1.2\pi`/`1.4\pi` vs `2\pi`), `chp_applications.tex:58` vs caption (10x vs 1.5 orders), 4000-core vs 1024-core figure range, `app_softwares.tex` omits QuasiEwald.jl | All statements about results. |

**Two remaining build warnings, both pre-existing and cosmetic:** 48 `destination with the same
identifier (name{figure.N.M})` and 49 hyperref `Token not allowed in a PDF string`. The first comes
from `ustthesis.cls`, which emits each float anchor twice; it does not affect links or output, and I
did not patch the class. (My original note attributing this to the two stray minipage labels was
wrong — those were unused labels, worth removing on their own, but not the cause.) The second comes
from math in captions/headings flowing into PDF bookmarks; `\texorpdfstring` would silence it.

---

## A. Must-fix — structural / renders wrong in the PDF

**A1. The abstract is missing from the thesis.**
`00III_abstract.tex` is never `\input`. `main.tex` §7 (lines 272–281) contains only comments —
there is no `\abstract ... \endabstract` block. Add:
```tex
\abstract
\input{00III_abstract.tex}
\endabstract
```

**A2. Four undefined references print as `??`.** All caused by an automated label-renaming pass that
appended the suffix twice in `\label` but once in `\ref`:

| `\ref` used | `\label` actually defined | Location of `\ref` |
|---|---|---|
| `eq::32_rbe` | `eq::32_rbe_rbe` | `chp_rbe2d.tex:81` |
| `eq::RBapp_qem` | `eq::RBapp_qem_qem` | `chp_quasiewald.tex:703` |
| `prop:unbiased_qem` | `prop:unbiased_qem_qem` | `chp_quasiewald.tex:736` |
| `lem::upper_bound_phiRB_qem` | **commented out** at `chp_quasiewald.tex:752` | `chp_quasiewald.tex:813` |

Also check the other doubled-suffix labels for the same problem: `eq::33_icm_icm`
(`chp_icmewald2d.tex:117`), `eq::hk_qem_qem` (669), `eq::xichi_qem_qem` (732).

**A3. `chp_quasiewald.tex:813` — the proof of Theorem `thm:unbiased_qem` cites a lemma that was
deleted from the thesis.** `lem::upper_bound_phiRB_qem` is commented out (lines 751–769), so the
key bound in the proof has no support. Either restore the lemma or supply the bound inline.

**A4. Revision-markup braces that split words** (`\rev{}` was replaced by ` {}` leaving a space):
- `chp_rbe2d.tex:108` — `the force exert {ed} on` → renders "the force exert ed on".
- `chp_rbe2d.tex:275` — `different batch size {s} $P$` → renders "batch size s P".

**A5. `chp_applications.tex:18` — "the ratio $H/H$ is fixed at $3$".** `$H/H$` = 1. Should be
`$L_z/H$` (cf. `chp_rbe2d.tex:289`, "as $L_z \geq 3H$").

**A6. `chp_icmewald2d.tex:482, 508` — a table referenced as a figure.**
`Fig.~\ref{tab:parameter_selection_results}` → `Table~\ref{...}`. Same class of error at
`chp_rbe2d.tex:128`: `Fig.~\ref{Alg::1}` → `Algorithm~\ref{Alg::1}`.

**A7. Duplicate bibliography entries — both keys cited, so both are printed with different numbers:**
- Julia: `Julia-2017` (`chp_quasiewald.tex:972`) and `Bezanson_Julia_A_fresh_2017` (`app_softwares.tex:1`).
- CellListMap.jl: `celllistmap` (`chp_quasiewald.tex:972`) and `MARTINEZ2022108452` (`app_softwares.tex:3`).

(9 further duplicate-title pairs exist in `references.bib` but only one key of each is cited, so they
are harmless dead weight: `Bagchi2020PNAS`/`bagchi2020surface`,
`barnes1986hierarchical`/`Barnes1986Nature`, `de1980simulation`/`de1980simulation1`,
`doi:10.1021/acs.jctc.3c01124`/`huang2024pmc`, `gan2024fast`/`gan2025fast`,
`grzybowski2000ewald`/`PhysRevB.61.6706`, `Hastings1970Bio`/`hastings1970monte`,
`jiang2008efficient`/`Jiang2008EfficientRO`, `LAMMPS`/`thompson2021lammps`,
`liang2023random`/`liang2023SISC`.)

**A8. `00II_acknowledgments.tex:61` — publications B and C are mapped to the wrong chapters.**
Publication B is the SISC *Random batch Ewald method for dielectrically confined Coulomb systems*
(= RBE2D = Chapter `chp_rbe2d`); publication C is *Fast algorithm for quasi-2D Coulomb systems*,
JCP **524**:113733 (= SOEwald2D = Chapter `chp_soewald2d`). The sentence currently says
Chapters `chp_soewald2d`, `chp_rbe2d`, `chp_quasiewald` come "from B, C, and E, respectively" —
should be **C, B, and E**.

---

## B. Mathematical errors in displayed equations

Each item was checked by independent re-derivation.

### Chapter 2 — Preliminaries

**B1. `:291` — wrong Ewald3D self-energy.** `\sqrt{\frac{\alpha}{\pi}}\sum q_i^2` should be
`\frac{\alpha}{\sqrt{\pi}}\sum q_i^2`. The correct form is used later at `:333` and `:465`, and
follows from `-\frac12\sum_i q_i\phi_{\rm self}^i` with `\phi_{\rm self}^i = 2\alpha q_i/\sqrt\pi`
(Eq. `eq::self`).

**B2. `:329` — RBE normalization constant double-corrects for `k=0`.**
`S = \sum_{\V k\neq\V 0} e^{-k^2/4\alpha^2} - 1` subtracts 1 from a sum that already excludes
`\V k = \V 0`. Drop the `-1` (as done correctly in `chp_rbe2d.tex:38`), or sum over all `\V k`.

**B3. `:333` — the RBE estimator is missing the factor `S/P`, so it is not unbiased.**
The text immediately after ("thus `U_{\ell,*}` is an unbiased estimator") is therefore false as
written. Should be `\frac{2\pi S}{V P}\sum_{\ell=1}^{P} \frac{|\rho(\V k_\ell)|^2}{k_\ell^2}` —
cf. `chp_soewald2d.tex:430` and `chp_rbe2d.tex:25`, which both have `S/P`.

**B4. `:252` — `\gamma_+^{(l)}` and `\gamma_-^{(l)}` are interchanged.**
With `z_+^{(1)} = 2H - z` (reflection in the *upper* interface, Eq. `eq:z_l`) the weight must be
`\gamma_{\rm u}`, but the stated `\gamma_+^{(1)} = \gamma_{\rm d}^{\lceil 1/2\rceil}\gamma_{\rm u}^{\lfloor 1/2\rfloor} = \gamma_{\rm d}`.
Checking `l=3`: the path u→d→u gives `\gamma_{\rm u}^2\gamma_{\rm d}`, the formula gives
`\gamma_{\rm d}^2\gamma_{\rm u}`. Correct definitions:
`\gamma_+^{(l)} = \gamma_{\rm u}^{\lceil l/2\rceil}\gamma_{\rm d}^{\lfloor l/2\rfloor}`,
`\gamma_-^{(l)} = \gamma_{\rm d}^{\lceil l/2\rceil}\gamma_{\rm u}^{\lfloor l/2\rfloor}`.
This propagates into `chp_icmewald2d.tex:132, 262`, `chp_rbe2d.tex:204–205`, `app_force.tex:96–103`,
so fix consistently. (Everything downstream is *internally* consistent with the current
definition, so the fix is a single swap.)

**B5. `:450` — `A_0 = \frac{2\pi}{L_xL_y}\sum_j q_j z \equiv 0`.** The stray `z` makes a constant
depend on `z`; should be `\sum_j q_j`.

**B6. `:444` — `\phi_s^{\bm 0}` line has three slips:** the `q_j` factor is missing (it is present on
the next line), `|z|` should be `|z-z_j|`, and `\xi^{\pm}(\bm h, z)` should be
`\xi^{\pm}(h, z-z_j)` (`\xi` takes the scalar `h`, per Eq. `eq::xi_def`).

**B7. `:516` — asymptotic expansion of `erfc` has the wrong variable and a wrong Pochhammer
definition.** `z^{-(2m+1)}` should be `r^{-(2m+1)}`; and `(x)_m = x(x-1)\cdots(x-m+1)` is the
*falling* factorial, whereas the expansion requires the *rising* one,
`(x)_m = x(x+1)\cdots(x+m-1) = \Gamma(x+m)/\Gamma(x)`. (For `m=2` the falling version gives `-1/4`
instead of the correct `3/4`.) Also `=` should be `\sim`.

**B8. `:524` — `|\bm h + \kappa| > k_c`** adds a 2-vector to a scalar. Should be
`\sqrt{h^2+\kappa^2} > k_c`.

**B9. `:531` — "the coordinate along `\cos\theta`"** contradicts the formula, which uses `\cos\varphi`.

**B10. `:465` — `e^{\i \V h\cdot\V r_{ij}}` should be `e^{\i \V h\cdot\V\rho_{ij}}`** (`\V h` is 2D).
Same in `chp_icmewald2d.tex:13, 57, 82, 120, 126, 185, 244, 286`.

**B11. `:602` — `\lim_{z\to\pm\infty}\phi_{\rm p-s} = \mp\frac{2\pi}{L_xL_y}(\ldots|z|\ldots)`** is
malformed: the RHS still contains `|z|` and `|z-L_z|`, which diverge. Either keep `\lim` on the RHS
with a single `-` sign (as in Eq. `eq::boundary2`), or evaluate it. Note Eq. `eq::phionwallzero`
uses `-`, not `\mp`.

**B12. `:673` — `U_* = \sum_i \phi_*` is missing `\tfrac12` and `q_i`.** Should be
`U_* = \tfrac12\sum_i q_i\phi_*`, which is what Eqs. `eq::Us_Ewald2D`/`eq::Ul_Ewald2D` actually use.

**B13. `:566` — `\mathscr E_{\phi_s}\approx\mathscr E_{\phi_\ell} \sim Q\sqrt{s/(\alpha V)}\,e^{-s^2}/s^2`
has `Q` where it should be `\sqrt{Q}`.** With `r_c=s/\alpha`, Theorem
`thm:ewald2d_phi_error` gives `\mathscr E_{\phi_s} = \sqrt{Q/V}\,\alpha^{-1/2}s^{-3/2}e^{-s^2}`.
The `Q` prefactor as written matches the *energy* error `\mathscr E_{U_s}` of Prop. `prop::2.12`,
not the potential error named in the equation.

**B14. `:133 / :215` — "absolutely convergent" is the wrong term.** The proof of Theorem `order`
relies on the odd dipole term cancelling by symmetry of `\mathcal S`, i.e. it is *conditional*
convergence: `\sum_{\bm m}|\bm{\mathcal M}|^{-2}` diverges logarithmically over a 2D lattice, so
absolute convergence fails. Condition (1) (symmetry of `\mathcal S`) would be unnecessary if the
series were absolutely convergent. Suggest "convergent" / "convergent in the shape-`\mathcal S`
sense". Same wording at `:386` for `\phi_s` is fine (that one really is absolutely convergent).

**B15. `:14` — `\Omega\in\Omega_{\rm c}`** should be `\Omega\subset\Omega_{\rm c}`.

**B16. `:229, :233` — `r \in \mathbb{R}^3` and `r\to\infty`** use italic `r` for the position
vector; should be `\V r` and `|\V r|\to\infty`.

**B17. `:232` — the `+`/`-` limits are swapped relative to the stated convention.** Line 238 defines
`+`/`-` as exterior/interior; line 231 correctly pairs `\eps_{\rm c}` with `|_-`, but line 232 pairs
`\eps_{\rm c}` with `|_+`.

**B18. `:249 vs :252` — `\V r_{\pm}^{\prime(l)}` (used in Eq. `eq:ICM`) vs `\V r_{\pm}^{(l)}`
(defined at 252).** Also the image positions are written with the *unprimed* `(x,y,z)`; they are
images of the *source*, so they should be `(x', y', z'^{(l)}_\pm)`.

**B19. `:262` — `|\gamma|\le 1` should be `<1`** for the image series to converge; `|\gamma|=1`
(the `\eps\to0` or `\eps\to\infty` limit) is exactly the marginal case.

**B20. `:26` (Chapter 1) — `\V m \in \mathbb Z^2\otimes\{0\}`** uses a tensor product for what is a
Cartesian product: `\mathbb Z^2\times\{0\}`.

**B21. Unit convention clash.** Eq. `eq::U_pair` (`:211`) uses Gaussian units (`q_iq_j/r`), while
Eq. `eq:U_direct` (`:223`) with Eq. `eq:Poisson_G` uses `G = 1/(4\pi\eps|r|)`, i.e. the two energy
definitions differ by `4\pi`. `chp_icmewald2d.tex:119` then drops the `1/4\pi` again. Pick one
convention and state it once.

### Chapter 3 — ICM-Ewald2D and error estimates

**B22. `:31` — Eq. `eq::Ga-integral` is missing the prefactor `h/(\pi\alpha)`.**
Verified two ways. (i) By substituting `\kappa = 2\alpha t` into Eq. `eq::integral`:
`\mathcal G_\alpha(h,z) = \frac{h}{\pi\alpha}\int \frac{e^{-h^2/(4\alpha^2)-t^2}}{h^2/(4\alpha^2)+t^2}e^{2\i\alpha zt}\,dt`.
(ii) At `z=0`: the stated RHS is `\frac{2\pi\alpha}{h}\,{\rm erfc}(h/2\alpha)` whereas
`\mathcal G_\alpha(h,0) = 2\,{\rm erfc}(h/2\alpha)` — ratio `\pi\alpha/h`.
Note the companion Eq. `eq::J0-integral` (`:36`) **is** correct (checked at `z=0`), and Eq.
`eq::UcFour`'s prefactor `2\pi/(L_xL_yL_z)` only comes out right *with* the missing factor
restored — so `eq::UcFour` is fine and only `eq::Ga-integral` needs the fix.

**B23. `:48` — index error in the structure factor.**
`\widetilde\rho^M_{\bm k} := \sum_{j=1}^N q_j[e^{-\i\bm k\cdot\bm r_i} + \ldots]` uses `\bm r_i`
inside a sum over `j`. Should be `\bm r_j`.

**B24. `:42 vs :45/:48` — `\bar\rho_{\bm k}^M` vs `\widetilde\rho_{\bm k}^M`.** Eq. `eq::UcFour`
uses the bar, the definition uses the tilde. (`chp_rbe2d.tex:171` then says
"`\rho^M_{\bm k}` the conjugate of `\bar\rho^M_{\bm k}`", which is circular.) Unify.

**B25. `:132` — `i` and `j` are interchanged in the second line of Eq. `eq::rGrrij`.**
From line 131, `|z_i - z_{j+}^{(l)}| = 2\lceil l/2\rceil H + (-1)^l z_j - z_i`, but line 132 writes
`+(-1)^l z_i - z_j`. Same swap in the second term. (Harmless for the final bound because the double
sum is symmetric, but wrong as a definition of `\mathcal G_M(z_i,z_j)`.)

**B26. `:155` — Eq. `eq:beta_ij` is not valid for odd `M`, despite "In both cases".**
The "fact" proved at `:158–163` is exactly the *even*-`M` bracket. For odd `M` the bracket contains
`e^{-hz_{ij}} + e^{+hz_{ij}}` with coefficient 1, and `e^{+hz_{ij}}` is bounded by `e^{hH}`, not by a
constant. Combined with the prefactor `e^{-(M+1)hH}` the odd case gives
`\lesssim 2(\gamma_{\rm u}\gamma_{\rm d})^{(M+1)/2}e^{-MhH}`, i.e. one power of `e^{-hH}` weaker
than claimed. Either restrict the lemma to even `M` or state the odd case separately.

**B27. `:170, :193, :271, :294` — the lattice-sum→integral conversions omit the mode density
`L_xL_y/(2\pi)^2`.** `\sum_{\bm h\neq\bm 0}F(h) \approx \frac{L_xL_y}{4\pi^2}\int 2\pi h F(h)\,dh`
— this is stated correctly in `chp_soewald2d.tex:313` (Eq. `eq::integral2`), but Chapter 3 uses
`\sum \to 2\pi\int h\,dh`. All the explicit prefactors in Eqs. `eq:U_gamma_bound` and the force bound
are therefore off by `L_xL_y/4\pi^2`. (The `M`-dependence, which is what the chapter actually uses,
is unaffected.)

**B28. `:171` — Eq. `eq:U_gamma_bound` drops the factor `1/(2\lfloor (M+1)/2\rfloor H)` from the
integral, and this invalidates the conclusion at `:204`.**
`\int_{h_0}^\infty 2\pi h\,\frac{e^{-\lambda h}}{h}dh = 2\pi e^{-\lambda h_0}/\lambda` with
`\lambda = 2\lfloor (M+1)/2\rfloor H`. Retaining it, the energy bound also carries
`\lfloor (M+1)/2\rfloor^{-1}`, so the energy and force errors have the **same** `M`-rate — and the
claim at `:204` that "the error in force converge[s] slightly faster than that of the energy by a
factor of `\lfloor (M+1)/2\rfloor^{-1}`" is an artifact of the dropped factor. (Corroborated by
`:307`, where the same technique applied to the ELC term correctly yields *identical* rates for
energy and force.) Eq. `eq:Uerr` (`:180`) should therefore also carry
`\lfloor (M+1)/2\rfloor^{-1}`.

**B29. `:194` — the `1 + \lambda h_0` factor has `2\pi` where it should be `4\pi`.**
`\lambda h_0 = 2\lfloor (M+1)/2\rfloor H \cdot 2\pi/\max\{L_x,L_y\} = 4\pi\lfloor\cdot\rfloor H/\max\{L_x,L_y\}`.
(The analogous factor `J_1` at `:295, :301` **is** correct.)

**B30. `:185, :187` — Eq. `eq::Ferr` has a spurious extra `l`-sum and a false "identity".**
`\sum_{l=M+1}^{\infty}\mathcal G_M(z_i,z_j)`: `\mathcal G_M` already contains the sum over
`l \ge M+1` and has no `l` dependence, so the outer sum diverges — delete it. And
"`\partial_{z_i}\sum_{ij}\mathcal G_M = 2h\sum_{ij}\mathcal G_M`" is not an identity: term by term
`|\partial_{z_i}| \le h\,(\cdot)`, with signs that may cancel. State it as an inequality and justify
the factor 2.

**B31. `:334 / :339` vs `:317/:326–327` — inconsistent exponential prefactor.**
Eq. `eq::Error_refo` and case (i) carry `e^{-2\pi(L_z-H)/\max\{L_x,L_y\}}`; cases (ii) and (iii)
carry `e^{-2\pi L_z/\max\{L_x,L_y\}}`. They differ by the `M`-independent factor
`e^{2\pi H/\max\{L_x,L_y\}}` and cannot both follow from Eq. `eq::Error_refo`.

**B32. `:457` — Eq. `eq::Lz_gamma_u_gamma_d_less1` is missing a factor 2.**
Inverting `\varepsilon \sim (g_{\rm u}+g_{\rm d}+2)e^{-2\pi L_z/\max\{L_x,L_y\}}` and expanding
`g = \gamma e^{2\pi H/\max L}` gives
`\log|\gamma_{\rm u}+\gamma_{\rm d}+ \mathbf{2}e^{-2\pi H/\max\{L_x,L_y\}}|`; the printed formula has
no `2`.

**B33. `:462` — Eq. `eq::Lz_gamma_u_gamma_d_greater1` should have `\frac{M}{2}\log|\gamma_{\rm u}\gamma_{\rm d}|`.**
From case (i) with `M` even,
`\varepsilon \sim 2|\gamma_{\rm u}\gamma_{\rm d}|^{M/2}e^{-2\pi(L_z-(M+1)H)/\max\{L_x,L_y\}}`, so the
`\log|\gamma_{\rm u}\gamma_{\rm d}|` term picks up `M/2`. (The two agree for `\gamma=1`, which is
why Table 3.1's `\gamma=1` rows still check out — I verified 32, 62, 91.)

**B34. `:447` — the sign of the `4\pi H/\max\{L_x,L_y\}` term in Eq. `eq::38` does not follow from a
direct inversion of Eq. `eq:Uerr`.** Setting `n=\lfloor (M+1)/2\rfloor` and `M = 2n-1` gives
numerator `2\log\varepsilon - \log|\gamma_{\rm u}\gamma_{\rm d}| + \frac{4\pi H}{\max\{L_x,L_y\}}`
(a `+`), whereas the text has `-`. The printed version is equivalent to `M = 2n+1`, i.e. one image
level more conservative — which is presumably intended, but the derivation should say so.

**B35. `:470` — "`where $s = r_c/\alpha$`" is inverted.** Everywhere else `r_c = s/\alpha`, i.e.
`s = \alpha r_c` (`:21`, Eq. `eq::rckc`).

**B36. `:480` — the claim `|g_{\rm u}g_{\rm d}| < 1` for `\gamma=0.6` is numerically false.**
With `L_x=L_y=10, H=1`: `|g_{\rm u}g_{\rm d}| = 0.36\,e^{4\pi/10} = 0.36\times3.514 = 1.265 > 1`.
(The earlier test at `:363–365` uses `H=0.5`, where `0.36\,e^{0.628}=0.675<1` — that one is fine.)
This matters because Step 2 then applies Case 1 (Eq. `eq::Lz_gamma_u_gamma_d_less1`), which requires
`|g_{\rm u}g_{\rm d}|<1`. Relatedly, Table 3.1's `\gamma=0.6` values of `L_z` (15, 30, 45) do **not**
follow either Case-1 or Case-2 formula; they follow the bare
`e^{-2\pi(L_z-H)/L_x} = \varepsilon` (giving 17, 30, 45 — close but not the printed 15). Worth
recomputing the table from the stated strategy.

**B37. `:68 / :38` — the trapezoidal remainder estimate ignores the image charges.**
`|U^M_{\rm trap}| \sim e^{-\alpha^2(L_z-H)^2}` and "`$L_z>H$ is a parameter`" hold for the
homogeneous case. `app_math.tex:158` (Eq. `eq::A.8`) gives the decay as
`e^{-\alpha^2(L_z-|z|)^2}` with `|z| = \max|z_i - z_{j\pm}^{(l)}| \approx (M+1)H` once images are
present — consistent with the requirement `L_z > (M+1)H` imposed later at `:277`. Restate as
`e^{-\alpha^2(L_z-(M+1)H)^2}` (or say explicitly that `H` here means the *extended* thickness).

**B38. `:348` — "the dimensionless *padding ratio* `$R$`, defined as `$P = (L_z-H)/L_x$`".**
Symbol declared as `R`, defined as `P`. The rest of the thesis then mixes both: `R` at
`:357, 358, 364, 369, 374, 391`, `P` at `:349, 384`, and `app_icmewald2d.tex` uses `R` in one half
of a caption and `P` in the other (`:27, :34, :54, :61`). Also collides with `R` = summation
truncation parameter (Ch. 2 Theorem `order`) and `P` = batch size. Pick one new symbol.

### Chapter 4 — SOEwald2D / RBSE2D

**B39. `:20` — Eq. `eq:inverse_laplace` is missing a factor 2 in the exponent.**
`\frac{1}{2\pi\i}\int_\Gamma e^z\sqrt{\pi/z}\,e^{-\sqrt z|x|}dz = e^{-x^2/\mathbf{4}}`
(from `\mathcal L\{(\pi t)^{-1/2}e^{-a^2/4t}\} = z^{-1/2}e^{-a\sqrt z}` at `t=1`). To get `e^{-x^2}`
the integrand needs `e^{-2\sqrt z |x|}`.

**B40. `:185` and `:413` — a leftover `k` and a wrong index variable inside displayed equations.**
`\left(2\alpha s_l e^{-hz_{ij}} - 2\mathbf{k}e^{-\alpha s_l z_{ij}}\right)` should have `2h`
(verified: `\frac{1}{\alpha s_l+h}-\frac{1}{\alpha s_l-h} = \frac{-2h}{\alpha^2s_l^2-h^2}`).
Both equations also write `\sum_{\ell=1}^{M}` with summand `w_l`, `s_l` — `\ell` vs `l` mismatch.
Same `\ell`/`l` mismatch at `app_force.tex:16, 23`.

**B41. `:191` — Eq. `eq::Ul0SOE` self term is missing `1/\sqrt\pi`.**
`\mathcal G^0_\alpha(0) = 1/(\alpha\sqrt\pi)`, so the term is
`-\frac{\pi}{L_xL_y}\cdot\frac{Q}{\alpha\sqrt\pi} = -\frac{\sqrt\pi Q}{\alpha L_xL_y}`, not
`-\frac{\pi Q}{\alpha L_xL_y}`. Also the equation writes `\frac12\sum_i q_i\phi^{\bm 0}_\ell` (the
*exact* potential) for the SOE-approximated energy — should be `\approx` or `\phi^{\bm 0}_{\ell,\rm SOE}`.

**B42. `:413, :436` — `e^{-h^2/(4\alpha)}` should be `e^{-h^2/(4\alpha^2)}`** (two instances).
Same slip at `chp_rbe2d.tex:141, 225` (`\bm k \sim e^{-k^2/(4\alpha)}`).

**B43. `:409/:413` — `\widetilde\varphi(\bm h)` omits the `\frac{Q\pi}{hL_xL_y}{\rm erfc}(h/2\alpha)`
self term of Eq. `eq::pairssum8`,** yet `:409` asserts
`\sum_{\bm h\neq\bm 0}U^{\bm h}_{\ell,\rm SOE} = \sum_{\bm h\neq\bm 0}\widetilde\varphi(\bm h)`.

**B44. `:390` — `\mu := \sum \frac{f(\bm h)}{g(\bm h)}\mathbf{h(\bm h)}`.** The third factor should be
`g(\bm h)`; as written `h` is also `|\bm h|`.

**B45. `:423` — Eq. `eq::Happ` has `\alpha` where it should be `\alpha^2`, twice.**
Poisson summation gives
`S = \frac{\alpha^2 L_xL_y}{\pi}\sum_{m_x,m_y}e^{-\alpha^2(m_x^2L_x^2+m_y^2L_y^2)} - 1`
(check: the `m=0` term `\alpha^2L_xL_y/\pi` matches
`\frac{L_xL_y}{4\pi^2}\int e^{-h^2/4\alpha^2}d^2h = \frac{\alpha^2L_xL_y}{\pi}`).
Also `:425` says `m_\xi = L_\xi k_\xi/2\pi` — should be `h_\xi`, and `m` is the Poisson-*dual*
index, not `L h/2\pi`.

**B46. `:574, :586` and `chp_quasiewald.tex:748` — the normalization constant is called `$H$`,
which is the box height.** An earlier rename `H \to S` was left incomplete:
`:574` "by the definition of normalization factor `$H$`"; `:586` `\frac{2\lambda_D^4Q^2\mathbf{H}}{\pi P L_xL_y}`;
`chp_quasiewald.tex:748` `\frac{1}{\mathbf{H}}\sum_{\bm h_2}\ldots`. All three should be `S`.

**B47. `:586, :594, :602` — 3D spherical volume element used for 2D lattice sums.**
`4\pi h^2\,dh` is used where the chapter's own Eq. `eq::integral2` prescribes
`\frac{L_xL_y}{4\pi^2}\,2\pi h\,dh`. Concretely, Eq. `eq::happrx` (`:602`) is internally
inconsistent: with `4\pi h^2dh` the integral gives `2L_xL_y\alpha^3/\sqrt\pi` (an `\alpha^3`),
whereas the stated result is `2\alpha^2L_xL_y/\sqrt\pi`; the correct 2D answer is
`\alpha^2L_xL_y/\pi` (which agrees with Eq. `eq::Happ` once B45 is fixed).

**B48. `:308 → :317` — Eq. (4.x) `\sum_{\bm h\neq\bm 0}\mathscr E^{\bm h}_{U_\ell,\rm SOE}\le
\frac{2\lambda_D^2\alpha^3 Q}{\sqrt\pi}\varepsilon` cannot be reproduced from Eq. `eq::47` via
Eq. `eq::integral2`.** The `L_xL_y` in `eq::47` disappears and `\alpha` jumps from `\alpha^1` to
`\alpha^3`; the 2D conversion yields `\frac{\alpha}{2\pi}\int_{h_0}^\infty \frac{e^{-h^2/4\alpha^2}}{h}dh`
(logarithmic, `\propto\alpha`). Note also that in Eq. `eq::51` the first bound has no `L_xL_y`
while the second does, although both are energies. Please re-derive.

**B49. `:519, :570` and `app_force.tex:50` — dimensionally inconsistent bounds.**
`\left[3+\frac{\alpha}{\sqrt\pi}\left(1+\frac{2\sqrt2\varepsilon}{h}\right)\right]` and
`(3\sqrt\pi+\alpha)` add a dimensionless number to `\alpha` (units 1/length). Contrast the
well-formed `(1+2\alpha^2H)` in Eq. `eq::51`. `app_force.tex:50`'s `(1+2\alpha)L_z` has the same
problem *and* uses `L_z`, which does not exist in SOEwald2D (should be `H`).

**B50. `:550–552` — Eq. `eq::estiZ` does not follow.** The middle line replaces
`e^{\pm hz}\int_{\pm z}^\infty e^{-\alpha^2t^2-ht}dt` by
`e^{\pm hz}\int_{-\infty}^{\infty}e^{-\alpha^2t^2-ht}dt`; but
`\int_{-\infty}^{\infty}e^{-\alpha^2t^2-ht}dt = \frac{\sqrt\pi}{\alpha}e^{h^2/4\alpha^2}`, which
*grows*, and the un-bounded `e^{\pm hz}` remains. The correct move is the one used in Eq.
`eq::118`: keep `\int_{\pm z}^\infty` and use `e^{\pm hz-ht}\le1`. The final bound
`(2+\frac{2\alpha}{\sqrt\pi})h` also does not match what that gives.

**B51. `:578–580` vs Prop. `prop:unbiased` (`:502`) — a factor 2 discrepancy.**
`\frac{1}{2P}\sum_{h_1}\sum_{h_2}e_1e_2(\varphi_1-\varphi_2)^2` equals the Proposition's
`\frac{S}{P}\sum e|\varphi - \bar\varphi|^2`; the proof uses `\frac1P`, so the final bounds at
`:570` are a factor 2 loose. Also `[\varphi_1-\varphi_2]^2` should be `|\varphi_1-\varphi_2|^2`
(`\varphi` is complex). Same `2S/P` appears at `chp_quasiewald.tex:819`.

**B52. `:627` — `\sqrt{\frac{2\beta}{\beta}} = \sqrt2`.** The Langevin noise amplitude is broken;
`\beta` is being used both for the friction and for the inverse temperature. Compare the correct
form at `chp_rbe2d.tex:105`: `\sqrt{2\gamma k_{\rm B}T}`. Also `:630` calls `\beta` "the reciprocal
characteristic time" while `:673` calls it "the relaxation time" — contradictory.

**B53. `:272 / :496` — index/argument mismatches.**
`:272` `\mathscr E_{U_\ell}(r_c,\alpha)` should be `(k_c,\alpha)` (per Prop. `prop::2.12`), and
`\sum_{\bm k\neq\bm 0}\mathscr E^{\bm h}_{\ldots}` sums over `\bm k` with an `\bm h`-indexed
summand. Same in Eq. `eq::xichi` (`:496`) and Eq. `eq::Exi` (`:502`, `e^{-k_1^2/(4\alpha^2)}` should
be `h_1`). Also `:526` "`By the definition of \widetilde\varphi^{\rm RB}(\bm k)`" → `\bm h`.
`:352` uses cutoff `h_c` where the rest of the thesis uses `k_c`.

**B54. `:374` — the spurious `H`.** `2\pi r_c^2\frac{N}{L_xL_y}N = \frac{2\pi s^2N^2}{\alpha^2L_xL_y}`,
but the printed result has `\ldots/(\alpha^2 L_xL_y H)`. (The `O(N^{3/2})` conclusion is right — I
verified it from the `H`-free form.)

**B55. Algorithm 1 (`:259`) computes `U_{\rm p-s}` in the total energy but the line that computes it
(`:256`) is commented out.**

**B56. `:676` — "`choosing $P=100$ is adequate`" but Section `subsec::RBSE2D` uses
`P=20,30,60,120` and the CPU section uses `P=120`.** `P=100` appears nowhere.

**B57. `:800` — wrong panel labels.** The MSDs are panels (b) and (c) of `fig:norm` (per its caption
at `:773`), not (a) and (b).

**B58. `:841` — the causality is inverted.** "the reduced number of interacting neighbors … allowing
a much smaller real space cutoff `$r_c$`" — the smaller `r_c` is what reduces the neighbor count.

### Chapter 5 — RBE2D

**B59. `:25, :179` — `k_\ell^2` should be `k_\eta^2`** (the batch index is `\eta`; `\ell` is the
long-range superscript). Two instances.

**B60. `:4` — the chapter introduction promises content that is not in the chapter.**
"First, we derive the accurate error estimates for the Ewald summation for dielectrically confined
Coulomb systems" — that is Chapter 3. This chapter has only §"Random batch Ewald2D method" and
§"Numerical results".

**B61. `:12` — `\ref{subsec::convergence}` points into Chapter 4.** That label is defined at
`chp_soewald2d.tex:614`; this chapter's subsection is labelled `sec::convergence` (`:97`). Fix the
reference (and the inconsistent `sec::`/`subsec::` prefixes).

**B62. `:132` — "a volume of `$4\pi r_c^3$`"** should be `\frac{4\pi}{3}r_c^3`
(cf. `chp_soewald2d.tex:349`, which is correct).

**B63. `:162, :179` — `\bm F_{\rm IBC}^M` uses an acronym that is never expanded before use.**
"IBC" (infinite boundary correction) is first spelled out in `chp_conclusion.tex:31`. Define it here
(and in `app_force.tex`).

**B64. `:167, :204–205` — `\sum_{\substack{l=1\\ \text{even}}}^{M}`** is self-contradictory (`l=1`
is odd). Use `l=2, l\text{ even}` or `\sum_{l\ \rm even}`. Two instances.

**B65. `:375` vs `:346/:372` — "`over $70\%$ when $\sim4000$ CPUs`" exceeds the stated range of the
figure** ("`$N_{\rm core}$ up to $1024$`"). Same tension at `chp_applications.tex:59` ("up to 4000
CPU cores"), where the referenced figure `fig:timedie` uses 64 and 1024 cores.

### Chapter 6 — QEM

**B66. `\alpha` has a different meaning in this chapter than in Chapters 2–5.**
Here the screening Gaussian is `\frac{\alpha}{\pi}e^{-\alpha\|\bm\rho-\bm\rho'\|^2}` (`:20`), so
`\alpha` has units of 1/length² and the decay factor is `e^{-h^2/(4\alpha)}`; in Chapters 2–5 it is
`e^{-\alpha^2r^2}` with `\alpha` in 1/length and `e^{-h^2/(4\alpha^2)}`. Both conventions are called
"the Ewald splitting parameter `\alpha`". This needs an explicit statement (or, better,
unification), since expressions are compared across chapters.

**B67. `:192` — Eq. `eq...hat_G_solution...Zeroth` has the wrong sign.**
`-\partial_z^2\hat G = \delta(z-z_0)` gives `\hat G(\bm 0,z;z_0) = -\frac{|z-z_0|}{2}` (+ linear).
Confirmed by the `h\to0` limit of the non-zero mode: `\frac{1}{2h}e^{-h|z-z_0|} =
\frac{1}{2h} - \frac{|z-z_0|}{2} + O(h)`. This sign propagates into `S_0` (`:537`) and Eq.
`eq...text...Energy`.

**B68. `:344` — Eq. `eq::U_s_truncated` is missing the factor `\tfrac12`** present in the definition
Eq. `eq::U_s` (`:326`).

**B69. `:581–582, :592` — `S_0` algebra.** Line 582 equals **twice** line 581
(`2\sum_i q_i\sum_{j\le i}q_j(z_i-z_j)` vs `\sum_i q_i\sum_{j\le i}q_j(z_i-z_j)`), yet they are
joined by `=`. Separately, all three lines omit the prefactor `\frac{1}{4L_xL_y}` from the
definition at `:537`. Correct:
`S_0 = \frac{1}{2L_xL_y}\sum_i q_i(z_iT_1(i)-T_2(i))`.

**B70. `:620` — `S_4(\bm h)` is written with `\exp(-h|z_i-z_j|)`; the definition at `:541` has
`\exp(+h|z_i-z_j|)`.** The subsequent lines (621–624) use the `+` version, so `:620` is the typo.

**B71. `:688` — `h_q(\bm h_\eta)` should be `\beta_q(\bm h_\eta)`** (and `h` is already `|\bm h|`).

**B72. `:776, :821` — `P^2` should be `P`.** `\frac{2S^2}{P}\left(\frac{q_i\lambda_D^2}{2L_xL_y}X\right)^2
= \frac{S^2q_i^2\lambda_D^4X^2}{2\mathbf{P}L_x^2L_y^2}`. As printed (`1/P^2`) it also contradicts the
`O(1/P)` claim at `:824–825`.

**B73. `:795–802` — `\nabla_{\bm r_i}S_q` appears to be missing a factor 2.** `S_q` sums over both
`i` and `j`, and the index `i` appears in both slots, so
`\nabla_{\bm r_k}S_q = 2\sum_{j\neq k}q_kq_j\nabla_{\bm r_k}(\cdot)`. Also `:804` "the choice of `+`
or `-` depends on the choice of particle" should read "on the sign of `z_i-z_j`".

**B74. `:806` — `g(z_{ij})` is used without definition and collides with the sampling density
`g(\bm h)`** of Eq. `eq::hk_qem_qem` (`:670`).

**B75. `:897` — Eq. `eq::int_I` has a sign error and drops a factor.**
`\frac{1}{1-\gamma_{\rm u}\gamma_{\rm d}e^{-2hH}} = \frac{-1}{e^{-2H(h-k_0)}-1}`, so the integrand
needs an overall minus. The `[1-e^{-h^2/(4\alpha)}]` factor from Eq. `eq::int_Gs` is also silently
dropped — say so if intentional.

**B76. `:925` — sign of the constant term.** For `|f_1|\to|J_0(bh)e^{-ah}|` as `h\to\infty` the
constant must be `-J_0(bk_0)e^{-ak_0}` (matching `f_1` in Eq. `eq::split`); the printed `+` gives
`-J_0(bh)e^{-ah}+2J_0(bk_0)e^{-ak_0}`.
(The splitting Eq. `eq::split` itself is exact — I verified `f_1+f_2` reproduces the integrand —
and the analytic singular integral at `:939–941` is correct as a principal value.)

**B77. `:933` — `k_f` in the limit, `h_f` in the result;** and the second `<` should be `=`
(`\int_{k_f}^\infty e^{-ah}dh = e^{-ak_f}/a`).

**B78. `:920` — `\frac{\partial\{J_0(bh)e^{-ah}\}}{\partial\{e^{-2H(h-k_0)}\}}`** is not standard
notation for L'Hôpital; write the ratio of `d/dh` derivatives.

**B79. `:948` — the display ends with `\;.,`** (period followed by comma).

**B80. `:242` — "for `$z\in[0,z_0]$`"** should be `z\in[0,H]`; the piecewise formula that follows
covers both `z>z_0` and `z<z_0`.

**B81. `:339, :375, :515` — `k` in the denominator where `h` is meant** (`\frac{\ldots}{k(1-\gamma_{\rm u}\gamma_{\rm d}e^{-2hH})}`).
Also `k` for `h` at `:209, :234, :235, :276, :435–438, :534, :918, :899`.

**B82. `:1134` — "`for $|\gamma|<1$ cases`" contradicts the paragraph,** which is about
`|\gamma|>1` (`\gamma` from `-10` to `+10`).

**B83. `:1153, :1158, :1159, :1166` — `\gamma_1,\gamma_2,L_z` are not this thesis's notation.**
Should be `\gamma_{\rm u},\gamma_{\rm d},H` (`L_z` means the *padded* box height elsewhere).
Also `\ln\gamma_1\gamma_2` needs parentheses.

**B84. `:1148` — `\eps` and `\eps'` are undefined** (presumably `\eps_{\rm c}` and `\eps_{\rm d}`).

**B85. `:986` vs `:994` — text says charges `+1` and `-1`; the caption says "a pair of *negatively*
charged particles".**

**B86. `:1097` — the caption says "symmetric electrolytes" for Case I, the *asymmetric* dielectric
case** (`\gamma_{\rm u}=0.95, \gamma_{\rm d}=-0.95`).

**B87. Self-energy is never discussed.** `U_\ell` (Eq. `eq::U_l`, `:327`) has no prime and so
includes the `i=j`, `\bm m=\bm 0` term, while `U_s` excludes it. For a dielectric system the
self-*polarization* energy is physical but the bare self-interaction is not; Eq. `eq::Green`
(Ch. 2) sets this up but the chapter never applies it. Please state explicitly what is kept.

### Appendices

**B88. `app_math.tex:24` — Lemma `lem::2dfourier` has the wrong argument on the LHS.**
`\widetilde f(\rho,z) = 2\pi\int_0^\infty J_0(h\rho)f(\rho,z)\rho\,d\rho` — the LHS should be
`\widetilde f(h,z)`. Also this lemma is about *radially symmetric* functions (a Hankel transform),
but Chapter 2 cites it (`:428, :613, :627`) as "the 2D Fourier transform (see Lemma …)", which is
not what it states.

**B89. `app_math.tex:69–70` — Eq. `eq::E.4` is missing the factor `\kappa^2` from
`\rho_> = -\kappa^2\phi`, which is exactly what produces the spurious `\lambda_D^2`.**
`|\mathscr G| \le C_f\int|\rho_>|d\bm r = C_f\kappa^2\int_a^\infty|\phi|4\pi r^2dr = q C_f`,
not `qC_f\lambda_D^2`. The stated result is also dimensionally wrong (charge × length² vs charge).
Note the fallback bound at `:74` is `|\mathscr G|\le C_sC_fq_i` — no `\lambda_D^2` — consistent with
this. Since every DH bound in Chapters 4 and 6 carries `\lambda_D^2`/`\lambda_D^4`, this needs
resolving before those bounds can be trusted.

**B90. `app_math.tex:65` — `\kappa=\sqrt{\beta q^2\rho}` is missing `4\pi`.**
From `:56`, `\Delta\phi = 4\pi\beta q^2\rho_r\phi`, so `\kappa = \sqrt{4\pi\beta q^2\rho_r}`.
(`\rho` should also be `\rho_r`.) The DH solution at `:60–63` is otherwise correct — I verified
continuity at `r=a`.

**B91. `app_math.tex:29–33` — the symbol changes mid-definition** from
`\bm{\mathcal\psi}` to `\bm{\mathcal S}_i`. Also `\bm{\mathcal{\psi}}` is invalid markup
(`\mathcal` applied to a Greek letter).

**B92. `app_math.tex:43` — "`where $Q$ represents the total charge of the system`" is wrong.**
`Q=\sum_i q_i^2` (Ch. 2 `:495`); the total charge is zero by neutrality. The following line
(`\delta\psi\approx\zeta Q/\sqrt N`) only works with `Q=\sum q_i^2`, confirming the mislabel.

**B93. `app_math.tex:110` — VACF should be `\langle v_z(t_0)v_z(t+t_0)\rangle`;** as written it is
not a function of lag alone (cf. the MSD definition two lines above).

**B94. `app_math.tex:114` — "The `$z$-$z$` component *subtracted from* the pressure tensor"** →
"the `$z$-$z$` component *of* the pressure tensor".

**B95. `app_math.tex:17` — `\mathcal H^2` is never defined** and is a third name for the reciprocal
lattice (`\mathcal K^2` in Ch. 2, `\fK^2` in Ch. 6). `\mathcal H` is also the Hamiltonian
(`chp_soewald2d.tex:670`) and the old name for the normalization constant.

**B96. `app_math.tex:98` — "The sampling procedure has a small error since
`\mathcal H(\bm m^{\rm new})\approx q(\bm m^{\rm new}|\bm m^{\rm old})`" is misleading.**
Metropolis is exact regardless; `\mathcal H\approx q` gives a high *acceptance rate*, not a small
error. (The proposal machinery itself checks out: I verified the variance
`(\alpha L_\xi)^2/2\pi^2` and the closed-form `erf` expressions.)

**B97. `app_math.tex:120` — §`app::trapezoidal` is never referenced from the text,**
although Chapter 3 (`:65`) is where it is needed.

**B98. `app_force.tex:9` — sign.** `U^{\bm h}_{\ell,\rm SOE}` has prefactor `+\frac{\pi}{L_xL_y}`
(Eq. `eq::pairssum8`), so `\grad_{\V r_i}U^{\bm h}_{\ell,\rm SOE} = +\frac{\pi q_i}{L_xL_y}[\ldots]`,
not `-`. (Line 10, for the `\bm 0` mode, **is** correct because `U^{\bm 0}_{\ell,\rm SOE}` carries a
minus.)

**B99. `app_force.tex:23, :31` — the `{\rm sgn}(z_i-z_j)` factor is missing.**
`z_{ij}=|z_i-z_j|`, so `\partial_{z_i} = {\rm sgn}(z_i-z_j)\,\partial_{z_{ij}}`. (The
`\partial_{z_{ij}}` derivatives themselves are correct — I verified both.)

**B100. `app_force.tex:103` — `\mathcal M_z^M` sums over `i` but the summand uses `z_j`,
`z_{j\pm}^{(l)}`.**

**B101. `app_force.tex:108` — the ion–wall force appears to have a spurious factor `z_i`.**
For uniform plates the force is `z`-independent: from Eq. `eq:spectial`,
`F_z = -q_i\partial_z\phi_{\rm p-s} = -2\pi q_i(\sigma_u-\sigma_d)`. The printed expression has
`-2\pi z_i q_i[\ldots]`, which looks like the *energy*. Please check.

**B102. `app_force.tex:108` — mixed subscript conventions in one formula:** `\gamma_{\text{top}}`
and `\gamma_d` for the two interfaces. Same file uses `\sigma_d/\sigma_u` (`:66`) where Ch. 2 uses
`\sigma_{\rm bot}/\sigma_{\rm top}`, and `\gamma_l^{\pm}` (`:82–90`) where the rest of the thesis
uses `\gamma_{\pm}^{(l)}`. `:93` also defines `\bm r_{i,\bm n} = \bm r_i + \bm\ell` with `\bm\ell`
undefined (should be `\V L_{\bm n}`).

**B103. `app_quasiewald.tex:81/:87` — the table caption and the column headers disagree.**
Caption says the surface charge densities are `(\sigma_u,\sigma_d)`; the header columns are
`\sigma_1` and `\sigma_2`. The caption also says "`($\gamma_u$ and $\gamma_d$) of the bottom and top
walls`" — reversed (`\gamma_u` is the upper wall), and the columns are ordered `\gamma_d, \gamma_u`.

**B104. `app_quasiewald.tex:93–94` — rows (C) and (D) of Table `Table::Dielec` are not charge
neutral,** which Eq. `eq::chargeneutrality` requires. Row (C), 2:1: ions
`23\times2 - 42 = +4`; slabs `(0.01-0.02)\times42.9^2 = -18.4`; total `-14.4`. Row (D), 3:1:
`16\times3 - 44 = +4`; total `-14.4`. Rows (A) and (B) do check out
(`109\times2=218`, `73\times3=219`, neutral walls). Please verify the tabulated numbers.

**B105. `app_quasiewald.tex:142` — panel (b) is used twice in the caption,** once for the clusters
and once for the `z`-force; the body text (`:157`) says the forces are in panels (c) and (d).

**B106. `app_quasiewald.tex:102, :110` — "`the ensemble average of the total potential energy
$\mathscr P(U)$`"**: `\mathscr P(U)` is the *distribution* of `U` (`chp_rbe2d.tex:330`), not its
average. Also "difference *on*" → "difference *in*". And `:102` reads "Data are shown for the same
systems used in producing Figs.~\ref{fig:3_1} and $P=50$ and $100$" — garbled.

---

## C. Leftovers from the journal papers

These render as broken or meaningless statements in a thesis:

| Location | Problem |
|---|---|
| `app_icmewald2d.tex:3` | "Figures **1 to 4** in the main text", "Eq. **(52)** of the main text" — hardcoded numbers |
| `app_icmewald2d.tex:10` | "as described in Eq. **(54)** in the main text" |
| `app_quasiewald.tex:155` | "as shown in **Fig. 1 (a) of the main text**" → `\ref{fig:MD}(a)` |
| `chp_icmewald2d.tex:386` | "summarized in **Section 2 of the SI**" → the appendix |
| `chp_rbe2d.tex:183–184` | "found in the **SM, Section S1**" and "**Eq. (SM1.3)**" → `\ref{app::surfacecharge}` / `\eqref{eq::F_real_detail}` (both exist and are currently orphaned) |
| `chp_applications.tex:59` | "in Fig.~\ref{fig:timedie} **in the SM**" → the figure is in Appendix D |
| `chp_soewald2d.tex:34` | "throughout **this paper**" |
| `chp_soewald2d.tex:222` | "described in **this article**" |
| `chp_quasiewald.tex:113, 547` | "Throughout **this paper**", "hereafter in **this paper**" |
| `app_math.tex:102` | "of interest in **this paper**" |
| `app_force.tex:71` | "In **this manuscript**" |

Also: `chp_soewald2d.tex:1` calls the method **RBSE2D** while the chapter title, §3, the
introduction (`chp_introduction.tex:113`) and the abstract all call it **SOEwald2D**;
`chp_soewald2d.tex:830` writes **SoEwald2D**; `:840` writes **RBSOEwald**. Four spellings.

---

## D. Notation collisions worth a systematic pass

These are single symbols carrying incompatible meanings, sometimes within one chapter:

- **`\varepsilon` / `\epsilon` / `\eps`** — error tolerance, LJ energy scale, and permittivity.
  Worst case `chp_rbe2d.tex`: `\varepsilon` = tolerance (`:160`), `\varepsilon_{\rm top/c/bot}` =
  permittivity (`:222`), `\varepsilon_{\rm LJ}` (`:262`). Also `\epsilon` vs `\varepsilon` used for
  the *same* tolerance in `chp_icmewald2d.tex` (Eq. `eq:error_icmelc` vs Table 3.1 and Figs. 3.6–3.7).
- **`P`** — batch size; padding ratio (`chp_icmewald2d.tex:349, 384`; `app_icmewald2d.tex:27,34,61`);
  and `\mathcal P` = downsampling factor vs `\mathscr P` = sampling density (`chp_rbe2d.tex:42–43`).
- **`M`** — number of SOE terms (Ch. 4); image truncation level (Chs. 3, 5); integral truncation
  point (`chp_quasiewald.tex:435–458`); also `\mathcal M_z` = dipole moment (`app_force.tex`) and
  `\bm{\mathcal M}` = lattice vector (Ch. 2).
- **`\gamma`** — dielectric reflection factor vs Langevin friction (`chp_rbe2d.tex:105`).
- **`\beta`** — SOE exponent (`chp_soewald2d.tex:201`), thermostat parameter (`:627`), relaxation
  time (`:673`), inverse temperature (`app_math.tex`), and `\beta_{1:4}` auxiliary functions
  (`chp_quasiewald.tex:529`).
- **`\sigma`** — surface charge density; LJ/ion diameter; the split charge densities
  `\sigma_s,\sigma_\ell` (Ch. 6); and the length unit in `Table::Dielec`.
- **`\eta`** — batch index; parallel efficiency (`chp_rbe2d.tex:369`); MSD/VACF
  (`app_math.tex:104–110`).
- **`S`** — RBE normalization constant vs the auxiliary sums `S_{0:4}` (`chp_quasiewald.tex:534`).
- **`H`** — box height; former name of the normalization constant (B46); `\mathcal H^2` reciprocal
  lattice; `\mathcal H` Hamiltonian; `\mathcal H(\bm m)` sampling density (`app_math.tex:79`).
- **Box height** is `H` in Chs. 2 (quasi-2D), 3, 5, 6 but `L_z` in the Ewald3D subsection
  (`chp_preliminaries.tex:293`) and in the charged-slab subsection (`:588`), where `L_z` also means
  the padded height. Please separate these.
- **Ion diameter / length unit** is `\sigma` (Ch. 4), `r_{\rm d}` (Ch. 5 `:310`), `\tau_0` (Ch. 6) —
  and `\tau` is the LJ *time* unit in Chs. 4–5 while `\tau_0` is a *length* in Ch. 6.
  LJ time unit is `\tau` (Chs. 4–5) and `t_0` (`app_quasiewald.tex:147`).
- **Reciprocal vector** is `\bm k` in Ch. 2 §2.2 and Ch. 5, `\bm h` in Ch. 2 §2.3, Chs. 4 and 6 —
  and the two are mixed *within* single proofs (`chp_preliminaries.tex:428–451, 608, 623, 671`).
- **Force notation** `\bm F^i` (Chs. 2, 3) vs `\bm F_i` (Chs. 4–6).
- **Structure factor** `\rho(\V k)=\sum q_ie^{-\i\V k\cdot\V r_i}` (Ch. 2 `:295`) vs
  `+\i` (Ch. 3 `:47`, Ch. 5 `:32`).
- **Error symbol** `\mathscr E` (Chs. 2, 4) vs `\mathcal E_r` (Chs. 3, 6).
- **Big-O** `O(\cdot)` vs `\fO(\cdot)` (`chp_quasiewald.tex:517, 566, 824`).
- **Imaginary unit** `\i`, `\m{i}`, `\mathrm{i}`, `{\T i}` (`chp_icmewald2d.tex:65`), plain `i`
  (`chp_preliminaries.tex:527`) all appear.
- **Ångström**: `\mathring{A}`, `\mathring{\text{A}}`, and the unused `\Ai` macro.
- **`\V{r_0}` / `\bm{\rho_j}`** bold the subscript; should be `\V{r}_0`, `\bm{\rho}_j`.
  Many instances in `chp_quasiewald.tex` and `chp_preliminaries.tex:169`.
- **Units set in math italic**: `$1fs$`, `$298K$`, `$10fs$`, `$10\mathring{A}$`, `$1\,g/cm^3$`,
  `($fs$ to $ns$)`, `$MSD_z$`, `$MSD_{xy}$` — use `\mathrm{}` or `\text{}`
  (`chp_applications.tex:16, 18, 40, 56`; `chp_rbe2d.tex:336`; `chp_soewald2d.tex:802`).

---

## E. Content / consistency issues

- **`chp_introduction.tex:40–41` — the same sentence twice.** "In recent years, various strategies
  have been developed to effectively simulate quasi-2D systems under dielectric confinement."
  followed by "\rev{Several approaches have been developed to effectively simulate quasi-2D systems
  under dielectric confinement.}" Delete one.
- **`chp_introduction.tex:37 and :54` repeat the YB/ELC correction discussion nearly verbatim,**
  with inconsistent dashes ("Yeh--Berkowitz" vs "Yeh-Berkowitz").
- **`chp_introduction.tex:4` — garbled clause:** "achieved through confinement, bulk-like, and are
  modeled as periodic in the transverse `$xy$` directions".
- **`chp_introduction.tex:42` — "A common approach is by combining the image charge method (ICM)…"**
  — combining it with *what*? Sentence is incomplete.
- **`chp_introduction.tex:84` — "kernel truncation methods (TKM)"** — the acronym does not match;
  Vico et al.'s method is the *truncated kernel method* (TKM).
- **`chp_introduction.tex:31/51` and `:31/55`** define FMM and FFT twice each; RBE is defined at
  `chp_introduction.tex:31`, `chp_preliminaries.tex:4` and `:312`.
- **`chp_preliminaries.tex:426` — `\eqref{eq::xi20}` points into Chapter 4.** `\xi^{\pm}` is defined
  *in this chapter* at `eq::xi_def` (`:400`), which is never referenced. Same problem at `:531`:
  `\eqref{eq::integral2}` is defined at `chp_soewald2d.tex:312`. Both are forward references from
  Chapter 2 into Chapter 4 — fix to the local labels (or move the definitions).
- **`chp_soewald2d.tex:705` — "aligns with the analysis presented in Section~\ref{sec::ewald2d}"** —
  the `N`-independence argument is in §`subsec::errSOEwald2D` (`:330–334`), not §`sec::ewald2d`.
- **`chp_soewald2d.tex:1 vs :6`** — `O(N)` vs `O(N^{7/5})` in adjacent sentences without saying
  which method each refers to.
- **`chp_soewald2d.tex:11/17/90`** — `f(x)` is defined as `e^{-x^2}` at `:17` and as
  `e^{-\alpha^2x^2}` at `:90`, but Eq. `eq::SOE1` uses `e^{-s_l|x|}` (no `\alpha`), matching the
  former. State once that `(w_l,s_l)` approximate `e^{-x^2}`.
- **`chp_applications.tex:58 vs the caption at :48`** — "more than 10 times faster" vs "about 1.5
  orders of magnitude faster" (≈32×) for the same comparison.
- **`chp_applications.tex:111` vs `chp_quasiewald.tex:1156`** — the fitted lattice-constant ratios
  `1.2\pi` and `1.4\pi` are called "consistent with our theoretical prediction", but the prediction
  is `\lambda k_0 = 2\pi`. The lattice constant and the field wavelength are different quantities
  (line 112 relates the former to the *second zero* of the induced-charge profile) — please
  qualify.
- **`chp_applications.tex:94` — "Anderson thermostat"** → **Andersen** (Hans C. Andersen).
- **`chp_applications.tex:104` — "likely-charged particles"** → "like-charged".
- **`chp_conclusion.tex:42` — "For negatively **charged** systems"** → negatively **confined**
  (negative permittivity). Changes the meaning.
- **`chp_conclusion.tex:31`** is where IBC is first expanded, after its use in Chapter 5 (B63).
- **`app_softwares.tex`** lists ExTinyMD.jl, EwaldSummations.jl, ChebParticleMesh.jl but **not
  QuasiEwald.jl**, which is what Chapter 6 cites (`\cite{QuasiEwald}`) and what
  `chp_soewald2d.tex:829` alludes to as "a self-developed package".
- **`app_quasiewald.tex:77/:128`** — described as an RBE2D-vs-PPPM comparison, but `fig:timedie`'s
  own caption says ICM-PPPM at `\gamma_{\rm u}=\gamma_{\rm d}=-0.9`, whereas
  `chp_applications.tex:39` sets `\gamma=-0.5` for the water system it is cited for.
- **`chp_icmewald2d.tex:123`** — "Note that the contribution of the `\bm h=\bm 0` mode is contained
  in the homogeneous term, not relevant in the error analysis here" is both ungrammatical and, as
  far as I can tell, not right: Eq. `eq::33_icm_icm`'s image series has its own nonzero `\bm h=\bm 0`
  contribution which is not part of the first (homogeneous) term. Please clarify.
- **Build warnings to clear:** 2× `Package amsmath Warning: Cannot use 'split' here` (likely the
  `\resizebox{...}{\begin{split}}` constructs in `chp_icmewald2d.tex:268, 291, 429` or
  `app_force.tex:79`), and a duplicate-hyperref-destination warning for every figure
  (`figure.2.1`…`figure.6.14`). The two stray `\label{fig:compare}` / `\label{fig:msd_norm}` inside
  caption-less minipages at `chp_soewald2d.tex:765, 770` should be deleted (they are unused and
  attach to the wrong counter).
- **`chp_preliminaries.tex:166, 478, 510` and `chp_soewald2d.tex:242`** — bare `\ref{app::…}`
  without "Appendix"/"Sec.", so they render as "in A.3". Add the word.
- **Cross-reference style** is `Fig.~` in Chs. 3, 5–7 but `Figure~` in Ch. 4 (`:701, 703, 712, 721,
  729, 740, 743`). Pick one.
- **`chp_soewald2d.tex:599`** uses `\ref` instead of `\eqref` for an equation.

---

## F. Language and typography (line-referenced)

### Confirmed misspellings / wrong words
| File:line | Fix |
|---|---|
| `chp_quasiewald.tex:1027` | simluation**s** → simulations |
| `chp_quasiewald.tex:1030` | shifted-**trunacted** → truncated |
| `chp_quasiewald.tex:901` | **principle** value → principal value |
| `chp_quasiewald.tex:342` | in **actually** calculations → in actual calculations |
| `chp_quasiewald.tex:1091` | **an** strongly → a strongly |
| `chp_quasiewald.tex:1159, 1168` | label `fig:k_wavelegth` → `wavelength` |
| `chp_rbe2d.tex:288` | **estimination** → estimates |
| `chp_conclusion.tex:13` | over-looked → overlooked |
| `chp_applications.tex:94` | Anderson → Andersen |
| `chp_applications.tex:104` | likely-charged → like-charged |
| `chp_icmewald2d.tex:70`, `chp_soewald2d.tex:32` | "Specially," → "Specifically," / "In particular," |
| `chp_soewald2d.tex:396` | stray `ee` after the colon: "random batch sampling:**ee**" |
| `chp_soewald2d.tex:275` | "The **remainder** two terms" → remaining |
| `chp_quasiewald.tex:1038` | "the **charged** distribution" → charge distribution |
| `chp_quasiewald.tex:1059` | "**simulations** time" → simulation time |
| `chp_soewald2d.tex:279, 283` etc., `chp_soewald2d.tex:512, 638`, `app_math.tex:47–48`, `app_force.tex:60` | `Debye-H$\ddot{\text{u}}$ckel` → `Debye--H\"uckel` (math-mode umlaut inside a name) |
| `chp_preliminaries.tex:642` | double closing paren: `(z_{i} - 0))` |
| `chp_icmewald2d.tex:417` | `$L_z=$$45$` — doubled `$$` |
| `chp_quasiewald.tex:948` | `\;.,` |

### "independent with" → "independent of"
`chp_icmewald2d.tex:336`; `chp_soewald2d.tex:740, 743`.

### "the force exert(s) on" → "exerted on"
`chp_soewald2d.tex:240, 735`; `chp_rbe2d.tex:22`.

### "As been mentioned/introduced" → "As mentioned/introduced"
`chp_preliminaries.tex:265, 313`; `chp_rbe2d.tex:17`; `chp_quasiewald.tex:1081`.

### Comma splices and run-ons
`chp_introduction.tex:11`; `chp_preliminaries.tex:123, 279, 346, 386`; `chp_icmewald2d.tex:312`;
`chp_soewald2d.tex:387`; `chp_rbe2d.tex:19`; `chp_quasiewald.tex:135, 968`; `app_math.tex` passim.

### Missing verbs / sentence fragments
| File:line | Text |
|---|---|
| `chp_preliminaries.tex:269` | "…technique~\cite{Ewald1921AnnPhys}**.** which is the most widely used…" — period should be a comma |
| `chp_preliminaries.tex:586` | "we present the Ewald2D formulation **with** such a situation can be well treated" — needs "with which" |
| `chp_soewald2d.tex:790` | "where the ions **confined** within the central region" — needs "are" |
| `chp_quasiewald.tex:869–877` | "Then by taking [Eq.] where `\varrho` represents the surface density." — no main verb |
| `chp_quasiewald.tex:1090` | "The spatial distribution of ions in `$z$` **as shown in** Fig. …" — needs "is" |
| `chp_quasiewald.tex:1108` | "when `$\gamma$` greater or smaller than `$0$`" — needs "is" |
| `chp_quasiewald.tex:1075` | "The blue and red dots **corresponding** CPU time of…" — needs "correspond to" |
| `chp_rbe2d.tex:160` | "let `\bm F` **denotes**" → denote |
| `chp_soewald2d.tex:676` | "is system-dependent, **relies on** performing…" — needs "and relies" |

### Subject–verb agreement
`chp_preliminaries.tex:476` ("truncation … are performed"), `:485` ("one only calculate"),
`:193` ("there exists two solutions"), `:610` ("each components");
`chp_soewald2d.tex:139` ("Achieving these properties are"), `:330` ("Theorems~\ref{}" for one),
`:340` ("various sorting algorithm are"), `:835` ("costs … scale as … which is");
`chp_icmewald2d.tex:204` ("the error in force converge"), `:216` ("increases" for two subjects);
`chp_quasiewald.tex:1043` ("the MSD stop"), `:1067` ("the same parameter … are applied";
"parameter are selected"), `:1152` ("field lines has"), `:130/:132` ("reads" for two objects);
`app_math.tex:4` ("Fourier transform" for two).

### Missing articles ("the")
`chp_preliminaries.tex:124, 215, 297, 299, 645`; `chp_soewald2d.tex:366, 642, 645, 664`;
`chp_quasiewald.tex:878`; `chp_applications.tex:104`. Recurring pattern:
"singularity of Coulomb kernel", "long-range nature of Coulomb interaction".

### Punctuation after displayed equations
Several displays end in `.` where the sentence continues with a lowercase "where"/"and", or in `,`
where a new sentence begins: `chp_preliminaries.tex:465/467`, `:281`;
`chp_quasiewald.tex:546/547`, `:692/694`, `:929/931`; `chp_icmewald2d.tex:295/299`.

### Dashes
- Two-author names need en-dashes: Yeh–Berkowitz, Debye–Hückel, Nosé–Hoover,
  Martyna–Tuckerman–Tobias–Klein, Fourier–Chebyshev, Gauss–Legendre/Laguerre/Hermite (Ch. 6 mixes
  `-` and `–` within one paragraph: `:488, 492, 498, 507`), Poisson–Boltzmann, Rotne–Prager–Yamakawa.
- Page ranges: `pp.~109-112` (`chp_preliminaries.tex:514`), `pp.~374-375`
  (`chp_soewald2d.tex:45`) → `--`. Also `pp.~688` for a single page → `p.~688`.
- `chp_conclusion.tex:61` — "phase separation **-** scenarios where" → em dash.

### Miscellaneous
- `chp_introduction.tex:97` — "Thesis outline**s**" → outline.
- `chp_introduction.tex:3` — `mazars2011long` is cited twice in one sentence.
- `chp_introduction.tex:12` — "sub-Angstrom" → sub-Ångström / sub-angstrom.
- `chp_preliminaries.tex:24` — "materials properties" → material properties.
- `chp_preliminaries.tex:129` — "with `$R\in\mathbb R$` **be** a truncation parameter" → being.
- `chp_preliminaries.tex:133` — `);and` → `); and`.
- `chp_preliminaries.tex:151` — `$\abs{\V m}=0, 1, 2\ldots, R$` — missing comma after `2`.
- `chp_preliminaries.tex:152, 176` — `z\to\infty` should be `z\to\pm\infty`.
- `chp_preliminaries.tex:162–164` and `chp_quasiewald.tex:359` — `\{\bm\rho\in[0,L_x]\times[0,L_y]\}`
  is set-builder abuse; write `:= [0,L_x]\times[0,L_y]`.
- `chp_preliminaries.tex:196` — needs brackets around the integrand and `|\nabla u|^2` rather than
  `(\nabla u)^2`; also state that `\Delta u = 0` is why the first integral vanishes.
- `chp_preliminaries.tex:304` — "Jin, {\it et al.}" — drop the comma before *et al.*; same at
  `chp_soewald2d.tex:26`.
- `chp_preliminaries.tex:552` — `k_c/2\alpha` → `k_c/(2\alpha)`; `chp_quasiewald.tex:899` —
  `\ln(\gamma_u\gamma_d)/2H` → `/(2H)`.
- `chp_preliminaries.tex:485` — `\bm{h }` (stray space) and `\mathcal M` should be `\V{\mathcal M}`.
- `chp_preliminaries.tex:580` — "catastrophic error cancellation" → "catastrophic cancellation".
- `chp_soewald2d.tex:337` — "irrelevant to the system size" → "independent of".
- `chp_soewald2d.tex:342, 829` — `\rho_s`/`\rho_\ell` collide with `\rho_{\bm k}` and `\bm\rho`;
  "a self-developed package developed in Julia language" is redundant and needs "the".
- `chp_soewald2d.tex:780` — "Lennard Jones" → Lennard-Jones (and LJ is already defined at `:645`).
- `chp_soewald2d.tex:791` — "contains 218 cations and anions" — ambiguous; say "218 of each".
- `chp_soewald2d.tex:792` — `$\Delta=0.001\tau$` → `$\Delta t$`.
- `chp_soewald2d.tex:793` — "relaxation time**s** `$0.1\tau$`"; "the **reduce** external temperature"
  → reduced.
- `chp_icmewald2d.tex:7` — "For dielectric confined systems, … can be extended to accommodate for
  dielectric-confined systems" — says it twice; also "accommodate for" → "accommodate".
- `chp_icmewald2d.tex:107` — "doubly-periodization" is not a word.
- `chp_icmewald2d.tex:210` — `s` is called "the Ewald splitting parameter", but that is `\alpha`
  (`chp_preliminaries.tex:279`); `s` is the accuracy parameter.
- `chp_icmewald2d.tex:404` — "in figure 2 of Ref~\cite{}" → "Fig. 2 of Ref.~\cite{}".
- `chp_icmewald2d.tex:444` — "(since there is no image charge reflection)" contradicts `M=1`;
  presumably "no *repeated* reflection".
- `chp_icmewald2d.tex:341–343` uses `\(...\)` inline math while the rest of the thesis uses `$...$`.
- `chp_rbe2d.tex:465` — the un-`\rev`'d markup leaves runs of 2–3 spaces in ~25 places in this file
  (harmless in output, but worth cleaning: `grep -n '  {' chapters/chp_rbe2d.tex`).
- `chp_quasiewald.tex:1106` — Unicode curly quotes `“pairing”` where the rest of the thesis uses
  ` ``…'' `.
- `chp_preliminaries.tex:102` — space before the comma: `\partial \Omega_{d}$ , respectively`.
- `chp_preliminaries.tex:32–34` — the tikz comments label the axes "y-axis"/"z-axis" but the nodes
  draw `$z$`/`$y$`; `:36` is missing `-- cycle` (the green twin at `:85` has it).
- `00II_acknowledgments.tex:3` — "Completing this thesis **is** a long journey" → has been.
- `00III_abstract.tex:5` — "quasi-Ewald **M**ethod" → method (lowercase, as elsewhere).
- Serial commas are inconsistent (the preamble sets `\creflastconjunction` to ", and~", suggesting
  Oxford style): `00III_abstract.tex:2`; `chp_preliminaries.tex:13, 573`;
  `chp_conclusion.tex:71`; `chp_introduction.tex:51` (spurious comma before "or").
</content>
</invoke>
