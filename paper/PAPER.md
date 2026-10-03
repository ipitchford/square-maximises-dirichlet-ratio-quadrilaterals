---
title: "The square maximises the ratio of the first two Dirichlet eigenvalues among quadrilaterals"
author: Anonymous
date: 2026-10-03 (version 0.1.0-candidate)
---

# Abstract

Let $\lambda_1(Q) < \lambda_2(Q)$ be the first two Dirichlet eigenvalues of the Laplacian on a bounded planar domain $Q$
and $\xi(Q) = \lambda_2(Q)/\lambda_1(Q)$. We prove that $\xi(Q) \le 5/2$ for every simple quadrilateral $Q$, convex or not, with equality only
for squares. This is the $k=4$ case of the polygonal Payne–Pólya–Weinberger conjecture stated by Siudeja and by Antunes and Freitas; the case $k=3$
(equilateral triangle, ratio $7/3$) was proved by Siudeja and by Arbon, Mannan, Psenka and Ragavan. The proof is computer-assisted.
(i) Two analytic lemmas remove near-triangles and thin shapes (width/diameter $\le 1/10$) and confine the problem to a compact parameter set.
(ii) A certified local step establishes the inequality on an explicit neighbourhood of the square, via an inertia criterion for a
$3\times3$ Rayleigh–Ritz pencil of transported eigenfunctions and a Kato–Temple bound expanded as Taylor models.
(iii) Certified box covers of the remaining parameter set bound $\lambda_2$ by Rayleigh–Ritz and $\lambda_1$ by a Lehmann–Goerisch inequality with an
exactly divergence-matched flux; away from the square, the test functions and the flux are transported from the box centre by a piecewise-affine
map and a contravariant Piola transform, so that all parameter dependence is explicit. (iv) The band $1/10 <$ width/diameter $\le 2/5$ is covered in
diameter-normalised charts by one-dimensional bounds: a Born–Oppenheimer lower bound for $\lambda_1$ certified by Sturm oscillation, and a weighted
one-dimensional Ritz upper bound for $\lambda_2$. All certificates are computed in ball arithmetic (Arb) or exact rational arithmetic; every cover is
replayed by a fail-closed checker. For convex quadrilaterals, coverage, exclusions and samples of certificates are re-checked by
separately written programs; for non-convex ones, coverage is.
(v) Non-convex quadrilaterals are handled in a separate chart centred on the reflex vertex. Analytic lemmas cover the degenerate and
unbounded parts of that shape space: shallow notches, thin wedges, flat hulls and thin stems, and long arms through
Dirichlet–Neumann bracketing and one-dimensional bounds with exponential tails. Certified covers of compact families, in charts
uniform in the arm lengths, cover the rest.

# 1. Introduction and statement

A *simple quadrilateral* is the open region bounded by a closed four-sided polygon without self-intersections whose four interior
angles all differ from $\pi$. (If an angle equals $\pi$ the region is a triangle, and $\xi\le 7/3 < 5/2$ by the triangle theorem of
Arbon, Mannan, Psenka and Ragavan (2022).)

**Theorem 1.** For every simple quadrilateral $Q \subset \mathbb R^2$, convex or not, $\lambda_2(Q)/\lambda_1(Q) \le 5/2$, with equality if and only if
$Q$ is a square.

*Proof.* A simple quadrilateral is either convex or has exactly one reflex vertex (a dart). For convex $Q$ this is Theorem 1$'$ (§4.6);
for darts, Theorem 5 (§5) gives the strict inequality. $\square$

**Context.** Payne, Pólya and Weinberger (1956) proved $\xi \le 3$ for all planar domains and conjectured that the disk is the
maximiser ($\xi = j_{1,1}^2/j_{0,1}^2 \approx 2.539$); Ashbaugh and Benguria (1991, 1992) proved this. Siudeja (2010, Conjecture 1.2;
preprint 2007) and, on numerical evidence, Antunes and Freitas (2008) conjectured that among $k$-gons the regular $k$-gon maximises $\xi$
(see also the survey of Laugesen and Siudeja 2017); Theorem 1 settles the case $k=4$. For triangles the conjecture was proved by
Siudeja (2010, acute triangles) and Arbon, Mannan, Psenka and Ragavan (2022, all triangles); Arbon (2022) gave a
second, computer-assisted proof of the triangle case and proved that the square is a strict local maximiser among quadrilaterals
in an explicit but tiny neighbourhood. Endo and Osting (2026) proved the Neumann analogue (the square maximises the first
nonzero Neumann eigenvalue among convex quadrilaterals of fixed area) by a certified box covering; our global step follows their
architecture, but the Dirichlet ratio needs a sharp *lower* bound on $\lambda_1$, which is the main technical difference.

**On the existing local result.** Arbon's Theorem 1.3 is stated for all quadrilaterals sufficiently close to the square, with the
explicit region $(t,s) \in [0,0.0015]^2$ in his vertex parametrisation. Two steps of the argument do not hold as written.
First, equation (3.16) there gives the first-order variation of $\xi$ for a simultaneous perturbation of two vertices as the sum of the
first-order variations for the two vertices separately; for the degenerate second eigenvalue this is not valid, because the first-order
shift under a sum of perturbations is the smallest eigenvalue of the *sum* of the two $2\times 2$ first-order matrices, not the sum
of the smallest eigenvalues. An exact calculation shows the difference. Moving one vertex of the unit square changes the boundary
along its two edges with a linear normal speed, and Hadamard's formula (Henrot 2006, Ch. 2; Henrot and Pierre 2018) gives a $2\times2$
first-order matrix $M_k$ for the double eigenvalue $5\pi^2$ (its eigenvalues are the one-sided derivatives of the double eigenvalue) and a first-order change $\dot\lambda_{1,k}$ of $\lambda_1$. For the isosceles trapezoid with parallel sides $1\pm t$ (direction
$b$ of §2.1) every vertex moves by $\pm t/2$ horizontally, and
$$M_1 = \begin{pmatrix}-\pi^2/2&-32/9\\-32/9&-2\pi^2\end{pmatrix},\quad M_2 = \begin{pmatrix}-\pi^2/2&32/9\\32/9&-2\pi^2\end{pmatrix},\quad
M_3 = \begin{pmatrix}\pi^2/2&32/9\\32/9&2\pi^2\end{pmatrix},\quad M_4 = \begin{pmatrix}\pi^2/2&-32/9\\-32/9&2\pi^2\end{pmatrix},$$
with $\dot\lambda_{1,k} = -\pi^2/2,-\pi^2/2,\pi^2/2,\pi^2/2$. The sums vanish, so the exact first-order term of $\xi$ along this family is $0$
(Lemma 3.2). Adding instead the separate first-order variations of $\xi$, as (3.16) does, gives $-\sqrt{16384+729\pi^4}/(18\pi^2)\approx-1.664$ per
unit $t$ when the four vertex moves are added, and $-3/2$ when the move is split into its two simultaneous two-vertex moves. For a single
two-vertex move (the two upper vertices moved towards each other) the true slope is $-3/4$, while the rule gives
$-\sqrt{16384+729\pi^4}/(36\pi^2)\approx-0.832$. (Exact arithmetic in SymPy; `referee_response/arbon_perturbation.py`. Summing the $M_k$
reproduces the matrices of Lemma 3.2 for all four directions.)
Second, the bound "$O(t^2) \le 295 t^2$" on the second-order remainder in the triangle case (and the analogous constant $336.972$ in the
quadrilateral case) is inferred from a crude global bound on $\xi$ over a large interval of $t$; a bound on the size of a remainder
term at one value of $t$ does not control its behaviour near $t=0$. The local statement itself is true (it follows from Theorem 2 below);
only the argument needs replacing.

**Notation.** $\xi = \lambda_2/\lambda_1$ throughout. In §§3–4, $\beta$ denotes a certified lower bound for $\lambda_2$ (the shift of the Temple and
Lehmann–Goerisch inequalities); in §5 and Supplement S, $\beta$ is the notch depth of a dart. $\rho$ is a Rayleigh quotient in §3 and
$\rho_A,\rho_C$ are arm lengths in §5; $\rho_0$ is the radius of the Taylor-model ball in §3 and $\mu_0$ the flux constant in §4.2.

# 2. Parametrisation and reductions

## 2.1 Shape parameters
We use the parametrisation of Endo and Osting. For $p = (a,b,c,d) \in \mathbb R^4$ let
$$v_1 = (0,0)+\tfrac12(a+d-b,\,-a+d+c),\quad v_2 = (1,0)+\tfrac12(-a+d+b,\,-a-d-c),$$
$$v_3 = (1,1)+\tfrac12(-a-d-b,\,a-d+c),\quad v_4 = (0,1)+\tfrac12(a-d+b,\,a+d-c),$$
and $Q_p$ the quadrilateral with these vertices in this cyclic order. $Q_0$ is the unit square; $a$ is a rectangle mode, $b,c$ are
isosceles-trapezoid modes, $d$ a rhombus mode; $|Q_p| = 1-a^2-d^2$; $Q_p$ is a nondegenerate convex quadrilateral iff the four cross
products $c_i(p) = (v_{i+1}-v_i)\times(v_{i+2}-v_i)$ are positive. Every convex quadrilateral is similar to some $Q_p$ (Endo–Osting,
Theorem 2.1), and $\xi$ is invariant under similarities, so it suffices to prove $\xi(Q_p) < 5/2$ for all $p \ne 0$ with all $c_i(p) > 0$.
The dihedral group $D_4$ acts on parameter space by the signed permutations
$\sigma_h : (a,b,c,d)\mapsto(a,b,-c,-d)$, $\sigma_v : (a,b,c,d)\mapsto(a,-b,c,-d)$, $\sigma_r:(a,b,c,d)\mapsto(-a,-c,-b,d)$,
and $\xi(Q_p)$ is invariant; a fundamental domain is $\{a\ge 0, b\ge 0, c \ge 0\}$.

## 2.2 Lemma A (near-triangle quadrilaterals)
For a convex quadrilateral with vertices $v_1,\dots,v_4$ write $c_i = (v_{i+1}-v_i)\times(v_{i+2}-v_i) > 0$ (indices mod 4), twice the area
of the triangle $v_iv_{i+1}v_{i+2}$.

**Lemma A.** If $c_i \le \eta\, c_{i+2}$ for some $i$, then $\xi(Q) \le \tfrac73(1+\eta)^2$. In particular $\xi(Q) < 5/2$ when $\eta \le 0.035$
(since $\tfrac73\cdot 1.035^2 = 2.4995\ldots$).

*Proof.* Put $T = \mathrm{conv}(v_{i+2},v_{i+3},v_i)$ and let $\ell$ be the line through $v_i$ and $v_{i+2}$. The distances of $v_{i+1}$ and
$v_{i+3}$ to $\ell$ are $\delta = c_i/|v_{i+2}-v_i|$ and $h = c_{i+2}/|v_{i+2}-v_i|$, so $\delta/h = c_i/c_{i+2} \le \eta$. Let $T'$ be the image of $T$
under the homothety with centre $v_{i+3}$ and ratio $1+\delta/h$; it is bounded by the lines $v_{i+3}v_i$, $v_{i+3}v_{i+2}$ and the line parallel
to $\ell$ at distance $\delta$ on the far side. Convexity of $Q$ at $v_i$ and $v_{i+2}$ means $(v_i - v_{i+3})\times(v_{i+1}-v_{i+3}) \ge 0$ and
$(v_{i+1}-v_{i+3})\times(v_{i+2}-v_{i+3})\ge0$, i.e. $v_{i+1}$ lies in the closed angle at $v_{i+3}$ spanned by $v_i$ and $v_{i+2}$; it lies on the far
side of $\ell$ at distance $\delta$, hence on the base of $T'$ (the segment of the line parallel to $\ell$ at distance $\delta$, cut by that angle).
So $v_{i+1}\in T'$ and $Q = \mathrm{conv}(T\cup\{v_{i+1}\}) \subset T'$. Domain monotonicity gives
$\lambda_2(Q)\le\lambda_2(T)$ and $\lambda_1(Q)\ge\lambda_1(T') = \lambda_1(T)/(1+\delta/h)^2$. Hence
$\xi(Q)\le(1+\delta/h)^2\,\xi(T)\le\tfrac73(1+\eta)^2$ by the triangle theorem $\xi(T)\le7/3$ for every triangle (Arbon, Mannan, Psenka and Ragavan 2022;
Siudeja 2010 for acute triangles). $\square$

## 2.3 Lemma B (thin convex domains)
**Lemma B.** Let $\Omega$ be a bounded convex planar domain with minimal width $w$ and diameter $D$, $s = w/D$. Then for every
$\theta\in(0,1)$
$$\xi(\Omega) \le \frac{4s^2}{\theta^2} + \frac{1}{(1-\theta)^2}.$$
With $\theta = 1/4$: $\xi(\Omega)\le 64s^2 + 16/9$, which is $\le 2.4178 < 5/2$ for $s \le 1/10$ (optimising $\theta$ extends this to $s\le 0.1065$).

*Proof.* Work with the compact convex set $K = \overline\Omega$ (same width and diameter). (i) $\Omega$ lies in an open strip of width $w$, so
$\lambda_1(\Omega)\ge\pi^2/w^2$.
(ii) Let $AB$ be a diameter of $K$; put $A = (0,0)$, $B = (D,0)$. Then $\Omega\subset\{0<x<D\}$ (a point with $x<0$ would be farther than $D$ from $B$),
and the vertical chords over $x = 0$ and $x = D$ reduce to $A$ and $B$. Write $\Omega = \{(x,y): \underline y(x) < y < \overline y(x)\}$ with $\overline y$ concave,
$\underline y$ convex, $\overline y(0)=\overline y(D)=\underline y(0)=\underline y(D)=0$, and let $H = \overline y(x^*)-\underline y(x^*)$ be the longest vertical chord.
(iii) $H\ge w$: the longest chord of a convex body in direction $u$ equals the radial function $\rho_{\Omega-\Omega}(u)$ of the difference body,
the width in direction $v$ equals its support function $h_{\Omega-\Omega}(v)$, and for a centrally symmetric convex body
$\min_u\rho = \min_v h$; hence every directional maximal chord is at least the minimal width.
(iv) By concavity of $\overline y$ and convexity of $\underline y$ (both vanishing at $0$ and $D$), on the window
$W = [(1-\theta)x^*,\, x^* + \theta(D-x^*)]$ of length $\theta D$ we have $\overline y \ge (1-\theta)\overline y(x^*)$ and $\underline y\le(1-\theta)\underline y(x^*)$. So the closed rectangle
$R = W\times[(1-\theta)\underline y(x^*),(1-\theta)\overline y(x^*)]$, of sides $\theta D$ and $(1-\theta)H \ge (1-\theta)w$, lies in $K$; its interior lies in
$\mathrm{int}\,K = \Omega$, and
$\lambda_2(\Omega)\le\lambda_2(R)\le\pi^2\big(4/(\theta D)^2 + 1/((1-\theta)w)^2\big)$.
Dividing (iv) by (i) gives the claim. $\square$

(Numerical remark, not used in any proof: in finite-element experiments the best inscribed rectangle fails at $s\approx0.107$ for every
shape family tried, so Lemma B is essentially sharp for this method; on 5066 samples the bound exceeds the computed ratio by at least
0.21 for $s\le0.1$.)

## 2.4 What remains
By §2.1–2.3 and the local theorem of §3 it suffices to prove $\xi(Q_p) < 5/2$ for $p$ in the fundamental domain $\{a,b,c\ge0\}$ with all
$c_i(p) > 0$, $c_i(p) > 0.035\,c_{i+2}(p)$ for all $i$, $\mathrm{width}/\mathrm{diam} > 0.1$, and $|p|_\infty > 1/20$. Since a convex domain contains the quadrilateral spanned by a diameter and the longest chord perpendicular to it, $|Q| \ge wD/2$, so
$D^2 \le 2|Q|/s < 20|Q|$ on the remaining set. Each diagonal is at most $D$, and the Endo–Osting identity
$|v_1-v_3|^2+|v_2-v_4|^2 = 2(1-d)^2+2a^2+2(1+d)^2+2a^2 = 4a^2+4(1+d^2)$ (direct computation from the vertex formulas) gives
$4a^2+4+4d^2 \le 2D^2 < 40(1-a^2-d^2)$, i.e. $a^2+d^2 < 9/11$, so $|a|,|d| < 0.9045$. Independently of thinness, convexity bounds $b$ and $c$:
with $S = 1-a^2-d^2 = |Q_p| > 0$, $X = b+c$, $Y = b-c$, the four conditions $c_i(p)>0$ read $|X(1-d)+aY| < S$ and $|Y(1+d)+aX| < S$;
inverting $\big(\begin{smallmatrix}1-d&a\\a&1+d\end{smallmatrix}\big)$ (determinant $S$) gives $|X| < 1+d+|a|$ and $|Y| < 1-d+|a|$, hence
$|b|,|c| < 1+|a| < 1.905$. The remaining set therefore lies in the box $a\in[0,0.905]$, $b,c\in[0,4.36]$, $|d|\le0.905$ used as the
root region of the covers in §4 (the $b,c$ range is generous; all endpoints are rounded outward in floating point). (Monte Carlo: the remaining set has volume $\approx0.29$ and extent $a\le0.83$, $b\le0.85$, $c\le1.27$, $|d|\le0.8$.)

# 3. The local step

## 3.1 Transported modes, the pencil, and the Kato–Temple bound
Let $\Phi_p : [-1,1]^2 \to Q_p$ be the bilinear map sending the corners to $v_1,\dots,v_4$, with Jacobian $J$ (affine in the reference
coordinates $(u,v)$, quadratic in $p$, $J = \tfrac14$ at $p=0$). For the square eigenfunctions
$g_{mn}(u,v) = 2\sin\frac{m\pi(u+1)}2\sin\frac{n\pi(v+1)}2$ the transported functions $w = g\circ\Phi_p^{-1}$ lie in $H^2\cap H^1_0(Q_p)$
and
$$K_{kl} = \int_{Q_p}\nabla w_k\cdot\nabla w_l = \int_{[-1,1]^2} \nabla g_k^{\!\top}\frac{B}{J}\nabla g_l,\qquad
M_{kl} = \int_{[-1,1]^2} g_kg_l\,J,\qquad LL_{kl} = \int_{[-1,1]^2} (Lg_k)(Lg_l)\,J,$$
where $B = \mathrm{adj}(D\Phi_p)\,\mathrm{adj}(D\Phi_p)^{\!\top}$ and $Lg = \frac1J\partial_\alpha\!\big(\frac{B^{\alpha\beta}}J\partial_\beta g\big)$
is the transported Laplacian. With the mode set $\{(m,n): 1\le m,n\le N\}$, $N=10$ (100 modes), we form three test functions
$$\tilde w_j = \Big(g_j + \sum_{k}\sum_{i=1}^4 C_{jki}\,p_i\, g_k\Big)\circ\Phi_p^{-1},\qquad j\in\{(1,1),(2,1),(1,2)\},$$
where the fixed rational numbers $C_{jki}$ approximate the first-order perturbation coefficients (they may be any fixed numbers; the
choice only affects the size of the certified neighbourhood). They are dyadic rationals with denominator $2^{40}$, rounded from
$-(\partial_iK_{kj}-\lambda_j\partial_iM_{kj})(0)/(\lambda_k-\lambda_j)$, with $C_{jki} = 0$ for $k = j$ and for $j,k$ both in the cluster; all 1 180
values are published (`local/replay_v2/corrections.json`). Let $\tilde K,\tilde M$ be the $3\times3$ Rayleigh–Ritz matrices and
$\rho = \tilde K_{11}/\tilde M_{11}$, $\varepsilon^2 = \widetilde{LL}_{11}/\tilde M_{11}-\rho^2 = \|(\Delta+\rho)\tilde w_1\|^2/\|\tilde w_1\|^2$.
For $\beta \le \lambda_2(Q_p)$ with $\beta > \rho$ (Temple's inequality; Temple 1928, Kato 1949) $\lambda_1(Q_p) \ge L := \rho - \varepsilon^2/(\beta-\rho)$; and
$\lambda_2(Q_p) \le \theta_2$, the second eigenvalue of the pencil $(\tilde K,\tilde M)$.

**Lemma 3.1 (inertia criterion).** Let $A = \tilde K - \tfrac52 L\tilde M$ and let $S = A_{cc} - A_{c1}A_{1c}/A_{11}$ be the Schur complement onto the
cluster block $c = \{(2,1),(1,2)\}$. If $A_{11} < 0$ and $\lambda_{\min}(S) < 0$ (i.e. $\det S < 0$ or $\operatorname{tr} S < 0$), then
$\theta_2 < \tfrac52 L$, hence $\xi(Q_p) < 5/2$.
*Proof.* Sylvester's law: the number of pencil eigenvalues $\theta_k \le \sigma$ equals the number of non-positive eigenvalues of
$\tilde K - \sigma\tilde M$. By Haynsworth's inertia additivity (Haynsworth 1968), $\mathrm{In}(A) = \mathrm{In}(A_{11}) + \mathrm{In}(S)$, so $A$ has at least two
negative eigenvalues; thus $\theta_2 < \tfrac52 L \le \tfrac52\lambda_1$ and $\lambda_2 \le \theta_2$. $\square$

**Lemma 3.2 (exact structure at the square; sympy, exact arithmetic).** With $\lambda_1 = 2\pi^2$, $\lambda_2 = 5\pi^2$:
(L0) $K(0) = \mathrm{diag}(\lambda_k)$, $M(0) = I$, $LL(0) = \mathrm{diag}(\lambda_k^2)$;
(L1) $\partial_i\lambda_1(0) = 0$ and $\partial_i\operatorname{tr}\big[K_{cc} - 5\pi^2 M_{cc}\big](0) = 0$ for $i \in\{a,b,c,d\}$; the first-order
cluster matrices are $\mathrm{diag}(6\pi^2,-6\pi^2)$ for $a$, $\big(\begin{smallmatrix}0&256/9\\256/9&0\end{smallmatrix}\big)$ for $d$, and $0$ for $b$ and $c$;
(L2) $\partial_i\big[LL_{11}-\lambda_1^2M_{11}-2\lambda_1(K_{11}-\lambda_1M_{11})\big](0) = 0$.
Consequently, along every ray $p = te$, $A_{cc}$, $A_{1c}$ and $\varepsilon^2$ vanish at $t=0$ exactly, $\varepsilon^2 = O(t^2)$, $S = O(t)$,
$\det S = O(t^2)$ and $\operatorname{tr} S = O(t^2)$; the exact first-order term of $\operatorname{tr}S$ vanishes.

**Corollary 3.3 (the corrected pencil).** (L0)–(L2) also hold for the corrected matrices $\tilde K,\tilde M,\widetilde{LL}$ of §3.1. Hence the following
Taylor coefficients of the corrected quantities are exactly zero: degree 0 of the off-diagonal entries of $\tilde K$ and $\tilde M$; degree 1 of $\rho$;
degrees 0 and 1 of $\varepsilon^2$; degree 0 of $A_{ij}$ for $(i,j)\ne(1,1)$; and degree 1 of $\operatorname{tr}S$.
*Proof.* At $p=0$ the functions $g_k$ are the orthonormal Dirichlet eigenfunctions of the square, so Green's identity gives
$K_{kj}(0) = \lambda_jM_{kj}(0) = \lambda_j\delta_{kj}$ and $LL_{kj}(0) = \lambda_j^2\delta_{kj}$ for all modes $k,j$. The corrections vanish at $p=0$, so
they do not change degree-0 coefficients, which gives the degree-0 statements (with $L(0)=\rho(0)=\lambda_1 = 2\pi^2$ and
$\lambda_2 - \tfrac52\cdot2\pi^2 = 0$). At degree 1 a correction $C_{jki}p_i$ enters $\tilde K_{jl}-\lambda\tilde M_{jl}$ only through
$C_{jki}\big(K_{kl}(0)-\lambda M_{kl}(0)\big) = C_{jki}(\lambda_l-\lambda)\delta_{kl}$, and $\widetilde{LL}_{11}$ similarly through
$LL_{k1}(0) = \lambda_1^2\delta_{k1}$. In $\rho$ ($\lambda = \lambda_1$, $l = 1$), in $\varepsilon^2$ (L2 form, $l=1$) and in $\operatorname{tr}[K_{cc}-5\pi^2M_{cc}]$
($\lambda = 5\pi^2 = \lambda_l$ for $l$ in the cluster) these contributions are multiplied by $\lambda_l - \lambda = 0$ or have $k = j$, where $C_{jki}=0$.
The products $A_{c1}A_{1c}$ in $S$ are $O(|p|^2)$. $\square$

## 3.2 Taylor models and the certificate
All entries $K_{kl}, M_{kl}, LL_{kl}$ are integrals over $[-1,1]^2$ of products of a rational function of $(u,v,p)$ (from $B$, $J$ and
$J^{-1},J^{-2},J^{-3}$) with trigonometric polynomials. Writing $4J = 1 + j_1$ with $|j_1| \le \hat g\,|p|_2$ on $[-1,1]^2\times\{|p|_2\le\rho_0\}$
(an explicit rational $\hat g = s_1 + s_2\rho_0$ from the exact coefficients of $j_1$), we expand $J^{-k}$ as a truncated binomial series with a
geometric tail, keep exact polynomial parts (rational coefficients, degree $\le m$ in $p$, exact in $(u,v)$) and bound all dropped terms by
$r\,|p|_2^{m+1}$ using $|u^av^b|\le1$ and $|p^\gamma|\le\sigma(\gamma)|p|_2^{|\gamma|}$ with $\sigma(\gamma) = \prod_i(\gamma_i/|\gamma|)^{\gamma_i/2}$, which is
submultiplicative. Integrals against the trigonometric factors reduce to exact one-dimensional tables evaluated in ball arithmetic (Arb).
The Taylor coefficients were validated against direct high-order quadrature (agreement $10^{-8}$ at $|p| = 0.06$).

**Schur complement as a Taylor model.** The $3\times3$ pencil, $\rho$, $\varepsilon^2$, the Kato–Temple bound $L$, $A = \tilde K - \tfrac52 L\tilde M$,
its Schur complement $S$, $\det S$ and $\operatorname{tr}S$ are all computed as multivariate Taylor models in $p$ (order $m$, remainder
$r|p|_2^{m+1}$). The coefficients listed in Corollary 3.3 are deleted. The justification is the exact identity of the corollary: a Taylor model with
ball coefficients asserts that the exact Taylor coefficients lie in the balls, so an exactly vanishing coefficient may be replaced by $0$.
Before each deletion the code also checks that the computed enclosure contains $0$; this is a consistency check only and is never used
as a justification. Only then are directions substituted. This ordering is essential: substituting
interval directions into the pencil first and forming the Schur complement afterwards inflates the enclosures roughly a hundredfold.
The shift $\beta$ is taken from a finite list of levels $36\le\beta_j\le49$; for each direction box and $s$-interval we use the largest level
not exceeding the circumscribed-rectangle lower bound for $\lambda_2$ (or $36 < 2\pi j_{0,1}^2 \le 2\pi j_{0,1}^2/|Q_p|$, the Krahn–Szegő
inequality: at fixed area $\lambda_2$ is smallest for two equal disks; Krahn 1926, Hong 1954, Pólya 1955).

**The cover.** Write $p = sE$ with $E$ on a face of the cube $[-1,1]^4$ and $s = |p|_\infty$. For each box of directions $E$ on a face we certify
$s\in(0, 1/20]$; since $|E|_2\le 2$, $|p|_2 \le 2s \le \rho_0 = 1/10$ (checked in exact rational arithmetic), so the Taylor models apply. On each
box and each of ten $s$-subintervals of length $1/200$ we verify $\det S/s^2 < 0$ or $\operatorname{tr}S/s^2 < 0$ by a second-order centred form in $E$ (exact value and gradient at the
box centre, interval Hessian over the box) plus the remainder $r\,s^{m-1}|E|^{m+1}$, together with $A_{11}<0$, $\tilde M\succ0$ (Gershgorin) and
$\beta>\rho$ (the last three, and $\tilde M\succ0$, are verified once on the whole ball $|p|_2\le\rho_0$ by Taylor-model range bounds, hence for
every box). When a Taylor model is contracted at an order lower than its own, the higher-degree terms are integrated exactly and folded
into the remainder. The shift level for each box and $s$-interval is chosen by comparing the levels with the exact rational lower endpoint
of the circumscribed-rectangle bound. For every box and $s$-interval the record `local/replay_v2/schur_cover2.json` stores the level, the
criterion used and the exact rational upper bound of $\det S/s^2$ or $\operatorname{tr}S/s^2$; all remainder radii are in
`local/replay_v2/remainders.json`.

**Theorem 2 (local step).** With $N = 10$ ($100$ transported modes), $m = 12$ (stiffness and mass), $m_{LL} = 14$ (Laplacian), $\rho_0 = 1/10$:
for every $p$ with $0 < |p|_\infty \le 1/20$ the criterion of Lemma 3.1 holds, hence $\xi(Q_p) < 5/2$. The certificate consists of $960$ direction
boxes on the eight faces of the cube, obtained from the $512$ initial boxes by at most four bisections ($1\,408$ box evaluations); of the
$9\,600$ pairs of a box and an $s$-subinterval, $8\,320$ are decided by $\det S$ and $1\,280$ by $\operatorname{tr}S$. The replay of 2 October 2026
(exact lemmas, $230$ s of exact assembly, Schur-complement Taylor models, the cover in $64$ s and the published records; $577$ s in total on one core,
peak memory $3.9$ GB) reproduced this with $0$ failures; its receipt `local/replay_v2/RECEIPT.json` has SHA-256 `9ebbaa46…`. The configuration
$N=8$, $m=10$, $\rho_0 = 2/25$ certifies the cube $|p|_\infty\le 1/25$ as a second, separate certificate. (Code: `local/engine2.py`, `assemble.py`, `localcert.py`, `mtm.py`, `schur_mtm.py`, `cover2.py`;
exact lemmas `lemma_firstorder.py`, `lemma_residual.py`; one-command replays `local/replay_local_v2.sh` (cube $1/20$) and
`local/replay_local.sh` (cube $1/25$); the first verifies the code against `local/MANIFEST.sha256` and ends with a receipt binding the hashes
of the code and of every output.)

# 4. The global step

By §2.4 it remains to treat the set $\mathcal R$, contained in the compact a-priori box of §2.4, of parameters $p$ in the fundamental domain $\{a,b,c\ge0\}$ with $Q_p$ convex,
$c_i(p) > 0.035\,c_{i+2}(p)$ for all $i$, $\mathrm{width}/\mathrm{diam} > 0.1$, and $|p|_\infty > 0.049$ (overlapping the cube $|p|_\infty\le 1/20$ of Theorem 2).

## 4.1 Test functions and the upper bound
On a box $B = p_0 + \prod_i[-r_i,r_i]$ we use three test functions $u_1,u_2,u_3\in H^1_0(Q_p)$ depending polynomially on $\delta = p-p_0$, and
the Rayleigh–Ritz value $\theta_2$ of the $3\times3$ pencil $K = (\nabla u_i,\nabla u_j)$, $M = (u_i,u_j)$ bounds $\lambda_2(Q_p)$ from above.
*Shell cover:* $u_k = b_p\,(P^0_k + \sum_i\delta_iP^i_k)$ with the edge bubble $b_p = \prod_{k=1}^4\ell_k$, $\ell_k(x) = (v_{k+1}-v_k)\times(x-v_k)$, fixed
polynomials $P^0_k$ (Rayleigh–Ritz eigenvectors at $p_0$ in a scaled Legendre basis of total degree $n=8$) and first-order corrections $P^i_k$
(finite differences of the aligned eigenvectors). *Transport cover:* $u_k(p) = \hat u_k(\cdot;\delta)\circ\Psi_p^{-1}$, where $\Psi_p$ is affine on the
two triangles of a diagonal split of $Q_{p_0}$ with linear parts $A_K(p) = M_K(p)M_K(p_0)^{-1}$ (certified $\det A_K>0$ on the box, so $\Psi_p$ is a
homeomorphism onto $Q_p$), and $\hat u_k$ are bubble polynomials on $Q_{p_0}$, optionally with first-order corrections. The diagonal is chosen so
that a vertex with angle near $\pi$ is an endpoint of it. Any such choice is admissible; it only affects tightness.

## 4.2 The lower bound: Lehmann–Goerisch with an exact flux
The Kato–Temple bound used in §3 fails for convex quadrilaterals with an angle near $\pi$ (for a $157^\circ$ angle the polynomial residual is
too large and the bound loses 40%). We use instead the Lehmann–Goerisch bound with $n=1$ (Lehmann 1949, 1950; Goerisch and Haunhorst
1985; Behnke and Goerisch 1994), in the flux form that Endo, Liu and Mariano (2026) also use for triangles.

For $u\in H^1_0(\Omega)$ and $w\in L^2(\Omega)^2$ write
$$A_0 = (\nabla u,\nabla u),\qquad A_1 = (u,u),\qquad A_2 = (w,w).$$
In the covers $u = u_1$ is the first test function, so $A_0 = K_{11}$ and $A_1 = M_{11}$.

**Lemma 4.1 (Lehmann–Goerisch, $n=1$).** Let $u\in H^1_0(\Omega)\setminus\{0\}$, $w\in L^2(\Omega)^2$ with $(\nabla f, w) = (f,u)$ for all
$f\in H^1_0(\Omega)$ (equivalently $w\in H(\mathrm{div})$, $\operatorname{div}w = -u$), and $0<\beta\le\lambda_2(\Omega)$ with $A_0<\beta A_1$. Then
$\beta A_2 - A_1 > 0$ and $\lambda_1(\Omega)\ge(\beta A_1 - A_0)/(\beta A_2 - A_1)$.

*Proof.* Let $(\varphi_k,\lambda_k)$ be the orthonormal Dirichlet eigenpairs and $c_k = (u,\varphi_k)$, so $A_1 = \sum c_k^2$ and $A_0 = \sum\lambda_kc_k^2$.
The field $w_* = \sum_k c_k\nabla\varphi_k/\lambda_k$ satisfies the same constraint, and $w - w_*$ is $L^2$-orthogonal to $\nabla H^1_0$, which contains
$w_*$; hence $A_2 = \|w\|^2\ge\|w_*\|^2 = \sum c_k^2/\lambda_k$. Therefore
$$\lambda_1(\beta A_2 - A_1) - (\beta A_1 - A_0) \ \ge\ \sum_k c_k^2\Big(\frac{\beta\lambda_1}{\lambda_k} - \lambda_1 - \beta + \lambda_k\Big)
= \sum_k c_k^2\,\frac{(\lambda_k-\lambda_1)(\lambda_k-\beta)}{\lambda_k}\ \ge\ 0,$$
since the $k=1$ term vanishes and $\lambda_k\ge\lambda_2\ge\beta$ for $k\ge2$. Finally $A_1 = (\nabla u, w)\le\sqrt{A_0A_2}$ gives
$\beta A_2 - A_1\ge A_1(\beta A_1 - A_0)/A_0 > 0$. $\square$

*The flux.* The two covers use two different constructions; both give $w\in H(\mathrm{div};Q_p)$ with $\operatorname{div}w = -u_1$ exactly for every
$p$ in the box.

(a) *Shell cover (global polynomial field).* Let $\mu_0\ne0$ be a fixed number (the floating-point Ritz value at $p_0$; any nonzero value is
admissible and no inequality is required of it). With $e = v_3(p)-v_1(p)$ the direction of a diagonal of $Q_p$, put $\tau = \varphi\,e$, where $\varphi$
is the polynomial with $e\cdot\nabla\varphi = -(\Delta u_1+\mu_0u_1)/\mu_0$ and $\varphi = 0$ on the line through $v_2(p)$ and $v_4(p)$ (computed exactly by
integration along lines parallel to $e$). Then $\operatorname{div}\tau = -u_1-\Delta u_1/\mu_0$, and
$w = \nabla u_1/\mu_0 + \tau + \operatorname{curl}\psi$ with a fixed polynomial $\psi$ (least-squares optimal at $p_0$, with first-order corrections) is a
global polynomial field, so no interface condition arises. Moreover $A_2 = (2\mu_0A_1 - A_0)/\mu_0^2 + \|\tau+\operatorname{curl}\psi\|^2$, which keeps the
large part $\nabla u_1/\mu_0$ out of the remainder chain.

(b) *Transport cover (piecewise Piola field).* Here $u_1(p) = \hat u_1\circ\Psi_p^{-1}$ is only piecewise polynomial and its normal derivative jumps on
the split diagonal, so $\nabla u_1/\mu_0$ is not in $H(\mathrm{div})$ and construction (a) is not used. On the reference quadrilateral $Q_{p_0}$, split
along a diagonal with direction $\hat e$ into triangles $K_1,K_2$, put $\hat\tau = \hat\varphi\,\hat e$ with $\hat e\cdot\nabla\hat\varphi = -\hat u_1$ and $\hat\varphi = 0$ on
the other diagonal, so that $\operatorname{div}\hat\tau = -\hat u_1$ (no Laplacian and no constant $\mu_0$), and let $\hat\psi$ be a single polynomial on
$Q_{p_0}$. On each triangle set $\hat w = d_K\,\hat\tau + \operatorname{curl}\hat\psi$ with $d_K = \det A_K > 0$ (certified on the box), and
$w = (A_K\hat w/d_K)\circ\Psi_p^{-1}$ (contravariant Piola map). Then $\operatorname{div}w = (\operatorname{div}\hat\tau)\circ\Psi_p^{-1} = -u_1(p)$ on each
triangle. Normal continuity across the image of the diagonal holds because (i) the two affine maps agree on the diagonal ($A_K\hat e$ is
the same for both triangles), (ii) $\hat\tau\cdot\hat n = 0$ there, so the normal component of $\hat w$ is the tangential derivative of $\hat\psi$, which
is single-valued because $\hat\psi$ is continuous, and (iii) the Piola map preserves normal flux, $(w\cdot n)\,ds = (\hat w\cdot\hat n)\,d\hat s$, so
different $d_K$ on the two sides do no harm. Equivalently $w = (\hat\varphi\circ\Psi_p^{-1})\,e(p) + \operatorname{curl}(\hat\psi\circ\Psi_p^{-1})$ with
$\hat\psi\circ\Psi_p^{-1}\in H^1(Q_p)$. (Specification in Appendix A, entries P3a–P3b.)

![The transported test functions and flux of §4.2(b) and §4.4: the box-centre quadrilateral $Q_{p_0}$ is split along a diagonal into
$K_1,K_2$ and mapped onto $Q_p$ by a map that is affine on each triangle; $\hat\tau$ is parallel to the diagonal, so the normal flux of the
Piola-transformed field is continuous across it.](figures/fig2_transport.pdf){width=82%}

The bound only needs the residual in a dual norm, which is why it survives near-straight angles (within 0.1% of $\lambda_1$ at $157^\circ$).

*The shift $\beta\le\lambda_2$.* Either a constant valid on the whole box, the maximum of: Krahn–Szegő $2\pi j_{0,1}^2/\max_B|Q_p|$; enclosing
rectangles in fixed directions, $\pi^2\min(4/A^2+1/B^2,\,1/A^2+4/B^2)$; and the Andrews–Clutterbuck gap $\lambda_2\ge\lambda_1+3\pi^2/D^2$
(Andrews and Clutterbuck 2011) with a lower bound for $\lambda_1$ from rectangles or Hersch's inequality $\lambda_1\ge\pi^2/(4r^2)$ ($r$ the
inradius; Hersch 1960), or a bootstrap that feeds a Lehmann–Goerisch value obtained with an
already valid shift back into the gap bound (both covers) — or (transport cover only) a Taylor model $\beta(p)\le\lambda_2(Q_p)$ valid
pointwise (Krahn–Szegő with the exact area $1-a^2-d^2$, rectangles with affine side bounds, the gap bound, and the bootstrap). The
Andrews–Clutterbuck gap and Hersch's inequality require convexity, which holds on every certified box. The shift $\beta$ and the flux
constant $\mu_0$ are unrelated. For elongated shapes, where no available $\beta$ exceeds $\lambda_1$, a *domain-monotonicity
mode* replaces $L$ by a scalar $L_0\le\lambda_1$ (enclosing rectangles); only 12 boxes of the two covers use it.

## 4.3 The criterion and its Taylor models
With $\tilde A := K(\beta A_2 - A_1) - \kappa(\beta A_1 - A_0)M = (\beta A_2 - A_1)(K - \kappa L M)$ (or $\tilde A = K - \kappa L_0M$ in the
domain-monotonicity mode), Lemma 3.1 applies verbatim: if $\tilde A_{11} < 0$ and $S' = -\tilde A_{11}\tilde A_{cc} + \tilde A_{c1}\tilde A_{1c}$ has a negative
eigenvalue (detected by $S'_{11}<0$, $S'_{22}<0$, $\operatorname{tr}S'<0$ or $\det S'<0$), then $\tilde A$ has two negative eigenvalues (Haynsworth),
and $\theta_2 < \kappa L\le\kappa\lambda_1$. The side conditions $\beta A_1 - A_0>0$ and $M\succ0$ are certified as well ($\beta A_2-A_1>0$ then follows).
The theorem uses $\kappa = 5/2$; every certificate is checked at $\kappa = 5/2$ first, with strict inequalities.

All quantities are exact integrals of polynomials: the triangles of the split are integrated with a Duffy–Gauss–Legendre product rule
with $m = n+6$ points, exact up to total degree $2n+10$, which covers every integrand; nodes and weights are Arb balls. In the shell cover the
base quantities are order-3 Taylor models in $t = \delta/r\in[-1,1]^4$ (vertices and nodes move linearly; fixed polynomials are expanded at the
moving nodes with a rigorous tail; every product folds the dropped coefficients into the remainder). In the transport cover all
$p$-dependence sits in $A_K(p)$ and $\delta$: $K$, $M$, $A_2$ are explicit rational functions of $p$ built from reference integrals on $Q_{p_0}$ that
are quadratic in $\delta$ and computed once per box in ball arithmetic; denominators are cleared with $\det A_1\det A_2>0$, and the resulting
polynomials are order-4 Taylor models with exact full products (coefficients of degree $\ge5$ folded into the remainder). Since the criterion
is homogeneous in $(K,M,A_2)$, the base models are first multiplied by a certified positive polynomial approximating $1/M_{11}$, which removes
the common dilation variation without changing any sign. Ranges over the box are bounded in ball arithmetic.

## 4.4 The cover
**Near-square shell (fixed test functions).** With fixed polynomial test functions (degree $n = 8$, first-order corrections in $p$),
order-3 Taylor models and the criterion above, the shell $0.049 < |p|_\infty \le 0.225$ of the fundamental domain is covered by
$8216$ certified boxes ($8209$ Lehmann–Goerisch, $7$ with a scalar domain-monotonicity bound $L_0\le\lambda_1$), with no failures,
in $1171$ s on $8$ threads; a separate replay script recomputes every certificate from the log and checks coverage (certified volume
equals the region's volume). The smallest certified margin is $5/2 - \xi_{\rm bound} = 2.36\times10^{-5}$.

**Theorem 3 (near-square region).** For every convex quadrilateral $Q_p$ with $0 < |p|_\infty\le 0.225$, $\xi(Q_p) < 5/2$.

**Remaining region (transported test functions).** Fixed test polynomials lose accuracy across a box when a vertex angle is close to
$\pi$ (the eigenfunction's corner behaviour moves with the vertex); the pilot predicted $\ge10^7$ boxes. We therefore transport the
box-centre test functions with the vertices, $u_k(p) = \hat u_k\circ\Psi_p^{-1}$ with $\Psi_p$ affine on the two triangles of a diagonal split
(the diagonal is chosen per box, as in §4.1; each of the two diagonals was used for about half of the certified boxes, 16 491 and 14 935
times), and the flux by construction (b) of §4.2, $w = A_K\hat w/\det A_K$ with $\hat w = \det A_K\,\hat\tau + \operatorname{curl}\hat\psi$, where
$\hat\tau = \hat\varphi\,\hat e$ and $\hat\varphi$ is minus the antiderivative of $\hat u_1$ along the diagonal direction (so $\operatorname{div}\hat\tau = -\hat u_1$ and
$\hat\tau\cdot\hat n = 0$ on the diagonal). Then $\operatorname{div} w = -u_1(p)$ exactly for every $p$, and
all $p$-dependence sits in the two $2\times2$ matrices $A_K(p)$, so every quantity is an explicit low-degree rational function of $p$ built from
reference integrals computed once per box.

The cover starts from the a-priori box of §2.4 cut at $a,b,c = 1/5$ and $d=\pm1/5$ (24 roots; the root $[0,\tfrac15]^3\times[-\tfrac15,\tfrac15]$
lies in the region of Theorem 3). A box is either certified, bisected, or discarded because one of the following holds at every point
of the box: E2 (Lemma A, $c_i\le0.035\,c_{i+2}$ for a fixed $i$), E4 ($c_i\le0$ for a fixed $i$: not convex), SHELL (inside the region of
Theorem 3), E3T (width/diam $\le 2/5$, decided by an exact sufficient test on one edge and one vertex pair), or a recursive combination
of these on sub-boxes (COMPOUND). The run produced 100 342 records forming 24 complete bisection trees with 50 183 leaves:

| leaves | count | share of the root volume |
|---|---|---|
| certified (Lehmann–Goerisch 31 421; domain monotonicity 5) | 31 426 | 0.45 % |
| E2 (Lemma A) | 13 933 | 2.50 % |
| E4 (not convex) | 138 | 4.04 % |
| E3T (width/diam $\le 2/5$; Lemma B$'$ of §4.5) | 1 568 | 0.36 % |
| COMPOUND (1 010 of them use E3T on some sub-box) | 3 117 | 92.62 % |
| SHELL | 1 | 0.01 % |

There were no failures and no errors. Every certificate holds at $\kappa=5/2$ with strict inequalities; on a ladder of trial values the
smallest certified $\kappa$ has median 2.31 and minimum 1.66, and about 1 % of the boxes certify only at $\kappa = 5/2$ itself.
The run took about 5 200 s on 8 threads of a shared machine (it was interrupted once by a system restart and resumed from its log; the
driver skips finished leaves and re-expands logged splits, so no box was certified twice or lost).

*Separately written checks.* (i) A separate exact-arithmetic program (`verify/xcover.py`) rebuilds the 24 roots and the bisection trees from the
log, checks that every leaf's box equals the reconstructed one (so the leaves tile the root box exactly: exact rational volume fractions
sum to 1), and re-decides every discard with its own exact polynomial bounds and recursive refinement: 0 problems. The same program
confirms the shell cover of Theorem 3 (8 216 certified leaves and the central cube; volume sum exactly 1). (ii) A clean replay of the frozen log, bound to its hash and to the hashes of the verifier sources, re-verified all 31 426 certificates and re-checked all 18 756 exclusion and compound leaves with 0 failures (§6). (iii) A second implementation of the transported certificate re-certifies a stratified random sample of certified boxes (§6).

## 4.5 Lemma B$'$: the thin band $1/10 < \mathrm{width}/\mathrm{diam} \le 2/5$

**Lemma B$'$.** Every convex quadrilateral $Q$ with $1/10 < s := \mathrm{width}(Q)/\mathrm{diam}(Q) \le 2/5$ satisfies $\xi(Q) < 5/2$.

Thin shapes defeat the two-dimensional certificates of §4.4 for a structural reason: $\lambda_1\approx\pi^2/w^2$ varies by $\pm 8r/w$ over a
parameter box of radius $r$, so any box-constant eigenvalue bound loses that factor twice. We therefore normalise by the diameter and
use one-dimensional bounds along it, which follow the $1/w^2$ scaling.

*Charts.* Normalise a diameter pair of vertices to $(0,0),(1,0)$ (all six vertex distances are then $\le 1$). If the pair is a diagonal
(chart I), the other vertices are $P_1 = (x_1,-t y_3)$ and $P_3 = (x_3, y_3)$ with $0<t\le1$, and by the reflections $y\mapsto-y$,
$x\mapsto1-x$ we may assume $x_3\le 1/2$; the triangle $P_0P_2P_3$ has longest side $1$, so $y_3\le s\le(1+t)y_3$, and the band forces
$(y_3,t,x_1,x_3)\in[1/20,2/5]\times[0,1]\times[0,1]\times[0,1/2]$. If the pair is a side (chart II), the other vertices are $P_2=(x_2,y_2)$,
$P_3 = (x_3,t y_2)$ with $0<t\le1$ after a reflection, $x_3<x_2$, and $s = y_2$ exactly, so $(y_2,t,x_2,x_3)\in[1/10,2/5]\times[0,1]^3$.
All cross products $c_i$ are multilinear in these coordinates, so their ranges over a box are attained at its corners.

*Lower bound for $\lambda_1$ (one-dimensional reduction).* If $h_Q(x)$ is the length of the slice $Q\cap\{x=\mathrm{const}\}$ and
$0\le V\le\pi^2/h_Q^2$, then $\lambda_1(Q)\ge\mu_1\big(-\tfrac{d^2}{dx^2}+V\big)$ with Dirichlet conditions on $[0,1]$: apply the one-dimensional
Poincaré inequality on each slice to $\varphi\in C_c^\infty(Q)$ and note that $v(x) = \|\varphi(x,\cdot)\|_{L^2}$ satisfies $|v'|\le\|\partial_x\varphi(x,\cdot)\|$.
For a parameter box we take a polyline $h_{\rm out}\ge h_Q$ valid for every $Q$ in the box and the cellwise-constant minorant
$V = \pi^2/\max_{\rm cell}h_{\rm out}^2$; $\mu_1 > E$ is certified by Sturm oscillation (the solution of $-v''+(V-E)v=0$, $v(0)=0$,
$v'(0)=1$ stays positive on $(0,1]$), using exact transfer matrices in ball arithmetic.

*Upper bound for $\lambda_2$.* For an inner domain $Q_{\rm in} = \{a(x)<y<b(x)\}$ contained in every $Q$ of the box ($a,b$ polylines) and
$u = f(x)\sin(\pi(y-a)/h)$, $h = b-a$, $f(0)=f(1)=0$:
$$\int|\nabla u|^2 = \int_0^1\Big[\tfrac{h}{2}f'^2 + \tfrac{h'}{2}ff' + \tfrac{\pi^2}{h}f^2\Big(\tfrac12+\tfrac{a'b'}{2}+h'^2\big(\tfrac16+\tfrac1{4\pi^2}\big)\Big)\Big]dx,\qquad
\int u^2 = \int_0^1\tfrac{h}{2}f^2\,dx .$$
With two piecewise-linear functions $f_1,f_2$ the $2\times2$ Gram matrices $N, D$ are enclosed exactly element by element (closed forms
with logarithms, or rigorous series when the element is nearly uniform).

*Certificate.* If $\mu_1>E$ is certified and $\tfrac52 E\,D - N\succ0$ (leading minors positive as balls), then
$\lambda_2(Q)\le\lambda_2(Q_{\rm in})<\tfrac52E<\tfrac52\lambda_1(Q)$ for every valid $Q$ in the box.

*Exclusions* (each decided exactly on the whole closed box): Lemma A; $s\le 1/10$; $s>2/5$ (using the exact formula
$w = \min_i\max(c_i,c_{i+3})/|v_{i+1}-v_i|$ for the width of a quadrilateral); and boxes containing no valid chart point.

*Repairs* after review are recorded in Appendix C; the numbers below are from the repaired code.

*The cover.* Chart I: 45 169 leaves (36 696 certified, 8 094 with $s>2/5$, 141 with $s\le1/10$, 110 Lemma A, 128 invalid); chart II:
8 048 leaves (7 572 certified, 161 Lemma A, 315 invalid); no box needed refinement below width $2^{-24}$; largest certified value of
$R_2/(\tfrac52E)$: $0.9999985$ and $0.99997$ (a box is certified as soon as the strict inequality holds). The replay of the re-certified logs recomputes every box from
its path, re-derives every exclusion and re-runs every certificate (583 s on three threads, 0 errors); it fails closed on unknown, duplicate, malformed or
pending records and binds each result to hashes of its log record and of the verifier sources; twenty deliberately corrupted logs were all
rejected. *Separately written checks:* a separately written exact program (`verify/xbandcover.py`) rebuilds both tilings from the paths (exact
volume), and re-derives all 8 949 exclusions with its own predicates (0 problems); a separately written implementation of both
one-dimensional bounds (`verify/xband.py`; validated against two-dimensional finite elements on 40 random band shapes, where it brackets
both eigenvalues with margins of 1–10 % and 0.07–6 %) re-certifies a random sample of 200 certified boxes per chart of the re-certified logs: all 400 certificates are reproduced (382 with its default mesh, the other 18, all with logged margins above $0.997$, after refining its own one-dimensional meshes). An adversarial review also found a defect in this separately written checker (its outer slice profile was too small for boxes touching $x=0$); it was repaired and the whole sample re-run before the numbers above were obtained. Every receipt records SHA-256 hashes of the checked logs and of the checker sources.

## 4.6 Convex quadrilaterals
**Theorem 1$'$.** For every convex quadrilateral $Q$, $\xi(Q)\le5/2$, with equality if and only if $Q$ is a square.

*Proof.* Let $Q$ be a convex quadrilateral that is not a square. By Endo–Osting (Theorem 2.1) $Q$ is similar to some $Q_p$ with all $c_i(p) > 0$,
and applying an element of $D_4$ we may take $p$ in $\{a,b,c\ge0\}$; $p\ne0$. If $c_i(p)\le0.035\,c_{i+2}(p)$ for some $i$, Lemma A gives
$\xi(Q)<5/2$. If $\mathrm{width}/\mathrm{diam}\le 1/10$, Lemma B gives $\xi(Q)\le 2.418$; if $1/10 < \mathrm{width}/\mathrm{diam}\le 2/5$, Lemma B$'$ gives
$\xi(Q)<5/2$. If $0<|p|_\infty\le 1/20$, Theorem 2 applies, and if $|p|_\infty\le 0.225$, Theorem 3. Otherwise $p$ lies in the a-priori box of
§2.4 (because width/diam $> 1/10$), hence in a leaf of the transported cover; that leaf is not a discard of type E2, E3T, E4 or SHELL
covering $p$ (those cases were treated above or are impossible for convex $Q_p$), so $p$ lies in a certified box and $\xi(Q_p) < 5/2$.
(For COMPOUND leaves the same reasoning applies to the sub-box containing $p$.) For the square $\xi = 5\pi^2/2\pi^2 = 5/2$. $\square$

# 5. Non-convex quadrilaterals (darts)

**Theorem 5.** For every simple non-convex quadrilateral $Q$, $\xi(Q) < 5/2$.

## 5.1 Setting
A dart has one reflex vertex $D$; call the opposite vertex $B$ and the others $A$, $C$ (boundary $A\to B\to C\to D\to A$). The diagonal $BD$
lies inside $Q$, $Q = \mathrm{int}(ABD\cup CBD\cup BD)$, its convex hull is the triangle $T = ABC$, and the notch is $ACD$. By a similarity,
$B=(0,0)$, $D=(1,0)$, $A = \rho_A(\cos b_A,\sin b_A)$, $C = \rho_C(\cos b_C,-\sin b_C)$ with $b_A,b_C\in(0,\pi)$ (the *BD chart*); the reflection in
$BD$ exchanges $A$ and $C$, so we may assume $b_A\le b_C$. The line $AC$ meets the ray $BD$ at abscissa
$f = \sin\theta_B/(\sin b_A/\rho_C + \sin b_C/\rho_A)$, $\theta_B = b_A+b_C$, and $Q$ is a dart with reflex vertex $D$ iff $\theta_B<\pi$ and $f>1$;
the notch depth $\beta = \mathrm{dist}(D,AC)/\mathrm{dist}(B,AC)$ equals $1-1/f$ (Figure 1).

![The BD chart of a dart (§5.1), drawn to scale for $b_A = 50^\circ$, $b_C = 60^\circ$, $\rho_A = 3$, $\rho_C = 2.6$ ($f\approx1.611$,
$\beta\approx0.379$).](figures/fig1_dart_chart.pdf){width=58%}

 Numerical evidence (finite elements, adversarial search; an
observation, not proved and not used here) suggests that the supremum of $\xi$ over darts is $7/3$, approached only along three degenerate curves: shallow notches over an equilateral hull ($\beta\to0$), and an
equilateral $CBD$ with a vanishing arm $ABD$ (a thin spike, or $A\to D$). The shape space is non-compact (arms of unbounded length, deep
notches $\beta\to1$), so the proof combines analytic lemmas for the degenerate and unbounded regions with certified covers of compact
families.

## 5.2 Analytic results (complete proofs in Supplement S)
* **Shallow notch** (Proposition 5.1): if $\beta\le0.0339$ then $\xi(Q)\le\tfrac73(1-\beta)^{-2}<5/2$ (homothety of the hull at $B$ and the
  triangle theorem).
* **Collapse collar** (Lemma 3.1): if $A$ is within $0.035$ times the inradius of $CBD$ from the segment $BD$ (or $A\leftrightarrow C$), then
  $\xi(Q)\le 2.4995$.
* **Thin wedge at $B$** (Lemmas 5.1, 5.3): if both $b_A$ and $b_C$ are acute and $\tan b_A + \tan b_C\le 1/10$, then $\xi\le2.4223$; if
  $\theta_B\le0.1772$ (which forces both angles to be acute), then $\xi\le2.4892$ (sector comparison, Bessel zeros certified in ball arithmetic).
  The acute hypothesis in Lemma 5.1 is essential: the tangent condition alone admits obtuse configurations with a negative tangent. The
  rigorous predicate for Lemma 5.1 checks $a_x>0$ and $c_x>0$ on the whole box, and no record of any final log uses Lemma 5.1; the
  decomposition of §5.4 uses Lemma 5.3. **Flat hull** (Lemma 5.2): hull angles at $A$ and $C$ at most
  $\arctan(1/10)$ give $\xi\le2.4356$. **Thin stem** (Lemma 12.1): if the vertical chords of $Q$ over $BD$ are at most $0.05\,|BD|$, then
  $\xi\le 2.4356$.
* **Long arms** (Lemmas 6.1, 9.1–9.5): cutting an arm by a segment and imposing a Neumann condition on the cut (Dirichlet–Neumann
  bracketing) bounds $\lambda_1$ below by the mixed eigenvalue of a bounded core and the chord bound $\pi^2/\ell^2$ of the cut-off tip, uniformly
  in the arm length; together with domain monotonicity for $\lambda_2$ this reduces arms of any length to compact families. A closed-form
  bound (Lemma 9.5) disposes of very long arms making an angle between $40^\circ$ and $140^\circ$ with $BD$.
* **One-dimensional bounds for dart families** (Lemmas 11.1–11.3): for any open set, $\lambda_1\ge\inf\sigma\big(-\tfrac{d^2}{dx^2} + \pi^2/\ell(x)^2\big)$,
  where $\ell(x)$ is the length of the *longest component* of the slice at $x$; slab-wise upper bounds of $\ell$ valid for whole families
  of darts (including arms of unbounded length, via half-strip containment and exponential tails), Sturm shooting in ball arithmetic
  (on unbounded intervals with a decaying-exponential start and an Allegretto–Piepenbrink argument; Allegretto 1974, Piepenbrink 1974), and a P2 Rayleigh–Ritz bound on an
  *inner dart* contained in every dart of the family give a certificate that is uniform in the arm lengths.

## 5.3 Certificates for boxes of darts
Three certificates are used, each valid for every dart of a parameter box. (a) *Finite elements* (Supplement S, §§8 and 10): a
Crouzeix–Raviart (1973) discretisation on an affine image of a fixed reference mesh with Liu's (2015) guaranteed lower bound
$\lambda_k\ge\lambda_{k,h}/(1+0.1893^2h^2\lambda_{k,h})$, a Temple step whose shift $\sigma\le\lambda_{2,h}$ is certified by one verified $LDL^T$
inertia count, transfer between parameter points by a strain bound, and a conforming P2 Rayleigh–Ritz bound for $\lambda_2$. (b) *BO2*:
the one-dimensional lower bound of §5.2 for $\lambda_1$ and a P2 Ritz bound for an inner dart. (c) *INC* (Lemma 11.4): an outer dart
containing every dart of the box and an inner dart contained in all of them, with a point finite-element lower bound for the outer
and a P2 upper bound for the inner one. All inequalities are decided in ball arithmetic; box geometry is exact.

## 5.4 The decomposition and the covers
Write $t=\tan b_A$. The dart condition in polar form (Lemma 13.1 of Supplement S) gives two facts used throughout: $1-\beta =
(\sin b_C/\rho_A + \sin b_A/\rho_C)/\sin\theta_B$, and every dart has $\rho_A>\sin b_C$ and $\rho_C>\sin b_A$.

**Theorem 5.4** (Theorem 13.2 of Supplement S). Every dart with $b_A\le b_C$ (so $b_A<\pi/2$) lies in at least one of the following sets,
and each set is covered by a result of §5.2 or by a certified run:

1. $\theta_B\le0.1772$: Lemma 5.3.
2. $\theta_B\ge174.29^\circ$: Lemma 5.2 (the hull angles at $A$ and $C$ sum to at most $5.71^\circ$).
3. $\beta\le0.0339$: Proposition 5.1.
4. every vertical chord of $Q$ over $BD$ has length at most $0.05$: Lemma 12.1.
5. *Collar*, $t<0.02$, every $\rho_A$: run SP_col for $\rho_C\in[e^{-4},e^3]$, together with the earlier finite-element run kf_collar for $C$ in
   the sector $38^\circ\le b_C\le100^\circ$, $0.5\le\rho_C\le1.8$; run SF_col for $\rho_C\ge20$; Lemma 12.1 for $\rho_C<e^{-4}$ (the stem chords
   are then below $0.02+e^{-4}<0.05$).
6. *Spike*, $0.02\le t\le3/10$ and $a_x\ge1$, every $\rho_A$: run SP_spk (with kf_spike in the same sector) for $\rho_C\in[e^{-4},e^3]$, and
   run SF_spk for $\rho_C\ge20$. For $\rho_C<e^{-4}$ there is no dart, since $\rho_C>\sin b_A>0.0199$.
7. All other darts. If $\rho_A<e^{-2.5}$ then $\sin b_C<\rho_A<0.083$, so $b_C<4.8^\circ$ (and $\theta_B\le2b_C$ is in set 1) or $b_C>175.2^\circ$
   (set 2). If $\rho_C<e^{-4}$ there is no dart. The polar box $\rho_A\in[e^{-2.5},e^3]$, $\rho_C\in[e^{-4},e^3]$ is run PP. Darts with an arm
   of length $\rho\ge20$ are runs PF, FP and FF; the polar and far grids overlap because $e^3>20$.

The certificates of the collar and spike runs hold for every arm length $|BA|$ at once. Lemmas 9.1–9.4 cut the arm $BA$ and put a Neumann condition on the
cut, and the BO2 bound treats the spike as a far arm. The far runs bound $\lambda_1$ from below with exponential tails (Lemma 11.3(ii)),
uniformly in $\rho\ge20$. Inside PP, a box with exactly one arm of $\rho\ge4$ is certified in the oblique-cut chart of Lemma 9.6 (FARL),
and two compact regions near the critical configurations are delegated to earlier finite-element runs in Cartesian BD coordinates.
These are `shallow_crit`, a box around shallow notches over an equilateral hull, and `main_cbox`, which takes $C$ in a box around the
equilateral apex with $t\ge0.3$ and $a_y\le1.5479$. Two further runs, `crit_collar` and `crit_spike_X1`, cover the spike charts for $C$ in the same box.
Every other box of a run is either certified or discarded by a re-decided predicate naming one of the sets 1–7, "not a dart"
($\theta_B\ge\pi$ on the box) or "mirror" ($b_C\le b_A$ on the box).

**Exactness of the grids.** A box is stored as floats $(p_0,r)$ and denotes the closed box $[p_0-r,p_0+r]$ in exact reals; every leaf,
skipped cell and chart box is checked to have exactly representable ends before any predicate or certificate uses it. Grid cells
whose nominal ends are not dyadic (the spike cells $0.02\cdot1.5^i$) are encoded outward. The processed roots of the earlier runs,
stored with nearest-rounded encodings, are represented by their outward dyadic snapping, which their leaves tile exactly. A separate
accounting program checks in rational arithmetic that the represented boxes of every family cover its parameter ranges. Decimal
thresholds enter every predicate as exact decimals compared with directed ball endpoints.

**Statistics.** The re-certified family logs (all leaves recomputed with the repaired verifier; 4 workers, 2.69 h) and the
earlier runs:

| run | region (Theorem 5.4) | leaves | discards | BO2 | INC | FEM |
|---|---|---:|---:|---:|---:|---:|
| PP | set 7, polar arms | 306 017 | 84 142 | 126 002 | 51 382 | 44 491 (27 861 FARL) |
| SP_col | set 5, $\rho_C\in[e^{-4},e^3]$ | 43 078 | 2 375 | 37 722 | — | 2 981 |
| SP_spk | set 6, $\rho_C\in[e^{-4},e^3]$ | 62 671 | 2 655 | 52 228 | — | 7 788 |
| PF | set 7, $\rho_C\ge20$ | 8 545 | 2 924 | 5 621 | — | — |
| FP | set 7, $\rho_A\ge20$ | 29 339 | 2 201 | 17 420 | — | 9 718 |
| FF | set 7, both arms $\ge20$ | 540 | 141 | 399 | — | — |
| SF_col / SF_spk | sets 5 / 6, $\rho_C\ge20$ | 393 / 996 | 23 / 91 | 370 / 905 | — | — |
| earlier runs (6 logs) | critical regions, KF sector | 11 270 | 1 042 | — | — | 10 228 |
| **total** | | **462 849** | | | | |

Every leaf of PP, SP_spk, PF, FP, FF and the far families re-certified at its original box. In SP_col, 232 finite-element leaves
whose $t$-range exceeds the collar chart's $t_1$ were subdivided and re-certified (+1 193 leaves); no leaf failed.

**Replay.** All dart logs, the earlier runs included, were replayed from scratch by the aggregate driver `replay_all.py` in proof mode.
The replay is fail-closed: unknown, duplicate, malformed or pending records are errors, the tiling of every root is checked
exactly, every discard is re-decided and every certificate is recomputed. Each result is bound to the SHA-256 of its log
record and to a binding hash. That hash covers every decision source, every cached chart, the manifest bytes, the family and the
runtime identity (interpreter, numerical libraries, floating-point environment). The driver runs the analytic gates (the certified
tables of Lemmas 5.3 and 9.5), the chart and accounting checks and the rigour lint before the replays. It validates every receipt record of every shard
against the bound log and writes an append-only, hash-chained receipt. Result: PASS. The replay reported here ran from a fresh extraction
of the frozen archive components, after a portability patch of the verifier (Appendix C). All 462 849 leaves of the 14 logs were
recomputed or re-decided, 2 shards per log, with no failure (5.40 h). The 33-record receipt has final digest
`bfe58c0b9ffe7759…` (receipt file SHA-256 `9055fcbce245ca81…`, pin `58bd96bf…`). It is verified by `verify/xreceipt.py`, an exact
state machine checked against the pin and the runtime anchors of the preceding structure-only run from the same extraction. The
mutation test of that checker rejects every forged variant of the receipt. An earlier replay of the same logs, before the portability
patch (pin `7e564f3f…`, 5.63 h), also passed; its receipt is kept in the archive.

**Separately written coverage check.** `verify/xdart.py`, written separately, samples darts and checks that each lies in a lemma region (re-decided
with exact or directed arithmetic) or in a box of a certified leaf of the logs (no enlargement). It uses 160 000 samples: uniform in each run's
parameters, concentrated near the boundaries of the lemma regions and of the grids (ten kinds of boundary), and 20 000 with
$\log\rho$ up to 40. On the pre-repair logs no sample was uncovered, and on the re-certified logs none was either (XDART PASS, 267 s; receipt with the SHA-256 of every log it read). This check is evidence of
correct bookkeeping, not part of the proof; the exact tiling and accounting are.

# 6. Verification, reproducibility and method statement

**Rigour conventions.** Floating point is used only to choose test functions, fluxes, correction coefficients and box splits; every inequality
that enters a proof is evaluated in ball arithmetic (Arb, Johansson 2017, in FLINT 3.6.0, through python-flint 0.9.0 and Arblib.jl 1.8.1 under
Julia 1.12.7) or with outward-rounded intervals, with exact rational box endpoints. Exact symbolic identities use SymPy 1.14.0 (Meurer et al.
2017); scikit-fem 12.0.2 (Gustafsson and McBain 2020) generates the connectivity of the dart reference meshes, whose node coordinates are
exact dyadic rationals. Every Taylor model carries an explicit remainder; truncating to a lower order folds the dropped terms into the remainder.

**Replays.** Every certificate log is replayed from scratch by a fail-closed script: unknown, duplicate, malformed or pending records are
errors; every box is recomputed from its log record (path or centre and radius); every exclusion is re-decided; every certificate is
recomputed; and every result is bound to SHA-256 hashes of its log record and of the verifier sources. Local step: `local/replay_local_v2.sh`,
all 960 direction boxes recomputed from the exact assembly, 0 failures (receipt `local/replay_v2/RECEIPT.json`).
Shell cover: exhaustive fail-closed replay of all 8 217 leaves (coverage exact; 8 216 certificates re-verified, 0 failures; 2 888 s on three threads). Transport cover: clean replay of all 50 182 non-shell leaves of the frozen log
(SHA-256 `cb187ebb…`): all 31 426 certificates re-verified and all 18 756 exclusion and compound leaves re-checked, 0 failures, smallest certified margin $5/2-\kappa^* = 5.05\times10^{-6}$ (7 481 s on three threads; receipt with result hash `c9fa7f8b…`). Lemma $B'$: 44 268 certificates and 8 949 exclusions, 0 errors. Darts: fresh fail-closed replay of all 14 logs, 462 849 leaves, PASS (§5.4).

**Fresh-extraction replay of the frozen archive.** After the component files of the archive (code, inputs, final logs, pins) were
frozen, they were extracted into a new directory and every component was replayed there, each run through the local job queue:
the local step (960 direction boxes, cube $1/20$; the exact correction coefficients, remainder radii and all 960 box records are reproduced byte for byte), Lemma $B'$ (44 268 certificates re-run, both charts tiled exactly), the shell (8 216 certificates, 8 217 leaves), the transport cover (31 426 certificates and all 18 756 exclusion and compound leaves) and the darts (all 462 849 leaves of the 14 logs, 2 shards per log); every replay passed with no failure. The release archive contains these component files byte for byte, which was checked against the manifest of the
replayed extraction. The receipts and logs of these runs are in the archive (`REPLAY_RECEIPT.md`).

**Separately written programs** (by the same producer, in Python with exact rationals and python-flint balls). (i) `verify/xcover.py` rebuilds the
bisection trees of the shell and transport covers from their logs, checks that every leaf box equals the reconstructed one (so the leaves
tile the root boxes exactly) and re-decides every discard with its own exact predicates and refinement: both covers pass (0 problems).
(ii) `verify/xcert.py` is a second implementation of the transported certificate (exact rational reference integrals, its own choice of
split diagonal, constant and pointwise Krahn–Szegő shifts, range bounding by Taylor shifts of exact polynomials, no first-order test-function
corrections). It treats a box by one certificate with adaptive range subdivision, then by subdivision of the box, with its own exact
exclusions (Lemma A, non-convexity) for sub-boxes. On a stratified random sample of 60 certified transport leaves it re-certifies all 60:
51 within a work budget of 40 units per box, and the other 9 with the current version of the checker and a budget of 200 units (at most
103 used; see Appendix C). (iii) `verify/xbandcover.py` and `verify/xband.py` re-check Lemma $B'$ (§4.5): exact tiling and all exclusions, and the certificates
of a random sample of 200 boxes per chart: all 400 are reproduced (382 with its default mesh, 18 after refining its own mesh; re-run
after a repair of the checker, Appendix C). (iv) `verify/xdart.py` checks that the dart decomposition is covered: 160 000 sampled
darts, including samples concentrated at every boundary of the lemma regions and grids and arms up to $e^{40}$, each lie in a lemma region or a
certified box of the final logs (0 uncovered). The dart certificates themselves have no second, independently written implementation; their
evidence is the fail-closed replay of §5.4. (v) `verify/xlocal.py` checks the published records of the local step in exact rational
arithmetic: the 960 direction boxes tile the eight faces exactly, the ten $s$-intervals of every box cover $(0,1/20]$, every shift level is
admissible by the circumscribed-rectangle or Krahn–Szegő bound, re-derived from the vertex formulas of §2.1, and every logged bound is
negative (PASS). It does not recompute the Taylor-model enclosures. It rejects each of five corrupted copies of the record (a flipped sign,
a missing or a duplicated box, a raised level, a missing $s$-interval). All receipts record the SHA-256 hashes of the inputs and of the
checker sources.

**Assurance by component.** The four kinds of evidence are kept separate: a written argument (reviewed adversarially, Appendix C); exact
coverage of the region; an exhaustive replay with the producing verifier; and a check by separately written code.

| component | written argument | exact coverage | exhaustive replay | separately written check |
|---|---|---|---|---|
| Lemmas A, B (§2) | complete proofs | — | — | exact arithmetic of the constants |
| Theorem 2, local step (§3) | Lemmas 3.1–3.2, Corollary 3.3 | 960 direction boxes, exact rational tiling of the faces, also checked by `xlocal.py` | all 960 boxes recomputed | the exact first-order structure (SymPy; also derived by the referee); tiling, shift levels and signs of the published records (`xlocal.py`); enclosures not reimplemented |
| Theorem 3, shell (§4.4) | §§4.1–4.3 | exact, also checked by `xcover.py` | all 8 216 certificates | coverage and exclusions (`xcover.py`); certificates not reimplemented |
| transport cover (§4.4) | §§4.1–4.3 | exact, also checked by `xcover.py` | all 31 426 certificates | coverage and exclusions (`xcover.py`); certificates of a sample (`xcert.py`) |
| Lemma $B'$, band (§4.5) | §4.5 | exact, also checked by `xbandcover.py` | all 44 268 certificates | coverage and exclusions (`xbandcover.py`); certificates of a sample (`xband.py`) |
| Theorem 5, darts (§5) | Supplement S | exact tiling and accounting in the replay | all 462 849 leaves | sampled coverage only (`xdart.py`) |

No component has been formally verified in a proof assistant, independently reproduced by an unaffiliated party, or refereed by a journal.

**Reviews.** Every part was reviewed adversarially by a second AI model prompted to refute it, and the complete draft was assessed in an
external referee report (2 October 2026). The findings, the repairs and the responses are recorded in Appendix C and in
`REFEREE_RESPONSE.md` of the archive.

**Method statement.** The problem selection, proofs, code and computations were produced by AI systems (Anthropic Claude, with an
OpenAI model used as an adversarial reviewer) operating autonomously; no human contributed mathematical content. The author line is
"Anonymous" accordingly.

# References
- Allegretto, W. (1974). On the equivalence of two types of oscillation for elliptic operators. *Pacific J. Math.* 55(2), 319–328. doi:10.2140/pjm.1974.55.319.
- Andrews, B., Clutterbuck, J. (2011). Proof of the fundamental gap conjecture. *J. Amer. Math. Soc.* 24(3), 899–916. doi:10.1090/S0894-0347-2011-00699-1.
- Antunes, P. R. S., Freitas, P. (2008). A numerical study of the spectral gap. *J. Phys. A: Math. Theor.* 41(5), 055201. doi:10.1088/1751-8113/41/5/055201.
- Arbon, R. (2022). Global and local bounds on the fundamental ratio of triangles and quadrilaterals. Preprint, arXiv:2207.05814.
- Arbon, R., Mannan, M., Psenka, M., Ragavan, S. (2022). A proof of the triangular Ashbaugh–Benguria–Payne–Pólya–Weinberger inequality. *J. Spectr. Theory* 12(2), 515–533. doi:10.4171/JST/409.
- Ashbaugh, M. S., Benguria, R. D. (1991). Proof of the Payne–Pólya–Weinberger conjecture. *Bull. Amer. Math. Soc. (N.S.)* 25(1), 19–29. doi:10.1090/S0273-0979-1991-16016-7.
- Ashbaugh, M. S., Benguria, R. D. (1992). A sharp bound for the ratio of the first two eigenvalues of Dirichlet Laplacians and extensions. *Ann. of Math. (2)* 135(3), 601–628. doi:10.2307/2946578.
- Behnke, H., Goerisch, F. (1994). Inclusions for eigenvalues of selfadjoint problems. In J. Herzberger (ed.), *Topics in Validated Computations*, Stud. Comput. Math. 5, North-Holland, Amsterdam, 277–322.
- Crouzeix, M., Raviart, P.-A. (1973). Conforming and nonconforming finite element methods for solving the stationary Stokes equations I. *Rev. Française Automat. Informat. Recherche Opérationnelle Sér. Rouge* 7(R-3), 33–75. doi:10.1051/m2an/197307R300331.
- Endo, R., Liu, X., Mariano, P. (2026). Sharp Dirichlet eigenvalue inequalities on triangles. Preprint, arXiv:2605.04331.
- Endo, R., Osting, B. (2026). Maximizing the fundamental Laplace–Neumann eigenvalue on quadrilaterals. Preprint, arXiv:2609.02784.
- Goerisch, F., Haunhorst, H. (1985). Eigenwertschranken für Eigenwertaufgaben mit partiellen Differentialgleichungen. *Z. Angew. Math. Mech.* 65(3), 129–135. doi:10.1002/zamm.19850650302.
- Gustafsson, T., McBain, G. D. (2020). scikit-fem: a Python package for finite element assembly. *J. Open Source Softw.* 5(52), 2369. doi:10.21105/joss.02369.
- Haynsworth, E. V. (1968). Determination of the inertia of a partitioned Hermitian matrix. *Linear Algebra Appl.* 1(1), 73–81. doi:10.1016/0024-3795(68)90050-5.
- Henrot, A. (2006). *Extremum Problems for Eigenvalues of Elliptic Operators.* Birkhäuser, Basel. doi:10.1007/3-7643-7706-2.
- Henrot, A., Pierre, M. (2018). *Shape Variation and Optimization: A Geometrical Analysis.* EMS Tracts in Mathematics 28, European Mathematical Society, Zürich. doi:10.4171/178.
- Hersch, J. (1960). Sur la fréquence fondamentale d'une membrane vibrante: évaluations par défaut et principe de maximum. *Z. Angew. Math. Phys.* 11(5), 387–413. doi:10.1007/BF01604498.
- Hong, I. (1954). On an inequality concerning the eigenvalue problem of membrane. *Kōdai Math. Sem. Rep.* 6(4), 113–114. doi:10.2996/kmj/1138843535.
- Johansson, F. (2017). Arb: efficient arbitrary-precision midpoint-radius interval arithmetic. *IEEE Trans. Comput.* 66(8), 1281–1292. doi:10.1109/TC.2017.2690633.
- Kato, T. (1949). On the upper and lower bounds of eigenvalues. *J. Phys. Soc. Japan* 4(4–6), 334–339. doi:10.1143/JPSJ.4.334.
- Krahn, E. (1926). Über Minimaleigenschaften der Kugel in drei und mehr Dimensionen. *Acta Comm. Univ. Tartu. (Dorpat.) A* 9, 1–44.
- Laugesen, R. S., Siudeja, B. A. (2017). Triangles and other special domains. In A. Henrot (ed.), *Shape Optimization and Spectral Theory*, De Gruyter Open, Warsaw, 149–200. doi:10.1515/9783110550887-006.
- Lehmann, N. J. (1949). Beiträge zur numerischen Lösung linearer Eigenwertprobleme. I. *Z. Angew. Math. Mech.* 29(11–12), 341–356. doi:10.1002/zamm.19502911005.
- Lehmann, N. J. (1950). Beiträge zur numerischen Lösung linearer Eigenwertprobleme. II. *Z. Angew. Math. Mech.* 30(1–2), 1–16. doi:10.1002/zamm.19500300101.
- Liu, X. (2015). A framework of verified eigenvalue bounds for self-adjoint differential operators. *Appl. Math. Comput.* 267, 341–355. doi:10.1016/j.amc.2015.03.048.
- Meurer, A., Smith, C. P., Paprocki, M., et al. (2017). SymPy: symbolic computing in Python. *PeerJ Comput. Sci.* 3, e103. doi:10.7717/peerj-cs.103.
- Payne, L. E., Pólya, G., Weinberger, H. F. (1956). On the ratio of consecutive eigenvalues. *J. Math. and Phys.* 35, 289–298. doi:10.1002/sapm1956351289.
- Pólya, G. (1955). On the characteristic frequencies of a symmetric membrane. *Math. Z.* 63, 331–337. doi:10.1007/BF01187944.
- Siudeja, B. (2010). Isoperimetric inequalities for eigenvalues of triangles. *Indiana Univ. Math. J.* 59(3), 1097–1120. doi:10.1512/iumj.2010.59.3744. Preprint arXiv:0707.3631 (2007).
- Temple, G. (1928). The theory of Rayleigh's principle as applied to continuous systems. *Proc. R. Soc. Lond. A* 119(782), 276–293. doi:10.1098/rspa.1928.0098.

*Software.* FLINT 3.6.0 (https://flintlib.org), python-flint 0.9.0 (https://github.com/flintlib/python-flint), Arblib.jl 1.8.1
(doi:10.5281/zenodo.22073429), Julia 1.12.7, SymPy 1.14.0, scikit-fem 12.0.2, NumPy 2.4.4 and SciPy 1.17.1 (floating point only, for choices
that the certificates then verify).

# Appendix A. Soundness specification of the certificate primitives

Each entry states what a primitive takes as input, the hypotheses the verifier establishes (in ball arithmetic or exact rational
arithmetic, on the whole closed parameter box unless stated otherwise), and the mathematical conclusion. File names refer to the
research archive; the replays of §6 recompute every primitive from the logged inputs.

**P1. Rayleigh–Ritz with the inertia criterion** (Lemma 3.1, §4.3; `local/schur_mtm.py`, `global2/criterion.jl`). *Input:* the Gram
matrices $K = ((\nabla u_i,\nabla u_j))$ and $M = ((u_i,u_j))$ of three functions $u_1,u_2,u_3\in H^1_0(Q_p)$; a number $L$ with $0<L\le\lambda_1(Q_p)$;
$\kappa = 5/2$ (exact); a factor $c>0$; $\tilde A = c\,(K-\kappa L M)$. *Hypotheses:* $M\succ0$; $\tilde A_{11}<0$; the matrix
$S' = -\tilde A_{11}\tilde A_{cc}+\tilde A_{c1}\tilde A_{1c}$ has a negative eigenvalue, established by a negative diagonal entry, a negative trace or a
negative determinant. *Conclusion:* $S' = (-\tilde A_{11})$ times the Schur complement of $\tilde A_{11}$, so by Haynsworth's inertia additivity $\tilde A$ has at
least two negative eigenvalues; $K-\kappa LM$ is negative definite on a two-dimensional subspace, and by the min–max principle
$\lambda_2(Q_p)<\kappa L\le\kappa\lambda_1(Q_p)$.

**P2. Temple's inequality** (§3.1; `local/schur_mtm.py`). *Input:* $\tilde w_1\in H^2\cap H^1_0(Q_p)$; $\rho = \tilde K_{11}/\tilde M_{11}$;
$\varepsilon^2 = \|(\Delta+\rho)\tilde w_1\|^2/\|\tilde w_1\|^2$; a shift $\beta$. *Hypotheses:* $\beta\le\lambda_2(Q_p)$ (P5) and $\beta>\rho$.
*Conclusion:* $\lambda_1(Q_p)\ge L := \rho-\varepsilon^2/(\beta-\rho)$.

**P3a. Shell flux** (§4.2(a); `global2/assemble.jl`). *Input:* box $B$; polynomial degree $n$; fixed polynomials (test function and its
first-order corrections, $\psi$ and its corrections); a constant $\mu_0\ne0$. *Hypotheses:* $Q_p$ is convex for $p\in B$; $u_1 = b_p\,P_1(\cdot;p)$ with the
edge bubble $b_p$; $\varphi$ is the polynomial with $e(p)\cdot\nabla\varphi = -(\Delta u_1+\mu_0u_1)/\mu_0$ and $\varphi = 0$ on the line through
$v_2(p),v_4(p)$. *Conclusion:* $u_1\in H^1_0\cap H^2(Q_p)$; $w = \nabla u_1/\mu_0+\varphi\,e(p)+\operatorname{curl}\psi$ is a polynomial vector field with
$\operatorname{div}w = -u_1$; $A_0 = K_{11}$, $A_1 = M_{11}$ and $A_2 = (2\mu_0A_1-A_0)/\mu_0^2+\|\varphi e+\operatorname{curl}\psi\|^2$ are enclosed by the
Taylor models on $B$ (the quadrature is exact for every integrand, up to total degree $2n+10$).

**P3b. Transported flux** (§4.2(b); `global2/transport1.jl`). *Input:* $Q_{p_0}$ split along a diagonal of direction $\hat e$ into $K_1,K_2$;
$A_K(p) = M_K(p)M_K(p_0)^{-1}$, $d_K = \det A_K$; the reference test function $\hat u_1(\cdot;\delta)$; one polynomial $\hat\psi$ on $Q_{p_0}$.
*Hypotheses:* $d_K>0$ on $B$ for both triangles; $\hat\tau = \hat\varphi\,\hat e$ with $\hat e\cdot\nabla\hat\varphi = -\hat u_1$ (so $\operatorname{div}\hat\tau = -\hat u_1$
and $\hat\tau\cdot\hat n = 0$ on the diagonal); the two affine maps agree on the diagonal. *Conclusion:* $\Psi_p$ is a piecewise-affine
homeomorphism of $Q_{p_0}$ onto $Q_p$ and $u_1(p) = \hat u_1\circ\Psi_p^{-1}\in H^1_0(Q_p)$; on each image triangle
$w = (A_K(d_K\hat\tau+\operatorname{curl}\hat\psi)/d_K)\circ\Psi_p^{-1}$ has normal trace equal to the tangential derivative of $\hat\psi\circ\Psi_p^{-1}$ from both
sides of the diagonal, hence $w\in H(\mathrm{div};Q_p)$ and $\operatorname{div}w = -u_1(p)$.

**P4. Lehmann–Goerisch** (Lemma 4.1; `global2/criterion.jl`). *Input:* $u = u_1$, $w$ from P3a or P3b, a shift $\beta$. *Hypotheses:*
$\beta\le\lambda_2(Q_p)$ (P5); $\beta A_1-A_0>0$ (certified); $\beta A_2-A_1>0$ (certified, and implied). *Conclusion:*
$\lambda_1(Q_p)\ge L = (\beta A_1-A_0)/(\beta A_2-A_1)$, used in P1 with $c = \beta A_2-A_1$. Where no admissible $\beta$ exists, the
domain-monotonicity mode uses $L = L_0$ from enclosing rectangles (12 boxes in total).

**P5. Shift $\beta\le\lambda_2$** (`local/cover2.py`, `global2/criterion.jl`, `global2/bounds_tm.jl`). *Input:* rigorous enclosures of the area, of the
sides of enclosing rectangles, of the diameter and of the inradius over the box; $j_{0,1}$ as a ball. *Hypotheses:* $Q_p$ is convex (needed by the
Andrews–Clutterbuck and Hersch bounds). *Conclusion:* each logged $\beta$ satisfies $\beta\le\lambda_2(Q_p)$, for all $p$ in the box when constant
or pointwise when a Taylor model, by one of: Krahn–Szegő $\lambda_2\ge2\pi j_{0,1}^2/|Q_p|$; domain monotonicity in an enclosing rectangle;
Andrews–Clutterbuck $\lambda_2\ge\lambda_1+3\pi^2/D^2$ combined with a lower bound for $\lambda_1$ (rectangles, Hersch, or a Lehmann–Goerisch
value obtained with an already valid shift). In the local step the levels are compared with the exact rational lower endpoint of the
circumscribed-rectangle bound, with the Krahn–Szegő floor $36<2\pi j_{0,1}^2\le2\pi j_{0,1}^2/|Q_p|$.

**P6. Taylor models** (`local/engine2.py`, `local/mtm.py`, `global2/bounds_tm.jl`, `darts/tm.py`). *Representation:* a polynomial in the box
variables with ball (or exact rational) coefficients and a remainder $r$, asserting that the function differs from a polynomial with
coefficients in the balls by at most $r\,|p|^{m+1}$ (local step) or by an element of $[-r,r]$ (covers). *Operations:* sums and products with
every coefficient of degree above the order integrated and folded into the remainder; reciprocals and inverse powers by truncated series
with a geometric tail; range bounds by centred forms. *Deletions:* a coefficient is deleted only when an exact identity proves that the
corresponding Taylor coefficient vanishes (Corollary 3.3); the accompanying check that its enclosure contains $0$ is a consistency
check only. *Conclusion:* every derived quantity is enclosed on the box.

**P7. Exclusion predicates** (`global2/cover_t.jl`, `verify/xcover.py`, `bandB/cover.jl`, `verify/xbandcover.py`, `darts/fam5.py`, `darts/regions.py`).
*Input:* a closed box with exact rational (convex covers, band) or exactly representable floating-point (darts) end points. *Hypotheses:*
the named condition holds at every point of the box, decided by exact polynomial bounds (multilinear cross products are evaluated at
the corners) or by directed ball arithmetic with decimal thresholds entered as exact decimals. *Conclusion:* the box is covered by the
named lemma or lies outside the region (not convex, not a dart, mirror image, outside a chart).

**P8. One-dimensional lower bound and Sturm certificate** (§4.5; `bandB/cert.jl`, `verify/xband.py`). *Input:* an outer profile $h_{\rm out}$,
a polyline with $h_{\rm out}\ge h_Q$ for every $Q$ of the box, where $h_Q$ is the slice length; nodes $0=\xi_0<\dots<\xi_K=1$ containing every
breakpoint of $h_{\rm out}$; $V\equiv V_k\le\pi^2/\max_{[\xi_k,\xi_{k+1}]}h_{\rm out}^2$ on cell $J_k$; a number $E$. *Hypotheses:* with $v$ the $C^1$ solution of
$-v''+(V-E)v = 0$, $v(0)=0$, $v'(0)=1$, propagated cell by cell with the closed-form transfer matrices
$\left(\begin{smallmatrix}\cosh\kappa\ell&\kappa^{-1}\sinh\kappa\ell\\ \kappa\sinh\kappa\ell&\cosh\kappa\ell\end{smallmatrix}\right)$, $\kappa = \sqrt{V_k-E}$, or
$\left(\begin{smallmatrix}\cos\omega\ell&\omega^{-1}\sin\omega\ell\\ -\omega\sin\omega\ell&\cos\omega\ell\end{smallmatrix}\right)$, $\omega = \sqrt{E-V_k}$, in outward-rounded ball
arithmetic: for every $k$, (a) the enclosure of $V_k-E$ does not contain $0$ (in `verify/xband.py` alternatively $|(V_k-E)\ell_k^2|<1$, with the
transfer entries from their power series and a rigorous remainder); (b) if $V_k<E$ then $\omega_k\ell_k<\pi$ is proved; (c) the enclosure of
$v(\xi_{k+1})$ lies in $(0,\infty)$. *Conclusion:* $v>0$ on $(0,1]$: on each cell a non-trivial solution is $\alpha e^{\kappa s}+\gamma e^{-\kappa s}$ (at most
one zero) or $R\cos(\omega s-\phi)$ (zeros $\pi/\omega>\ell_k$ apart), every zero is simple, so positivity at both ends of a cell excludes an
interior zero; on the first cell $v = \sinh(\kappa s)/\kappa$ or $\sin(\omega s)/\omega$. By the Wronskian argument
$\mu_1(-d^2/dx^2+V)>E$, and by the slice Poincaré inequality $\lambda_1(Q)\ge\mu_1>E$ for every $Q$ of the box. Positivity inside cells is a
consequence of (a)–(c), not a separate numerical test.

**P9. One-dimensional Ritz upper bound** (§4.5; `bandB/cert.jl`, `verify/xband.py`). *Input:* polylines $a<b$ defining
$Q_{\rm in} = \{a(x)<y<b(x)\}$; two piecewise-linear $f_1,f_2$ with $f(0)=f(1)=0$; $E$ from P8. *Hypotheses:* $Q_{\rm in}\subset Q$ for every $Q$ of the
box (checked with the inner profile of the box); the $2\times2$ Gram matrices $N,D$ of the energy and mass forms of §4.5 are enclosed element by
element (closed forms with logarithms, or series with rigorous remainders); $\tfrac52E\,D-N\succ0$ (leading minors positive). *Conclusion:*
$\lambda_2(Q)\le\lambda_2(Q_{\rm in})<\tfrac52E<\tfrac52\lambda_1(Q)$.

**P10. Finite-element certificate for darts** (Supplement S §§8, 10; `darts/gcert.py`, `darts/gsigma.py`). *Input:* a chart, that is triangular
pieces that are affine images of exact dyadic reference meshes with vertices polynomial in the box parameters; a reference level $k$;
fixed test vectors. *Hypotheses:* the triangulation is conforming and resolves the Dirichlet and Neumann parts of the boundary; the
piece determinants have a constant sign on the box; Liu's bound $\lambda_k\ge\lambda_{k,h}/(1+0.1893^2h^2\lambda_{k,h})$ is applied with $h$ an upper
bound of the element diameters over the box; Temple's inequality for the discrete pencil uses a shift $\sigma\le\lambda_{2,h}$ certified on the
whole box by one verified $LDL^\top$ inertia count of a Loewner minorant (strain transfer, Supplement 1.10); every residual and quadratic form is
enclosed with explicit error terms; the P2 Rayleigh–Ritz matrix $K_2-\tfrac52L_1M_2$ is negative definite (Schur test) with $M_2\succ0$.
*Conclusion:* $\lambda_2(Q)<\tfrac52\lambda_1(Q)$ for every dart of the box. Neumann cuts (P13) add tail constants $c$, entering as
$\lambda_1\ge\min(L_1,c)$.

**P11. BO2** (Supplement S Lemmas 11.1–11.3; `darts/bo1d.py`, `darts/bo2.py`). *Input:* arm specifications (polar boxes or far arms $\rho\ge R$), a
direction $u$, an inner dart. *Hypotheses:* slab-wise upper bounds of the longest slice component $\ell(x)$ valid for every dart of the family
(convex hulls of the annular sectors, half-strip containment for far arms); on unbounded intervals a decaying-exponential start and an
Allegretto–Piepenbrink positivity argument; the inner dart is contained in every dart of the box (Lemma 11.2(d)). *Conclusion:*
$\lambda_1\ge\inf\sigma(-d^2/dx^2+\pi^2/\ell^2)>E$ and $\lambda_2\le\lambda_2(\text{inner dart})<\tfrac52E$, uniformly in the arm lengths.

**P12. INC** (Supplement S Lemma 11.4; `darts/inc.py`). *Input:* a polar box; an outer and an inner dart. *Hypotheses:* the outer apexes satisfy
the angle conditions of Lemma 11.4(b), so the outer dart contains every dart of the box, and the inner apexes those for the inner dart
(`inc.verify_inclusion`, exact); a finite-element lower bound $E$ for $\lambda_1$ of the outer dart at a point (P10 at a single parameter); a P2
upper bound for $\lambda_2$ of the inner dart with $\lambda_2<\tfrac52E$. *Conclusion:* by domain monotonicity
$\lambda_1(Q)\ge\lambda_1(Q_{\rm out})\ge E$ and $\lambda_2(Q)\le\lambda_2(Q_{\rm in})<\tfrac52E$ for every dart $Q$ of the box.

**P13. Neumann bracketing and cuts** (Supplement S 1.2, 1.3, Lemmas 9.1–9.4, 9.6; `darts/charts.py`, `darts/fam5.py`). *Input:* a cut $\Sigma$, a
finite union of segments, splitting the dart into a chart domain and tips. *Hypotheses:* for the oblique cut of Lemma 9.6 the conditions
(H1)–(H3) of that lemma, verified in ball arithmetic for every parameter of the box; every chord of a tip component, in the stated direction,
ends on the Dirichlet boundary and has length at most $\ell$. *Conclusion:* $\lambda_1(Q)\ge\min(\mu_1(\text{chart; Neumann on }\Sigma),\pi^2/\ell^2)$,
uniformly in the arm length; $\lambda_2$ is bounded above on a domain contained in $Q$ by domain monotonicity.

**P14. Exact boxes and grids** (§5.4; `darts/nanguard.py`, `darts/grid.py`, `darts/accounting.py`, `global2`, `bandB`). *Hypotheses:* every box denotes
the closed box $[p_0-r,p_0+r]$ in exact reals; convex and band boxes are dyadic rationals reconstructed from bisection paths; dart leaves, skipped
cells and chart boxes have exactly representable end points (checked); processed legacy roots are represented by outward dyadic snapping.
*Conclusion:* the leaves of each root tile it exactly, and the represented roots of every family cover its parameter ranges (checked in
rational arithmetic).

**P15. Replay semantics** (§6). A log proves the statement for its region if and only if its fail-closed replay passes: every root is complete,
the leaves tile the roots exactly, every discard is re-decided by P7, every certificate is recomputed from its record, and unknown,
duplicate, malformed or pending records are errors.

# Appendix B. Dependencies of the theorems

For each result: the written argument, the certificate primitives of Appendix A, the logs and records, the replay, and the separately
written check. File names refer to the research archive.

**Lemmas A and B** (§§2.2–2.3). Written proofs; no certificates. The constants were checked in exact arithmetic, also in the
referee's audit.

**Theorem 2** (local step, $0<\lvert p\rvert_\infty\le1/20$; §3, Lemmas 3.1, 3.2 and Corollary 3.3).

- Primitives: P1, P2, P5, P6.
- Records: `local/replay_v2/` (exact correction coefficients, remainder radii, the 960 direction boxes with exact bounds, receipt).
- Replay: `local/replay_local_v2.sh`.
- Separately written checks: Lemma 3.2 by SymPy; the first-order matrices were also derived by the referee and in
  `referee_response/arbon_perturbation.py`; `verify/xlocal.py` checks tiling, shift levels and signs of the published records.

**Theorem 3** (shell, $0.049<\lvert p\rvert_\infty\le0.225$; §§4.1–4.4).

- Primitives: P1, P3a, P4, P5, P6, P7.
- Log: `global2/runs/shell0225_cert.jsonl`.
- Replay: fail-closed replay of all 8 217 leaves.
- Separately written check: `verify/xcover.py`, exact tiling and every exclusion.

**Transport cover** (the remaining convex region; §§4.1–4.4).

- Primitives: P1, P3b, P4, P5, P6, P7.
- Log: `global2/runs/transport_cert.jsonl`.
- Replay: clean hash-bound replay of all 50 182 non-shell leaves.
- Separately written checks: `verify/xcover.py` (tiling, exclusions); `verify/xcert.py` (certificates of a stratified sample).

**Lemma $B'$** ($1/10<s\le2/5$; §4.5).

- Primitives: P7, P8, P9.
- Logs: `bandB/runs/recert_chart1.jsonl` and `recert_chart2.jsonl`.
- Replay: fail-closed replay.
- Separately written checks: `verify/xbandcover.py` (tiling, exclusions); `verify/xband.py` (certificates of a sample).

**Theorem 1$'$** (convex quadrilaterals; §4.6). Assembles Lemmas A, B, $B'$, Theorems 2, 3 and the transport cover with Endo–Osting's
Theorem 2.1 and the $D_4$ symmetry.

**Theorem 5** (darts; §5 and Supplement S).

- Primitives: P7, P10–P14.
- Logs: `darts/certs4/` (8 families), `darts/certs/` and `darts/certs2/` (6 earlier runs).
- Replay: aggregate fail-closed replay of all 14 logs; the receipt is verified by `verify/xreceipt.py`.
- Separately written check: `verify/xdart.py` (coverage, by sampling). The certificates have no second implementation.

**Theorem 1** (§1). Theorem 1$'$ and Theorem 5.

# Appendix C. Assurance record

This appendix records the review history that the main text summarises. None of the findings below changed a mathematical statement
of the paper; several changed code, and every affected log was re-certified or replayed afterwards.

**Adversarial reviews.** A second model (GPT-5.6 "Sol", OpenAI), prompted to refute, reviewed every part. (a) Lemmas A and B and the
local step: two defects in the certificate code (an endpoint rounded to $0.039999999999$; Taylor terms dropped when lowering the order
without bounding them), both corrected before the replay reported here; Lemmas A and B, the Kato–Temple and inertia arguments and the
analysis of Arbon's local argument were confirmed. (b) The global covers: the mathematics of both certificates (conformity of the
transported test functions, exact divergence and normal continuity of the Piola flux, the shifts $\beta$, the inertia criterion, the
Taylor-model algebra, coverage and the a-priori box) passed; the review found assurance defects (a replay cache keyed only by path,
status handling that did not fail closed — also in our separately written coverage checker — and the then-unfinished
recertification by the separately written checker), all repaired before the receipts above. (c) Lemma $B'$: a fatal implementation defect (floating-point products formed
before conversion to balls in the Gram matrices), an unproved floating-point remainder radius, and a bug in our separately written checker
(§4.5); all repaired and the cover re-certified. (d) Darts, in five rounds. Round 1 refuted the first package. It found nearest-float
conversions in rigorous decisions (the finite-element Schur test, strain factors, arctangent and cotangent endpoints of spike and FARL boxes, the FARL tail and
the BO2 projection range), 81 uncovered slivers between skipped spike grid cells, an unguarded NaN path in a Taylor-model range,
a replay not bound to its inputs, gaps in Lemmas 5.1, 5.2, 9.3 and 9.4, and a missing lemma for oblique cuts of long arms (now Lemma 9.6).
All were repaired and every family leaf re-certified. Round 2 cleared the numerics and Lemma 9.6, and rounds 2 and 3 found defects
only in the replay harness: shard-count and quick-mode validation, time-of-check snapshots, the floating-point environment, the
aggregate receipt and reuse of earlier replays. These were repaired and tested by probes that corrupt each input in turn. Rounds 4 and 5 found no defect in the certificates and none in the kernel. They found gaps in the provenance of the runs: inputs that a concurrent honest edit could change between the start and end checks, an interpreter that could drift between logs, and a receipt validator that accepted a run with no terminal record. These were closed operationally. Every pinned file and the directories holding them were made immutable for the duration of the runs. The reviewed pin digest became a launch precondition, and the runtime identity was anchored to an earlier run. A separately written receipt verifier (`verify/xreceipt.py`, an exact state machine over the 33 receipt records) checks the completed run. It was itself mutation-tested: 27 re-chained forgeries were rejected. Before the final replay, a recomputation of a sample containing every certificate method found a defect introduced by the round-4 changes. A local variable in the replay script shadowed the module used by the end-of-run checks, so every shard containing finite-element leaves crashed after its last leaf. The certificates themselves were unaffected. The variable was renamed, the pin regenerated (it differs only in that file), and the sample re-run on shards containing 234, 174, 106 and 16 finite-element leaves of PP, SP_col, SP_spk and FP, with no failure.
The structural identities of Lemma 3.2 are verified symbolically (sympy), and the Taylor
coefficients of the local step agree with direct high-order quadrature.

*Lemma $B'$.* An adversarial review found that the first implementation formed products of nodal values in floating point before converting them to balls, so the Gram enclosures did not contain the exact forms (errors near $10^{-15}$ against radii near $10^{-34}$), and that one series remainder radius was computed without directed rounding. Both were repaired (all inputs promoted to balls before any arithmetic; the remainder computed with directed upper bounds), new tests check strict containment of every element integral against exact rational references (the old assembly fails them on all 72 test meshes), and every leaf was re-certified with the repaired code: no leaf changed status or needed further subdivision.

*The separately written transported certificate.* While building `verify/xcert.py` we found and fixed two defects of our own in it: a fixed
choice of split diagonal, degenerate for near-triangles, and a criterion demanded at non-convex points of boxes that straddle the convexity
boundary. The draft's account of its sample ("52 of 60, 8 over budget, none fails") was wrong: the records gave 51 verified, 8 over
the work budget and 1 that had failed with the version of the checker before the second fix and had not been retried. All 9 were re-run
with the current version and a budget of 200 work units, and all 9 were verified (`referee_response/xcert9.log`).

*Portability and wording (2 October 2026).* While freezing the archive we found that nine pinned files of the dart verifier, and four
checking and test files, contained the producer's absolute package path. Run from any other location, the verifier would import
modules from that path; its runtime binding check is fail-closed, so a copy of the archive could not have been replayed at all. The
constants were replaced by paths derived from each file's location, and the launcher now takes its interpreter from `PY`. The patch
has 78 changed lines in 13 files (`darts/review/PORTABILITY_PATCH_2026-10-02.diff`), and no predicate, certificate, tolerance or
log changed. The verifier was re-pinned, the rigour lint and all 82 infrastructure probes passed, and the dart replay of §5.4 was then
run from the fresh extraction. While typesetting Supplement S, a token-level comparison with its source found the following, all
corrected in `darts/LEMMAS.md` (`darts/review/LEMMAS_WORDING_2026-10-02.diff`):

- the Weyl step of §1.9 is strict;
- Lemma 11.3(iv) certifies $\mu\ge E$, not $\mu>E$, and only $\mu\ge E$ is used;
- the reflection in the proof of Lemma 5.2 applies to $Q_R$;
- the parallel sides in Lemma 9.6(c) have lengths $L\sqrt{1+u^2}$ and $(L-k_0u)\sqrt{1+u^2}$;
- in Lemma 11.1 it is $v$, not $v^2$, that is smooth where positive;
- the second sector eigenvalue in Lemma 5.3(b) is $j^2_{\nu,2}/x^{*2}$;
- §7 now states that outside $[40^\circ,140^\circ]$ the long-arm face is covered by the one-dimensional families alone;
- a citation year was wrong.

No statement used by a certificate changed.

*External referee report (2 October 2026).* The report recommended major revisions. Its requests and their dispositions are recorded in
`REFEREE_RESPONSE.md`; in summary: the complete supplementary proofs (Supplement S), the certificate inputs, logs, receipts and checker sources
are published in the archive; the angular hypothesis of Lemma 5.1 is stated in the main text (it was in the supplement and in the predicate,
and no final log uses the lemma); Appendix A specifies every certificate primitive; Corollary 3.3 states the exact identities behind every
deleted Taylor coefficient; the assurance table above separates the kinds of evidence; the comparison with Arbon's (3.16) is an exact
calculation; the bibliography is completed and checked. During the response we also found and corrected two defects of our own
(see the response): a floating-point remainder in a rarely used branch of the separately written band checker, and a nearest-rounded bound
in the choice of shift levels of the local step; both were made rigorous and the affected checks were re-run.
