---
title: "Supplement S. Non-convex quadrilaterals (darts): lemmas and proofs"
author: "Anonymous"
date: "2026-10-02"
---

*This supplement is a typeset rendering of `darts/LEMMAS.md` in the research archive; that file is the reviewed source and governs in case of any discrepancy.*

# Darts: analytic lemmas and the rigour of the box certificate \newline (2026-10-01)

Everything here is about the Dirichlet Laplacian. $\xi(\Omega) = \lambda_2(\Omega)/\lambda_1(\Omega)$. The target is $\xi(Q) < 5/2$ for every dart.

## 0. Setting, charts, symmetries

A **dart** is a simple quadrilateral with one reflex vertex. Label the reflex vertex $D$ and the opposite vertex $B$, and write
the remaining two vertices $A$, $C$ so that the boundary is $A\to B\to C\to D\to A$. The diagonal $BD$ lies inside $Q$, and
$$Q = \operatorname{int}(ABD \cup CBD \cup BD),\qquad T := ABC\ \text{(convex hull)},\qquad N := ACD\ \text{(notch)},\qquad Q = \operatorname{int}(T) \setminus N.$$

**$BD$ chart.** By a similarity, $B = (0,0)$, $D = (1,0)$, $A = (a_x, a_y)$ with $a_y > 0$ and $C = (c_x, c_y)$ with $c_y < 0$. Every
quadrilateral of the form $\operatorname{int}(ABD \cup CBD \cup BD)$ with $A$, $C$ on opposite sides of $BD$ is simple; it is a dart with reflex vertex $D$ iff
$$f_D := (1 - c_x)\,a_y + c_y\,(a_x - 1) < 0 \qquad \text{(cross product } (D-C)\times(A-D) < 0\text{)}.$$

A dart whose reflex vertex is $B$ is the same shape as one with reflex vertex $D$ after relabelling, so the darts are exactly the
points $(A, C)$ of the chart with $f_D < 0$. The reflection in the line $BD$ gives the symmetry $(a_x,a_y,c_x,c_y) \mapsto (c_x,-c_y,a_x,-a_y)$.

Angles: $b_A := \angle ABD$, $b_C := \angle CBD$ (so $\theta_B = \angle ABC = b_A + b_C$), $\beta := \operatorname{dist}(D, AC)/\operatorname{dist}(B, AC)$ (notch depth).
The shape space of darts is the open 4-simplex $\{(\angle BAD, b_A, b_C, \angle BCD) > 0,\ \text{sum} < \pi\}$; its faces are: the shallow face
(reflex angle $\to \pi$, $D \to AC$, $\beta \to 0$), the faces $b_A \to 0$ or $b_C \to 0$ (one triangle becomes a sliver along $BD$: collapse $A \to D$ or a
thin spike), and the faces $\angle BAD \to 0$ or $\angle BCD \to 0$ (an arm of infinite length).

## 1. Tools (all standard; proofs given or referenced)

**1.1 Domain monotonicity.** $\Omega \subset \Omega'$ (open) $\Rightarrow$ $\lambda_k(\Omega) \ge \lambda_k(\Omega')$.

**1.2 Mixed problems.** For a bounded Lipschitz domain $\Omega$ with $\partial\Omega = \Gamma_D \cup \Gamma_N$ ($\Gamma_N$ relatively open, a finite union of segments)
let $V(\Omega)$ be the $H^1$-closure of $\{u \in C^1(\bar\Omega) : u = 0 \text{ near } \Gamma_D\}$ and $\mu_k$ the min–max values of $\int|\nabla u|^2/\int u^2$ on $V$.

*Neumann bracketing:* if $\Omega \setminus \Sigma = \Omega_1 \cup \Omega_2$ (disjoint open, $\Sigma$ a finite union of segments) then $\lambda_1(\Omega) \ge \min \mu_1(\Omega_i)$, where on $\Omega_i$ the
cut $\Sigma \cap \partial\Omega_i$ is Neumann and $\partial\Omega \cap \partial\Omega_i$ Dirichlet (restriction maps $H^1_0(\Omega)$ into $V(\Omega_1) \oplus V(\Omega_2)$).

*Monotonicity for mixed problems:* if $\Omega \subset \Omega'$, $\partial\Omega \cap \Omega' \subset \bar\Gamma_D$ and $\Gamma_N \subset \Gamma'_N$, then zero extension maps $V(\Omega)$ isometrically into
$V(\Omega')$, so $\mu_1(\Omega) \ge \mu_1(\Omega')$.

**1.3 Chord bound.** Let $e$ be a unit vector. If every connected component of $\Omega \cap (x + \mathbb R e)$ is a segment of length $\le \ell$ whose two
endpoints lie in $\bar\Gamma_D$, then $\mu_1(\Omega) \ge \pi^2/\ell^2$. *Proof.* For $u \in C^1$ vanishing near $\Gamma_D$, on each component $I$ the 1-D Dirichlet
Poincaré inequality gives $\int_I |\partial_e u|^2 \ge (\pi/|I|)^2 \int_I u^2 \ge (\pi/\ell)^2 \int_I u^2$; integrate over the transverse variable (Fubini) and use
density. $\square$ (Pure Dirichlet case: $\Gamma_N = \emptyset$.)

**1.4 Parallelogram.** For $P = \{(x,y) : 0 < x < L,\ sx < y < sx + h\}$, $\lambda_2(P) \le \pi^2(1+s^2)/h^2 + 4\pi^2/L^2$.

*Proof.* $u_k = \sin(k\pi x/L)\sin(\pi(y - sx)/h)$, $k = 1, 2$. With $\xi = x/L$, $\eta = (y-sx)/h$: $\partial_x u = (k\pi/L)\cos(k\pi\xi)\sin(\pi\eta) - s(\pi/h)\sin(k\pi\xi)\cos(\pi\eta)$,
$\partial_y u = (\pi/h)\sin(k\pi\xi)\cos(\pi\eta)$. Every mixed product contains $\int_0^1 \sin(\pi\eta)\cos(\pi\eta)\,d\eta = 0$ or $\int_0^1 \cos(\pi\xi)\cos(2\pi\xi)\,d\xi = 0$ or
$\int_0^1 \sin(\pi\xi)\sin(2\pi\xi)\,d\xi = 0$, so $u_1$, $u_2$ are orthogonal in $L^2$ and in energy, with Rayleigh quotients $\pi^2(1+s^2)/h^2 + k^2\pi^2/L^2$. $\square$

**1.5 Parallelogram in a convex region with an apex.** Let $K = \{(x,y) : 0 < x < X,\ \beta(x) < y < \tau(x)\}$ with $\tau$ concave, $\beta$ convex,
$\tau(0) = \beta(0) = 0$, and let $c(x) = \tau(x) - \beta(x)$ attain its maximum $\ell$ at $x^* \in (0, X]$. Then for every $\theta \in (0,1)$, $K$ contains the
parallelogram $\{(1-\theta)x^* < x < x^*,\ |y - sx| < (1-\theta)\ell/2\}$ with $s = (\tau(x^*) + \beta(x^*))/(2x^*)$. *Proof.* For $0 \le x \le x^*$, concavity and
$\tau(0) = 0$ give $\tau(x) \ge (x/x^*)\tau(x^*)$; convexity and $\beta(0) = 0$ give $\beta(x) \le (x/x^*)\beta(x^*)$. The vertical section of $K$ at $x$ therefore
contains the interval centred at $sx$ of length $(x/x^*)\ell \ge (1-\theta)\ell$. $\square$

**1.6 Temple's inequality for a matrix pencil.** $K$ symmetric, $M$ symmetric positive definite, eigenvalues $\lambda_1 \le \lambda_2 \le \dots$ of $Kv = \lambda Mv$.
For $x \ne 0$ let $\rho = x^{\mathsf T}Kx/x^{\mathsf T}Mx$ and $\varepsilon^2 = (Kx-\rho Mx)^{\mathsf T}M^{-1}(Kx-\rho Mx)/x^{\mathsf T}Mx$. If $\sigma \le \lambda_2$ and $\rho < \sigma$ then $\lambda_1 \ge \rho - \varepsilon^2/(\sigma - \rho)$.

*Proof.* Expand $x = \sum c_k\varphi_k$ ($M$-orthonormal eigenvectors), $b = x^{\mathsf T}Mx$. Then $\sum(\lambda_k-\rho)(\lambda_k-\sigma)c_k^2 = \varepsilon^2 b$ (the cross term vanishes since
$\sum(\lambda_k-\rho)c_k^2 = 0$), and $\sum(\lambda_k-\sigma)(\lambda_k-\lambda_1)c_k^2 \ge 0$ termwise ($\lambda_k \ge \sigma$ for $k \ge 2$). Subtracting, $\varepsilon^2 b + (\rho-\lambda_1)(\rho-\sigma)b \ge 0$. $\square$

**1.7 Liu's lower bound for Crouzeix–Raviart eigenvalues** (X. Liu, Appl. Math. Comput. 267 (2015); interpolation constant from
X. Liu, F. Kikuchi, J. Math. Sci. Univ. Tokyo 17 (2010)). Let a mixed problem (1.2) on a polygon be discretised by
Crouzeix–Raviart elements on a conforming triangulation that resolves $\Gamma_D$ (DOFs on $\Gamma_D$-edges set to 0, DOFs on $\Gamma_N$-edges free),
max element diameter $h$, discrete eigenvalues $\lambda_{k,h}$ of $(K, M)$. Then $\lambda_k \ge \lambda_{k,h}/(1 + C_h^2\lambda_{k,h})$, $C_h = 0.1893\,h$.

*Proof (reproduced, so that the only external input is the constant).* Let $\Pi$ be the CR interpolation (edge means). For $u \in V$,
$\Pi u$ is admissible (edge means vanish on $\Gamma_D$), and on each element $K$, $\int_K \nabla(u-\Pi u)\cdot\nabla v_h = \nabla v_h\cdot\int_{\partial K}(u-\Pi u)\,n = 0$ for every
CR function $v_h$; hence $a_h(u - \Pi u, v_h) = 0$ (broken energy) and $a(u,u) = a_h(\Pi u,\Pi u) + a_h(u-\Pi u,u-\Pi u)$. On every triangle,
$u - \Pi u$ has zero mean on the three edges, so $\|u - \Pi u\|_{L^2(K)} \le 0.1893\,h_K\,\|\nabla(u - \Pi u)\|_{L^2(K)}$ (Liu–Kikuchi), and summing,
$\|u - \Pi u\| \le C_h\,\|u - \Pi u\|_{a_h}$. Let $E_k$ be the span of the first $k$ exact eigenfunctions. If $C_h^2\lambda_k \ge 1$ the claim is trivial
($x/(1+C^2x) < 1/C^2$). Otherwise $\Pi$ is injective on $E_k$ ($\Pi u = 0$ would give $\|u\| \le C_h\|u\|_a \le C_h\sqrt{\lambda_k}\,\|u\| < \|u\|$), and for $u \in E_k$
with $\|u\| = 1$, $A := \|u\|_a^2 \le \lambda_k$, $t := \|u-\Pi u\|_{a_h}$: $\mathrm{RQ}(\Pi u) \le (A - t^2)/(1 - C_h t)^2 \le A/(1 - C_h^2A)$ (maximise over $t$; the maximum is
at $t = C_hA$). Min–max on $\Pi(E_k)$: $\lambda_{k,h} \le \lambda_k/(1 - C_h^2\lambda_k)$, i.e. $\lambda_k \ge \lambda_{k,h}/(1 + C_h^2\lambda_{k,h})$. $\square$

*Check of the constant (numerical, not part of the proof):* the sharp constant $\sup \|w\|/(h\|\nabla w\|)$ over $w \in H^1(K)$ with zero means on
all three edges, computed by P3 FEM with the three linear constraints (two refinement levels agree to 5 digits), is 0.18918 for the
equilateral triangle, 0.1885 for a needle, 0.1681 for the right isosceles, 0.162–0.175 for other shapes — consistent with 0.1893 as
the supremum over triangles.

**1.8 Rayleigh–Ritz.** For linearly independent $u_1, u_2 \in H^1_0(\Omega)$, $\lambda_2(\Omega) \le$ largest eigenvalue of the $2\times2$ pencil
$((\nabla u_i,\nabla u_j), (u_i,u_j))$. Continuous P2 finite element functions vanishing on $\partial\Omega$ are in $H^1_0(\Omega)$.

**1.9 Verified inertia.** If the float $LDL^{\mathsf T}$ factorisation (no pivoting, band $b$) of a symmetric float matrix $A_{\mathrm{fl}}$ has exactly one
negative pivot and no zero pivot, then $\hat L\hat D\hat L^{\mathsf T} = A_{\mathrm{fl}} + \Delta A$ with $|\Delta A| \le \gamma_{b+1}|\hat L||\hat D||\hat L^{\mathsf T}|$ (Higham, ASNA, Thm 9.3 and the remark on
banded matrices; symmetric elimination is Doolittle with $U = DL^{\mathsf T}$), and by Sylvester $A_{\mathrm{fl}} + \Delta A$ has exactly one negative and no zero
eigenvalue. With $A_{\mathrm{ex}} = K - \sigma'M$ (exact) and $\|A_{\mathrm{fl}} - A_{\mathrm{ex}}\|_2 \le e_A$, Weyl gives $\lambda_2(A_{\mathrm{ex}}) > -(e_A + \|\Delta A\|_2)$ (strict, since $\lambda_2(A_{\mathrm{fl}} + \Delta A) > 0$); since $M$ is diagonal
positive (CR), $\lambda_2(K - \sigma M) \ge \lambda_2(A_{\mathrm{ex}}) + (\sigma'-\sigma)\min_i M_{ii} > 0$ for $\sigma := \sigma' - (e_A + \|\Delta A\|_2)/\min M_{ii}$. Hence the pencil has at most one
eigenvalue $\le \sigma$, i.e. $\lambda_{2,h} > \sigma$. ($\|\cdot\|_2 \le$ max row sum for the symmetric bound matrices.) Implemented in `ldl.py`, `sigma.py`.

**1.10 Strain transfer between parameter points.** In the $BD$ chart each triangle $X \in \{A, C\}$ is the image of one fixed reference mesh
under $F_X = \bigl[\begin{smallmatrix}1 & x\\ 0 & y\end{smallmatrix}\bigr]$ ($X = (x, y)$); for CR (or P2) vectors $v$, $a_p(v) = \sum_X \int\hat\nabla\hat v^{\mathsf T} G_X(p)\,\hat\nabla\hat v$, $m_p(v) = \sum_X J_X(p)\,\hat m_X(v)$, with
$G_X = |\det F_X|\,F_X^{-1}F_X^{-\mathsf T}$ ($\det G_X = 1$) and $J_X = |y|$. Since $v^{\mathsf T}G_X(p)v/v^{\mathsf T}G_X(q)v \ge (\eta/\eta_q)\,\lambda_{\min}(TT^{\mathsf T})$,
$T = F_X(q)F_X(p)^{-1} = \bigl[\begin{smallmatrix}1 & (x_q - x)/\eta\\ 0 & \eta_q/\eta\end{smallmatrix}\bigr]$ ($\eta = |y|$), min–max gives
$$\lambda_{k,h}(p) \ge \lambda_{k,h}(q)\cdot\min_X\bigl[(\eta_X/\eta_{X,q})\,\lambda_{\min}(T_XT_X^{\mathsf T})\bigr]\big/\max_X\bigl(J_X(p)/J_X(q)\bigr),\qquad \lambda_{\min} = 1 + \bigl(a - \sqrt{a^2 + 4b^2}\bigr)/2,$$
$a = b^2 + c^2 - 1$, $b = (x_q - x)/\eta$, $c = \eta_q/\eta$, evaluated in interval arithmetic over a box (monotone in $a$ and $b^2$).

## 2. Prop 5.1 (shallow notch) — from paper §5, restated

If $\beta \le 0.0339$ then $\xi(Q) \le (7/3)(1 - \beta)^{-2} < 5/2$. (Homothety at $B$: $T_s := B + (1-\beta)(T - B) \subset Q \subset T$; triangle theorem.)
In the chart, $\beta = |n\cdot(D - A)|/|n\cdot(B - A)|$ with $n = (a_y - c_y, c_x - a_x)$ normal to $AC$; on a box this is enclosed by interval
arithmetic (`cover.beta_max`), and boxes with $\beta \le 0.0339$ throughout are discarded.

## 3. Lemma 3.1 (collapse collar: the reflex vertex close to a neighbouring vertex)

Let $\delta := \operatorname{dist}(A, \text{segment } BD)$ and $r :=$ inradius of the triangle $CBD$. Then $\xi(Q) \le (7/3)(1 + \delta/r)^2$, which is $< 5/2$ when $\delta \le 0.035\,r$.
(Symmetrically with $A$ and $C$ exchanged.)

*Proof.* $\operatorname{dist}(\cdot, \text{segment } BD)$ is convex, so its maximum on the triangle $ABD$ is attained at a vertex, i.e. it is $\delta$; hence
$ABD \subset \bar N_\delta(BD) \subset \bar N_\delta(CBD)$. With $I$ the incentre of $CBD$, each side line of $CBD$ is at distance $r$ from $I$, so
$\bar N_\delta(CBD) \subset I + (1 + \delta/r)(CBD - I) =: H$. Thus $Q \subset \operatorname{int} H$ and $\lambda_1(Q) \ge \lambda_1(CBD)/(1+\delta/r)^2$; $Q \supset CBD$ gives
$\lambda_2(Q) \le \lambda_2(CBD) \le (7/3)\lambda_1(CBD)$ (triangle theorem: Siudeja 2010; Arbon–Mannan–Psenka–Ragavan 2022). Finally $(7/3)\cdot 1.035^2 = 2.4995\ldots$ $\square$

*Coverage.* In the $BD$ chart, $A = D + \rho(\cos\varphi, \sin\varphi)$: $\delta \le \rho$. For $C$ equilateral ($|BD| = 1$, $r = 0.2887$) this is $\rho \le 0.0101$.
This collar is a neighbourhood of the face $b_A = 0$ *only where $A$ is close to $D$*; where the sliver $ABD$ is long ($A$ far along the ray
beyond $D$: a thin spike) the lemma does not apply (points of $ABD$ near $A$ are far from $CBD$).

## 4. Lemma 4.1 ($BD$-chord bound) — valid for every dart

$\lambda_1(Q) \ge \pi^2/|BD|^2$.

*Proof.* Write points as $P = B + a\,e_A + c\,e_C$ ($e_A$, $e_C$ unit vectors along $BA$, $BC$; $a, c \ge 0$ on $T$) and $D = B + d_A e_A + d_C e_C$,
$d_A, d_C > 0$. (i) If $P \in T$ with $a \ge d_A$ and $c \ge d_C$ then $P \in N$: indeed $P - D \in \operatorname{cone}(e_A, e_C)$, and $e_A \in \operatorname{cone}(A-D, C-D)$ because
$e_A = \mu(A-D) + \nu(C-D)$ with $\nu = \mu d_C/(|BC| - d_C) \ge 0$ and $\mu\cdot|BA||BC|(1 - d_A/|BA| - d_C/|BC|)/(|BC| - d_C) = 1$, $\mu > 0$ ($D$ is
interior to $T$); similarly $e_C$. So $P \in (D + \operatorname{cone}(A-D, C-D)) \cap T =$ triangle $ACD = N$ (the cone meets $T$ exactly in $ACD$ because $D$
and $B$ lie on the same side of $AC$). Hence $Q \subset W := \{a > 0,\ c > 0,\ (a < d_A \text{ or } c < d_C)\}$.

\(ii) A line parallel to $BD$, parametrised as $(a, c) = (a_0 + \tau d_A, c_0 + \tau d_C)$, meets $W$ exactly for
$\tau_{\min} < \tau < \tau_{\min} + 1$, $\tau_{\min} := \max(-a_0/d_A, -c_0/d_C)$: the conditions $a, c > 0$ give $\tau > \tau_{\min}$ and "$a < d_A$ or $c < d_C$" gives
$\tau < 1 + \tau_{\min}$. So every chord of $Q$ parallel to $BD$ lies in a segment of length $|BD|$; apply 1.3. $\square$

## 5. Lemma 5.1 (thin wedge at $B$) — a complete analytic treatment of the corner $\theta_B \to 0$

If **both angles are acute**, $b_A, b_C \in (0, \pi/2)$ (equivalently $a_x > 0$ and $c_x > 0$), and $\tan b_A + \tan b_C \le 1/10$, then
$\xi(Q) \le 2.4223 < 5/2$. In particular this holds if $\theta_B = b_A + b_C \le 0.0996$ rad ($5.7^\circ$): then both angles are acute and, $\tan$ being
superadditive on $[0, \pi/2)$, $\tan b_A + \tan b_C \le \tan\theta_B \le \tan 0.0996 = 0.09993\ldots < 1/10$.

(The acute-angle hypothesis is essential: $\tan b_A + \tan b_C \le 1/10$ alone admits obtuse angles with negative tangents, e.g.
$b_A = 10^\circ$, $b_C = 160^\circ$, $\rho_A = \rho_C = 10$, a dart with $c_x < 0$ to which the proof below does not apply. The rigorous predicate
`regions.P_lemma51` checks $a_x > 0$, $c_x > 0$, $a_y \ge 0 \ge c_y$ and the tangent sum on the whole box; the family runs use Lemma 5.3,
whose hypothesis $\theta_B \le 0.1772$ implies acute angles.)

(The proof does not use the reflex condition; it holds for every quadrilateral $\operatorname{int}(ABD \cup CBD \cup BD)$ of the chart with these angles.)

*Proof.* $BD$ chart; $b_A, b_C < \pi/2$, so $A$ and $C$ have positive abscissae and $Q \subset \{x > 0\}$. Put $\tau_B := \tan b_A + \tan b_C$.

\(i) *Stem.* For $0 < x < 1$ the vertical section of $ABD$ is $[0, h_A(x)]$, where $h_A$ is linear (if $a_x \ge 1$) or a tent with apex at $a_x$
(if $a_x < 1$); in both cases $h_A$ is concave on $[0,1]$, $h_A(0) = 0$ and $h_A(x) \le x\tan b_A$ ($ABD$ lies in the angle $\angle DBA$). Likewise $CBD$
gives $[-h_C(x), 0]$. Hence for $0 < x < 1$ the section of $Q$ is the single interval $(-h_C(x), h_A(x))$, of length $c(x) := h_A + h_C \le x\tau_B$,
and the stem $S := \{0 < x < 1,\ -h_C(x) < y < h_A(x)\} \subset Q$ is convex.

\(ii) *Arms.* For $x > 1$ the section of $ABD$ is empty if $a_x \le 1$, and otherwise is $[y_{DA}(x), y_{BA}(x)]$, of length decreasing linearly
from $\tan b_A$ ($= h_A(1)$) at $x = 1$ to 0 at $x = a_x$; it lies above the $x$-axis. Likewise $CBD$ below. The two pieces are disjoint (the open
ray beyond $D$ is not in $Q$), so every vertical chord component for $x \ge 1$ has length $\le \max(h_A(1), h_C(1)) \le c(1^-)$.

\(iii) Hence every vertical chord of $Q$ has length $\le \ell := \max_{(0,1]} c = c(x^*)$, and $\lambda_1(Q) \ge \pi^2/\ell^2$ by 1.3.

\(iv) By 1.5 applied to $S$ ($\tau = h_A$, $\beta = -h_C$), with $\theta = 1/4$, $S$ contains a parallelogram with vertical sides $(3/4)\ell$, length $x^*/4$ and
shear $s$, $|s| \le \max(\tan b_A, \tan b_C)/2 \le \tau_B/2$ (because $h_A(x^*) \le x^*\tan b_A$ and $h_C(x^*) \le x^*\tan b_C$). By 1.4,
$\lambda_2(Q) \le \pi^2(1 + s^2)(4/(3\ell))^2 + 4\pi^2(4/x^*)^2$.

\(v) With $\ell \le x^*\tau_B$: $\xi(Q) \le (16/9)(1 + \tau_B^2/4) + 64\tau_B^2 \le (16/9)(1.0025) + 0.64 = 2.42222\ldots$ $\square$

## 5a. Lemma 5.3 (thin wedge at $B$, sector version) — extends Lemma 5.1 to $\theta_B \le 0.1772$ rad ($10.15^\circ$)

If $\theta_B = b_A + b_C \le 0.1772$ (radians), then $\xi(Q) < 5/2$.

*Proof.* Use the notation of the proof of Lemma 5.1 ($BD$ chart, stem $S$, maximiser $x^*$ of the vertical chord length $c$, $\ell = c(x^*)$).
Let $\varphi_{\mathrm{hi}} := \arctan(h_A(x^*)/x^*) \le b_A$ and $\varphi_{\mathrm{lo}} := \arctan(h_C(x^*)/x^*) \le b_C$, $\alpha := \varphi_{\mathrm{hi}} + \varphi_{\mathrm{lo}} \le \theta_B$,
$T := \ell/x^* = \tan\varphi_{\mathrm{hi}} + \tan\varphi_{\mathrm{lo}} \le \tan\alpha$ ($\tan$ is superadditive on $[0, \pi/2)$). The triangle $W^* = \{0 < x < x^*,\ -x\tan\varphi_{\mathrm{lo}} < y < x\tan\varphi_{\mathrm{hi}}\}$ lies in $S \subset Q$ (proof of 5.1)
and contains the circular sector $\Sigma(\alpha, x^*) = \{0 < r < x^*,\ -\varphi_{\mathrm{lo}} < \theta < \varphi_{\mathrm{hi}}\}$ (a point of the sector has abscissa $r\cos\theta < x^*$).

\(a) If $\alpha \le 0.0997$ then $T \le \tan\alpha \le \tan 0.0997 = 0.1000316611\ldots$ (note: slightly above $1/10$) and the parallelogram argument of
Lemma 5.1 applies verbatim with $\tau_B$ replaced by $T$ (it used only $\ell \le x^*T$ and $|s| \le T/2$), giving
$\xi \le (16/9)(1 + T^2/4) + 64T^2 \le 2.4226304 < 5/2$.

\(b) If $\alpha \in [\alpha_k, \alpha_{k+1}]$ with $\alpha_k = 0.0997 + 0.0025k$ (exact decimals), $k = 0, \dots, 30$ (31 intervals; $\alpha_{31} = 0.1772$), then
$\Sigma(\alpha, x^*)$ contains a sector of angle $\alpha_k$
with the same apex and radius, whose second Dirichlet eigenvalue is $j^2_{\nu,2}/x^{*2}$ with $\nu = \pi/\alpha_k$ (eigenfunctions
$J_{m\nu}(jr/x^*)\sin(m\nu(\theta - \theta_0))$; the second eigenvalue is $j^2_{\nu,2}/x^{*2}$ because $j_{\nu,2} < 2\nu < j_{2\nu,1}$). Hence $\lambda_2(Q) \le j^2_{\nu,2}/x^{*2}$, while
$\lambda_1(Q) \ge \pi^2/\ell^2 \ge \pi^2/(x^*\tan\alpha_{k+1})^2$ (chord bound, as in 5.1). So $\xi(Q) \le (j_{\nu,2}\tan\alpha_{k+1}/\pi)^2$. For each $k$, `lemma_sector.py`
certifies in ball arithmetic (128 bits) an upper bound $X_k > j_{\nu,2}$ by the sign pattern $J_\nu(Y_k) < 0 < J_\nu(X_k)$ with $0 < Y_k < X_k$
checked explicitly ($J_\nu > 0$ on $(0, j_{\nu,1})$), checks $X_k < 2\nu$, and checks $(X_k\tan\alpha_{k+1}/\pi)^2 < 5/2$; the largest value is
$2.489242112\ldots < 2.4893$ ($k = 30$). Log: `lemma_sector.log`. $\square$

## 5b. Lemma 5.2 (flat hull triangle) — the corner $\theta_B \to \pi$

Let $\theta_A = \angle BAC$, $\theta_C = \angle BCA$ be the angles of the hull $T$ at $A$ and $C$. If **both are at most $\arctan(1/10) = 5.7106^\circ$** (so both are acute
and $\tau := \max(\tan\theta_A, \tan\theta_C) \le 1/10$), then $\xi(Q) \le (16/9)(1 + \tau^2) + 64\tau^2 \le 2.4356 < 5/2$.

(A bound on the tangents alone would admit an obtuse angle with a negative tangent; the hypothesis is on the angles. The runs use
it through $\theta_B \ge 174.29^\circ$: $\theta_B$ is the hull angle at $B$, so $\theta_A + \theta_C = \pi - \theta_B \le 5.71^\circ < \arctan(1/10)$, and both angles qualify.)

*Proof.* Place $C = (0,0)$, $A = (1,0)$, $B$ above the $x$-axis. Both base angles are acute, so $B$ projects into $(0,1)$, and $D \in \operatorname{int} T$
gives $0 < d_x < 1$. With $t(x)$ the upper boundary of $T$ (tent $C\to B\to A$) and $n(x)$ the upper boundary of $N = CDA$ (tent $C\to D\to A$),
$Q = \{0 < x < 1,\ n(x) < y < t(x)\}$, so every vertical chord of $Q$ is a single interval of length $c(x) = t(x) - n(x) \le \ell := \max c$,
and $\lambda_1(Q) \ge \pi^2/\ell^2$ (1.3). $Q_L := Q \cap \{x < d_x\} = \{0 < x < d_x,\ (d_y/d_x)x < y < t(x)\}$ is convex with apex $C$ (lower boundary
linear, upper concave, both 0 at $x = 0$); likewise $Q_R := Q \cap \{x > d_x\}$ with apex $A$. Let $x^*$ be a maximiser of $c$. If $x^* \le d_x$ apply
1.5 to $Q_L$ (whose apex $C$ is at the origin): $Q$ contains a parallelogram with vertical sides $(3/4)\ell$, length $x^*/4$ and
shear $|s| \le \max(t(x^*), n(x^*))/x^* \le \tan\theta_C$ ($\theta = 1/4$); here $\ell \le t(x^*) \le x^*\tan\theta_C$. If $x^* \ge d_x$ apply 1.5 to $Q_R$, reflected by $x \mapsto 1 - x$ so that its apex $A$ is at the origin, with $x^*$ replaced by
$1 - x^*$ and $\theta_C$ by $\theta_A$. In both cases 1.4 gives $\lambda_2(Q) \le \pi^2(1+\tau^2)(4/(3\ell))^2 + 4\pi^2(4/L)^2$ with $\ell/L \le \tau$, $L$ the distance from $x^*$ to the
apex used. Divide. $\square$

## 6. Lemma 6.1 (tails of long arms)

Let $S$ be a segment cutting $Q$ into $Q_1$ and $Q_2$, where $Q_2$ is a triangle with $S$ as a side and its other two sides on $\partial Q$. Then the mixed
problem on $Q_2$ (Neumann on $S$) has $\mu_1 \ge \pi^2/|S|^2$, hence $\lambda_1(Q) \ge \min(\mu_1(Q_1; \text{Neumann on } S), \pi^2/|S|^2)$.

*Proof.* Chords of $Q_2$ parallel to $S$ end on the two other sides (Dirichlet) and are shorter than $S$ (similar triangles); apply 1.3 and
1.2. $\square$

Consequences used in the design of the (FEM) far-arm charts, superseded in the final cover by Lemma 11.3 (see §13): (a) for a cut $S$ parallel to $BD$ inside the arm
$ABD$, $|S| \le |BD|$ (so the tail bound is at least the $BD$-chord bound); (b) the truncated core $Q_1$ of a dart with a long arm is contained,
with the same cut lines, in the truncated core of the limiting bent waveguide $W$ of Lemma 4.1 (walls = rays from $B$ and from $D$ parallel
to $BA$, $BC$), so by 1.2 its mixed eigenvalue is bounded below by that of the waveguide core, a two-parameter compact family.

## 7. The faces of the shape simplex: what is regular for the box certificate

The box certificate (§8) uses the $BD$ chart with each triangle meshed as an affine image of a fixed reference mesh.

* **Shallow face ($\beta \to 0$, reflex angle $\to \pi$):** regular — nothing degenerates in the chart, and the certificate is valid for convex
  quadrilaterals as well, so boxes may straddle the face. Prop 5.1 removes $\beta \le 0.0339$.
* **Faces $b_A \to 0$ / $b_C \to 0$ (one triangle degenerates to a sliver along $BD$):** NOT regular for the $BD$-chart certificate: $G_X$ has
  entries $\propto 1/|y|$, the reference elements in the sliver degenerate, and box sizes must shrink in proportion to $a_y$. The part of the
  face with the apex close to $D$ (collapse) is covered by Lemma 3.1. The part with a long thin spike ($A$ far along a ray close to $BD$;
  one of the critical curves, $CBD$ equilateral, $\xi \to 7/3$) is made regular by Lemmas 9.1, 9.2 and 9.4: the certificates there use the
  charts SPIKE ($C\,B\,E\,D$ with a Neumann cut $DE$, for $\tan b_A \le t_1$) and SPIKE2 (with the stub of Lemma 9.4, for $t_1 \le \tan b_A \le t_2$),
  whose parameters are $(t, C)$ only — uniform in $|A|$, including $|A| \to \infty$. Implemented (`charts.spike_problem`,
  `charts.spike2_problem`) and run for $C \in \mathrm{CBOX}$ (Theorem S) and for $C$ in the KF sector (certs2/kf\_\*); all other $C$ are covered by
  the families SP_col/SP_spk/SF_col/SF_spk with the BO2 certificate (Lemma 11.3, spike arms as 'far' arms with $R = 0$) (§13).
* **Faces $\angle BAD \to 0$ / $\angle BCD \to 0$ (infinitely long arm), including the deep notch $\beta \to 1$:** not compact in the $BD$ chart; Lemma 9.3
  (+ Lemma 6.1) reduces them to compact families in ($C$, direction of $BA$, $L$). Charts implemented (`charts.far_problem`,
  `deep_problem`, `far2_problem`) and test boxes certified, but too weak for steep arms (PILOT §5); not used. The final cover handles
  this face with the BO2 certificate with tails (Lemma 11.3 (ii), families PF/FP/FF, $\rho \ge 20$; §13). Lemma 9.5 closes the
  part of this face with $b_X \in [40^\circ, 140^\circ]$ analytically (for $|BX| \ge L_k|BD|$); outside $[40^\circ, 140^\circ]$ Lemma 9.5 does not apply (its bound cannot reach $5/2$ at all when $b_X \le 39.23^\circ$ or $b_X \ge 140.77^\circ$), and there the face is covered by the BO2 families alone.
* **Corner $\theta_B \to 0$:** Lemmas 5.1 and 5.3 (complete for $\theta_B \le 0.1772$). **Corner $\theta_B \to \pi$ (flat hull):** Lemma 5.2.

## 8. The box certificate (cert.py) — statement and rigour

For a box $P = p_0 + r\cdot[-1,1]^4$ in the chart ($a_y > 0 > c_y$ on $P$), a reference mesh level $k$ ($\kappa = 2$ grading toward $D$) and a number
$\sigma_2$, `BoxCert.certify` returns ok = True only if the following chain is verified in ball arithmetic (python-flint arb, 106 bits):

1. Exact reference data: all reference CR and P2 matrices are exact rationals rounded once to doubles (`refexact.py`);
   every float product used later carries a rigorous error bound ($\gamma_n$-type bounds + the rounding of the stored entries).
2. Fixed test vectors (any vectors are admissible): CR $x(t) = x_0 + \sum t_i x_i + \sum t_it_j x_{ij}$ (eigenvector at $p_0$ and its scaled first
   and second parameter derivatives), P2 $y_a(t) = y_a + \sum t_i y_{a,i}$, $a = 1, 2$.
3. The $p$-dependence enters only through $G_A$, $G_C$ (3 entries each) and $J_A$, $J_C$: $G_X = (1/\eta)\bigl[\begin{smallmatrix}\eta^2+x^2 & -x\\ -x & 1\end{smallmatrix}\bigr]$, $J_X = \eta$.
   $\rho(p)$ is computed in the cancelled form $\lambda_0 + \sum_X \eta_X^{-1}[x^{\mathsf T}P_X x - \lambda_0\eta_X^2\,x^{\mathsf T}\hat M_X x]/b$ with exact polynomial brackets; all scalar
   functions are Taylor models in $t$ of order 3 (`tm.ORDER`; order 2 in early probes only) with rigorous remainders (every
   truncated term is folded into the remainder).
4. The residual $Kx - \rho Mx$ is scaled on each DOF set ($A$-interior by $\eta_A$, $C$-interior by $\eta_C$, $BD$ by $\eta_A\eta_C$) so that its coefficients are
   exact polynomials; its Taylor coefficient vectors are formed with explicit error vectors; $\varepsilon^2$ is bounded by
   $(\sum_\gamma\|R_\gamma\|_w + \|\mathrm{err}\|_w)^2$ per set (weights = inverse diagonal mass, rescaled by sup-bounds of $\eta$-powers over $P$).
5. Temple (1.6) with $\sigma_2 \le \lambda_{2,h}$ on $P$ from 1.9 at a centre + 1.10 over $P$; requires $\sigma_2 > \rho$ on $P$.
6. Liu (1.7) with $h = \sup$ over $P$ of the element diameters (`hmax_box`, monotone bound in $x$ and $\eta$).
7. Rayleigh–Ritz (1.8) with the transported P2 vectors: $N = K_2 - (5/2)L_1M_2$ is shown negative definite on $P$ via $N_{11} < 0$ and
   $N_{22} + N_{12}^2/|N_{11}| < 0$ (Schur complement), with $M_2$ positive definite.

Then for every $p \in P$: $\lambda_2(Q(p)) \le$ largest pencil eigenvalue $< (5/2)L_1(p) \le (5/2)\lambda_1(Q(p))$, i.e. $\xi(Q(p)) < 5/2$.

## 9. Reduction lemmas for the non-compact faces (spikes and long arms)

**Lemma 9.1 (spike head).** In the $BD$ chart let $0 < b_A < \pi/2$ and $a_x \ge 1$ ($A$ beyond the abscissa of $D$). Put $t := \tan b_A$,
$E := (1, t)$ (the point of $BA$ above $D$), $H :=$ triangle $BDE$, $P := \operatorname{int}(CBD \cup H \cup BD)$ (a quadrilateral $C, B, E, D$). Then
$$\lambda_1(Q) \ge \min\bigl(\mu_1(P; \text{Neumann on } DE, \text{Dirichlet elsewhere}), \pi^2/t^2\bigr),\qquad \lambda_2(Q) \le \lambda_2(P)\ \ \text{(Dirichlet)}.$$
In particular $\xi(Q) \le \lambda_2(P)/\min(\mu_1(P), \pi^2/t^2)$ depends only on $(C, t)$, not on $|A|$.

*Proof.* $DE$ cuts $Q$ into $P$ and the triangle $DEA$ (if $a_x > 1$); the chords of $DEA$ parallel to $DE$ end on $EA \subset BA$ and $AD$ and have
length $\le |DE| = t$ (Lemma 6.1), and $P \cup DEA \cup DE = Q$, so Neumann bracketing (1.2) gives the first claim. $H \subset ABD$ because $a_x \ge 1$
($H = \{0 < x < 1,\ 0 < y < xt\}$ and $BA$ lies above it for $x \le a_x$), so $P \subset Q$ and 1.1 gives the second. $\square$

**Lemma 9.2 (monotonicity of the spike head in $t$).** For $0 < t \le t_1$: $\mu_1(P(t)) \ge \mu_1(P(t_1))$ and $\lambda_2(P(t)) \le \lambda_2(CBD)$.

*Proof.* $P(t) \subset P(t_1)$, $\partial P(t) \cap P(t_1) =$ the segment $BE(t) \subset$ Dirichlet part, Neumann part $DE(t) \subset DE(t_1)$: 1.2 (mixed monotonicity).
$P(t) \supset CBD$: 1.1. $\square$

Remark ($a_x < 1$). If $A = a(1, t)$ with $0 < a < 1$ and $t \le t_1$, then $A \in \bar H(t_1) = \operatorname{tri}(B, D, (1, t_1))$, so $ABD \subset \bar H(t_1)$, $Q \subset P(t_1)$ and
$\lambda_1(Q) \ge \lambda_1(P(t_1)) \ge \mu_1(P(t_1); \text{Neumann on } DE)$, while $\lambda_2(Q) \le \lambda_2(CBD)$ still holds. Hence the collar certificate below covers
EVERY $a > 0$.

Consequence: the spike collar $\{a_x \ge 1,\ \tan b_A \le t_1\}$ is covered as soon as $\lambda_2(CBD) < (5/2)\cdot\min(\mu_1(P(C, t_1)), \pi^2/t_1^2)$ — a
two-parameter ($C$) certificate; and $\{a_x \ge 1,\ \tan b_A \in [t_1, t_2]\}$ needs a three-parameter ($C$, $t$) certificate. Neither depends on $|A|$.

**Lemma 9.3 (long arm, cut parallel to $BD$; design only — the final cover uses the oblique cut of Lemma 9.6).** Let
$e_A = (A - B)/|A - B|$ and $S_A := \{B + s\,e_A + \tau(D - B) : s > 0,\ 0 < \tau < 1\}$ (the half-strip spanned by $e_A$ and $BD$). Then $ABD \subset \bar S_A$,
hence $Q \subset \operatorname{int}(CBD \cup \bar S_A)$. For $L > 0$ let $\Sigma_L := \{B + L\,e_A + \tau(D - B) : 0 \le \tau \le 1\}$ (a cut parallel to $BD$, $|\Sigma_L| = |BD| = 1$) and
suppose (i) $0 < L < |A - B|$ and (ii) the closed triangle $CBD$ lies in the open half plane on $B$'s side of the line through $\Sigma_L$
(i.e. all its points have $s < L$ in the coordinates $(s, \tau)$ of the strip). Then the line of $\Sigma_L$ crosses the segments $BA$ and $DA$ and not
$BD$, so $\Sigma_L \cap Q$ is a single segment cutting off the far part of $ABD$, and
$$\lambda_1(Q) \ge \min\bigl(\mu_1(\mathrm{core}_L; \text{Neumann on } \Sigma_L), \pi^2\bigr),\qquad \mathrm{core}_L := \operatorname{int}\bigl(CBD \cup (\bar S_A \cap \{s \le L\})\bigr),$$
which depends on ($C$, direction of $e_A$, $L$) but not on $|A|$; and for fixed $C$ and direction $e_A$, $\lambda_2(Q)$ is non-increasing in $|A|$
($ABD(R) \subset ABD(R')$ for $R < R'$ because $A(R)$ lies on the segment $B\,A(R')$).

*Proof.* $ABD = \{B + \lambda(A-B) + \mu(D-B) : \lambda, \mu \ge 0,\ \lambda+\mu \le 1\} \subset \bar S_A$. The far part $ABD \cap \{s > L\}$ is a triangle whose chords parallel to $\Sigma_L$
end on $BA$ and $DA$ and are shorter than $|BD|$ (Lemma 6.1); the near part of $Q$ is contained in $\mathrm{core}_L$ with $Q$'s piece of $\Sigma_L$ inside the
Neumann segment $\Sigma_L$ and the rest of $\partial(Q \cap \{s < L\})$ inside $\mathrm{core}_L$ Dirichlet (1.2, mixed monotonicity). $\square$

(Both arms long — the deep notch $\beta \to 1$ — is the same construction applied to both triangles; the core then depends only on the two
directions, i.e. on $(b_A, b_C)$, and $L$.)

These lemmas make every non-compact face of the shape simplex compact for the box certificate, except that the certificates
themselves (a chart with the pieces $CBD$, $H$ or $\mathrm{core}_L$, and Neumann DOFs on the cuts) have to be run. (Superseded for long arms by
Lemma 11.3; the spike-head charts are used by the KF runs and as the FEM fallback of the SP families.)

**Lemma 9.4 (long spike with a stub; used by `charts.spike2_problem`).** $BD$ chart, $t := \tan b_A > 0$, $A = a_x(1, t)$ with $a_x \ge 1$.
Put $w := 3t$, $E := (1, t)$, $D' := (1 + w, tw)$, $E' := (1 + w, t(1 + w))$ and $P' := \operatorname{int}(CBD \cup BDE \cup DD'E'E \cup \text{common edges})$, the
stub $DD'E'E$ being the parallelogram between the line through $D$ of slope $t$ and $BA$. Then
$$\lambda_1(Q) \ge \min\bigl(\mu_1(P'; \text{Neumann on } D'E', \text{Dirichlet elsewhere}), \pi^2/t^2\bigr),$$
and for every $X_0 \in [1, a_x]$: $\lambda_2(Q) \le \lambda_2(Q_{X_0})$, $Q_{X_0}$ the dart with the same $C$ and $A$ replaced by $X_0(1, t)$.

*Proof.* $ABD = BDE \cup (\text{needle } DEA)$. The line $DA$ has slope $t\,a_x/(a_x - 1) > t$, so for $1 \le x \le \min(a_x, 1 + w)$ the needle lies
above the line $DD'$ and below $BA$, i.e. inside the stub. If $a_x \le 1 + w$ then $Q \subset P'$ with $\partial Q \cap P' \subset \partial Q$ (Dirichlet), and 1.2 gives
$\lambda_1(Q) \ge \mu_1(P')$. Otherwise cut the needle by $\sigma = \{1 + w\} \times [y_{DA}(1 + w), t(1 + w)] \subset D'E'$ ($y_{DA}(1 + w) \ge tw$); the far part has vertical
chords ending on $DA$ and $BA$ of length $\le t$ (Lemma 6.1), the near part lies in $P'$ with its cut inside the Neumann segment $D'E'$ (1.2).
For $\lambda_2$: $X_0(1, t)$ lies on the segment from $B$ to $A$, so $ABD(X_0) \subset ABD$ and $Q_{X_0} \subset Q$ (1.1).

*The endpoint $a_x = 1$.* Then $A = E$, $ABD = BDE$ and $Q = \operatorname{int}(CBD \cup BDE \cup BD) \subset P'$; the boundary of $Q$ inside $P'$ is the segment $DE$ (a
common edge of $BDE$ and of the stub), on which the functions of $H^1_0(Q)$ vanish: extending them by zero gives functions in the form
domain of the mixed problem on $P'$ (Dirichlet except on $D'E'$), so $\lambda_1(Q) \ge \mu_1(P') \ge \min(\mu_1(P'), \pi^2/t^2)$ (1.2); and $X_0 = 1$ is the only
admissible $X_0$, $Q_{X_0} = Q$. $\square$

Note: in the $BD$ chart every dart with $c_x < 1$ has $a_x > 1$ ($f_D < 0$ gives $a_x - 1 > (1 - c_x)a_y/|c_y|$).

**Lemma 9.6 (oblique cut of a long arm: the FARL chart; `charts.farl_problem`, validity `fam5.farl_setup`).** $BD$ chart,
$B = (0, 0)$, $D = (1, 0)$. Let the arm apex be $A = a\,(u, 1)$ with $u = \cot b_A \ge 0$ and $a \ge R > 0$, and let $C$ lie in the open lower half
plane. Fix rationals $L > 0$ and $k_0 > 0$ and put
$$\begin{aligned} Q_1 &:= L\,(u, 1)\quad \text{(on the ray } BA\text{)},\\ Q_2 &:= D + (L - k_0u)(u, 1)\quad \text{(on the parallel to } BA \text{ through } D\text{)},\\ \ell &:= \text{the line } Q_1Q_2,\\ P &:= \operatorname{int}\bigl(CBD \cup \operatorname{tri}(B, D, Q_2) \cup \operatorname{tri}(B, Q_2, Q_1) \cup \text{common edges}\bigr)\\ &\qquad \text{(the chart domain: pieces } C,\ \mathrm{S1} = BDQ_2,\ \mathrm{S2} = BQ_2Q_1\text{)}. \end{aligned}$$
Assume (exactly the conditions verified by `farl_setup`, in ball arithmetic, for every $u$ of the box):

(H1) $R > L > k_0u$ (so $0 < L - k_0u$, i.e. $Q_2$ lies beyond $D$ on its ray, and $Q_1$ lies strictly between $B$ and $A$ on $BA$);

(H2) $B$, $D$ and the closed triangle $CBD$ lie strictly on one side of $\ell$ (checked at $B$, $D$ and the corners of the Cartesian box that
contains $C$; the triangle is the convex hull of $B$, $D$, $C$);

(H3) the point $R(u, 1)$ lies strictly on the other side of $\ell$.

Then
$$\lambda_1(Q) \ge \min\bigl(\mu_1(P; \text{Neumann on } Q_1Q_2, \text{Dirichlet elsewhere}), \pi^2/|Q_1Q_2|^2\bigr),\qquad |Q_1Q_2|^2 = (1 - k_0u^2)^2 + (k_0u)^2,$$
and $\lambda_2(Q) \le \lambda_2(Q_R)$ with $Q_R$ the dart with $A$ replaced by $R(u, 1)$.

*Proof.* (a) The ray $B + s(u, 1)$, $s \ge 0$, meets $\ell$ only at $Q_1$ ($s = L$); by (H3) the points with $s > L$, in particular $A$ ($s = a \ge R > L$),
lie on the far side of $\ell$, while $B$ lies on the near side (H2). Hence $\ell$ crosses $BA$ at $Q_1$ and, as $D$ is on the near side and $A$ on the
far side, crosses $DA$ at one point $X$; it does not meet $BD$ (H2). The closed strip between the line $BA$ and its parallel through $D$
contains the segment $DA$, and $\ell$ meets this strip in the segment $[Q_1, Q_2]$; so $X \in [Q_1, Q_2]$.

\(b) The tip $T := ABD \cap \{\text{far side of } \ell\}$ is the triangle $A\,Q_1X$; its sides $AQ_1 \subset BA$ and $AX \subset DA$ lie on $\partial Q$ (the triangle $CBD$ is on the
near side by (H2), so it does not meet $T$), and its third side is $\sigma := [Q_1, X]$. By Lemma 6.1 the mixed problem on $T$ (Neumann on $\sigma$)
has $\mu_1 \ge \pi^2/|\sigma|^2 \ge \pi^2/|Q_1Q_2|^2$, since $\sigma \subset [Q_1, Q_2]$; and $|Q_1Q_2|^2 = |(1 + (L - k_0u - L)u, -k_0u)|^2 = (1 - k_0u^2)^2 + (k_0u)^2$.

\(c) The near part $Q' := Q \cap \{\text{near side of } \ell\} = \operatorname{int}(CBD \cup \text{polygon } B\,D\,X\,Q_1 \cup \dots)$ is contained in $P$: the polygon $BDXQ_1$ lies in the
convex quadrilateral $B\,D\,Q_2\,Q_1$ ($X \in [Q_1, Q_2]$; $BQ_1$ and $DQ_2$ are parallel, of positive lengths $L\sqrt{1 + u^2}$ and $(L - k_0u)\sqrt{1 + u^2}$), which is the union of
the pieces S1 and S2. The boundary of $Q'$ consists of $\sigma \subset Q_1Q_2$ (the Neumann edge of $P$) and of parts of $\partial Q$; the part of $\partial Q'$ inside
$P$ (the segment $DX$ of $DA$) carries the Dirichlet condition. Extension by zero maps the form domain of the mixed problem on $Q'$
(functions vanishing on $\partial Q' \setminus \sigma$) into that of the mixed problem on $P$, so $\mu_1(Q'; \text{Neumann on } \sigma) \ge \mu_1(P; \text{Neumann on } Q_1Q_2)$ (1.2).

\(d) Neumann bracketing along $\sigma$ (1.2): $\lambda_1(Q) \ge \min(\mu_1(Q'; \text{N on } \sigma), \mu_1(T; \text{N on } \sigma))$, and (b), (c) give the claim.

\(e) $R(u, 1)$ lies on the segment $BA$ ($a \ge R$), so $\operatorname{tri}(B, D, R(u, 1)) \subset ABD$, $Q_R \subset Q$ and $\lambda_2(Q) \le \lambda_2(Q_R)$ (1.1). $\square$

Use. The FEM certificate of a FARL box bounds $\mu_1(P; \text{Neumann on } Q_1Q_2)$ from below over $u$ in the $u$-interval and $C$ in the Cartesian box
(CR chart: pieces $C$, S1, S2 with the Neumann cut), uses the tail constant $\pi^2/\max|Q_1Q_2|^2$ (`charts.cut_tail`: the maximum of the
convex quadratic $g(s) = (1 - k_0s)^2 + k_0^2s$ over the range of $s = u^2$, rounded upward, then $\pi^2/\cdot$ rounded downward), and bounds
$\lambda_2(Q_R)$ from above (P2 chart: the $BD$ chart with $A = R(u, 1)$). The $u$-interval is $[\cot b_2, \cot b_1]$ of the arm's angle interval
$[b_1, b_2] \subset (0, \pi/2]$ (directed: the lower end rounded down, the upper end rounded up), $R$ is a dyadic lower bound of
$\rho_{\min}\sin b_1$ (so $a = \rho\sin b \ge R$ for every $\rho \ge \rho_{\min}$), and the chart box $(u, C)$ is an exact-endpoint float box containing these
ranges (`nanguard.enclose`). For PP boxes with one long arm ($\rho \ge 4$ on the box), $\rho_{\min}$ is that arm's smallest radius on the box
(the superset family with $\rho \ge \rho_{\min}$); for FP/PF, $\rho_{\min} = \mathrm{RF} = 20$; a long $C$ arm is mirrored first (reflection in $BD$).

**Lemma 9.5 (very long steep arm; closes part of the long-arm face and of the deep notch $\beta \to 1$).** $BD$ chart, $|BD| = 1$,
$L := |BA|$, $b := b_A \in (0, \pi)$. For every $\ell \in (0, L - 1)$:
$$\xi(Q) \le (L + 1)^2 / \bigl(\sin^2 b\cdot(L - 1 - \ell)^2\bigr) + 4/\ell^2.$$

*Proof.* $\lambda_1(Q) \ge \pi^2$ (Lemma 4.1). Use coordinates $s$ along $e_A = (A - B)/L$ and $\tau$ along the unit normal pointing to $D$'s side, so
$B = (0,0)$, $A = (L,0)$, $D = (\cos b, \sin b)$. Put $s_0 := \max(\cos b, 0) \le 1$. For $s_0 \le s \le L$ the section of the triangle $ABD$ at abscissa $s$
is $\{0 \le \tau \le \tau_{DA}(s)\}$ with $\tau_{DA}(s) = \sin b\,(L - s)/(L - \cos b)$: the edge $BA$ is $\tau = 0$, the edge $BD$ lies in $\{s \le s_0\}$, and the third
edge is $DA$. $\tau_{DA}$ is decreasing, so the rectangle $R := (s_0, s_0 + \ell) \times (0, w)$, $w := \sin b\,(L - s_0 - \ell)/(L - \cos b)$, lies in
$ABD \subset Q$, and $w \ge \sin b\,(L - 1 - \ell)/(L + 1)$. By 1.1, $\lambda_2(Q) \le \lambda_2(R) \le \pi^2/w^2 + 4\pi^2/\ell^2$ (the (1,2) mode
$\sin(\pi\tau/w)\sin(2\pi(s - s_0)/\ell)$ and the (1,1) mode span a 2-dimensional space). Divide. $\square$

The right-hand side decreases in $L$ for fixed $\ell$ and depends on $b$ only through $\sin b$. `lemma_longarm.py` (arb, 128 bits; log
`lemma_longarm.log`) certifies, for each row, a value $\ell_k$ with $\mathrm{RHS}(L_k, \ell_k, b_k) < 5/2$, hence:
**$\xi(Q) < 5/2$ whenever $b_A \in [b_k, 180^\circ - b_k]$ and $|BA| \ge L_k\cdot|BD|$** (and the same with $A \leftrightarrow C$), for
$(b_k, L_k) = (40^\circ, 1318)$, $(41^\circ, 418)$, $(42^\circ, 230)$, $(43^\circ, 155)$, $(45^\circ, 92)$, $(47.5^\circ, 61)$, $(50^\circ, 46)$, $(55^\circ, 31)$, $(60^\circ, 25)$, $(65^\circ, 21)$,
$(70^\circ, 19)$, $(80^\circ, 17)$, $(90^\circ, 16)$. The bound can never reach $5/2$ when $\sin^2 b_A \le 2/5$ ($b_A \le 39.23^\circ$ or $b_A \ge 140.77^\circ$).

*Deep notch.* In the hull placement $D - B = (1 - \beta)(F - B)$ ($F$ = foot of the cevian $BD$ on $AC$), so $|BD| = (1 - \beta)|BF|$,
$b_A = \angle ABF$, $b_C = \angle FBC$. Hence, for fixed $T$ and $F$, the limit $\beta \to 1$ is closed by Lemma 9.5 as soon as one of $\angle ABF$, $\angle FBC$ lies in
$[40^\circ, 140^\circ]$ — for instance whenever $80^\circ \le \theta_B \le 140^\circ$ — namely for $1 - \beta \le |BX|/(L_k|BF|)$ ($X$ the steep arm). It is NOT closed when
both $\angle ABF$, $\angle FBC$ are $< 40^\circ$ (e.g. $T$ equilateral with $F$ the midpoint: $30^\circ/30^\circ$) or one of them exceeds $140^\circ$: there the chord bound
$\lambda_1 \ge \pi^2/|BD|^2$ is too weak (the true $\lambda_1$ is near the waveguide threshold $\pi^2/(|BD|^2\max(\sin^2 b_A, \sin^2 b_C))$), and a numerical
certificate of the bent-waveguide core (Lemma 9.3; `charts.deep_problem`) would be needed.

## 10. The general certificate (gcert.py) — used for all runs

Charts are unions of triangular pieces $q$, each the affine image of an exact dyadic reference mesh (`chart.py`); the vertex
positions are polynomials in the box parameters (`charts.py`), so per piece $G_q = N_q/d_q$ with $N_q$, $d_q$ exact polynomials
($d_q = |\det E_q|$, sign checked on the box), $J_q = d_q/|\det R_q|$. Steps 1–7 of §8 hold verbatim with: the CR chart may carry Neumann
edges (cuts; Liu's bound holds for mixed problems, 1.7) and its own geometry; the P2 chart (Dirichlet everywhere) may be a different
domain contained in every $Q$ of the box (Lemmas 9.1, 9.4: $CBD$, or the dart $Q_{X_0}$); extra Neumann-bracketing constants $c$ (tails)
enter as $\lambda_1 \ge \min(L_1(p), c)$, so the Schur test is run for each of $L_1$ and the $c$'s; the residual is scaled on each DOF group by the
product of the $d_q$ of the pieces containing it, which makes all stiffness coefficients exact polynomials; $\sigma_2$ comes from ONE verified
$LDL^{\mathsf T}$ per box of the piecewise Loewner minorant $\sum_q \mu_q K_q(p_c) - s\sum_q j_q M_q(p_c)$ (`gsigma.sigma_box`): for every $v$,
$v^{\mathsf T}K_q(p)v \ge \mu_q\,v^{\mathsf T}K_q(p_c)v$ with $\mu_q = \min_{\mathrm{box}} (J_q(p)/J_q(p_c))\,\lambda_{\min}(TT^{\mathsf T})$, $T = E_q(p_c)E_q(p)^{-1}$ (LEMMAS 1.10), and
$M_q(p) \le j_q M_q(p_c)$, $j_q = \max_{\mathrm{box}} J_q(p)/J_q(p_c)$; hence $K(p) - sM(p) \ge A_{\mathrm{low}}$ in the Loewner order and $\lambda_2(K(p) - sM(p)) \ge \lambda_2(A_{\mathrm{low}})$.

## 11. One-dimensional (Born–Oppenheimer) bounds for families of darts (`bo1d.py`, `bo2.py`, `p2fix.py`, `cover5.py`)

**Lemma 11.1 (1-D reduction with the longest slice component).** Let $\Omega \subset \mathbb R^2$ be open (bounded or not), $u$ a unit vector,
$x = \langle\cdot, u\rangle$, and $\ell(x)$ the length of the LONGEST connected component of the slice $\Omega \cap \{\langle\cdot,u\rangle = x\}$ ($\ell(x) = 0$ if the slice is
empty). Let $V$ be measurable with $0 \le V(x) \le \pi^2/\ell(x)^2$ wherever the slice is non-empty. Then for every $\varphi \in C_c^\infty(\Omega)$
$$\int_\Omega |\nabla\varphi|^2 \ge \int (v'^2 + V v^2)\,dx,\qquad v(x) := \Bigl(\int \varphi(x, y)^2\,dy\Bigr)^{1/2}\qquad (y = \langle\cdot, u^\perp\rangle),$$
and $v \in H^1$ with compact support in the projection of $\Omega$. Consequently $\lambda_1(\Omega) \ge \inf\sigma(-d^2/dx^2 + V)$ on that projection
(Dirichlet at finite ends).

*Proof.* Fix $x$. The slice is a union of disjoint open intervals $J$ (its components); $\varphi(x, \cdot)$ has compact support in the slice,
so on each $J$ it vanishes near both ends and the 1-D Poincaré inequality gives $\int_J |\partial_y\varphi|^2 \ge (\pi/|J|)^2 \int_J \varphi^2 \ge (\pi/\ell(x))^2 \int_J \varphi^2$;
summing over $J$, $\int |\partial_y\varphi|^2\,dy \ge V(x)\,v(x)^2$. For the $x$-derivative: $v^2$ is smooth, $v$ is smooth where positive, and
$d(v^2)/dx = 2\int \varphi\,\varphi_x\,dy$ (differentiation under the integral is legitimate because $\varphi$ has compact support in $\Omega$, so no
boundary terms arise), hence $|v'| \le (\int \varphi_x^2\,dy)^{1/2}$ by Cauchy–Schwarz where $v > 0$, and $v' = 0$ a.e. on $\{v = 0\}$. Integrate
$\int|\nabla\varphi|^2 = \iint (\varphi_x^2 + \varphi_y^2) \ge \int (v'^2 + V v^2)$. Density of $C_c^\infty(\Omega)$ in $H^1_0(\Omega)$ gives the eigenvalue statement. $\square$

**Lemma 11.2 (geometry of an arm).** Let $X$ be a point in the open upper half plane (lower for $C$), $b = \angle XBD \in (0, \pi)$,
$e(b) = (\cos b, \pm\sin b)$ and $T_X = \operatorname{tri}(B, D, X)$, $B = (0, 0)$, $D = (1, 0)$.

\(a) $T_X \subset \bar S(b) := \{B + s\,e(b) + \tau(D - B) : s \ge 0,\ 0 \le \tau \le 1\}$ (the half-strip spanned by $e(b)$ and $BD$; its width is $w(b) = \sin b$).

\(b) A line $L = \{\langle\cdot, u\rangle = x\}$ ($u$ a unit vector) not meeting the segment $BD$ meets $T_X$ in a segment (or not at all) whose
endpoints lie on $BX$ and $DX$; if $\langle e(b), u\rangle \ne 0$, its length is $\le w(b)/|\langle e(b), u\rangle|$.

\(c) The angle $\angle XDB = \operatorname{atan2}(\rho\sin b, 1 - \rho\cos b)$ ($X = \rho\,e(b)$) is strictly increasing in $\rho$.

\(d) If $X'$ satisfies $\angle X'BD \le \angle XBD$ and $\angle X'DB \le \angle XDB$ (both on the same side of $BD$), then $T_{X'} \subset T_X$.

*Proof.* (a) $T_X = \operatorname{conv}\{B, D, X\}$ and $X = B + \rho\,e(b)$, so every point is $B + \lambda\rho\,e(b) + \mu(D - B)$ with $\lambda, \mu \ge 0$, $\lambda + \mu \le 1$:
$s = \lambda\rho \ge 0$, $\tau = \mu \in [0, 1]$. (b) $L \cap T_X$ is convex, and its endpoints lie on $\partial T_X = BD \cup BX \cup DX$; $L$ misses $BD$, so they lie on
$BX$ and $DX$, hence in $\bar S(b)$ and in the full strip between the parallel lines $B + \mathbb R e(b)$ and $D + \mathbb R e(b)$, whose width (the distance
of $D$ from the line $B + \mathbb R e(b)$) is $|e(b) \times (D - B)| = \sin b$. The line $L$ has direction $u^\perp$; a line with unit direction $d$ crosses a
strip of width $w$ and unit normal $n$ in a segment of length $w/|\langle d, n\rangle|$; with $d = u^\perp$ and $n = e(b)^\perp$ (rotation by $90^\circ$ preserves
inner products) $\langle d, n\rangle = \langle u, e(b)\rangle$. (c) $d/d\rho\,[\rho\sin b/(1 - \rho\cos b)] = \sin b/(1 - \rho\cos b)^2 > 0$ on each branch, and the
angle is continuous through $\rho\cos b = 1$ (where it equals $\pi/2$). (d) A point $P$ above $BD$ lies in $T_X$ iff $\angle PBD \le \angle XBD$ and
$\angle PDB \le \angle XDB$; apply this to $P = X'$ and use convexity. $\square$

**Lemma 11.3 (the BO2 box certificate).** Let $\mathcal F$ be a family of darts given by arm specifications (`bo2.Arm`):
'polar' ($X = \rho\,e(b)$, $b \in [b_1, b_2]$, $\rho \in [r_1, r_2]$) or 'far' ($b \in [b_1, b_2]$, $\rho \ge R \ge 0$), one for $A$ and one for $C$. Fix a unit $u$.

\(i) Compact part: the union over $\mathcal F$ of $T_X$ is contained in $\operatorname{conv}(P)$ with
$$P = \{B, D, r_1e(b_1), r_1e(b_2), r_2e(b_1), r_2e(b_2), (r_2/\cos((b_2-b_1)/2))\,e((b_1+b_2)/2)\}$$
(the annular sector lies in the convex hull of its four corners and of the intersection
point of the tangents at the ends of its outer arc).

\(ii) Far arms: if $\langle e(b), u\rangle > 0$ on $[b_1, b_2]$, then for any $x_{\mathrm{cut}} \ge \max(\langle B,u\rangle, \langle D,u\rangle)$ the union over $\mathcal F$ of
$T_X \cap \{\langle\cdot,u\rangle \le x_{\mathrm{cut}}\}$ lies in $\operatorname{conv}\{B, D$, the four points where the rays
$B + \mathbb R_+e(b_i)$, $D + \mathbb R_+e(b_i)$ ($i = 1, 2$) meet $\{\langle\cdot,u\rangle = x_{\mathrm{cut}}\}\}$ (by (a): $\bar S(b) \cap \{x \le x_{\mathrm{cut}}\}$ is the quadrilateral with these
vertices for $b = b_i$, and the crossing points move monotonically along the cut line as $b$ varies); and for $x > x_{\mathrm{cut}}$ the
slice of each dart consists only of far-arm pieces, each of length $\le c := \max_{b\in[b_1,b_2]} \sin b/\langle e(b), u\rangle$ by (b).

\(iii) For a slab $S_k = \{x_k \le \langle\cdot,u\rangle \le x_{k+1}\}$, let $\mathrm{ext}_X(S_k)$ be the $y$-extent of $\operatorname{conv}(P_X) \cap S_k$, computed from the points of
$P_X$ inside $S_k$ and from the intersections of ALL segments between points of $P_X$ with the two slab lines (a superset of the
vertices of $\operatorname{conv}(P_X) \cap S_k$), and, if $S_k$ misses the projection of $BD$, also $\le c_X$ by (b). If $S_k$ meets the projection of $BD$
let $H_k = \mathrm{ext}_A + \mathrm{ext}_C$, otherwise $H_k = \max(\mathrm{ext}_A, \mathrm{ext}_C)$ (slices not meeting $BD$ split into an $A$-part and a $C$-part lying on
different sides of the line $BD$, which are different components). Then every slice of every dart of $\mathcal F$ with $x \in S_k$ has all
components of length $\le H_k$.

\(iv) With $V = \pi^2/H_k^2$ on $S_k$ and $V = \pi^2/c^2$ beyond a cut on a side towards which some far arm extends ($V_\infty = \pi^2/(\max c)^2$),
Lemma 11.1 gives $\lambda_1(Q) \ge \mu$ for every $Q \in \mathcal F$, where $\mu$ is the bottom of the spectrum of $-d^2/dx^2 + V$. $\mu \ge E$ is certified by
Sturm shooting in ball arithmetic (`bo2._no_zero_arb`): the solution with $v(x_L) = 0$, $v'(x_L) = 1$ (or, at an unbounded left
end, $v = 1$, $v' = \sqrt{V_L - E}$: the decaying exponential) stays positive across all cells (hyperbolic cells: positivity at
both ends; oscillatory cells: phase condition), and at an unbounded right end $v' + \sqrt{V_R - E}\,v > 0$ (then $v$ stays positive
for all larger $x$). A positive solution on the whole line/half-line implies $E \le \inf\sigma$ (Allegretto–Piepenbrink), since the
operator is a limit of Dirichlet problems on bounded intervals for which Sturm comparison applies.

\(v) Inner dart: by Lemma 11.2(c)–(d), $X_{\mathrm{in}}$ with $\angle X_{\mathrm{in}}BD \le b_1$ and $\angle X_{\mathrm{in}}DB \le \min_{\mathcal F}\angle XDB$ (the minimum over $b$ at the smallest
radius, enclosed by interval arithmetic over subintervals) satisfies $T_{X_{\mathrm{in}}} \subset T_X$ for all $X$ of the spec; for spikes
($A = a(1, t)$, $a \ge 1$, $t \ge t_1$) the point $(1, t_1)$ has this property (Lemma 9.4). $Q_{\mathrm{in}} = \operatorname{int}(T_{A_{\mathrm{in}}} \cup T_{C_{\mathrm{in}}} \cup BD)$ (or $T_{C_{\mathrm{in}}}$
alone) is contained in every $Q \in \mathcal F$, so $\lambda_2(Q) \le \lambda_2(Q_{\mathrm{in}})$ (domain monotonicity), and `p2fix.p2_cert` bounds $\lambda_2(Q_{\mathrm{in}})$ by conforming
P2 Rayleigh–Ritz with rigorous $2\times2$ forms (§8 step 7 with zero box radius).

**Conclusion:** if the Schur test certifies $\lambda_2(Q_{\mathrm{in}}) < (5/2)\cdot E$, then $\xi(Q) < 5/2$ for every $Q \in \mathcal F$. $\square$

Remarks. (1) The family may contain non-darts (convex quadrilaterals); the certificate holds for them too and is harmless.
(2) Far arms make the certificate uniform in the arm length: the far family $\{\rho \ge R\}$ is one box in $\rho$. (3) Validation:
`test_bo1d.py` checks the slab bounds of (iii) against brute-force slices of sampled darts (35\ 370 samples, no violation).
(4) Arm replacement: a polar arm spec ($b \in [b_1, b_2]$, $\rho \in [r_1, r_2]$) may be replaced in (i)–(iv) by the far spec ($b \in [b_1, b_2]$,
$\rho \ge r_1$): the far family contains the polar one, so its lower bound $E$ is a lower bound for the polar family as well (the inner
dart of (v) is still taken from the polar family). BO2 records with key `arms` $\in \{$farA, farC, farAC$\}$ use this replacement for the
named arm(s) (`bo2.arm_variant`); the cover tries it for arms with $r_2 \ge 4$, where the chord and tail bounds of the far arm are
sharper than the slab extents of the hull. (5) Any P2 mesh level $k$ may be used in (v) (the record stores $k$; the cover tries
$k = 4, 5, 6$, since the P2 bound of long inner darts converges slowly in $k$); likewise any slab count $K$ and slicing direction $\theta$
(stored in the record).

**Lemma 11.4 (the INC box certificate: outer and inner darts).** Let $\mathcal B$ be a polar $\times$ polar box of darts: $b_A \in [b_1^A, b_2^A]$,
$\rho_A \in [r_1^A, r_2^A]$, $b_C \in [b_1^C, b_2^C]$, $\rho_C \in [r_1^C, r_2^C]$, all angles in $(0, \pi)$. For $X = \rho\,e(b)$ write
$d(X) = \angle XDB = \operatorname{atan2}(\rho\sin b, 1 - \rho\cos b)$ (for $C$ the mirror image; the formulas are the same).

\(a) $\partial d/\partial\rho = \sin b/|X - D|^2 > 0$ and $\partial d/\partial b = \rho(\cos b - \rho)/|X - D|^2$. Hence over a polar arm box, $\max d = \max_{b \in [b_1, b_2]} d(b, r_2)$,
attained at $b_1$, at $b_2$ or at $b^* = \arccos r_2$ (only if $r_2 < 1$ and $b^* \in [b_1, b_2]$), and $\min d = \min(d(b_1, r_1), d(b_2, r_1))$.

\(b) Outer apex: if $X^+$ lies in the open upper half plane with $\angle X^+BD \ge b_2$, $\angle X^+DB \ge \max d$ and $\angle X^+BD + \angle X^+DB < \pi$, then
$T_X \subset T_{X^+}$ for every $X$ of the arm box. Inner apex: if $\angle X^-BD \le b_1$ and $\angle X^-DB \le \min d$, then $T_{X^-} \subset T_X$ for every $X$.

\(c) Consequently $Q^- := \operatorname{int}(T_{A^-} \cup T_{C^-} \cup BD) \subset Q \subset Q^+ := \operatorname{int}(T_{A^+} \cup T_{C^+} \cup BD)$ for every $Q$ of $\mathcal B$, and
$\lambda_1(Q) \ge \lambda_1(Q^+)$, $\lambda_2(Q) \le \lambda_2(Q^-)$ (Dirichlet domain monotonicity).

\(d) If $L_1 \le \lambda_1(Q^+)$ is certified by the general certificate (§10: CR + Temple + Liu, $\sigma_2$ from the verified $LDL^{\mathsf T}$) on a Cartesian
$BD$-chart box containing the point $Q^+$ (in practice the box $Q^+ \pm 2^{-30}$), $E$ is a float $\le$ the exact lower end of the $L_1$ enclosure,
and `p2fix.p2_cert` certifies $\lambda_2(Q^-) < (5/2)E$ (§8 step 7 at zero box radius), then $\xi(Q) < 5/2$ for every $Q$ of $\mathcal B$.

*Proof.* (a) With $x = 1 - \rho\cos b$, $y = \rho\sin b$, $d(\operatorname{atan2}(y, x)) = (x\,dy - y\,dx)/(x^2 + y^2)$. For $\rho$: $dy = \sin b$, $dx = -\cos b$, so
$x\,dy - y\,dx = (1 - \rho\cos b)\sin b + \rho\sin b\cos b = \sin b > 0$. For $b$: $dy = \rho\cos b$, $dx = \rho\sin b$, so
$x\,dy - y\,dx = \rho\cos b - \rho^2 = \rho(\cos b - \rho)$. For fixed $\rho$ the sign of $\partial d/\partial b$ is that of $\cos b - \rho$: if $\rho \ge 1$, $d$ is decreasing on $(0, \pi)$; if $\rho < 1$
it increases up to $b^* = \arccos\rho$ and decreases after. So on $[b_1, b_2]$ the maximum of $d(\cdot, \rho)$ is at $b_1$, $b_2$ or $b^*$, and the minimum
(of a function increasing then decreasing) at $b_1$ or $b_2$. Since $d$ increases in $\rho$, the box maximum is attained at $\rho = r_2$ and the box
minimum at $\rho = r_1$. (b) The closed triangle $T_{X^+}$ on $BD$ is the set of points $P$ of the closed upper half plane with
$\angle PBD \le \angle X^+BD$ and $\angle PDB \le \angle X^+DB$ (intersection of its angular sectors at $B$ and at $D$). Every $X$ of the arm box has $\angle XBD = b \le b_2$ and
$\angle XDB \le \max d$, so $X \in T_{X^+}$, and $T_X = \operatorname{conv}\{B, D, X\} \subset T_{X^+}$ by convexity. The inner statement is Lemma 11.2(d).
(c) Unions and interiors are monotone; Dirichlet eigenvalues decrease under domain inclusion. (d)
$\lambda_2(Q) \le \lambda_2(Q^-) < (5/2)E \le (5/2)L_1 \le (5/2)\lambda_1(Q^+) \le (5/2)\lambda_1(Q)$. $\square$

Implementation (`inc.py`): `d_range` evaluates (a) in ball arithmetic (if $r_2$'s enclosure straddles 1, the maximum is enclosed by
evaluating $d$ over the whole interval $[b_1, b_2]$); the apexes are floats constructed from the target angles with a small outward
(outer) or inward (inner) pad and re-checked in ball arithmetic (`check_outer`, `check_inner`); apex radius $\le 8$. The replay
(`replay3.py`) re-checks the inclusion for the recorded apexes, recomputes $L_1$ at the recorded $Q^+$ (must be $\ge$ the recorded $E$) and
re-runs the P2 test of $Q^-$ with the recorded $E$. Use: family PP, where the 1-D bound of Lemma 11.3 is about $3\times$ off for 'S-shaped'
short darts ($A$ short and near $D$, $C$ beyond $D$) and the box FEM certificate needs very small boxes (its strain factors $\mu_q$ penalise
the whole stiffness of a thin piece); the INC loss depends only on the spread of the four apex angles over the box.

## 12. A lemma for darts that are thin along $BD$

**Lemma 12.1 (thin stem).** $BD$ chart ($B = (0,0)$, $D = (1,0)$). For $0 \le x \le 1$ let $h_A(x)$, $h_C(x)$ be the lengths of the vertical
slices of the triangles $ABD$, $CBD$ at abscissa $x$, and $c^* := \max_{0\le x\le 1} (h_A(x) + h_C(x))$. Then for every $\theta \in (0, 1)$
$$\xi(Q) \le (1 + 4c^{*2})/(1 - \theta)^2 + 16\,c^{*2}/\theta^2;\qquad \text{in particular } \xi(Q) \le 2.4356 < 5/2 \text{ if } c^* \le 0.05\ (\theta = 1/4).$$

*Proof.* (i) Chords. For $0 < x < 1$ the vertical slice of $Q$ is the single interval $(-h_C(x), h_A(x))$ (it crosses $BD$). For $x > 1$
the slice consists of the part in $ABD$ (between $DA$ and $BA$, non-empty only if $a_x > 1$, of length decreasing linearly from $h_A(1)$
to 0) and the part in $CBD$ (likewise $\le h_C(1)$), which are different components (the ray beyond $D$ is not in $Q$); so every component
has length $\le \max(h_A(1), h_C(1)) \le c^*$. For $x < 0$ likewise $\le \max(h_A(0), h_C(0)) \le c^*$. Every vertical chord ends on $\partial Q$, so
$\lambda_1(Q) \ge \pi^2/c^{*2}$ by 1.3.

\(ii) Parallelogram. $h_A$ is concave on $[0, 1]$ (the upper boundary of $ABD$ over $[0, 1]$ is the lower envelope of the lines $BA$ and $DA$,
or one of them), likewise $h_C$; let $x^*$ maximise $c = h_A + h_C$ on $[0, 1]$ and let $e \in \{0, 1\}$ be the endpoint farther from $x^*$
($|e - x^*| \ge 1/2$). For $x = \lambda x^* + (1 - \lambda)e$, $\lambda \in [1 - \theta, 1]$, concavity gives $h_A(x) \ge \lambda h_A(x^*) + (1-\lambda)h_A(e)$ and the same for $h_C$,
so the slice of $Q$ at $x$ contains the interval centred at $m(x) = [\lambda(h_A - h_C)(x^*) + (1 - \lambda)(h_A - h_C)(e)]/2$ (affine in $x$) of length
$\lambda c^* + (1 - \lambda)c(e) \ge (1 - \theta)c^*$. Hence $Q$ contains the parallelogram $\{x \text{ between } x^* \text{ and } x^* + \theta(e - x^*),\ |y - m(x)| < (1 - \theta)c^*/2\}$ of
horizontal length $L = \theta|e - x^*| \ge \theta/2$, vertical side $(1 - \theta)c^*$ and shear $|m'| \le (|h_A - h_C|(x^*) + |h_A - h_C|(e))/(2|e - x^*|) \le 2c^*$.
By 1.4 and 1.1, $\lambda_2(Q) \le \pi^2(1 + 4c^{*2})/((1 - \theta)c^*)^2 + 16\pi^2/\theta^2$. Divide by (i). For $c^* \le 0.05$, $\theta = 1/4$:
$(1.01)/(0.5625) + 16\cdot 0.0025/0.0625 = 1.7956 + 0.64 = 2.4356$. $\square$

*Use.* Any box whose darts all have $c^* \le 0.05$ is discarded ('Lemma 12.1'); $c^*$ is bounded on a box by the vertical-slicing slab
extents of Lemma 11.3(iii) over $[0, 1]$ (rigorous upper bound of $h_A + h_C$). In particular every dart with $\tan b_A + \rho_C \le 0.05$
($h_A \le x\tan b_A$, $h_C \le |c_y| \le \rho_C$), e.g. the corner $\{\tan b_A \le 0.02,\ |C| \le e^{-4} = 0.0183\}$ left by the spike/collar grids.

## 13. The global decomposition (every dart is covered)

Notation: $BD$ chart, $b_A = \angle ABD$, $b_C = \angle CBD \in (0, \pi)$, $\theta_B = b_A + b_C$, $\rho_A = |BA|$, $\rho_C = |BC|$ ($|BD| = 1$), $t = \tan b_A$,
$A = \rho_A\,e(b_A)$, $C = \rho_C\,\bar e(b_C)$. $\beta$ is the notch depth.

**Lemma 13.1 (dart condition in polar form).** The line $AC$ meets the $x$-axis at abscissa $f = \sin\theta_B/(\sin b_A/\rho_C + \sin b_C/\rho_A)$,
and $Q$ is a dart with reflex vertex $D$ iff $\theta_B < \pi$ and $f > 1$; then $1 - \beta = 1/f$, i.e.
$$1 - \beta = (\sin b_C/\rho_A + \sin b_A/\rho_C)/\sin\theta_B.$$

*Proof.* $f = (a_x|c_y| + c_x a_y)/(a_y + |c_y|)$ (the zero of the affine interpolation of $y$ along $AC$); insert the polar
coordinates: numerator $\rho_A\rho_C\sin(b_A + b_C)$, denominator $\rho_A\sin b_A + \rho_C\sin b_C$. $D = (1, 0)$ is interior to the hull iff it
lies strictly between $B$ and the point $F = (f, 0)$ of $AC$ on the ray $BD$, i.e. $f > 1$; and $D = (1 - \beta)F$ by similar triangles. $\square$

Consequences used as discards: (a) if $\sin b_C/\rho_A + \sin b_A/\rho_C - 0.9661\sin\theta_B \ge 0$ on a box, every dart of the box has
$\beta \le 0.0339$ (Prop 5.1) or is not a dart (for $b_A, b_C \in [0, \pi]$ the left side is non-increasing in $\rho_A$ and in $\rho_C$, so on a polar box
it suffices to check it at the largest radii, over a subdivision of the two angle intervals: `fam5.PP` predicate; for the
spike / collar families, where $\rho_A$ is unbounded, the term $\sin b_C/\rho_A \ge 0$ is dropped and $\sin b_A/\rho_C - 0.9661\sin\theta_B \ge 0$ is
checked at the largest $\rho_C$: `fam5.SP` predicate); (b) a dart has $\rho_A > \sin b_C$ and $\rho_C > \sin b_A$ (both terms are $< \sin\theta_B \le 1$).

**Theorem 13.2 (cover).** Every dart $Q$ satisfies $\xi(Q) < 5/2$ provided the runs listed in README §2 are complete with FAIL 0
(each verified by its fail-closed replay). *Proof.* By the reflection in $BD$ we may assume $b_A \le b_C$ (the reflection swaps the
roles of $A$ and $C$ and preserves $\xi$); then $b_A < \pi/2$ (as $b_A + b_C < \pi$). Every such dart is in one of the following sets, each
covered by a lemma or by a run (runs cover their parameter grids exactly; their discards are re-checked predicates naming one of
these lemmas or one of the other sets):

(1) $\theta_B \le 0.1772$: Lemma 5.3.
(2) $\theta_B \ge 174.29^\circ$: the hull angles at $A$ and $C$ sum to $\le 5.71^\circ$, Lemma 5.2.
(3) $\beta \le 0.0339$: Prop 5.1.
(4) max vertical stem chord $\le 0.05$: Lemma 12.1.
(5) $t \le t_1 := \mathrm{float}(0.02) = 0.0200000000000000004163\ldots$ ('collar', every $\rho_A$; covers $t < 0.02$): $C$ polar with
    $\rho_C \in [e^{-4}, e^3]$: run SP_col (BO2 and FEM collar chart, Lemmas 9.1–9.2, 11.3)
    and, for $C$ in the sector $[38^\circ, 100^\circ] \times [0.5, 1.8]$, the FEM run certs2/kf_collar; $\rho_C \ge 20$: run SF_col; $\rho_C < e^{-4}$:
    $c^* \le t + \rho_C < 0.02 + 0.0184 < 0.05$, Lemma 12.1. (The polar and far grids overlap: $e^3 > 20$.)
(6) $t_1 \le t \le 3/10$ and $a_x \ge 1$ ('spike', every $\rho_A$): $\rho_C \in [e^{-4}, e^3]$: run SP_spk (BO2 and FEM spike chart, Lemma 9.4) and
    certs2/kf_spike for $C$ in the KF sector; $\rho_C \ge 20$: run SF_spk; $\rho_C < e^{-4}$: no dart (Lemma 13.1(b): $\rho_C > \sin b_A \ge 0.0199$).
(7) all other darts with $b_A \le b_C$ ($t \ge 0.3$, or $0.02 \le t \le 0.3$ with $a_x < 1$): $\rho_A < e^{-2.5}$: then $\sin b_C < \rho_A < 0.083$ by
    13.1(b), so $b_C < 4.8^\circ$ (and $\theta_B \le 2b_C < 0.1772$: set (1)) or $b_C > 175.2^\circ$ (set (2)); $\rho_C < e^{-4}$: no dart ($\rho_C > \sin b_A \ge 0.0199$);
    $\rho_A \in [e^{-2.5}, e^3]$ and $\rho_C \in [e^{-4}, e^3]$: run PP (BO2, INC (Lemma 11.4) and FEM: $BD$ chart, or for a box with exactly one long arm ($\rho \ge 4$ on the box) the FARL chart (Lemma 9.6) of
    the far families applied to the superset family in which that arm is far with $\rho \ge$ its smallest radius on the box; the boxes inside the shallow box of
    Theorem H or in the region of Theorem M are discarded and covered by certs/shallow_crit, certs/main_cbox); $\rho_A \ge 20$ and/or
    $\rho_C \ge 20$ (the other polar): runs FP, PF, FF (BO2 with tails, Lemma 11.3 (ii)).

In every run the remaining discards are: $\theta_B \ge \pi$ (not a dart), $b_C < b_A$ on the box (the mirror image is in the domain), and the
sets (1)–(7) themselves. $\square$

*Represented boxes (exactness of the grid; adversarial review 2026-10-02).* A box is stored as floats $(p_0, r)$ and denotes the closed
box $[p_0 - r, p_0 + r]$ in exact reals. Every cell (root or skipped) of the cover5 family manifests has exact float ends (the float
expressions $p_0 \mp r$ are exact; spike $t$-cells, whose nominal ends $0.02\cdot 1.5^i$ are not dyadic, are encoded outward by
`run_fam.exact_cell`: the cell $[0.02, 0.03]$ becomes $[t_1, 0.030000000000000002]$, the last cell $[0.2278125\ldots, 0.3]$ becomes
$[0.22781249999999997, 0.30000000000000004] \supseteq [\cdot, 3/10]$, etc.). The processed roots of the legacy and KF runs (certs/, certs2/) are
stored with nearest-rounded encodings that are in general NOT exact-endpoint boxes (4 of 4 crit_collar, 52 of 52 crit_spike_X1,
63 of 75 main_cbox, 256 of 256 shallow_crit and 378 of 882 kf_spike roots are not). No predicate or certificate ever uses a
processed root's own encoding: every processed root, legacy or family, is represented by its outward dyadic snapping
(`grid.dyadic_root`, exact rational arithmetic), which its leaves tile exactly and which `accounting.py` checks to contain the
root's nominal cell. Every leaf box, every skipped cell and every chart box has exact float ends, checked by `nanguard.exact_box`
before a predicate or certificate forms the float expressions $p_0 \mp r$, which are then the true ends. `accounting.py` checks, in rational arithmetic, that the represented boxes of each
family cover its parameter ranges (contiguous cells per axis; the col and spk $t$-cells together cover $[0, 3/10]$), that the KF runs
cover the $t$-ranges $[0, t_1]$ (collar) and $[1310/65536, 19661/65536] \supseteq [t_1, 0.30000000000000004]$ (spike, snapped slabs) used by the
SP predicate, and that the legacy manifests equal their generators. Decimal thresholds (0.0339, 0.05, 0.1772, $174.29^\circ$, 0.3, the
KF sector $38^\circ$, $100^\circ$, 0.5, 1.8, …) are exact decimals in every predicate (`nanguard.dec`, compared with directed ball endpoints).
