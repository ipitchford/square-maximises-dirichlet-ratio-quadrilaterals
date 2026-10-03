"""Figures 1-2 of the paper (exact geometry from the paper's formulas; matplotlib mathtext, no TeX dependency)."""
import math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Arc, FancyArrowPatch
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})
INK, ACC, GREY = "#1b1b1b", "#b5481d", "#8a8a8a"

def fig_dart(out):
    bA, bC, rA, rC = math.radians(50), math.radians(60), 3.0, 2.6
    B, D = np.array([0.0, 0.0]), np.array([1.0, 0.0])
    A = rA*np.array([math.cos(bA), math.sin(bA)]); C = rC*np.array([math.cos(bC), -math.sin(bC)])
    thB = bA + bC; f = math.sin(thB)/(math.sin(bA)/rC + math.sin(bC)/rA); F = np.array([f, 0.0])
    assert thB < math.pi and f > 1, (thB, f)   # the drawn quadrilateral must be a dart (Lemma 13.1)
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.add_patch(Polygon([A, B, C], closed=True, fill=False, ls=(0, (4, 3)), lw=0.9, ec=GREY))
    ax.add_patch(Polygon([A, B, C, D], closed=True, fc="#f3e1d6", ec=INK, lw=1.4))
    ax.add_patch(Polygon([A, D, C], closed=True, fc="white", ec="none", alpha=0.0))
    ax.plot([B[0], D[0]], [B[1], D[1]], color=INK, lw=0.8, ls=":")
    ax.plot([D[0], F[0]], [0, 0], color=GREY, lw=0.8, ls=":")
    for P, lab, off in ((A, "$A$", (0.05, 0.05)), (B, "$B$", (-0.22, -0.05)), (C, "$C$", (0.06, -0.12)), (D, "$D$", (-0.2, -0.24)), (F, "$F$", (0.06, 0.06))):
        ax.plot(*P, "o", ms=3, color=INK); ax.text(P[0]+off[0], P[1]+off[1], lab, color=INK)
    ax.add_patch(Arc(B, 0.9, 0.9, theta1=0, theta2=math.degrees(bA), color=ACC, lw=1.1)); ax.text(0.5, 0.14, "$b_A$", color=ACC)
    ax.add_patch(Arc(B, 0.65, 0.65, theta1=-math.degrees(bC), theta2=0, color=ACC, lw=1.1)); ax.text(0.36, -0.3, "$b_C$", color=ACC)
    ax.text(*(A/2 + np.array([-0.32, 0.06])), r"$\rho_A$", color=INK); ax.text(*(C/2 + np.array([-0.36, -0.02])), r"$\rho_C$", color=INK)
    ax.text(1.24, 0.42, "notch", color=GREY, fontsize=8.5)
    ax.text(-1.0, 2.55, "dart $Q = ABCD$\nreflex vertex $D$\nhull $T = ABC$ (dashed)", fontsize=8.5, color=INK, va="top")
    ax.text(-1.0, -2.55, r"$B=(0,0)$, $D=(1,0)$; $Q$ is a dart iff $\theta_B=b_A+b_C<\pi$ and $f>1$;" "\n" r"notch depth $\beta = 1-1/f$, $F=(f,0)$ on $AC$", fontsize=8.5, color=INK, va="top")
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.05, 2.45); ax.set_ylim(-3.0, 2.7)
    fig.savefig(out + ".pdf", bbox_inches="tight"); fig.savefig(out + ".png", dpi=220, bbox_inches="tight"); plt.close(fig)
    return dict(f=f, beta=1-1/f)

def fig_transport(out):
    v0 = [np.array(v) for v in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]]
    def Q(a, b, c, d):
        return [np.array([0, 0]) + 0.5*np.array([a+d-b, -a+d+c]), np.array([1, 0]) + 0.5*np.array([-a+d+b, -a-d-c]),
                np.array([1, 1]) + 0.5*np.array([-a-d-b, a-d+c]), np.array([0, 1]) + 0.5*np.array([a-d+b, a+d-c])]
    P0 = Q(0.15, 0.10, 0.20, 0.05); P1 = Q(0.22, 0.16, 0.26, 0.12)
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.0))
    for ax, P, ttl in ((axes[0], P0, r"reference $Q_{p_0}$ (box centre)"), (axes[1], P1, r"$Q_p = \Psi_p(Q_{p_0})$, $p$ in the box")):
        ax.add_patch(Polygon([P[0], P[1], P[2]], closed=True, fc="#f3e1d6", ec=INK, lw=1.2))
        ax.add_patch(Polygon([P[0], P[2], P[3]], closed=True, fc="#e3ecf2", ec=INK, lw=1.2))
        ax.plot([P[0][0], P[2][0]], [P[0][1], P[2][1]], color=ACC, lw=1.6)
        c1 = (P[0]+P[1]+P[2])/3; c2 = (P[0]+P[2]+P[3])/3
        ax.text(*c1, "$K_1$", ha="center"); ax.text(*c2, "$K_2$", ha="center")
        mid = (P[0]+P[2])/2; e = (P[2]-P[0])/np.linalg.norm(P[2]-P[0])
        ax.annotate("", xy=mid + 0.22*e, xytext=mid - 0.02*e, arrowprops=dict(arrowstyle="->", color=ACC, lw=1.2))
        ax.set_title(ttl, fontsize=9); ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-0.35, 1.35); ax.set_ylim(-0.35, 1.35)
    axes[0].text(0.62, 0.43, r"$\hat\tau\parallel\hat e$", color=ACC, fontsize=9)
    fig.text(0.5, 0.55, r"$\Psi_p$: affine on each $K$," "\n" r"$A_K = M_K(p)M_K(p_0)^{-1}$", ha="center", fontsize=8.5)
    fig.patches.append(FancyArrowPatch((0.44, 0.47), (0.56, 0.47), transform=fig.transFigure, arrowstyle="->", mutation_scale=12, color=INK))
    fig.text(0.5, 0.04, r"$w = (A_K\,\hat w/\det A_K)\circ\Psi_p^{-1}$, $\hat w = \det A_K\,\hat\tau + \mathrm{curl}\,\hat\psi$: normal flux continuous across the diagonal",
             ha="center", fontsize=8.5)
    fig.savefig(out + ".pdf", bbox_inches="tight"); fig.savefig(out + ".png", dpi=220, bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    import os; here = os.path.dirname(os.path.abspath(__file__))
    print(fig_dart(os.path.join(here, "fig1_dart_chart"))); fig_transport(os.path.join(here, "fig2_transport")); print("figures written")
