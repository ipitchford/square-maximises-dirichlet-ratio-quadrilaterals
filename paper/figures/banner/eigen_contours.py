"""Level curves of second Dirichlet eigenfunctions for the Evidence Press banner (equal-area shapes, unit area).
Square: u2 = sin(2 pi x) sin(pi y) (exact). Disk: u2 = J1(j11 r) cos(theta) (exact). Equilateral triangle: the eigenfunction of the
double second eigenvalue that is odd under reflection in the altitude, computed with P2 finite elements (scikit-fem; refined mesh).
Exact ratios: square 5/2, equilateral triangle 7/3, disk j11^2/j01^2. The triangle's computed ratio is printed as a check.
Output: eigen_contours.json (normalised coordinates: each shape centred at the origin, unit area; y up)."""
import json, numpy as np
from scipy.special import jv, jn_zeros
import contourpy
import skfem as fe
from skfem.models.poisson import laplace, mass
from scipy.sparse.linalg import eigsh

LEVELS = [-0.8, -0.55, -0.3, 0.3, 0.55, 0.8]
def contours(F, X, Y, mask, levels=LEVELS):
    out = {}
    G = np.ma.masked_array(F, mask=~mask)
    gen = contourpy.contour_generator(X, Y, G)
    for lv in levels:
        paths = [np.round(c[::3], 4).tolist() for c in gen.lines(lv) if len(c) > 8]
        out[str(lv)] = paths
    return out

res = {}
n = 601
# square of unit area centred at 0
s = 1.0; X, Y = np.meshgrid(np.linspace(-s/2, s/2, n), np.linspace(-s/2, s/2, n))
F = -np.sin(2*np.pi*(X + s/2)) * np.sin(np.pi*(Y + s/2))      # sign chosen so that the negative lobe is at x < 0, as for the disk and the triangle
res["square"] = {"outline": [[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]], "nodal": [[0, -.5], [0, .5]], "levels": contours(F, X, Y, np.ones_like(F, bool)),
                 "ratio": "5/2", "ratio_value": 2.5}
# disk of unit area
R = 1/np.sqrt(np.pi); j11 = jn_zeros(1, 1)[0]; j01 = jn_zeros(0, 1)[0]
X, Y = np.meshgrid(np.linspace(-R, R, n), np.linspace(-R, R, n)); r = np.hypot(X, Y)
F = jv(1, j11*r/R) * np.where(r > 0, X/np.maximum(r, 1e-300), 0); F /= np.nanmax(np.abs(np.where(r <= R, F, 0)))
res["disk"] = {"radius": R, "nodal": [[0, -R], [0, R]], "levels": contours(F, X, Y, r <= R), "ratio": "j11^2/j01^2", "ratio_value": (j11/j01)**2}
# equilateral triangle of unit area, apex up, centroid at 0
a = np.sqrt(4/np.sqrt(3)); h = a*np.sqrt(3)/2
P = np.array([[-a/2, -h/3], [a/2, -h/3], [0, 2*h/3]])
m = fe.MeshTri(P.T, np.array([[0, 1, 2]]).T).refined(6)
basis = fe.Basis(m, fe.ElementTriP2())
K = fe.asm(laplace, basis); M = fe.asm(mass, basis)
D = basis.get_dofs().all()
I = basis.complement_dofs(D)
w, V = eigsh(K[I][:, I], k=4, M=M[I][:, I], sigma=0)
order = np.argsort(w); w = w[order]; V = V[:, order]
print("triangle eigenvalues*area:", w[:3], "ratio", w[1]/w[0], "exact 7/3 =", 7/3, "double:", w[2]/w[1], flush=True)
# choose the combination in the 2nd eigenspace that is odd under x -> -x (vertex values; the refined mesh is mirror-symmetric)
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); from matplotlib.tri import Triangulation; import matplotlib.pyplot as plt
nv = m.p.shape[1]
u2 = np.zeros(K.shape[0]); u3 = np.zeros(K.shape[0]); u2[I] = V[:, 1]; u3[I] = V[:, 2]
z2, z3 = u2[:nv], u3[:nv]                       # P2 dofs on the vertices come first in skfem's numbering
mirror = cKDTree(m.p.T).query(np.c_[-m.p[0], m.p[1]])[1]
A = np.c_[z2 + z2[mirror], z3 + z3[mirror]]     # even parts
_, _, vt = np.linalg.svd(A, full_matrices=False); c = vt[-1]
z = c[0]*z2 + c[1]*z3; z /= np.max(np.abs(z))
if np.mean(z[m.p[0] < 0]) > 0: z = -z          # negative lobe on the left (as for the square and the disk)
print("odd-mode residual even part:", np.max(np.abs(z + z[mirror])) / 2, flush=True)
tri = Triangulation(m.p[0], m.p[1], m.t.T)
cs = plt.tricontour(tri, z, levels=sorted(LEVELS))
tri_levels = {}
for lv, segs in zip(sorted(LEVELS), cs.allsegs):
    tri_levels[str(lv)] = [np.round(sg[::3], 4).tolist() for sg in segs if len(sg) > 8]
res["triangle"] = {"outline": P.round(5).tolist(), "nodal": [[0, -h/3], [0, 2*h/3]], "levels": tri_levels,
                   "ratio": "7/3", "ratio_value": 7/3, "computed_ratio": float(w[1]/w[0])}
for k in ("square", "disk"):
    res[k]["sign"] = "negative lobe x<0"
json.dump(res, open("eigen_contours.json", "w"))
print("written; ratios:", {k: res[k]["ratio_value"] for k in res})
