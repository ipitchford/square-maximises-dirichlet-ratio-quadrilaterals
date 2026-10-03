"""Conversion between the canonical (thA, thC, alpha, beta) dart parametrisation and the BD chart."""
import numpy as np
from dartgeom import dart

def to_bd(A, B, C, D):
    """similarity sending B->(0,0), D->(1,0); returns A', C' with A' in the upper half plane (reflect if needed)"""
    z = lambda P: complex(P[0], P[1])
    b, d = z(B), z(D)
    f = lambda P: (z(P) - b)/(d - b)
    a, c = f(A), f(C)
    if a.imag < 0: a, c = a.conjugate(), c.conjugate()
    return np.array([a.real, a.imag]), np.array([c.real, c.imag])

def from_angles(thA, thC, al, be):
    return to_bd(*dart(thA, thC, al, be))

def fD(A, C):
    return (1 - C[0])*A[1] + C[1]*(A[0] - 1)

if __name__ == '__main__':
    for p in [(60, 60, 0.35, 0.3), (60, 60, 0.05, 0.9), (60, 60, 0.475, 0.05), (25, 95, 0.405, 0.55), (45, 75, 0.675, 0.25),
              (90, 45, 0.3, 0.5), (30, 30, 0.2, 0.6)]:
        A, C = from_angles(np.radians(p[0]), np.radians(p[1]), p[2], p[3])
        print(p, 'A', np.round(A, 4), 'C', np.round(C, 4), 'fD', round(fD(A, C), 4), '|A|', round(np.linalg.norm(A), 3), '|C|', round(np.linalg.norm(C), 3),
              'bA', round(np.degrees(np.arctan2(A[1], A[0])), 2), 'bC', round(np.degrees(np.arctan2(-C[1], C[0])), 2))
