"""Dart geometry. Canonical placement: C=(0,0), A=(1,0), B=(bx,by), by>0 (A->B->C ccw);
reflex vertex D = alpha*A + beta*B + gamma*C (barycentric, all >0).  Dart Q = ABD u BCD (boundary A->B->C->D->A).
Parameters (thA, thC, alpha, beta): angles of the hull triangle T=ABC at A and C; gamma = 1-alpha-beta.
Symmetry: x -> 1-x swaps A<->C, thA<->thC, alpha<->gamma."""
import numpy as np

def B_from_angles(thA, thC):
    h = 1.0/(1.0/np.tan(thA) + 1.0/np.tan(thC))
    return np.array([h/np.tan(thC), h])

def dart(thA, thC, al, be):
    C = np.array([0.0, 0.0]); A = np.array([1.0, 0.0]); B = B_from_angles(thA, thC)
    ga = 1.0 - al - be
    D = al*A + be*B + ga*C
    return A, B, C, D

def dart_from_B(bx, by, al, be):
    C = np.array([0.0, 0.0]); A = np.array([1.0, 0.0]); B = np.array([bx, by], float)
    ga = 1.0 - al - be
    return A, B, C, al*A + be*B + ga*C

def area2(P, Q, R):
    return (Q[0]-P[0])*(R[1]-P[1]) - (Q[1]-P[1])*(R[0]-P[0])

def angles_T(A, B, C):
    def ang(P, Q, R):
        u = Q - P; v = R - P
        return np.arccos(np.clip(np.dot(u, v)/np.linalg.norm(u)/np.linalg.norm(v), -1, 1))
    return ang(A, B, C), ang(B, C, A), ang(C, A, B)
