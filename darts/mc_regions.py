"""Monte Carlo volumes of the cover regions in root-normalised coordinates (for the pilot extrapolation).
MAIN : (log|A|, b_A, log|C|, b_C), root = 0.1^4 (half-widths 5% of |A|, |C| ~ 0.1 in log and angle)
SHORT: same coordinates as MAIN
SPIKE: (log t, log|C|, b_C), root = 0.2*0.1*0.1
COLLAR: (log|C|, b_C), root = 0.1*0.1"""
import numpy as np
rng = np.random.default_rng(5)
R = 3.0; T1 = 0.02; T2 = 0.3; X0 = 1.5
N = 4_000_000
uA = rng.uniform(np.log(0.5), np.log(R), N); bA = rng.uniform(0, np.pi, N)
uC = rng.uniform(np.log(0.05), np.log(R), N); bC = rng.uniform(0, np.pi, N)
A = np.exp(uA)[:, None]*np.stack([np.cos(bA), np.sin(bA)], 1); C = np.exp(uC)[:, None]*np.stack([np.cos(bC), -np.sin(bC)], 1)
fD = (1 - C[:, 0])*A[:, 1] + C[:, 1]*(A[:, 0] - 1)
n1 = A[:, 1] - C[:, 1]; n2 = C[:, 0] - A[:, 0]
beta = np.abs(n1*(1 - A[:, 0]) + n2*(0 - A[:, 1]))/np.abs(n1*(0 - A[:, 0]) + n2*(0 - A[:, 1]))
tA = np.where(bA < np.pi/2, np.tan(np.minimum(bA, 1.5707)), np.inf); tC = np.where(bC < np.pi/2, np.tan(np.minimum(bC, 1.5707)), np.inf)
thin = tA + tC <= 0.1
dom = (fD < 0) & (bA <= bC) & (beta > 0.0339) & ~thin
vol = (np.log(R) - np.log(0.5))*np.pi*(np.log(R) - np.log(0.05))*np.pi
main = dom & ((tA >= T2) | (bA >= np.pi/2))
short = dom & (tA >= T1) & (tA < T2) & (A[:, 0] <= X0)
print(f"MAIN  log-polar volume {main.mean()*vol:.3f}  -> roots (0.1^4) {main.mean()*vol/1e-4:.0f}")
print(f"SHORT log-polar volume {short.mean()*vol:.3f}  -> roots (0.1^4) {short.mean()*vol/1e-4:.0f}")
# SPIKE / COLLAR in (log t, log|C|, b_C) and (log|C|, b_C): sample C and t
M = 2_000_000
lt = rng.uniform(np.log(T1), np.log(T2), M); uC = rng.uniform(np.log(0.05), np.log(R), M); bC = rng.uniform(0, np.pi, M)
t = np.exp(lt); C = np.exp(uC)[:, None]*np.stack([np.cos(bC), -np.sin(bC)], 1)
Ax = X0; Ay = X0*t
fD = (1 - C[:, 0])*Ay + C[:, 1]*(Ax - 1)
ok = (fD < 0) & (np.arctan(t) <= bC) & (bC <= np.radians(170))
v3 = (np.log(T2) - np.log(T1))*(np.log(R) - np.log(0.05))*np.pi
print(f"SPIKE volume (log t, log|C|, b_C) {ok.mean()*v3:.3f} -> roots (0.2*0.1*0.1) {ok.mean()*v3/2e-3:.0f}")
okc = (bC >= np.arctan(T1)) & (bC <= np.radians(170)) & (C[:, 1] < 0)
v2 = (np.log(R) - np.log(0.05))*np.pi
print(f"COLLAR volume (log|C|, b_C) {okc.mean()*v2:.3f} -> roots (0.1^2) {okc.mean()*v2/1e-2:.0f}")
