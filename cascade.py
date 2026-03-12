import numpy as np
import qutip as qt

# -------------------------
#
# https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.129.193604/Supplement_Impact_TPE_new.pdf
#
# -------------------------
# 1) Basis: |G>, |XH>, |XV>, |B>
# -------------------------
G, XH, XV, B = [qt.basis(4, i) for i in range(4)]
PG  = G*G.dag()
PXH = XH*XH.dag()
PXV = XV*XV.dag()
PB  = B*B.dag()

# Transition operators sigma_H/V = |G><X| + |X><B|
sigma_H = (G*XH.dag()) + (XH*B.dag())
sigma_V = (G*XV.dag()) + (XV*B.dag())

# -------------------------
# 2) Atomic units:
# -------------------------
Eh = 27211.  # Eh = 27211 meV
Th = 2.42e-5 # Th = 2.42 x 10^-5 ps

# -------------------------
# and Parameters:
# -------------------------
hbar = 1.0              # atomic units
delta = 0.001/Eh        # fine structure splitting (meV)
EB = 4.0/Eh             # biexciton binding energy (meV)
gammaX = 0.005*Th       # exciton decay rate (1/Th)
gammaB = 0.010*Th       # biexciton decay rate (1/Th)

# Two-photon resonance: Delta_XL = EB/2
Delta_XL = EB/2

# Laser pulse
tL = 100.0/Th           # (1/Th)
FWHM = 20./Th          # (1/Th). There are problems below 12.33 - dont know why
Theta = np.sqrt(EB*FWHM/np.sqrt(2.*np.pi*np.log(2.)))*np.pi
#Theta = 20.0            # (2\pi)
print(Theta)

def Omega(t):
    arg = 4*np.log(2)*((t-tL)/FWHM)**2
    return np.sqrt(4*np.log(2)/np.pi) * Theta/FWHM * np.exp(-arg)

# Laser polarization coefficients (real)
alphaH = 1.0
alphaV = 0.0
sigma_L = alphaH*sigma_H + alphaV*sigma_V

# -------------------------
# 3) Hamiltonian H(t) (Eq. A1) :contentReference[oaicite:12]{index=12}
# -------------------------
H0 = (Delta_XL + delta/2)*PXH + (Delta_XL - delta/2)*PXV + (2*Delta_XL - EB)*PB

# Interaction: -(ħ/2) Ω(t) (σ_L + σ_L^†)
Hint_op = -(hbar/2)*(sigma_L + sigma_L.dag())

# QuTiP time-dependent Hamiltonian format: [H0, [Hint_op, f(t)]]
H = [H0, [Hint_op, lambda t, args: Omega(t)]]

# -------------------------
# 4) Collapse operators (Eq. A6 with Lindblad terms) :contentReference[oaicite:13]{index=13}
# -------------------------
c_ops = []
# X -> G at rate gammaX
c_ops += [np.sqrt(gammaX) * (G*XH.dag()), np.sqrt(gammaX) * (G*XV.dag())]
# B -> X at rate gammaB/2 per polarization
c_ops += [np.sqrt(gammaB/2) * (XH*B.dag()), np.sqrt(gammaB/2) * (XV*B.dag())]

# -------------------------
# 5) Time grids for integration over t and tau (Eq. A10b) :contentReference[oaicite:14]{index=14}
# -------------------------
# Choose windows "long enough" to cover pulse + full decay.
multiplier = 5
tmax = multiplier*FWHM + multiplier/max(gammaX, gammaB)
taumax = multiplier/max(gammaX, gammaB)

Nt = 200
Ntau = 200
tlist = np.linspace(0, tmax, Nt)
taulist = np.linspace(0, taumax, Ntau)

# Initial state: ground state 
rho0 = PG

# -------------------------
# 6) Build integrated G^(2) tensor in {H,V}x{H,V} indices
#    G^(2)_{jk,lm} = ∫dt ∫dτ  < σj†(t) σk†(t+τ) σm(t+τ) σl(t) >
# -------------------------
pol = {"H": sigma_H, "V": sigma_V}

# We'll fill a 4x4 matrix in basis |HH>,|HV>,|VH>,|VV|
basis_2p = [("H","H"), ("H","V"), ("V","H"), ("V","V")]
G2_int = np.zeros((4,4), dtype=complex)

# Precompute state evolution rho(t) so QRT can use it internally
# (QuTiP correlation functions will solve as needed, but this keeps things consistent.)
_ = qt.mesolve(H, rho0, tlist, c_ops, e_ops=[])

# Use QRT: correlation_2op_2t expects <A(t) B(t+τ)> with QRT.
# We need four-operator structure, so we use: A = σj† ; B = σk†σm ; and condition with σl at time t.
# The most direct in QuTiP is to use correlation_2op_2t with "state prep" by applying σl on rho(t).
# We'll implement Eq. (A11)-style by looping over t and using mesolve for the τ-propagation.

def prop_to_t_states():
    # Returns list of rho(t) for each t in tlist
    sol = qt.mesolve(H, rho0, tlist, c_ops, e_ops=[])
    return sol.states

states_t = prop_to_t_states()

# plot biexciton population
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(tlist*Th, Omega(tlist)/Th, label='$\Omega(t)$')
ax.plot(tlist*Th, [s_t[0,0].real for s_t in states_t], label='$\\rho_{GG}(t)$')
ax.plot(tlist*Th, [s_t[3,3].real for s_t in states_t], label='$\\rho_{BB}(t)$')
ax.plot(tlist*Th, [np.abs(s_t[3,0]) for s_t in states_t], label='$\\rho_{BG}(t)$')
#ax.set_xlim([0.,200.])
ax.legend()
fig.savefig('biexciton.png')   # save the figure to file
plt.close(fig) 


# Do double integral by Riemann sums
dt = tlist[1] - tlist[0]
dtau = taulist[1] - taulist[0]

for a,(j,k) in enumerate(basis_2p):
    for b,(l,m) in enumerate(basis_2p):
        acc = 0.0 + 0.0j
        sig_j = pol[j]
        sig_k = pol[k]
        sig_l = pol[l]
        sig_m = pol[m]

        for it, t in enumerate(tlist):
            rho_t = states_t[it]

            # Apply first detection operators at time t:  σ_l ρ(t) σ_j†
            rho_cond = sig_l * rho_t * sig_j.dag()

            # Propagate from t to t+τ: use mesolve with initial rho_cond
            # with a time list in absolute time: [t, t+τ] -> we can propagate across taulist by shifting
            t_abs_list = t + taulist

            sol_tau = qt.mesolve(H, rho_cond, t_abs_list, c_ops, e_ops=[])

            # For each τ: take Tr[ σ_k† σ_m  ρ(t+τ | cond) ]
            # (matches Eq. A11 structure)
            vals = [ (sig_k.dag()*sig_m*st).tr() for st in sol_tau.states ]
            acc += np.sum(vals) * dtau
            breakpoint()
        G2_int[a,b] = acc * dt

# Normalize to get rho_2p (Eq. A10a)
rho_2p = G2_int / np.trace(G2_int)

print("rho_2p =\n", rho_2p)

# -------------------------
# 7) Concurrence (Eqs. A12-A13) :contentReference[oaicite:19]{index=19}
# -------------------------
T = np.array([[0,0,0,-1],
              [0,0,1, 0],
              [0,1,0, 0],
              [-1,0,0,0]], dtype=complex)

M = rho_2p @ T @ rho_2p.conj() @ T
evals = np.sort(np.real(np.linalg.eigvals(M)))[::-1]
C = max(0.0, np.sqrt(evals[0]) - np.sqrt(evals[1]) - np.sqrt(evals[2]) - np.sqrt(evals[3]))
print("Concurrence =", C)
