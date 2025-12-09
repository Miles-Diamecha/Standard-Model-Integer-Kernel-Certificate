#!/usr/bin/env python3
"""
sm_integer_kernel_certifier.py

Reproducibility script for the Standard Model one-loop integer certificate.

Functions:
- Computes the Smith normal form (SNF) of the 2×3 decoupling matrix ΔW_EM.
- Extracts the unique primitive integer kernel χ_EM = (-10, -18, 1).
- Applies the unimodular transport matrix M to obtain χ = (16, 13, 2) in the
  gauge–log (αs, α2, α) basis.
- Optionally loads pins.json to evaluate α(MZ), α2(MZ), Ω, α_G^pp, closure, and LOO α_s*.

This script contains:
- No geometry, no gravitational modeling, and no extra assumptions.
- Only exact-integer SNF algebra + Standard Model coupling pins.

For the accompanying manuscript: this script certifies all integer claims and
numeric evaluations reported in the text.
"""


import json
from sympy import Matrix, ilcm, igcd, ZZ

# ---------- SNF / kernel block ----------

def compute_snf_kernel():
    A = Matrix([[8, 8, 224],
                [0, 1,  18]])  # ΔW_EM in (SU3, SU2, EM)

    U = D = V = None

    # 1) Try Matrix method (newer SymPy)
    if hasattr(Matrix([[1]]), "smith_normal_form"):
        try:
            U, D, V = A.smith_normal_form()  # U*A*V = D
        except Exception:
            U = D = V = None

    # 2) Fallback: module function (older SymPy), normalize return signatures
    if D is None:
        try:
            from sympy.matrices.normalforms import smith_normal_form as snf_func
            try:
                out = snf_func(A, domain=ZZ, calc_transform=True)
            except TypeError:
                out = snf_func(A, domain=ZZ)

            # Normalize various return signatures: (D,U,V), (U,D,V), or (D,(U,V))
            if isinstance(out, tuple):
                if len(out) == 3:  # could be (D,U,V) or (U,D,V) or (U,V,D)
                    for Dm, Um, Vm in [(out[0], out[1], out[2]),
                                       (out[1], out[0], out[2]),
                                       (out[2], out[0], out[1])]:
                        try:
                            if Um*A*Vm == Dm:
                                D, U, V = Dm, Um, Vm
                                break
                        except Exception:
                            pass
                elif len(out) == 2 and isinstance(out[1], tuple) and len(out[1]) == 2:
                    D, (U, V) = out
            else:
                D = out  # D only
        except Exception:
            pass

    m, n = A.shape
    if D is not None:
        assert D.shape == (m, n)
        # rank = number of nonzero diagonal entries
        r = sum(1 for i in range(min(m, n)) if D[i, i] != 0)
        assert r == 2, f"Expected rank 2; got {r}"
        # columns beyond rank must be all zeros (here: the 3rd column)
        for j in range(r, n):
            assert all(D[i, j] == 0 for i in range(m)), "Trailing column not zero in D"
    else:
        r = 2  # expected for this A; continue without D/U/V assertions

    # --- Kernel from SNF if V present (preferred) ---
    # SNF gives a right-integer kernel ker_Z(A) via the last column of V.
    chiZ_snf = None
    if V is not None and D is not None:
        chiZ_snf = V[:, -1]  # last column spans ker_Z(A) since n - r = 1
        if chiZ_snf[-1] < 0:
            chiZ_snf = -chiZ_snf

    # --- Fallback: rational nullspace → integerize → primitive ---
    chiQ = A.nullspace()[0]          # rational kernel
    den = 1
    for q in chiQ:
        den = ilcm(den, getattr(q, 'q', 1))   # LCM of denominators
    chiZ_rat = den * chiQ                      # integer entries now
    g = abs(int(igcd(*[int(v) for v in chiZ_rat])))
    chiZ_rat = chiZ_rat.applyfunc(lambda v: v // g)  # primitive integer vector
    if chiZ_rat[-1] < 0:
        chiZ_rat = -chiZ_rat

    # Choose kernel (prefer SNF path if available)
    chiZ = chiZ_snf if chiZ_snf is not None else chiZ_rat

    # Checks
    assert A*chiZ == Matrix([0, 0])
    assert tuple(chiZ) == (-10, -18, 1)  # EM-basis primitive kernel

    # Unimodular transport to (αs, α2, α)
    M = Matrix([[-5, -3, -2],
                [ 2,  1,  1],
                [ 2,  1,  0]])
    assert M.det() in (1, -1)
    chi_gauge = M.T * chiZ
    assert tuple(chi_gauge) == (16, 13, 2)

    # Invariant factors (if D is available)
    diag_list = None
    if D is not None:
        diag_list = [int(D[i, i]) for i in range(min(m, n)) if D[i, i] != 0]

    return {
        "DeltaW_EM": [[int(x) for x in row] for row in A.tolist()],
        "snf_diag": diag_list,          # expected [1, 8]
        "chi_EM": tuple(int(v) for v in chiZ),
        "chi_gauge": tuple(int(v) for v in chi_gauge)
    }

# ---------- Pins / Omega block ----------

def load_pins(path="pins.json"):
    with open(path, "r") as f:
        return json.load(f)

def alpha2(alpha_em, sin2w):
    return alpha_em / sin2w

def omega(a_s, a2, a):
    return (a_s**16) * (a2**13) * (a**2)

def alpha_G_pp(G, mp, hb, c):
    return G * (mp**2) / (hb * c)

def loo_alpha_s_star(aGpp, a2, a):
    return (aGpp / (a2**13 * a**2))**(1.0 / 16.0)

def compute_pins_results(pins_json="pins.json"):
    P = load_pins(pins_json)
    Pp = P["pins"]

    aem = 1.0 / float(Pp["inv_alpha_MZ"])
    s2w = float(Pp["sin2_thetaW_MZ"])
    a_s = float(Pp["alpha_s_MZ"])

    a2 = alpha2(aem, s2w)
    aGpp = alpha_G_pp(
        float(Pp["G_N_SI"]),
        float(Pp["m_p_SI_kg"]),
        float(Pp["hbar_SI_Js"]),
        float(Pp["c_SI_mps"])
    )

    Om = omega(a_s, a2, aem)
    Z_G     = aGpp / Om      # Z_G       = alphaG_pp / Omega
    Z_G_inv = Om / aGpp      # Z_G^{-1}  = Omega / alphaG_pp
    a_s_star = loo_alpha_s_star(aGpp, a2, aem)

    return {
        "alpha_em_MZ": aem,
        "alpha2_MZ": a2,
        "alpha_s_MZ": a_s,
        "Omega": Om,
        "alpha_G_pp": aGpp,
        "Z_G": Z_G,
        "Z_G_inverse": Z_G_inv,
        "alpha_s_star_MZ": a_s_star
    }


# ---------- Main ----------

if __name__ == "__main__":
    snf_data = compute_snf_kernel()
    pins_data = compute_pins_results("pins.json")

    out = {
        "snf": snf_data,
        "pins_results": pins_data
    }

    with open("results.json", "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    # Readable stdout
    s = []
    if snf_data["snf_diag"] is not None:
        s.append(f"SNF invariant factors: {snf_data['snf_diag']}  (expected [1, 8])")
    else:
        s.append("SNF transform matrices not available; used rational nullspace path.")
    s.append(f"Primitive kernel in (SU3,SU2,EM): {snf_data['chi_EM']}")
    s.append(f"Transported kernel in (alpha_s,alpha_2,alpha): {snf_data['chi_gauge']}")
    s.append("")
    s.append(f"alpha(MZ)      = {pins_data['alpha_em_MZ']:.9f}")
    s.append(f"alpha2(MZ)     = {pins_data['alpha2_MZ']:.9f}")
    s.append(f"alpha_s(MZ)    = {pins_data['alpha_s_MZ']:.9f}")
    s.append(f"Omega          = {pins_data['Omega']:.12e}")
    s.append(f"alphaG_pp      = {pins_data['alpha_G_pp']:.12e}")
    s.append(f"Z_G        = alphaG_pp/Omega      = {pins_data['Z_G']:.8f}")
    s.append(f"Z_G^-1     = Omega/alphaG_pp      = {pins_data['Z_G_inverse']:.8f}")
    s.append(f"alpha_s* (LOO) = {pins_data['alpha_s_star_MZ']:.9f}")


    stdout_str = "\n".join(s) + "\n"
    print(stdout_str, end="")
    with open("stdout.txt", "w") as f:
        f.write(stdout_str)
