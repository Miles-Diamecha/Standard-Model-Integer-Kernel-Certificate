Standard Model One-Loop Integer Kernel and Ω Evaluation
Reproducibility Package for “Integer Structure of the Standard Model One-Loop Decoupling Matrix”

Overview
--------
This repository provides a fully reproducible computation of:

1. The Smith normal form (SNF) of the Standard Model one-loop decoupling matrix ΔW_EM.
2. The unique primitive integer kernel χ_EM = (-10, -18, 1).
3. Its unimodular transport to the gauge–log basis χ = (16, 13, 2).
4. Standard Model coupling evaluations at μ = MZ: α(MZ), α2(MZ), αs(MZ).
5. The derived composite quantity Ω = αs^16 α2^13 α^2.
6. The dimensionless proton–proton gravitational coupling αG_pp = G_N m_p^2 / (ħ c).
7. The closure ratio Z_G = αG_pp / Ω and its inverse Z_G^{-1} = Ω / αG_pp.
8. The leave-one-out strong coupling αs*(MZ).

All integer algebra uses exact arithmetic, and all SM quantities are computed directly from pins.json.
No geometric modeling, beyond-SM physics, or additional assumptions are used.

Repository Contents
-------------------
    sm_integer_kernel_certifier.py   # Main reproducibility script
    pins.json                        # Standard Model and SI input constants
    results.json                     # Auto-generated machine-readable output
    stdout.txt                       # Auto-generated human-readable summary

How to Run
----------
Requires Python ≥ 3.8 and SymPy ≥ 1.9.

    python3 sm_integer_kernel_certifier.py

This produces results.json and stdout.txt.

Matrix Definition
-----------------
The Standard Model one-loop decoupling matrix:

    ΔW_EM = [[8, 8, 224],
             [0, 1,  18]]

Columns correspond to (SU3, SU2, EM). All entries arise from integerized one-loop weights.

SNF and Integer Kernel
----------------------
The script computes:

- SNF invariant factors: (1, 8)
- χ_EM = (-10, -18, 1)
- χ = (16, 13, 2) via unimodular transport M^T

All steps are certified internally in the script.

SM Couplings and Ω
------------------
From pins.json, the script computes:

- α(MZ), α2(MZ), αs(MZ)
- Ω = αs^16 α2^13 α^2
- αG_pp = G_N m_p^2 / (ħ c)
- Z_G and Z_G^{-1}
- αs*(MZ)

Expected Output
---------------
SNF invariant factors: [1, 8]
Primitive kernel in (SU3,SU2,EM): (-10, -18, 1)
Transported kernel in (alpha_s,alpha_2,alpha): (16, 13, 2)

alpha(MZ)      = 0.007815248
alpha2(MZ)     = 0.033789820
alpha_s(MZ)    = 0.118000000
Omega          = 6.459725598437e-39
alphaG_pp      = 5.906149417424e-39
Z_G            = 0.91430345
Z_G^-1         = 1.09372878
alpha_s* (LOO) = 0.117341100

Citing
------
If used academically, cite the Zenodo DOI for this repository and the associated manuscript.

License
-------
MIT License.
