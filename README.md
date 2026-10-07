<div align="center">

# ⚛️ QuBindLab

### Validating VQE on H₂, the first step toward quantum drug-target binding energies

*A reproducible, benchmarked Variational Quantum Eigensolver pipeline: real molecular Hamiltonian → UCCSD ansatz → exact-baseline validation → dissociation curve.*

![Python](https://img.shields.io/badge/Python-3.9--3.12-3776AB?logo=python&logoColor=white)
![Qiskit](https://img.shields.io/badge/Qiskit-1.x-6929C4)
![Qiskit Nature](https://img.shields.io/badge/qiskit--nature-0.7.x-6929C4)
![PySCF](https://img.shields.io/badge/PySCF-chemistry-green)
![License](https://img.shields.io/badge/License-MIT-yellow)
![Status](https://img.shields.io/badge/Status-Research%20Prototype-orange)

</div>

---

## 📌 Overview

Drug discovery hinges on **binding energy**:

```
ΔE_bind = E(complex) − E(ligand) − E(target)
```

Each term is a many-electron ground-state energy, which is exponentially expensive to compute exactly on classical hardware. **QuBindLab** builds and validates the quantum-computing machinery for this problem on the smallest meaningful system, the **H₂ molecule**, where every result can be checked against an exact classical solution.

This project replaces an earlier toy prototype (hard-coded Hamiltonian, circular "ground-state" reference) with a scientifically grounded pipeline.

## ✨ Features

- **Real Hamiltonian** from first principles (PySCF, STO-3G basis) via `qiskit-nature`
- **Parity mapping + two-qubit reduction**: 4 spin orbitals → **2 qubits**
- **UCCSD ansatz** initialized from the **Hartree-Fock** state
- **Exact baseline** (classical diagonalization in the correct symmetry sector), so no circular validation
- **Multi-restart COBYLA** with seeded randomness for reproducibility
- **Dissociation curve scan** over H–H distance (0.4 to 2.5 Å)
- **Error analysis** against **chemical accuracy (1.6 mHa ≈ 1 kcal/mol)**
- Publication-ready 300-dpi plots

## 🧠 How it works

```
 Geometry (d) ─► PySCF integrals ─► Second-quantized H ─► Parity mapper (2 qubits)
                                                               │
                         ┌─────────────────────────────────────┤
                         ▼                                     ▼
                Exact diagonalization                HF state + UCCSD ansatz
                  (reference energy)                          │
                         │                         Statevector estimator ⟨H⟩
                         │                                     │
                         │                          COBYLA optimizer (θ)
                         │                                     │
                         └──────────► compare ◄────────────────┘
                                         │
                            Dissociation curve + error plot
```

The **variational principle** guarantees `E_VQE ≥ E_exact`. A VQE value below the exact one signals a bug, so it doubles as a built-in sanity check.

## 🚀 Quick Start

### 1. Clone

```bash
git clone https://github.com/RahulPatil-Tech/QuBindLab.git
cd QuBindLab
```

### 2. Create an environment

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
```

### 3. Install dependencies

```bash
pip install qiskit qiskit-nature qiskit-algorithms pyscf scipy matplotlib
```

### 4. Run

```bash
python -u h2_vqe_scan.py
```

### Sample console output

```
d=0.40 Å | exact=... | VQE=... | err=... Ha
...
d=0.74 Å | exact=-1.137... | VQE=-1.137... | err=~1e-07 Ha
...
Equilibrium (scan minimum): ~0.74 Å (experimental ≈ 0.74 Å)
```

> Exact digits depend on your library versions and scan grid.

## 📊 Output

Running the script generates `h2_dissociation.png` with two panels:

| Panel | Shows |
|---|---|
| **Top** | Exact (FCI/STO-3G) curve vs VQE points over bond distance |
| **Bottom** | \|VQE − exact\| on a log scale with the 1.6 mHa chemical-accuracy line |

**Expected physics**

- Energy minimum near **0.74 Å**, around **−1.137 Ha**
- Curve flattens toward two separated H atoms at large distance
- VQE error stays below chemical accuracy at every geometry (noiseless simulation)

## 🔬 Worked example (H₂ at 0.735 Å)

| Quantity | Value (Ha) |
|---|---|
| Nuclear repulsion | ≈ 0.7200 |
| Hartree-Fock total energy | ≈ −1.117 |
| VQE / FCI total energy | ≈ −1.137 |
| Correlation energy | ≈ −0.020 (≈ 12.7 kcal/mol) |

Correlation energy is more than 12× larger than chemical accuracy, which is why Hartree-Fock alone is not enough and VQE with UCCSD is useful. A full hand derivation is in [`H2_VQE_Documentation.md`](H2_VQE_Documentation.md).

## 📁 Project Structure

```
QuBindLab/
├── h2_vqe_scan.py              # Main pipeline: scan, VQE, exact baseline, plots
├── H2_VQE_Documentation.md     # Full theory, derivations, worked calculation
├── h2_dissociation.png         # Generated output plot
└── README.md
```

## ⚙️ Configuration

Edit these in `h2_vqe_scan.py`:

| Parameter | Location | Default | Purpose |
|---|---|---|---|
| Distance grid | `np.linspace(0.4, 2.5, 15)` | 15 points | Scan range (Å) |
| Basis set | `PySCFDriver(basis=...)` | `sto3g` | Larger basis = more qubits |
| Restarts | `vqe_energy(restarts=...)` | 2 | Guard against local minima |
| Max iterations | `minimize(options={"maxiter": ...})` | 300 | COBYLA budget |
| Seed | `vqe_energy(seed=...)` | 0 | Reproducibility |

## ⚠️ Limitations

- **Noiseless statevector simulation**: no hardware or shot noise
- **Minimal basis (STO-3G)**: qualitatively right, quantitatively approximate
- **2 qubits only**: validates the method, shows no quantum advantage (classical solves H₂ instantly)
- **Gas-phase electronic energy**: does not capture solvation, entropy or protein dynamics

## 🗺️ Roadmap

- [x] Real H₂ Hamiltonian + UCCSD + exact baseline
- [x] Dissociation curve and chemical-accuracy error analysis
- [ ] Noisy simulation (Qiskit Aer noise model, shot-based estimator)
- [ ] Error mitigation (ZNE)
- [ ] Optimizer comparison (COBYLA vs SPSA vs L-BFGS-B)
- [ ] Ansatz comparison (UCCSD vs EfficientSU2 vs ADAPT-VQE)
- [ ] LiH / H₂O with active-space selection
- [ ] Water-dimer interaction energy (first real "binding" calculation)
- [ ] Ligand + binding-site fragments from PDB/SDF (RDKit, Biopython)
- [ ] Docking-guided pose scan (AutoDock Vina) with quantum refinement
- [ ] Excited states (SSVQE / qEOM)
- [ ] Unit tests and CI

## 🛠️ Troubleshooting

| Problem | Fix |
|---|---|
| `ModuleNotFoundError: qiskit_nature` | Activate your venv, or `python3 -m pip install qiskit-nature` |
| `externally-managed-environment` | Use a virtual environment (see Quick Start) |
| PySCF fails to install | Use Python 3.9 to 3.12 |
| API/import errors | Pin `qiskit-nature==0.7.*` with Qiskit 1.x |
| VQE below exact energy | Bug: check mapper/problem consistency and nuclear-repulsion handling |

## 📚 References

1. Peruzzo et al., *A variational eigenvalue solver on a photonic quantum processor*, Nat. Commun. 5, 4213 (2014)
2. O'Malley et al., *Scalable Quantum Simulation of Molecular Energies*, Phys. Rev. X 6, 031007 (2016)
3. Kandala et al., *Hardware-efficient VQE for small molecules and quantum magnets*, Nature 549, 242 (2017)
4. McArdle et al., *Quantum computational chemistry*, Rev. Mod. Phys. 92, 015003 (2020)
5. [Qiskit Nature docs](https://qiskit-community.github.io/qiskit-nature/) · [PySCF docs](https://pyscf.org/)

## 🤝 Contributing

Issues and pull requests are welcome. Good first contributions: unit tests (VQE vs exact on H₂), a noise-model benchmark, or a new molecule.

---

<div align="center">

<sub>Built with ⚛️ by **[Rahul Patil](https://github.com/RahulPatil-Tech)** • Released under the [MIT License](LICENSE)</sub>

</div>
