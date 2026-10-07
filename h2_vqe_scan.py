"""
H2 dissociation curve: VQE (UCCSD) vs exact diagonalization.
 
Install:
    pip install qiskit qiskit-nature qiskit-algorithms pyscf scipy matplotlib
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from qiskit.primitives import StatevectorEstimator
from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_nature.units import DistanceUnit


def build_problem(distance: float):
    driver = PySCFDriver(
        atom=f"H 0 0 0; H 0 0 {distance}",
        basis="sto3g",
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM,
    )
    return driver.run()


def exact_energy(problem, mapper) -> float:
    """Classical reference (total energy incl. nuclear repulsion, in Hartree)."""
    solver = GroundStateEigensolver(mapper, NumPyMinimumEigensolver())
    return float(solver.solve(problem).total_energies[0])


def vqe_energy(problem, mapper, seed: int = 0, restarts: int = 1):
    """VQE with UCCSD ansatz. Returns (total_energy, history)."""
    hamiltonian = mapper.map(problem.hamiltonian.second_q_op())
    nuc = problem.hamiltonian.nuclear_repulsion_energy

    hf = HartreeFock(problem.num_spatial_orbitals,
                     problem.num_particles, mapper)
    ansatz = UCCSD(
        problem.num_spatial_orbitals,
        problem.num_particles,
        mapper,
        initial_state=hf,
    )

    estimator = StatevectorEstimator()
    rng = np.random.default_rng(seed)
    best_e, best_hist = np.inf, []

    for r in range(restarts):
        history = []

        def cost(params):
            pub = (ansatz, hamiltonian, params)
            e = float(estimator.run([pub]).result()[0].data.evs)
            history.append(e + nuc)
            return e

        # first restart starts at Hartree-Fock (all zeros), others are random
        x0 = np.zeros(
            ansatz.num_parameters) if r == 0 else rng.uniform(-0.1, 0.1, ansatz.num_parameters)
        res = minimize(cost, x0, method="COBYLA", options={"maxiter": 300})
        if res.fun + nuc < best_e:
            best_e, best_hist = res.fun + nuc, history

    return best_e, best_hist


def main():
    distances = np.linspace(0.4, 2.5, 15)  # Angstrom
    exact, vqe = [], []

    for d in distances:
        problem = build_problem(d)
        mapper = ParityMapper(num_particles=problem.num_particles)
        e_exact = exact_energy(problem, mapper)
        e_vqe, _ = vqe_energy(problem, mapper, restarts=2)
        exact.append(e_exact)
        vqe.append(e_vqe)
        print(
            f"d={d:.2f} Å | exact={e_exact:.6f} | VQE={e_vqe:.6f} | err={abs(e_vqe - e_exact):.2e} Ha")

    exact, vqe = np.array(exact), np.array(vqe)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    ax1.plot(distances, exact, "k-", label="Exact (FCI in STO-3G)")
    ax1.plot(distances, vqe, "o", color="tab:blue", label="VQE (UCCSD)")
    ax1.set_ylabel("Total energy (Hartree)")
    ax1.set_title("H2 dissociation curve")
    ax1.legend()
    ax1.grid(True, linestyle=":", alpha=0.6)

    ax2.semilogy(distances, np.abs(vqe - exact) + 1e-12, "o-", color="tab:red")
    ax2.axhline(1.6e-3, color="gray", linestyle="--",
                label="Chemical accuracy (1.6 mHa)")
    ax2.set_xlabel("H–H distance (Å)")
    ax2.set_ylabel("|VQE − exact| (Hartree)")
    ax2.legend()
    ax2.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig("h2_dissociation.png", dpi=300)
    plt.show()

    eq = distances[np.argmin(exact)]
    print(f"\nEquilibrium (scan minimum): ~{eq:.2f} Å (experimental ≈ 0.74 Å)")


if __name__ == "__main__":
    main()
