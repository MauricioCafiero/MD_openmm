"""All tunable parameters in one place."""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Config:
    # --- Force fields (AMBER family) ---
    protein_ff: str = "amber14/protein.ff14SB.xml"
    water_ff: str = "amber14/tip3p.xml"
    # amber14/tip3p.xml already bundles the ion set (NA, CL, K, MG, CA, ...),
    # so no separate ion file is needed for neutralize/ionicStrength. Set this
    # only if you switch to a water model whose file lacks ions (e.g. opc/tip3pfb).
    ion_ff: str | None = None
    ligand_ff: str = "gaff-2.11"  # GAFFTemplateGenerator forcefield string

    # --- Solvation ---
    water_model: str = "tip3p"  # tip3p | tip3pfb | opc (must match water_ff)
    # nm around solute. A dimer/elongated solute whose longest axis nearly fills
    # the box will straddle periodic boundaries during dynamics, so the *unwrapped*
    # DCD shows it split across faces (a visualization artifact, not real motion).
    # >=1.5 nm keeps the solute clear of the walls for typical systems.
    padding: float = 1.5  # nm around solute
    ionic_strength: float = 0.15  # molar
    neutralize: bool = True
    ph: float = 7.0

    # --- Bound cofactors to KEEP and parameterize with GAFF2 ---
    # Auto-discovery (default on): every non-water / non-free-ion hetero residue is
    # kept as a GAFF2 cofactor EXCEPT a resname whose coordinates overlap the docked
    # ligand (a crystal ligand -- dropped from all copies, so both monomers of a
    # symmetric dimer lose it). SMILES is resolved built-in (KNOWN_COFACTORS) ->
    # PubChem by the PDB HETNAM name; set --cofactor RES:SMILES to override.
    # `cofactors` force-keeps those resnames (value None = auto-resolve the SMILES),
    # overriding the overlap-drop. --no-auto-cofactors reverts to explicit-only.
    cofactors: dict[str, str | None] = field(default_factory=dict)
    auto_cofactors: bool = True
    cofactor_clash_threshold: float = 5.0  # A: a hetero resname within this of the
    # docked ligand (per-resname min) is treated as a crystal ligand and dropped.

    # --- Structural metals (Zn2+, Fe2+/3+, Mg2+, Ca2+, ...) ---
    # A bound metal (within metal_contact_threshold of a protein atom) is kept in the
    # topology as a charged LJ particle typed by the water FF and locked into its
    # crystal coordination site with harmonic metal-ligand bonds + ligand-metal-ligand
    # angles (a bonded model: the coordination polyhedron is preserved, the metal
    # cannot drift). Coordinating Cys -> CYM (thiolate) and His -> the tautomer with
    # the coordinating N free, so the site protonation is correct. Force constants are
    # generic (not QM-derived) and ligand charges are the FF defaults -- for
    # publication-quality metal-center dynamics use MCPB.py (QM-derived params).
    # --no-keep-metals drops all metals instead. Metal clusters / heme / unsupported
    # resnames raise (they need custom templates).
    keep_metals: bool = True
    metal_contact_threshold: float = 2.8  # A: protein atoms within this of a metal are ligands
    metal_bond_k: float = 50000.0  # kJ/mol/nm^2  coordination bond force constant
    metal_angle_k: float = 200.0   # kJ/mol/rad^2 coordination angle force constant

    # --- System ---
    nonbonded_cutoff: float = 1.0  # nm
    constraints: str = "HBonds"  # None | HBonds | AllBonds | HAngles
    rigid_water: bool = True

    # --- Dynamics ---
    temperature: float = 300.0  # K
    # Linear temperature ramp over the *production* phase (simulated heating /
    # annealing): if set, the Langevin thermostat (and barostat, if NPT) target is
    # ramped from `temperature` to `temperature_ramp_end` across the production
    # steps. None = constant `temperature` throughout. Equilibration always runs
    # at the base `temperature`.
    temperature_ramp_end: float | None = None
    pressure: float = 1.0  # atm (NPT); set <= 0 to disable barostat (NVT)
    barostat_interval: int = 25  # steps between volume moves (MonteCarloBarostat)
    friction: float = 1.0  # 1/ps
    timestep: float = 2.0  # fs

    # --- Protocol ---
    # Minimize to a force tolerance (not a fixed iteration count): a solvated
    # ~150k-particle system needs far more than ~1000 L-BFGS steps to converge,
    # and leftover water clashes blow up the first 2 fs dynamics steps. 0 = no
    # iteration cap (run until `minimize_tolerance` is met).
    minimize_max_iter: int = 0
    minimize_tolerance: float = 10.0  # kJ/mol/nm force tolerance
    equilibrate_steps: int = 5000
    production_steps: int = 5_000_000  # 10 ns at 2 fs
    report_interval: int = 500  # steps between StateDataReporter rows
    traj_interval: int = 500  # steps between trajectory frames
    checkpoint_interval: int = 2500

    # --- Equilibration restraints ---
    restrain_protein: bool = True
    restraint_weight: float = 100.0  # kJ/mol/nm^2 on protein heavy atoms

    # --- Platform ---
    platform: str = "auto"  # auto | OpenCL | CPU | Reference
    # Apple's deprecated OpenCL GPU device lacks double precision, so only
    # "single" works here; "mixed"/"double" fall back to CPU.
    opencl_precision: str = "single"  # single | mixed | double