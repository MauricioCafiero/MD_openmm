# openmm-md

An OpenMM molecular-dynamics pipeline for **two kinds of solvated system**, both
using the AMBER force-field family:

- **Protein/ligand complex** — ff14SB protein + TIP3P water + GAFF2 ligand (and any
  bound cofactors, auto-parameterized from structure). Built and tested on a
  SULT1A3 (PDB **2A3R**) dimer + PAP cofactor + docked ligand.
- **Single small molecule** — the molecule alone in a TIP3P + salt box, GAFF2-typed.
  Useful for cofactor/fragment conformational sampling, ligand hydration, or a fast
  end-to-end check that needs no protein input.

The five stages — **prep → build → run → analyze** (+ optional protein prep and
cofactor handling) — are each a separate `omd` subcommand so you can inspect
intermediates. All tunable parameters live in one dataclass (`src/openmm_md/config.py`).

---

## What each stage does

### Protein prep — `omd prep-protein` (`prepare_protein.py`)

Repairs a raw protein PDB with **PDBFixer**:

1. `findMissingResidues()` — finds gaps in the chain.
2. By default (`--add-missing-residues` off) it **clears** `missingResidues`, so only
   missing *atoms* are patched, not whole invented terminal loops (what you want for a
   structured binding-site protein). Pass `--add-missing-residues` to model loops too.
3. `findMissingAtoms()` → `addMissingAtoms()` — patches missing heavy atoms.
4. `addMissingHydrogens(pH=ph)` — adds H at the target pH (default 7.0).
5. Writes a clean PDB (keepIds=True preserves chain/residue numbering).

`build` can also take a raw PDB directly (it calls `addHydrogens` itself), so
`prep-protein` is optional — use it when the structure has missing heavy atoms or
you want pH-controlled protonation before building.

### Ligand prep — `omd prep-ligand` (`prepare_ligand.py`)

Two input modes:

- **SDF passthrough** (`--sdf`): a pre-docked pose is supplied. We read it without
  sanitizing, clear `NoImplicit` on each atom (some exporters — notably OpenBabel —
  write MDL valence/H-count fields that set `NoImplicit` *and* omit the explicit H,
  which makes RDKit perceive radicals instead of implicit H and rejects the mol in
  OpenFF), recompute the property cache, then sanitize. The 3D docked pose is preserved.
- **SMILES → SDF** (`--smiles`): `AddHs` → embed 5 conformers (ETKDGv3, fixed seed) →
  MMFF-optimize each (threaded) → keep the lowest-energy conformer → write SDF.

`--site "x y z"` (Å) translates a SMILES-built ligand so its centroid sits at a
coordinate (used when there's no pre-docked pose). Returns the chosen conformer's
MMFF energy.

### Build — `omd build` (complex) / `omd build-mol` (single molecule) (`build_system.py`)

Both paths share one tail (`_solvate_and_create`): register GAFF templates →
(optionally) add H → solvate → `createSystem` (PME) → write `system.xml` +
`complex.pdb` + `ligand.smi`.

**Complex (`build`):**

1. Read the protein PDB (`app.PDBFile`), grab its Å positions.
2. **Docked ligand**: load the SDF, optionally translate to `--site`, make H explicit
   (coords from the pose), convert to an OpenFF `Molecule` + OpenMM `Topology` + nm
   positions. The ligand residue is stamped **`LIG`** (OpenFF leaves small molecules
   `UNK`, which mdtraj misreads as an unknown *amino acid* — naming it `LIG` keeps it
   out of the protein selection in analysis).
3. **Discover cofactors vs. crystal ligands** (`discover_cofactors`; on by default,
   `--no-auto-cofactors` to opt out). Parse the PDB **HETNAM** records for each hetero
   resname's chemical name. For every hetero residue that is not water, not a standard
   free ion (NA/CL/K/MG/CA/ZN/…, handled by `addSolvent`/neutralize), and not protein,
   decide per **resname**:
   - compute the minimum heavy-atom distance from any residue of that resname to the
     docked ligand;
   - if that minimum is below `cofactor_clash_threshold` (default **5.0 Å**), the
     resname is a **crystal ligand** the docked pose displaces → **DROP all** residues
     of that resname (so a ligand present in both monomers of a symmetric dimer is
     dropped from both, even though the docked ligand only overlaps one site);
   - otherwise **KEEP** it as a cofactor and resolve a SMILES: built-in
     (`KNOWN_COFACTORS`: **A3P**/**PAP** = adenosine 3',5'-diphosphate, neutral
     phosphate-protonated form) → **PubChem** by the HETNAM name (`pubchempy`).
     If neither resolves, the build errors and asks for `--cofactor RES:SMILES`.
   A discovery table (resname | HETNAM name | min dist | KEEP/DROP) is printed.
   `--cofactor RES:SMILES` (or bare `--cofactor RES` to force-keep with an
   auto-resolved SMILES) overrides the overlap-drop for that resname.
4. **Parameterize kept cofactors**: for each kept cofactor residue, extract its
   heavy-atom graph + CONECT bonds + PDB coords as an RDKit mol, match it against the
   resolved SMILES template (`AssignBondOrdersFromTemplate`) to recover bond orders
   and formal charges while keeping the bound 3D pose, add H, and re-add it as a
   standalone residue stamped with its PDB resname. **Why re-add:** the
   `GAFFTemplateGenerator` keys registered molecules by full molecular formula
   *including H*, but a PDB residue carries no H — so a heavy-atom-only residue can
   never match a hydrogen-full registered molecule. Giving the cofactor its own
   H-bearing residue (exactly like the docked ligand) makes the formulas agree. One
   unique `OFFMol` per distinct SMILES is registered.
5. **Remove other hetero**: water, free ions, and dropped crystal ligands are deleted
   from the protein topology.
6. **Structural metals** (`discover_hetero` runs in the same pass as cofactors; on by
   default, `--no-keep-metals` to drop all metals instead). A monoatomic metal
   (`SUPPORTED_METAL_RESNAMES`: ZN/FE/FE2/MG/CA/MN/CU/NI/CO/CD/HG) within
   `metal_contact_threshold` (default **2.8 Å**) of a protein heavy atom is a **bound
   structural metal** — kept (not dropped) and **locked into its crystal coordination
   site** with a bonded model:
   - the metal stays as a **+2 charged LJ particle typed by the water FF** (`tip3p.xml`
     carries the ion templates), at its crystal position;
   - **harmonic metal–ligand bonds** (k = `metal_bond_k` = 50000 kJ/mol/nm²) at the
     crystal metal–ligand distances, plus **ligand–metal–ligand angles** (k =
     `metal_angle_k` = 200 kJ/mol/rad²) at the crystal geometry, so the coordination
     polyhedron is preserved and the metal cannot drift;
   - **nonbonded exclusions** on the 1-2 (metal–ligand) and 1-3 (ligand–ligand) pairs
     (chargeProd = 0, epsilon = 0) — without these, the +2/−1 electrostatics collapse
     the site;
   - **re-protonation of the coordinating residues**: Cys → **CYM** (thiolate S⁻, no
     HG) and His → the tautomer with the coordinating N free (NE2-coordinating → **HID**,
     ND1-coordinating → **HIE**). Protonation is applied *after* `addHydrogens` (force the
     His tautomer via `variants=`, then strip HG + rename Cys→CYM) — pre-renaming breaks
     `addHydrogens`, which can't template a H-deficient CYM.
   A metal with no protein contact within the threshold is a free ion (left to
   `addSolvent`/neutralize). **Unsupported cases raise** (they need custom templates):
   multi-atom metal clusters, heme, and any metal resname not in
   `SUPPORTED_METAL_RESNAMES`. The bonded model uses **generic force constants, not
   QM-derived ones**, and ligand charges are the FF defaults — good enough to keep a
   structural site intact through dynamics; for publication-quality metal-center
   dynamics use **MCPB.py** (QM-derived parameters).
7. Build a `Modeller` from the cleaned protein topology, re-add the cofactor(s), then
   add the ligand.
8. Force field = `amber14/protein.ff14SB.xml` + `amber14/tip3p.xml` (tip3p.xml also
   bundles the ion set, so no separate ion file). Register a `GAFFTemplateGenerator`
   for the ligand + cofactor OFFMols.
9. `addHydrogens(ff, pH)` (adds protein H, with metal-coordinating His tautomers forced
   and coordinating Cys left for the post-`addHydrogens` CYM deprotonation), then
   `addSolvent` (TIP3P, `padding` nm, 0.15 M NaCl, neutralize).
10. `createSystem` — PME, `nonbonded_cutoff` (1.0 nm), `constraints=HBonds`,
    `rigidWater=True`. If metals are present, the metal–ligand bonds/angles + nonbonded
    exclusions are added to the `System` here. Write `system.xml`, `complex.pdb`,
    `ligand.smi`.

**Single molecule (`build-mol`):**

1. Load the ligand SDF (or build it from `--smiles` first), optionally translate to
   `--site`, make H explicit, convert to OFFMol + topology + nm positions, stamped
   **`MOL`**.
2. `Modeller` from just the ligand topology (no protein). `addSolvent` defines the box.
3. Force field = **only `amber14/tip3p.xml`** — no protein FF is loaded; the solute is
   typed entirely by the GAFF template generator, water/ions come from the water FF.
   `addHydrogens` is **skipped** (the SDF already carries explicit H).
4. Same `addSolvent` + `createSystem` + write outputs as the complex path.

### Run — `omd run` (`dynamics.py`)

`minimize → equilibrate (restraints on) → release restraints → production`, with
platform probing and a hybrid-platform path for this Mac.

1. Deserialize `system.xml`, read positions from the topology PDB.
2. Add an NPT `MonteCarloBarostat` (1 atm, 300 K, every 25 steps) if `pressure > 0`
   (set `pressure <= 0` for NVT).
3. **Positional restraints** on protein heavy atoms (`PROTEIN_RES`, non-H) via a
   `CustomExternalForce` with a **global** `k` parameter (so it can be zeroed for
   production without rebuilding the system). **Auto-disable:** if there are no
   protein heavy atoms (a single-molecule system, or `restrain_protein=False`), no
   restraint force is added and the hybrid path is skipped — so `run` works unchanged
   on a protein-free system with the default config.
4. **Platform probing** (`get_platform("auto")`): try OpenCL → CPU → Reference, each
   verified by actually creating a 1-particle context (a platform can be "available"
   yet fail on a real PME system). Precision is applied per-context, not as a default
   (setting an unsupported Precision default poisons the OpenCL singleton here).
5. **Hybrid platform (this Mac).** Apple's deprecated OpenCL GPU device is
   single-precision only, and single precision cannot stably run the *restrained*
   minimize+equilibration of a large solvated PME system (L-BFGS sticks at a
   catastrophic energy; the first restrained steps NaN). CPU (double) handles it, and
   OpenCL single is stable for production once restraints are released. So when
   restraints are on and the production platform is single-precision OpenCL, `run`
   **equilibrates on CPU (double)**, then transfers the equilibrated state (positions,
   velocities, periodic box) to a fresh **OpenCL** context and produces with `k=0`.
   You'll see `equilibrating on CPU, producing on OpenCL` in the log. With no
   restraints (single molecule), it runs one platform throughout.
6. **Minimize** to a **force tolerance** (`minimize_tolerance`, default 10 kJ/mol/nm),
   not a fixed iteration count — a solvated ~150k-particle system needs far more than
   ~1000 L-BFGS steps, and leftover water clashes blow up the first 2 fs dynamics
   steps. `minimize_max_iter=0` means no iteration cap.
7. `setVelocitiesToTemperature`, equilibrate `equilibrate_steps` (default 5000) with
   restraints on.
8. Production: zero `k`, attach reporters — `StateDataReporter` (energy.csv, every
   `report_interval`), `DCReporter` (traj.dcd, every `traj_interval`),
   `CheckpointReporter` (checkpoint.chk, every `checkpoint_interval`). Run `steps`
   (default 5,000,000 = 10 ns at 2 fs).
9. Write `final.pdb`. Returns the traj path.

### Analyze — `omd analyze` (`analyze.py`)

Loads the DCD against the topology PDB, then branches on whether a protein is present.

**Both paths:** reassemble molecules across PBC (`image_molecules`), keep the solute
only (drop water/ions), center it, and write a **wrapped, centered, solute-only**
viewing trajectory: `traj_wrapped.xtc` + a matching first-frame `traj_wrapped.pdb`
topology. (The raw `traj.dcd` is *unwrapped* — a solute that drifts across a periodic
face paints itself split across opposite box walls; that's a visualization artifact,
not real motion. Load the wrapped xtc with the wrapped pdb in PyMOL/VMD.)

**Protein path:** protein-CA RMSD (aligned on CA); ligand RMSD in the protein-aligned
frame (align on protein heavy, measure the `LIG` atoms); per-residue CA RMSF; energy
plot. Writes `rmsd.csv` (`protein_ca_rmsd_nm`, `ligand_rmsd_nm`), `rmsf.csv`,
`analysis.png`.

**Single-molecule path** (no protein CA): whole-solute heavy-atom RMSD (aligned on the
heavy atoms); per-atom RMSF over the solute heavy atoms; **radius of gyration** of the
solute (mass-weighted, a compactness/conformational proxy — computed on a solute-only
slice since Rg over the full box is meaningless); energy plot. Writes `rmsd.csv`
(`solute_rmsd_nm`, `solute_rg_nm`), `rmsf.csv`, `analysis.png`.

---

## Environment setup (conda-forge via Miniforge)

The OpenFF/GAFF ligand tooling is **not cleanly installable from PyPI** (only yanked,
broken-metadata uploads exist), so this project uses a conda-forge environment.
Miniforge is installed at `~/miniforge3` with base auto-activation **off** (activate
explicitly).

```sh
# one-time: create the env from environment.yml (installs all deps + this package)
conda env create -f environment.yml
conda activate openmm-md
```

`environment.yml` pins `python=3.11` + `openmm`, `openff-toolkit`, `openmmforcefields`,
`rdkit`, `pdbfixer`, `mdtraj`, `numpy/scipy/pandas/matplotlib`, and installs this
package (`-e .`) via pip **without** pulling PyPI deps (deps come from conda).

Re-install the package after editing code (deps already come from conda, so `--no-deps`
keeps pip from touching them):

```sh
pip install -e . --no-deps
```

### Fixes you need on this Apple Silicon mac

- **libomp is keg-only.** Homebrew installs `libomp` but it is NOT on the dyld default
  path (`ctypes.util.find_library('omp')` returns `None`). Add to `~/.zshrc`:
  ```sh
  export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/opt/libomp/lib:/lib:/usr/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
  export OMP_NUM_THREADS=4
  ```
  This makes libomp *loadable by name*; only code built against `-fopenmp` gets a
  speedup (most pip/conda wheels are not OMP-linked, so don't expect speedups unless
  you compile something against OpenMP).
- **No MPS/Metal platform.** OpenMM on Apple Silicon exposes `Reference`, `CPU`, and
  `OpenCL` only. The OpenCL GPU device is **single-precision only** (`mixed`/`double`
  fall back to CPU) — see the hybrid-platform note under *Run* above. This is why the
  production run equilibrates on CPU and produces on OpenCL.
- **No separate ion XML — ions live in the water file.** This conda `openmm` build
  ships `amber14/protein.ff14SB.xml` and the water models, but **no `amber14/*-ions.xml`
  files**. An earlier config tried to glob `amber14/*-ions.xml` for
  `neutralize=True` + 0.15 M salt and hit
  `ValueError: Could not locate file "amber14/*-ions.xml"`. The fix is **not** to fetch
  anything: the ion templates (NA, CL, K, MG, CA, …) are **already inside
  `amber14/tip3p.xml`** alongside HOH. So the build loads `protein.ff14SB.xml` +
  `tip3p.xml` only, and `config.ion_ff = None`. Set `ion_ff` only if you switch to a
  water model whose XML lacks ions (e.g. `amber14/opc.xml` for OPC) — and even then the
  file ships with openmm, no download needed.

Disk: the conda env is ~1.5–2 GB (openmm + openff-toolkit + rdkit + mdtraj +
scipy/matplotlib + their native libs). Miniforge base is ~500 MB.

---

## Running from the command line

### Protein/ligand complex

```sh
# 1. (optional) protein prep
omd prep-protein --pdb data/protein/3pose.pdb --out outputs/protein_prep.pdb

# 2. ligand prep (one of):
omd prep-ligand --sdf data/ligands/pose.sdf --out outputs/ligand.sdf        # pre-docked passthrough
omd prep-ligand --smiles "CC(=O)Nc1ccc(O)cc1" --out outputs/ligand.sdf      # build from SMILES

# 3. build solvated system. Cofactors are auto-discovered: every non-water/non-ion
#    hetero residue is kept + GAFF2-parameterized EXCEPT a resname overlapping the
#    docked ligand (a crystal ligand, dropped from all copies). SMILES resolved
#    built-in (A3P/PAP) -> PubChem by the PDB HETNAM name. Bound structural metals
#    (Zn2+/Fe2+/...) are kept and locked in their crystal coordination site (bonded
#    model). --cofactor RES:SMILES overrides; --no-auto-cofactors / --no-keep-metals
#    opt out of either.
omd build --protein outputs/protein_prep.pdb --ligand outputs/ligand.sdf \
          --site "12.3  4.5  6.7" --out-dir work
#    (offline / full control: add --no-auto-cofactors --cofactor A3P)

# 4. run MD (positions read from the topology PDB; --steps in MD steps)
omd run --system work/system.xml --topology work/complex.pdb \
        --steps 5000000 --out-dir work

# 5. analyze (RMSD/RMSF/energy plots + wrapped viewing trajectory)
omd analyze --traj work/traj.dcd --topology work/complex.pdb --out-dir work
```

### Single small molecule

```sh
omd build-mol --smiles "c1ccccc1O" --out-dir workmol --padding 1.5     # from SMILES
# or: omd build-mol --ligand data/ligands/pose.sdf --out-dir workmol   # pre-made SDF
omd run --system workmol/system.xml --topology workmol/complex.pdb \
        --steps 5000000 --out-dir workmol
omd analyze --traj workmol/traj.dcd --topology workmol/complex.pdb --out-dir workmol
```

`run` auto-detects the protein-free system and disables restraints + the hybrid path,
so the single-molecule path needs no special flags. **Box size for NPT:** a single
small molecule in a tight box can NPT-crash with `The periodic box size has decreased
to less than twice the nonbonded cutoff` — with `nonbonded_cutoff=1.0` nm the box must
stay above 2.0 nm/side, and a ~230-water box (padding 1.0) is small enough that NPT
pressure noise lets the barostat overshoot. At the **default `padding=1.5`** the box is
~3.2 nm and NPT is stable (settles to ~0.99 g/mL). Use the default padding for NPT
single-molecule runs; only drop to padding 1.0 with `pressure<=0` (NVT) for tiny checks.

### Long runs: launch fully detached

Multi-hour background jobs on this laptop die ~30–50 min in if not **fully detached**
from the Claude session. `nohup` alone is not enough (it only blocks SIGHUP, not the
harness's reaped-task SIGTERM). Put `caffeinate -w $$ &` at the top of the run script
(holds off macOS idle sleep) and launch it detached with a one-line
`subprocess.Popen(..., start_new_session=True)` (macOS has no `setsid` binary; this
calls `setsid()` in the child → new session, reparented to launchd):

```sh
~/miniforge3/bin/python -c "import subprocess; subprocess.Popen(['zsh','./run.sh'], stdout=open('out.log','ab'), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True, cwd='$(pwd)')"
```

Verify with `ps -ax -o pid,ppid,sess,command | grep run.sh` — expect PPID 1. Monitor by
tailing the log / `energy.csv`; do not sleep-poll.

---

## Artifacts

After `build`/`run`/`analyze`, the output dir contains:

| file | produced by | what it is |
|---|---|---|
| `system.xml` | build | serialized OpenMM `System` |
| `complex.pdb` | build | solvated topology + starting coords (topology for traj) |
| `ligand.smi` | build | docked ligand / molecule SMILES record |
| `energy.csv` | run | step/time/energy/temperature/volume/density |
| `traj.dcd` | run | **unwrapped** production trajectory |
| `checkpoint.chk` | run | restart checkpoint |
| `final.pdb` | run | final coordinates |
| `rmsd.csv`, `rmsf.csv`, `analysis.png` | analyze | RMSD/RMSF tables + plots |
| `traj_wrapped.xtc` | analyze | wrapped, centered, **solute-only** trajectory |
| `traj_wrapped.pdb` | analyze | matching topology (first frame) for the xtc |

## Viewing the trajectory in PyMOL

`traj.dcd` is **unwrapped** — a solute that drifts across a periodic face appears split
across opposite box walls (a dimer looks like it "separates", a ligand looks like it
"flies out"). This is a **visualization artifact**, not real dynamics. `analyze` writes
a wrapped, centered, **solute-only** trajectory (`traj_wrapped.xtc` + `traj_wrapped.pdb`).
Load that pair — never the unwrapped `traj.dcd` against the full `complex.pdb`.

### Quick start (SULT1A3 run)

From the repo root, in PyMOL's command line (or a `.pml` script):

```pymol
# 1. load the solute-only topology, then the trajectory into the same object
load outputs/prodtest/traj_wrapped.pdb, complex
load outputs/prodtest/traj_wrapped.xtc, complex

# 2. representations: protein cartoon, ligand + cofactors as sticks
hide everything, all
show cartoon, polymer
color grey80, polymer
show sticks, resn LIG           # docked ligand (stamped LIG by the build)
color cyan, resn LIG
show sticks, resn A3P            # PAP cofactors (one per monomer)
color magenta, resn A3P
show spheres, resn ZN and name ZN   # structural metal (zinc-finger run only)
color orange, resn ZN

# 3. play / scrub
mplay                           # loop the trajectory
mstop                           # stop
```

### Tips

- **Do not load `traj_wrapped.xtc` against the full `complex.pdb`** (152k atoms incl.
  solvent) — a solute-only xtc needs the solute-only `traj_wrapped.pdb` topology or
  mdtraj/PyMOL reject the atom-count mismatch. The `.pdb` and `.xtc` from `analyze` are a
  matched pair; load them into the **same** object name (`complex` above).
- **Residue names the build stamps:** **`LIG`** (docked ligand), the cofactor's PDB
  resname (e.g. **`A3P`**), **`ZN`** (kept metal). Waters/ions are already stripped from
  the wrapped trajectory. Use these in selections: `select lig, resn LIG`.
- **Measure a distance over time** (e.g. Zn–thiolate in the zinc-finger run):
  ```pymol
  distance zn_sg, resn ZN and name ZN, resn CYM and name SG
  ```
  ~2.3 Å throughout confirms the metal stayed coordinated (the `test_zinc_finger_md`
  assertion checks the same thing on the last frame).
- **Export a movie:**
  ```pymol
  set ray_trace_frames, 0        # 1 = higher quality, much slower
  mpng frames                    # writes frames/*.png; stitch with ffmpeg
  ```
- **Save the session** for later: `save outputs/prodtest/view.pse`.

---

## Example: SULT1A3 dimer + PAP + ligand (~1 ns production test)

Inputs (kept in git): `data/protein/sult1a3_2A3R.pdb` (the dimer, with bound **A3P** =
PAP cofactors and a crystal `LDP` = L-dopamine ligand) and `data/ligands/sult1a3_2A3R_c0.sdf`
(an OpenBabel-exported docked pose).

**Build** (`omd build --protein … --ligand …` — cofactors auto-discovered, no
`--cofactor` needed):
- Discovery table printed: `A3P` (ADENOSINE-3'-5'-DIPHOSPHATE) min 9.17 Å → **KEEP**
  (built-in SMILES); `LDP` (L-DOPAMINE) min 3.37 Å → **DROP** (overlaps the docked
  ligand; dropped from both monomers by the per-resname rule).
- 2 A3P cofactors kept + GAFF2-parameterized (one per monomer); the crystal `LDP`
  (both copies) and waters were stripped.
- **175,674 particles, 55,314 waters.**

**Run** (`omd run --steps 500000`, 1.0 ns at 2 fs): hybrid path — equilibrated on CPU
(double), produced on OpenCL (single). Wall-clock ~4h18m (08:24:57 → 12:43:12 BST).
Stable throughout — no NaN in 1000 energy rows.

| metric | final value |
|---|---|
| steps / time | 500,000 / 1000.0 ps (1.0 ns) |
| final temperature | 300.53 K |
| final density | 1.0082 g/mL |
| protein CA RMSD (vs frame 0) | 0.364 nm |
| ligand RMSD (vs frame 0, protein-aligned) | 0.273 nm |
| total-energy NaN | none |

**Analyze** wrote `analysis.png`, `rmsd.csv`, `rmsf.csv`, and the wrapped/centered
solute-only viewing trajectory: **`outputs/prodtest/traj_wrapped.xtc`** (1000 frames,
9420 solute atoms) + `traj_wrapped.pdb` topology. Load that pair in PyMOL/VMD (not the
unwrapped `traj.dcd`, which shows the dimer split across periodic faces — a
visualization artifact).

---

## Cofactor robustness

The SULT1A3 build above is the worked case of the cofactor fix: the pipeline parses the
PDB **HETNAM** records, and for each non-water/non-ion hetero resname decides per
**resname** (not per residue) whether it's a bound cofactor to keep or a crystal ligand
the docked pose displaces. A resname whose minimum heavy-atom distance to the docked
ligand is below `cofactor_clash_threshold` (default **5.0 Å**) is dropped from **all**
copies — so L-dopamine (`LDP`) is removed from **both** monomers of the symmetric dimer
even though the docked ligand only overlaps one site — while `A3P`/PAP (9.17 Å away) is
kept and GAFF2-parameterized. SMILES resolve built-in (`KNOWN_COFACTORS`) → PubChem by
the HETNAM name; `--cofactor RES:SMILES` overrides and `--no-auto-cofactors` opts out.
Mechanism details are in the Build section (steps 3–4); the discovery table for
SULT1A3 is shown in the example above.

---

## Example: 1ZNF zinc finger + ligand (metal robustness)

A minimal end-to-end demo of the **structural-metal** robustness fix on a single small
system: **PDB 1ZNF** (a Cys2His2 zinc-finger domain) with a **bound Zn²⁺** and a small
proxy ligand (phenol) placed beside it. 1ZNF is **apo** (no bound small molecule), so
the ligand here is a stand-in placed at a defined site — the point is to exercise the
structural-metal bonded model and confirm the metal stays coordinated through dynamics.

Inputs (kept in git): `data/protein/1znf_zinc_finger.pdb` (1ZNF model 1, all H
stripped — 212 heavy atoms, one ZN at chain A resSeq 27, coordinated by CYS3/SG,
CYS6/SG, HIS19/NE2, HIS23/NE2). Build + run + analyze in one shot:

```sh
# build a phenol ligand, place it just outside the finger, build the solvated complex
# (Zn2+ auto-kept + locked into its Cys2His2 site), run 300 steps on CPU, analyze.
omd prep-ligand --smiles "c1ccccc1O" --site "17 -0.2 -0.2" --out outputs/zf_lig.sdf
omd build --protein data/protein/1znf_zinc_finger.pdb --ligand outputs/zf_lig.sdf \
          --site "17 -0.2 -0.2" --out-dir outputs/zf --padding 1.0
omd run --system outputs/zf/system.xml --topology outputs/zf/complex.pdb \
        --steps 300 --platform CPU --out-dir outputs/zf
omd analyze --traj outputs/zf/traj.dcd --topology outputs/zf/complex.pdb --out-dir outputs/zf
```

**Build** — `discover_hetero` finds the Zn²⁺ within 2.8 Å of four protein ligands and
classifies it KEEP (bonded model). The build log:

```
[build] hetero discovery: ZN  Zn  4-coordinate  KEEP (bonded model; typed by water FF)
        ligands: CYS3/SG, CYS6/SG, HIS19/NE2, HIS23/NE2
[build] protonation: A/3 -> CYM, A/6 -> CYM (metal-coordinating)
[build] protonation: A/19 -> HID, A/23 -> HID (metal-coordinating)
[build] deprotonated 2 coordinating CYS -> CYM (thiolate)
[build] metal coordination: 4 bonds + 6 angles + 10 nonbonded exclusions
[build] 8162 particles, 2569 waters
```

**Run + analyze** (300-step CPU smoke): protein CA RMSD **0.050 nm**, ligand RMSD
**0.186 nm**. The robustness check — Zn–ligand distances in the **last frame** vs.
crystal:

| ligand | crystal (Å) | last frame (Å) | range over run (Å) |
|---|---|---|---|
| CYM3/SG  | 2.29 | 2.42 | 2.20–2.49 |
| CYM6/SG  | 2.31 | 2.29 | 2.28–2.40 |
| HIS19/NE2 | 1.95 | 2.19 | 2.05–2.20 |
| HIS23/NE2 | 1.99 | 2.12 | 2.07–2.17 |

The Zn²⁺ stays in its Cys2His2 site (all four distances 2.05–2.49 Å throughout) — the
bonded model holds the metal through dynamics, no drift, no NaN. This is the
`test_zinc_finger_md` smoke test (always on; artifacts in `outputs/smoke_znmd/`).

---

## Tests

```sh
conda activate openmm-md
python tests/smoke.py             # ligand gen + platform probe + single-molecule
                                 # solvated run (no inputs needed)
# full protein/ligand pipeline once you have a protein + ligand:
OMD_PROTEIN_PDB=data/protein/x.pdb OMD_LIGAND_SDF=data/ligands/pose.sdf python tests/smoke.py
# or build the ligand from SMILES and place it at a site centroid:
OMD_PROTEIN_PDB=data/protein/x.pdb OMD_LIGAND_SMILES="CC(=O)Nc1ccc(O)cc1" \
OMD_SITE="12.3 4.5 6.7" python tests/smoke.py
# keep+GAFF a cofactor (built-in name) in the full pipeline:
OMD_COFACTORS=A3P OMD_PROTEIN_PDB=data/protein/x.pdb OMD_LIGAND_SDF=data/ligands/pose.sdf python tests/smoke.py
```

The always-on smoke tests need no inputs:
- **ligand generation** from SMILES (ethanol → SDF);
- a **platform probe** (auto);
- a **single-molecule** solvated run (paracetamol from SMILES → solvate → short NVT MD
  → analyze);
- **auto-cofactor discovery** (offline, on the SULT1A3 fixture: A3P kept, LDP dropped
  from both monomers by the per-resname overlap rule);
- **structural metal** (build-only on the 1ZNF fixture: Zn²⁺ kept +2, 4 bonds + 6
  angles + 10 exclusions, 2 Cys→CYM + 2 His→HID);
- **zinc-finger + ligand MD** (full build → 300-step run → analyze on 1ZNF + a placed
  phenol, asserting the Zn²⁺ stays coordinated ~2.3 Å through dynamics).

The full protein/ligand pipeline runs when you set `OMD_PROTEIN_PDB`
+ `OMD_LIGAND_SDF`/`OMD_LIGAND_SMILES`. Artifacts persist to `outputs/smoke/`,
`outputs/smoke_mol/`, `outputs/smoke_metal/`, `outputs/smoke_znmd/` (override with
`OMD_OUT_DIR` / `OMD_MOL_OUT_DIR` / `OMD_METAL_OUT_DIR` / `OMD_ZNMD_OUT_DIR`).

## Layout

```
openmm/
├── README.md
├── .gitignore
├── environment.yml          # conda-forge env (the supported install path)
├── pyproject.toml           # package metadata + `omd` console script
├── data/                    # INPUTS (kept in git)
│   ├── protein/             #   raw protein PDBs
│   └── ligands/             #   ligand SDFs or a SMILES list
├── outputs/                 # GENERATED (git-ignored)
├── src/openmm_md/           # CODE
│   ├── config.py            #   all tunable parameters (one dataclass)
│   ├── prepare_protein.py   #   PDBFixer pipeline
│   ├── prepare_ligand.py    #   SMILES->conformers->MMFF->SDF, SDF passthrough, centroid placement
│   ├── build_system.py      #   build_system (complex) + build_mol_system (single molecule)
│   ├── dynamics.py          #   platform probing, hybrid CPU-equil/OpenCL-produce, minimize, reporters
│   ├── analyze.py           #   mdtraj RMSD/RMSF/Rg + energy plots + wrapped trajectory
│   └── cli.py               #   argparse entrypoint (omd)
└── tests/
    └── smoke.py             # ligand+platform+single-molecule always; full pipeline with OMD_* env vars
```

## Known caveats / next steps

- **No restart-from-checkpoint.** `omd run` always starts a fresh minimize+equilibration;
  `checkpoint.chk` is written but not consumed. A mid-run kill means restarting from
  scratch (the partial production is lost). Add a `--restart` path that loads the
  checkpoint and continues production if you want resumable long runs (the Rotaxanes
  project's guarded/resumable stages are a model for this).
- **PAP protonation.** The A3P/PAP template is the neutral form; at physiological pH
  the phosphates are largely deprotonated. Good enough to build a system; refine if you
  need the deprotonated state.
- **Box padding.** Default `padding` is 1.5 nm (raised from 1.0) so an elongated dimer
  doesn't straddle periodic faces during dynamics. For single molecules, 1.5 nm also
  keeps the NPT box safely above 2× the nonbonded cutoff.
- **Under-thermalized smoke run.** The smoke config equilibrates only 100–200 steps,
  leaving the system under-thermalized. Fine for a smoke test; a real run should use a
  few ps of staged NVT→NPT equilibration (raise `cfg.equilibrate_steps`, default 5000).