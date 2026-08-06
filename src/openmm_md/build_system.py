"""Build the solvated protein/ligand system.

Protein (ff14SB) and small molecules (GAFF2) are typed in one ForceField via
openmmforcefields' GAFFTemplateGenerator, registered with the same OpenMM
ForceField that holds the AMBER protein force field.

Two kinds of small molecule get GAFF templates:
  * the docked ligand, read from an SDF (with bond orders / charges intact);
  * any bound cofactors kept by auto-discovery (see ``discover_cofactors``).
    Each cofactor residue's 3D heavy-atom graph is read from the protein PDB
    (CONECT bonds) and matched to a SMILES template to recover bond orders /
    formal charges while keeping the bound coordinates; H are added and the
    cofactor is re-added to the topology as a standalone residue (mirroring the
    docked ligand path). Waters, free ions, and crystal ligands are removed.

Cofactor discovery is robust to whatever hetero residues a PDB happens to carry:
every non-water / non-free-ion hetero residue is kept as a GAFF2 cofactor EXCEPT
a resname whose coordinates overlap the docked ligand (a crystal ligand -- e.g.
L-DOPAMINE in SULT1A3, which the docked ligand displaces). The SMILES is resolved
built-in (``KNOWN_COFACTORS``) -> PubChem by the PDB HETNAM name; ``--cofactor
RES:SMILES`` overrides. ``--no-auto-cofactors`` reverts to explicit-only.

Re-adding the cofactor (rather than leaving it in the protein topology) is
required because the GAFFTemplateGenerator keys registered molecules by their
full molecular formula *including H* (template_generators.py:154), while a PDB
residue carries no H -- so a heavy-atom-only residue can never match a
hydrogen-full registered molecule. Giving the cofactor its own residue with
explicit H (exactly like the docked ligand) makes the formulas agree.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from openmm import app, unit, XmlSerializer
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit import Molecule as OFFMol
from openmmforcefields.generators import GAFFTemplateGenerator

from .config import Config
from .prepare_ligand import load_ligand_rdkit, translate_to
from .dynamics import PROTEIN_RES

CONSTRAINTS = {
    "None": None,
    "HBonds": app.HBonds,
    "AllBonds": app.AllBonds,
    "HAngles": app.HAngles,
}

# Known cofactor residue-name -> SMILES, for template matching from a PDB.
# A3P is PAP (adenosine 3',5'-diphosphate). NOTE: the template encodes the
# *neutral* (phosphate-protonated) form (charge 0); at physiological pH the
# phosphates are largely deprotonated. Good enough to build a system; refine if
# you need the deprotonated state.
KNOWN_COFACTORS = {
    "A3P": "C1=NC2=C(N=CN2[C@@H]3O[C@H](COP(O)(=O)O)[C@@H](OP(O)(=O)O)[C@H]3O)C(=N1)N",
    "PAP": "C1=NC2=C(N=CN2[C@@H]3O[C@H](COP(O)(=O)O)[C@@H](OP(O)(=O)O)[C@H]3O)C(=N1)N",
}

# Standard free ions, handled by addSolvent/neutralize rather than kept + GAFF-typed.
# A bound *structural* metal (e.g. a catalytic ZN) is dropped here -- GAFF cannot
# type metals, so keeping it would fail template generation. That is a known
# limitation; override with --cofactor only for molecular (non-metal) cofactors.
FREE_ION_RESNAMES = {"NA", "CL", "K", "MG", "CA", "ZN", "FE", "CU", "MN", "BR", "F", "I", "LI", "RB", "CS", "SR", "BA"}
WATER_RESNAMES = {"HOH", "WAT", "TIP3", "SOL"}

# Free solvent ions (mobile, always dropped -- addSolvent re-adds them): alkali + halide.
SOLVENT_ION_RESNAMES = {"NA", "CL", "K", "BR", "F", "I", "LI", "RB", "CS"}
# Structural metals handled by the bonded model (monoatomic, typeable by amber14/tip3p.xml
# as a charged LJ particle, + geometric keep/drop + coordination restraints). Resname ==
# the tip3p.xml template name.
SUPPORTED_METAL_RESNAMES = {"ZN", "FE", "FE2", "MG", "CA", "MN", "CU", "NI", "CO", "CD", "HG"}
# Element symbols used to spot a metal inside a multi-atom residue (a metal cluster /
# cofactor like FES/SF4/HEM that this pipeline cannot parameterize -> error).
METAL_ELEMENTS = {
    "ZN", "FE", "MG", "CA", "MN", "CU", "NI", "CO", "CD", "HG", "MO", "W", "V", "CR",
    "AL", "AU", "AG", "PT", "PD", "RU", "OS", "IR", "RH", "RE", "U", "PB",
}
# Protonation fix for metal-coordinating residues (resname, ligating-atom name) -> FF
# residue name with the coordinating atom in its metal-binding form:
#   CYS/SG -> CYM (thiolate S-, no HG); HIS/ND1 -> HIE (H on NE2, ND1 free);
#   HIS/NE2 -> HID (H on ND1, NE2 free). ff14SB carries all these templates.
LIGAND_RENAMES = {("CYS", "SG"): "CYM", ("HIS", "ND1"): "HIE", ("HIS", "NE2"): "HID"}


def parse_hetnam(pdb_path):
    """Parse PDB HETNAM records -> {resname: full chemical name}.

    Handles multi-line continuation (no hetID on continuation lines). Columns:
    hetID at 12-14 (0-indexed 11:14), name at 16-70 (0-indexed 15:70), matching
    the wwPDB fixed-column layout as written by the SULT1A3 PDB.
    """
    het = {}
    last = None
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("HETNAM"):
                continue
            hetid = line[11:14].strip()
            name = line[15:70].strip()
            if hetid:
                het[hetid] = (het[hetid] + " " + name).strip() if hetid in het else name
                last = hetid
            elif last:
                het[last] = (het[last] + " " + name).strip()
    return het


def smiles_from_pubchem(name):
    """Resolve a chemical name to a SMILES via PubChem (pubchempy). None on any failure."""
    try:
        import pubchempy as pcp
    except Exception:
        return None
    try:
        res = pcp.get_compounds(name, "name")
        return res[0].smiles if res else None
    except Exception:
        return None


def resolve_cofactor_smiles(resname, hetnam_name):
    """Resolve a cofactor SMILES: built-in KNOWN_COFACTORS -> PubChem by HETNAM name."""
    if resname in KNOWN_COFACTORS:
        return KNOWN_COFACTORS[resname]
    if hetnam_name:
        smi = smiles_from_pubchem(hetnam_name)
        if smi:
            return smi
    raise ValueError(
        f"could not resolve SMILES for cofactor {resname!r} "
        f"(HETNAM={hetnam_name!r}); pass --cofactor {resname}:SMILES explicitly"
    )


def _ligand_heavy_coords_ang(rdmol):
    """Heavy-atom coords (angstrom) of the docked ligand, for the overlap test."""
    conf = rdmol.GetConformer(0)
    coords = []
    for i, a in enumerate(rdmol.GetAtoms()):
        if a.GetAtomicNum() == 1:
            continue
        p = conf.GetAtomPosition(i)
        coords.append([p.x, p.y, p.z])
    return np.array(coords, dtype=float)


def _residue_min_dist_to_ligand(residue, positions_ang, lig_coords_ang):
    """Min heavy-atom distance (angstrom) from a PDB residue to the docked ligand."""
    if len(lig_coords_ang) == 0:
        return float("inf")
    pts = []
    for a in residue.atoms():
        if a.element is None or a.element.symbol == "H":
            continue
        p = positions_ang[a.index]
        pts.append([float(p[0]), float(p[1]), float(p[2])])
    if not pts:
        return float("inf")
    pts = np.array(pts)
    d = np.linalg.norm(pts[:, None, :] - lig_coords_ang[None, :, :], axis=-1)
    return float(d.min())


def _find_metal_ligands(metal_pos_ang, topology, positions_ang, threshold):
    """Protein heavy atoms within `threshold` (angstrom) of a metal -> its coordinating
    ligands (the first coordination shell). Each entry records the atom's identity
    (chain, resSeq, resname, name) + position so coordination restraints can be added
    after solvation and the coordinating residues can be re-protonated."""
    lig = []
    for res in topology.residues():
        if res.name not in PROTEIN_RES:
            continue
        for a in res.atoms():
            if a.element is None or a.element.symbol == "H":
                continue
            p = positions_ang[a.index]
            pos = np.array([float(p[0]), float(p[1]), float(p[2])])
            d = float(np.linalg.norm(pos - metal_pos_ang))
            if d <= threshold:
                lig.append({
                    "chain": res.chain.id, "resSeq": res.id, "resname": res.name,
                    "name": a.name, "pos_ang": pos.tolist(), "dist": d,
                })
    return lig


def _add_coordination_restraints(system, topology, metals, cfg):
    """Lock each kept metal into its crystal coordination site: harmonic metal-ligand
    bonds + ligand-metal-ligand angles at the crystal geometry, plus the nonbonded
    exclusions a bonded model requires (1-2 metal-ligand and 1-3 ligand-ligand pairs
    are zeroed -- otherwise the +2/-1 electrostatics collapse the site onto the metal).

    Atom indices are re-located in the *built* topology (post addHydrogens/addSolvent)
    by (chain, resSeq, name); equilibrium distances/angles come from the crystal
    positions recorded at discovery time.
    """
    from openmm import HarmonicBondForce, HarmonicAngleForce

    def _idx(chain, resseq, name):
        for r in topology.residues():
            if r.chain.id == chain and r.id == resseq:
                for a in r.atoms():
                    if a.name == name:
                        return a.index
        raise ValueError(f"coordination: atom {name} (res {resseq}, chain {chain}) not in built topology")

    def _metal_idx(chain, resseq, resname):
        for r in topology.residues():
            if r.chain.id == chain and r.id == resseq and r.name == resname:
                return list(r.atoms())[0].index
        raise ValueError(f"coordination: metal {resname} (res {resseq}, chain {chain}) not in built topology")

    nb = next(f for f in system.getForces() if f.__class__.__name__ == "NonbondedForce")
    existing = set()
    for i in range(nb.getNumExceptions()):
        p1, p2, _, _, _ = nb.getExceptionParameters(i)
        existing.add(frozenset((p1, p2)))

    def _exclude(i, j):
        key = frozenset((i, j))
        if key not in existing:
            nb.addException(i, j, 0.0, 1.0, 0.0)  # chargeProd=0, epsilon=0 -> no nonbonded
            existing.add(key)
            return True
        return False

    bond_force = HarmonicBondForce()
    angle_force = HarmonicAngleForce()
    n_bonds = n_angles = n_excl = 0
    for m in metals:
        mi = _metal_idx(m["chain"], m["resSeq"], m["resname"])
        mpos = np.array(m["pos_ang"]) / 10.0  # A -> nm
        ligs = []
        for lig in m["ligands"]:
            li = _idx(lig["chain"], lig["resSeq"], lig["name"])
            lpos = np.array(lig["pos_ang"]) / 10.0
            r0 = float(np.linalg.norm(mpos - lpos))
            bond_force.addBond(mi, li, r0, cfg.metal_bond_k)
            n_bonds += 1
            if _exclude(mi, li):
                n_excl += 1
            ligs.append((li, lpos))
        for a in range(len(ligs)):
            for b in range(a + 1, len(ligs)):
                li, lpos_i = ligs[a]
                lj, lpos_j = ligs[b]
                v1, v2 = lpos_i - mpos, lpos_j - mpos
                n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
                cos = float(np.dot(v1, v2) / (n1 * n2)) if n1 and n2 else 0.0
                cos = max(-1.0, min(1.0, cos))
                angle_force.addAngle(li, mi, lj, float(np.arccos(cos)), cfg.metal_angle_k)
                n_angles += 1
                if _exclude(li, lj):
                    n_excl += 1
    system.addForce(bond_force)
    system.addForce(angle_force)
    print(f"[build] metal coordination: {n_bonds} bonds + {n_angles} angles + "
          f"{n_excl} nonbonded exclusions (k_bond={cfg.metal_bond_k} kJ/mol/nm^2, "
          f"k_angle={cfg.metal_angle_k} kJ/mol/rad^2)")


def discover_hetero(topology, positions_ang, lig_coords_ang, cfg, hetnam):
    """Classify every hetero residue (non-protein, non-water) as KEEP or DROP, and
    collect coordination info for any bound structural metals. Returns
    (keep_cofactors, drop, metals, renames):
      * keep_cofactors: {resname: smiles} -- organic cofactors to GAFF2-parameterize
        (per-resname overlap rule: a resname overlapping the docked ligand is dropped
        from all its residues).
      * drop: {resname: reason} -- for reporting only.
      * metals: list of dicts for bound metals to keep + lock (bonded model).
      * renames: {(chain, resSeq): new_name} -- coordinating CYS->CYM / HIS->HID|HIE.

    Solvent ions (Na/K/Cl/...) are dropped. A bound structural metal (Zn/Fe/Mg/Ca/...
    within `metal_contact_threshold` of a protein atom) is kept as a charged LJ
    particle typed by the water FF and locked into its crystal site; a free metal
    (no protein contacts) is dropped. Multi-atom metal centers (Fe-S clusters, heme,
    ...) and unsupported metal resnames raise -- they need custom templates (MCPB.py).
    """
    explicit = cfg.cofactors  # resname -> SMILES-or-None (force-keep organic cofactor)

    # pass 1: organic hetero per-resname min-d to the docked ligand + collect metals
    order = []
    min_d = {}
    metal_residues = []
    for res in topology.residues():
        rn = res.name
        if rn in PROTEIN_RES or rn in WATER_RESNAMES or rn in SOLVENT_ION_RESNAMES:
            continue
        atoms = list(res.atoms())
        if len(atoms) == 1:
            metal_residues.append(res)  # monoatomic ion (metal or otherwise)
            continue
        if any(a.element and a.element.symbol in METAL_ELEMENTS for a in atoms):
            raise ValueError(
                f"metal center/cofactor {rn!r} (HETNAM={hetnam.get(rn, '?')!r}) is "
                f"multi-atom and not supported by the auto pipeline (needs custom "
                f"templates, e.g. MCPB.py). Drop it or supply a pre-parameterized PDB."
            )
        if rn not in min_d:
            order.append(rn)
            min_d[rn] = float("inf")
        d = _residue_min_dist_to_ligand(res, positions_ang, lig_coords_ang)
        if d < min_d[rn]:
            min_d[rn] = d

    # pass 2: organic cofactors (per-resname overlap rule)
    keep, drop = {}, {}
    print("[build] hetero discovery (resname | HETNAM | min dist to ligand | action):")
    for rn in order:
        name = hetnam.get(rn, "?")
        d = min_d[rn]
        if rn in explicit:
            smi = explicit[rn]
            if smi is None:
                smi = resolve_cofactor_smiles(rn, name)
            keep[rn] = smi
            print(f"  {rn:4s} {name:32s} {d:6.2f} A  KEEP (explicit --cofactor)")
        elif cfg.auto_cofactors:
            if d < cfg.cofactor_clash_threshold:
                drop[rn] = f"overlaps docked ligand (min {d:.2f} A < {cfg.cofactor_clash_threshold} A)"
                print(f"  {rn:4s} {name:32s} {d:6.2f} A  DROP ({drop[rn]})")
            else:
                keep[rn] = resolve_cofactor_smiles(rn, name)
                print(f"  {rn:4s} {name:32s} {d:6.2f} A  KEEP -> GAFF2 ({keep[rn][:40]})")
        else:
            drop[rn] = "auto-cofactors disabled; not in --cofactor"
            print(f"  {rn:4s} {name:32s} {d:6.2f} A  DROP ({drop[rn]})")

    # pass 3: structural metals (per-residue geometric keep+lock / drop)
    metals = []
    renames = {}
    for res in metal_residues:
        rn = res.name
        atom = list(res.atoms())[0]
        el = atom.element.symbol if atom.element else "?"
        if rn not in SUPPORTED_METAL_RESNAMES:
            if el in METAL_ELEMENTS:
                raise ValueError(
                    f"unsupported metal {rn!r} (element {el}); supported resnames: "
                    f"{sorted(SUPPORTED_METAL_RESNAMES)}. Rename it or use --no-keep-metals."
                )
            drop[rn] = "unrecognized monoatomic residue (not a known ion/metal)"
            print(f"  {rn:4s} {el:3s}      --  DROP ({drop[rn]})")
            continue
        p = positions_ang[atom.index]
        mpos = np.array([float(p[0]), float(p[1]), float(p[2])])
        ligands = _find_metal_ligands(mpos, topology, positions_ang, cfg.metal_contact_threshold)
        if cfg.keep_metals and ligands:
            metals.append({
                "resname": rn, "chain": res.chain.id, "resSeq": res.id,
                "pos_ang": mpos.tolist(), "ligands": ligands,
            })
            for lig in ligands:
                key = (lig["resname"], lig["name"])
                if key in LIGAND_RENAMES:
                    renames[(lig["chain"], lig["resSeq"])] = LIGAND_RENAMES[key]
            lig_summary = ", ".join(f"{lg['resname']}{lg['resSeq']}/{lg['name']}" for lg in ligands)
            print(f"  {rn:4s} {el:3s} {len(ligands)}-coordinate  KEEP (bonded model; "
                  f"typed by water FF)  ligands: {lig_summary}")
        else:
            reason = ("--no-keep-metals" if not cfg.keep_metals
                      else f"no protein contacts within {cfg.metal_contact_threshold} A")
            drop[rn] = f"free metal ion ({reason})"
            print(f"  {rn:4s} {el:3s} 0-coordinate  DROP ({drop[rn]})")
    return keep, drop, metals, renames


def _explicit_h_rdmol_to_openmm(rdmol, residue_name="UNK"):
    """RDKit Mol (already bearing explicit H + a conformer) -> (OFFMol, Topology, nm positions).

    The topology comes from OpenFF (so it matches the GAFF template generator's
    atom perception); coordinates come from the RDKit conformer (angstrom) and
    are converted to nanometers for OpenMM. ``residue_name`` is stamped onto the
    topology residue -- OpenFF leaves it as "UNK", which mdtraj then misreads as
    an unknown *amino acid* (protein), so we name the ligand "LIG" and cofactors
    by their PDB resname to keep them out of the protein selection.
    """
    offmol = OFFMol.from_rdkit(rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    omm_topology = offmol.to_topology().to_openmm()
    for res in omm_topology.residues():
        res.name = residue_name

    conf = rdmol.GetConformer(0)
    coords = np.array(
        [[p.x, p.y, p.z] for p in (conf.GetAtomPosition(i) for i in range(rdmol.GetNumAtoms()))],
        dtype=float,
    )
    positions = unit.Quantity(coords / 10.0, unit.nanometer)  # A -> nm
    return offmol, omm_topology, positions


def _rdkit_to_openmm(rdmol, residue_name="LIG"):
    """Docked-ligand path: make H explicit (coords from the docked pose) then convert."""
    rdmol = Chem.AddHs(rdmol, addCoords=True)
    return _explicit_h_rdmol_to_openmm(rdmol, residue_name=residue_name)


def _residue_to_rdmol(topology, positions_ang, residue, template_smiles):
    """Extract a PDB residue as an RDKit Mol (with H) via template matching.

    Builds an RDKit mol from the residue's atoms + CONECT bonds + PDB coords
    (single bonds, no H), then matches it against ``template_smiles`` to recover
    bond orders and formal charges. Hydrogens are added with 3D coords placed
    off the existing geometry, preserving the bound heavy-atom pose.
    """
    atoms = list(residue.atoms())
    idxmap = {a: i for i, a in enumerate(atoms)}
    rmol = Chem.RWMol()
    conf = Chem.Conformer(len(atoms))
    for i, a in enumerate(atoms):
        rmol.AddAtom(Chem.Atom(a.element.symbol))
        p = positions_ang[a.index]
        conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))
    rmol.AddConformer(conf, assignId=True)
    for b in topology.bonds():
        if b[0] in idxmap and b[1] in idxmap:
            rmol.AddBond(idxmap[b[0]], idxmap[b[1]], Chem.BondType.SINGLE)

    template = Chem.MolFromSmiles(template_smiles)
    if template is None:
        raise ValueError(f"could not parse cofactor SMILES: {template_smiles!r}")
    matched = AllChem.AssignBondOrdersFromTemplate(template, rmol.GetMol())
    matched = Chem.AddHs(matched, addCoords=True)
    return matched


def _solvate_and_create(
    modeller, gaff_offmols, cfg: Config, ff_files, add_hydrogens: bool,
    out_dir: Path, ligand_smiles: str | None, metals=None,
):
    """Common tail: register GAFF templates, (optionally) add H, solvate, create System.

    ``ff_files`` is the ForceField file list. A protein-bearing system passes the
    protein FF + water FF (the protein needs ff14SB; tip3p.xml also provides the
    ions); a single-molecule system passes only the water FF -- the solute is
    typed by the GAFFTemplateGenerator, so no protein FF is loaded.
    ``add_hydrogens`` adds protein H (the protein path needs this; the molecule
    path's ligand already carries explicit H from the SDF, so it's skipped).
    ``metals`` (protein path only) are bound structural metals kept in the topology;
    after createSystem they are locked into their crystal coordination site with
    harmonic bonds/angles + nonbonded exclusions (bonded model).
    """
    ff = app.ForceField(*ff_files)
    if gaff_offmols:
        gaff = GAFFTemplateGenerator(molecules=gaff_offmols, forcefield=cfg.ligand_ff)
        ff.registerTemplateGenerator(gaff.generator)

    if add_hydrogens:
        # Force the tautomer of any metal-coordinating HIS (ND1 ligand -> HIE so ND1
        # is free; NE2 ligand -> HID so NE2 is free) via addHydrogens' variants arg.
        his_variants = {}
        cys_deprotonate = set()
        for m in metals or []:
            for lig in m["ligands"]:
                if lig["resname"] == "HIS" and lig["name"] == "ND1":
                    his_variants[(lig["chain"], lig["resSeq"])] = "HIE"
                elif lig["resname"] == "HIS" and lig["name"] == "NE2":
                    his_variants[(lig["chain"], lig["resSeq"])] = "HID"
                elif lig["resname"] == "CYS" and lig["name"] == "SG":
                    cys_deprotonate.add((lig["chain"], lig["resSeq"]))
        variants = None
        if his_variants:
            variants = [None] * modeller.topology.getNumResidues()
            for i, r in enumerate(modeller.topology.residues()):
                v = his_variants.get((r.chain.id, r.id))
                if v:
                    variants[i] = v
        modeller.addHydrogens(ff, pH=cfg.ph, variants=variants)
        # Deprotonate coordinating CYS -> CYM (thiolate): addHydrogens added HG (CYS
        # has no thiolate variant), so remove it and rename. ff14SB's CYM template
        # then types the S- with the metal-binding charge.
        if cys_deprotonate:
            rm = []
            for r in modeller.topology.residues():
                if (r.chain.id, r.id) in cys_deprotonate and r.name == "CYS":
                    r.name = "CYM"
                    rm.extend(a for a in r.atoms() if a.name == "HG")
            if rm:
                modeller.delete(rm)
                print(f"[build] deprotonated {len(rm)} coordinating CYS -> CYM (thiolate)")
    modeller.addSolvent(
        ff,
        model=cfg.water_model,
        padding=cfg.padding * unit.nanometer,
        ionicStrength=cfg.ionic_strength * unit.molar,
        neutralize=cfg.neutralize,
    )

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=cfg.nonbonded_cutoff * unit.nanometer,
        constraints=CONSTRAINTS[cfg.constraints],
        rigidWater=cfg.rigid_water,
    )
    if metals:
        _add_coordination_restraints(system, modeller.topology, metals, cfg)

    with open(out_dir / "complex.pdb", "w") as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)
    with open(out_dir / "system.xml", "w") as f:
        f.write(XmlSerializer.serialize(system))
    if ligand_smiles is not None:
        (out_dir / "ligand.smi").write_text(ligand_smiles + "\n")

    n = system.getNumParticles()
    n_wat = sum(1 for r in modeller.topology.residues() if r.name in ("HOH", "WAT"))
    print(f"[build] {n} particles, {n_wat} waters -> {out_dir/'system.xml'}, "
          f"{out_dir/'complex.pdb'}")
    return out_dir / "system.xml", out_dir / "complex.pdb"


def build_system(
    protein_pdb: Path,
    ligand_sdf: Path,
    out_dir: Path,
    cfg: Config | None = None,
    site_xyz_ang=None,
):
    """Combine protein + ligand (+ kept cofactors), solvate, write system.xml + complex.pdb."""
    cfg = cfg or Config()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb = app.PDBFile(str(protein_pdb))
    pos_ang = pdb.positions.value_in_unit(unit.angstrom)

    # --- docked ligand (from SDF, with bond orders intact) ---
    rdmol = load_ligand_rdkit(ligand_sdf)
    if site_xyz_ang is not None:
        rdmol = translate_to(rdmol, site_xyz_ang)
        print(f"[build] placed ligand centroid at {list(site_xyz_ang)} A")
    lig_heavy_ang = _ligand_heavy_coords_ang(rdmol)
    lig_offmol, lig_topology, lig_positions = _rdkit_to_openmm(rdmol)

    # --- discover cofactors, crystal ligands to drop, and bound metals ---
    # HETNAM gives the chemical name for each hetero resname, used to resolve a
    # SMILES via PubChem when the resname isn't a built-in. Auto mode keeps every
    # non-water/non-ion organic hetero residue as a GAFF2 cofactor EXCEPT resnames
    # whose coordinates overlap the docked ligand (crystal ligands). Bound metals
    # are kept as typed ions + locked into their crystal site (bonded model).
    # --no-auto-cofactors reverts to "only keep explicit --cofactor resnames".
    need_hetnam = cfg.auto_cofactors or any(v is None for v in cfg.cofactors.values())
    hetnam = parse_hetnam(protein_pdb) if need_hetnam else {}
    keep_cofactors, _drop, metals, renames = discover_hetero(
        pdb.topology, pos_ang, lig_heavy_ang, cfg, hetnam
    )
    kept_metal_ids = {(m["chain"], m["resSeq"]) for m in metals}
    for (chain, resseq), new_name in renames.items():
        print(f"[build] protonation: {chain}/{resseq} -> {new_name} (metal-coordinating)")

    # Each kept cofactor residue is extracted with its own bound coords, H-added,
    # and re-added as a standalone residue (so its H-bearing formula matches the
    # registered molecule). One unique OFFMol per distinct SMILES is registered.
    # Kept metals stay in the protein topology (typed by the water FF); everything
    # else non-protein (waters, solvent ions, dropped metals, crystal ligands) goes.
    cofactor_adds = []          # (topology, positions) per cofactor residue
    cofactor_offmols = []       # unique OFFMols for GAFF registration
    seen_smiles = set()
    delete_atoms = []
    for res in pdb.topology.residues():
        if res.name in keep_cofactors:
            matched = _residue_to_rdmol(pdb.topology, pos_ang, res, keep_cofactors[res.name])
            offmol, top, pos = _explicit_h_rdmol_to_openmm(matched, residue_name=res.name)
            cofactor_adds.append((top, pos))
            smi = offmol.to_smiles()
            if smi not in seen_smiles:
                seen_smiles.add(smi)
                cofactor_offmols.append(offmol)
                print(f"[build] cofactor {res.name} -> GAFF2 ({smi[:60]})")
            delete_atoms.extend(res.atoms())  # remove from protein topology; re-added below
        elif (res.chain.id, res.id) in kept_metal_ids:
            pass  # keep the bound metal in the topology (typed by water FF, locked later)
        elif res.name not in PROTEIN_RES:
            delete_atoms.extend(res.atoms())  # HOH, solvent ions, dropped metals, crystal ligands, ...

    modeller = app.Modeller(pdb.topology, pdb.positions)
    if delete_atoms:
        modeller.delete(delete_atoms)
        print(f"[build] removed {len(delete_atoms)} non-protein/non-cofactor atoms "
              f"(waters, crystal ligands); re-adding {len(cofactor_adds)} cofactor(s), "
              f"keeping {len(metals)} bound metal(s)")
    for top, pos in cofactor_adds:
        modeller.add(top, pos)
    modeller.add(lig_topology, lig_positions)

    ff_files = [cfg.protein_ff, cfg.water_ff]
    if cfg.ion_ff:
        ff_files.append(cfg.ion_ff)
    return _solvate_and_create(
        modeller, [lig_offmol, *cofactor_offmols], cfg, ff_files,
        add_hydrogens=True, out_dir=out_dir, ligand_smiles=lig_offmol.to_smiles(),
        metals=metals,
    )


def build_mol_system(
    ligand_sdf: Path,
    out_dir: Path,
    cfg: Config | None = None,
    site_xyz_ang=None,
):
    """Solvate a single small molecule (no protein) in TIP3P + ions, GAFF2 for the solute.

    The molecule is the only solute, placed at the box center (or at ``site_xyz_ang``
    Å if given). No protein force field is loaded -- the solute is typed entirely
    by the GAFFTemplateGenerator and the water/ions come from the water FF. The
    ligand SDF already carries explicit H, so ``addHydrogens`` is skipped.
    """
    cfg = cfg or Config()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rdmol = load_ligand_rdkit(ligand_sdf)
    if site_xyz_ang is not None:
        rdmol = translate_to(rdmol, site_xyz_ang)
        print(f"[build-mol] placed molecule centroid at {list(site_xyz_ang)} A")
    lig_offmol, lig_topology, lig_positions = _rdkit_to_openmm(rdmol, residue_name="MOL")

    # Modeller from just the ligand topology; addSolvent will define the box.
    modeller = app.Modeller(lig_topology, lig_positions)
    print(f"[build-mol] solute {lig_offmol.to_smiles()} "
          f"({lig_offmol.n_atoms} atoms) -> GAFF2, solvating with {cfg.water_model}")
    return _solvate_and_create(
        modeller, [lig_offmol], cfg, [cfg.water_ff],
        add_hydrogens=False, out_dir=out_dir, ligand_smiles=lig_offmol.to_smiles(),
    )