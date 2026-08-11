"""Command-line entrypoint: ``omd <subcommand> ...``.

Stages can be run independently so you can inspect intermediates between steps.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Put the env's bin/ on PATH before any openff/openmmforcefields import: those
# build the toolkit registry at import time, and GAFF AM1-BCC needs the
# antechamber/sqm binaries to resolve (else AmberToolsToolkitWrapper is skipped).
_env_bin = os.path.dirname(sys.executable)
if _env_bin and _env_bin not in os.environ.get("PATH", "").split(os.pathsep):
    os.environ["PATH"] = _env_bin + os.pathsep + os.environ.get("PATH", "")

from .config import Config
from .prepare_protein import prepare_protein
from .prepare_ligand import smiles_to_sdf, load_ligand_rdkit, write_sdf
from .build_system import build_system, build_mol_system, build_multimol_system, KNOWN_COFACTORS
from .dynamics import run as run_dynamics
from .analyze import analyze
from . import mmgbsa


def _site(s):
    if s is None:
        return None
    xyz = [float(v) for v in s.replace(",", " ").split()]
    if len(xyz) != 3:
        raise argparse.ArgumentTypeError("--site expects 3 numbers 'x y z' (angstrom)")
    return xyz


def main(argv=None):
    parser = argparse.ArgumentParser(prog="omd", description="OpenMM protein/ligand MD pipeline")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prep-protein", help="PDBFixer: repair + hydrogens")
    p.add_argument("--pdb", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("outputs/protein_prep.pdb"))
    p.add_argument("--ph", type=float, default=7.0)
    p.add_argument("--add-missing-residues", action="store_true")

    p = sub.add_parser("prep-ligand", help="SMILES->SDF or SDF passthrough")
    p.add_argument("--smiles")
    p.add_argument("--sdf", type=Path, help="input SDF (passthrough: sanitize + re-save)")
    p.add_argument("--out", type=Path, default=Path("outputs/ligand.sdf"))
    p.add_argument("--num-confs", type=int, default=5)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--site", type=_site, help="place SMILES-ligand centroid at 'x y z' (angstrom)")

    p = sub.add_parser("build", help="solvated protein/ligand system")
    p.add_argument("--protein", type=Path, required=True)
    p.add_argument("--ligand", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, default=Path("outputs"))
    p.add_argument("--site", type=_site, help="place ligand centroid at 'x y z' (angstrom)")
    p.add_argument("--padding", type=float, default=Config().padding)
    p.add_argument(
        "--cofactor", action="append", default=[], metavar="RES[:SMILES]",
        help="force-keep a bound cofactor resname + parameterize with GAFF2: a built-in "
             "name (e.g. A3P), RES:SMILES (custom), or bare RES (SMILES auto-resolved "
             "from the PDB HETNAM name via PubChem). Overrides the auto overlap-drop "
             "for that resname. Repeatable.",
    )
    p.add_argument(
        "--no-auto-cofactors", action="store_true",
        help="disable auto cofactor discovery; only keep explicit --cofactor resnames "
             "(offline / full control).",
    )
    p.add_argument(
        "--no-keep-metals", action="store_true",
        help="drop all structural metal ions (Zn2+/Fe2+/...) instead of keeping bound "
             "metals locked in their crystal coordination site (bonded model).",
    )

    p = sub.add_parser("build-mol", help="solvate a single small molecule (no protein)")
    p.add_argument("--ligand", type=Path, help="input ligand SDF (or use --smiles to build one)")
    p.add_argument("--smiles", help="build the ligand from this SMILES, then solvate")
    p.add_argument("--out-dir", type=Path, default=Path("outputs"))
    p.add_argument("--site", type=_site, help="place molecule centroid at 'x y z' (angstrom)")
    p.add_argument("--padding", type=float, default=Config().padding)
    p.add_argument("--num-confs", type=int, default=5)

    p = sub.add_parser("build-multimol",
                       help="solvate N covalently-separate molecules (no protein); "
                            "a rotaxane rod+wheel SDF -> two non-bonded GAFF2 solutes")
    p.add_argument("--sdf", type=Path, required=True,
                   help="multi-fragment SDF (one mol block per fragment)")
    p.add_argument("--out-dir", type=Path, default=Path("outputs"))
    p.add_argument("--padding", type=float, default=Config().padding)
    p.add_argument("--resnames", nargs="+", default=None,
                   help="residue name per fragment (default ROD WHL for 2, else FRG0 FRG1..)")

    p = sub.add_parser("run", help="minimize + equilibrate + production MD")
    p.add_argument("--system", type=Path, required=True)
    p.add_argument("--topology", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, default=Path("outputs"))
    p.add_argument("--steps", type=int, default=Config().production_steps)
    p.add_argument("--platform", default=Config().platform)
    p.add_argument("--pressure", type=float, default=Config().pressure,
                   help="pressure in atm (NPT); set <= 0 for NVT. Default 1.0.")
    p.add_argument("--ramp-end", type=float, default=None,
                   help="linearly ramp the thermostat temperature from the base "
                        "temperature to this value (K) over the production steps "
                        "(heating / annealing). Default: constant temperature.")

    p = sub.add_parser("analyze", help="RMSD/RMSF + energy plots")
    p.add_argument("--traj", type=Path, required=True)
    p.add_argument("--topology", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, default=Path("outputs"))

    p = sub.add_parser("mmgbsa", help="MM/GBSA binding energy (ParmEd prmtop + MMPBSA.py)")
    p.add_argument("--protein", type=Path, required=True,
                   help="original protein PDB (needs CONNECT for cofactor extraction)")
    p.add_argument("--ligand", type=Path, required=True, help="docked ligand SDF")
    p.add_argument("--traj", type=Path, required=True,
                   help="solute-only trajectory (e.g. traj_wrapped.xtc)")
    p.add_argument("--topology", type=Path, required=True,
                   help="solute-only topology matching --traj (e.g. traj_wrapped.pdb); "
                        "the prmtop is built from this (production protein H/tautomers, "
                        "trajectory atom order)")
    p.add_argument("--out-dir", type=Path, default=Path("outputs/mmgbsa"))
    p.add_argument("--no-auto-cofactors", action="store_true",
                   help="disable auto cofactor discovery; only explicit --cofactor resnames")
    p.add_argument("--cofactor", action="append", default=[], metavar="RES[:SMILES]",
                   help="force-keep a cofactor resname (RES:SMILES or bare RES)")
    p.add_argument("--run", action="store_true",
                   help="also run the MMPBSA.py energy evaluation (CPU-heavy). "
                        "Default: prep only (prmtop + traj.nc + mmgbsa.in).")

    args = parser.parse_args(argv)
    cfg = Config()

    if args.cmd == "prep-protein":
        prepare_protein(args.pdb, args.out, ph=args.ph,
                        add_missing_residues=args.add_missing_residues)

    elif args.cmd == "prep-ligand":
        if args.smiles and args.sdf:
            parser.error("provide either --smiles or --sdf, not both")
        if args.smiles:
            out, _ = smiles_to_sdf(args.smiles, args.out, num_confs=args.num_confs)
            if args.site is not None:
                from .prepare_ligand import translate_to
                mol = load_ligand_rdkit(out)
                translate_to(mol, args.site)
                write_sdf(mol, out)
                print(f"[prep-ligand] centroid placed at {args.site} A -> {out}")
        elif args.sdf:
            mol = load_ligand_rdkit(args.sdf)
            if args.site is not None:
                from .prepare_ligand import translate_to
                translate_to(mol, args.site)
            write_sdf(mol, args.out)
            print(f"[prep-ligand] SDF passthrough -> {args.out}")
        else:
            parser.error("provide --smiles or --sdf")

    elif args.cmd == "build":
        cfg.padding = args.padding
        cfg.auto_cofactors = not args.no_auto_cofactors
        cfg.keep_metals = not args.no_keep_metals
        cofactors = {}
        for c in args.cofactor:
            if ":" in c:
                name, smi = c.split(":", 1)
                cofactors[name] = smi
            elif c in KNOWN_COFACTORS:
                cofactors[c] = KNOWN_COFACTORS[c]
            else:
                cofactors[c] = None  # force-keep; SMILES auto-resolved at build
        cfg.cofactors = cofactors
        build_system(args.protein, args.ligand, args.out_dir, cfg=cfg,
                     site_xyz_ang=args.site)

    elif args.cmd == "build-mol":
        if args.smiles and args.ligand:
            parser.error("provide either --smiles or --ligand, not both")
        if args.smiles:
            cfg.padding = args.padding
            lig = args.out_dir / "ligand.sdf"
            smiles_to_sdf(args.smiles, lig, num_confs=args.num_confs)
            if args.site is not None:
                from .prepare_ligand import load_ligand_rdkit, translate_to, write_sdf
                m = load_ligand_rdkit(lig); translate_to(m, args.site); write_sdf(m, lig)
        elif args.ligand:
            cfg.padding = args.padding
            lig = args.ligand
        else:
            parser.error("provide --smiles or --ligand")
        build_mol_system(lig, args.out_dir, cfg=cfg, site_xyz_ang=args.site)

    elif args.cmd == "build-multimol":
        cfg.padding = args.padding
        build_multimol_system(args.sdf, args.out_dir, cfg=cfg,
                             residue_names=args.resnames)

    elif args.cmd == "run":
        cfg.platform = args.platform
        cfg.pressure = args.pressure
        if args.ramp_end is not None:
            cfg.temperature_ramp_end = args.ramp_end
        run_dynamics(args.system, args.topology, args.out_dir, cfg=cfg, steps=args.steps)

    elif args.cmd == "analyze":
        analyze(args.traj, args.topology, args.out_dir)

    elif args.cmd == "mmgbsa":
        cfg.auto_cofactors = not args.no_auto_cofactors
        for c in args.cofactor:
            name, smi = (c.split(":", 1) + [None])[:2] if ":" in c else (c, None)
            cfg.cofactors[name] = smi
        prepared = mmgbsa.prep(args.protein, args.ligand, args.traj, args.topology,
                               args.out_dir, cfg=cfg)
        if args.run:
            mmgbsa.run(args.out_dir)
        else:
            print(f"[mmgbsa] prep done. Run the (CPU-heavy) evaluation with:\n"
                  f"  omd mmgbsa --protein {args.protein} --ligand {args.ligand} "
                  f"--traj {args.traj} --topology {args.topology} --out-dir {args.out_dir} --run")


if __name__ == "__main__":
    main()