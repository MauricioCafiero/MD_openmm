"""Command-line entrypoint: ``omd <subcommand> ...``.

Stages can be run independently so you can inspect intermediates between steps.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from .config import Config
from .prepare_protein import prepare_protein
from .prepare_ligand import smiles_to_sdf, load_ligand_rdkit, write_sdf
from .build_system import build_system, build_mol_system, KNOWN_COFACTORS
from .dynamics import run as run_dynamics
from .analyze import analyze


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

    p = sub.add_parser("run", help="minimize + equilibrate + production MD")
    p.add_argument("--system", type=Path, required=True)
    p.add_argument("--topology", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, default=Path("outputs"))
    p.add_argument("--steps", type=int, default=Config().production_steps)
    p.add_argument("--platform", default=Config().platform)

    p = sub.add_parser("analyze", help="RMSD/RMSF + energy plots")
    p.add_argument("--traj", type=Path, required=True)
    p.add_argument("--topology", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, default=Path("outputs"))

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

    elif args.cmd == "run":
        cfg.platform = args.platform
        run_dynamics(args.system, args.topology, args.out_dir, cfg=cfg, steps=args.steps)

    elif args.cmd == "analyze":
        analyze(args.traj, args.topology, args.out_dir)


if __name__ == "__main__":
    main()