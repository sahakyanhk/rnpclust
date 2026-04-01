#!/usr/bin/env python3
"""
Extract interface residues from a PDB file based on inter-chain distance.

Usage:
  # Extract everything within 5A of chain A (all other chains near A)
  python extract_interface.py --input complex.pdb --chain A --cutoff 5

  # Extract interface between chains A and B (residues of A near B + residues of B near A)
  python extract_interface.py --input complex.pdb --chain A,B --cutoff 5

  # Custom output path
  python extract_interface.py --input complex.pdb --chain A,B --cutoff 5 -o interface.pdb
"""

import argparse
import sys

from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch


class InterfaceSelect(Select):
    """Select only residues in the interface set."""

    def __init__(self, interface_residues):
        self._residues = interface_residues

    def accept_residue(self, residue):
        return residue in self._residues


def get_interface_residues(structure, chains, cutoff):
    """
    Find interface residues using KDTree neighbor search.

    If one chain is given: find residues from all other chains within
    cutoff of that chain, plus residues of that chain near any other.

    If two chains are given: find residues of chain A near chain B
    and residues of chain B near chain A.
    """
    model = structure[0]
    chain_ids = {c.id for c in model.get_chains()}

    for c in chains:
        if c not in chain_ids:
            print(
                f"ERROR: chain '{c}' not found in structure. "
                f"Available chains: {', '.join(sorted(chain_ids))}",
                file=sys.stderr,
            )
            sys.exit(1)

    if len(chains) == 1:
        query_chain = chains[0]
        partner_chains = chain_ids - {query_chain}
        if not partner_chains:
            print("ERROR: only one chain in structure, nothing to compare", file=sys.stderr)
            sys.exit(1)

        query_atoms = list(model[query_chain].get_atoms())
        partner_atoms = [a for cid in partner_chains for a in model[cid].get_atoms()]

    else:
        chain_a, chain_b = chains
        query_atoms = list(model[chain_a].get_atoms())
        partner_atoms = list(model[chain_b].get_atoms())

    # Find query residues near partner
    ns_partner = NeighborSearch(partner_atoms)
    query_interface = set()
    for atom in query_atoms:
        neighbors = ns_partner.search(atom.get_vector().get_array(), cutoff, level='R')
        if neighbors:
            query_interface.add(atom.get_parent())

    # Find partner residues near query
    ns_query = NeighborSearch(query_atoms)
    partner_interface = set()
    for atom in partner_atoms:
        neighbors = ns_query.search(atom.get_vector().get_array(), cutoff, level='R')
        if neighbors:
            partner_interface.add(atom.get_parent())

    return query_interface | partner_interface


def main():
    p = argparse.ArgumentParser(
        description="Extract interface residues from a PDB file"
    )
    p.add_argument("--input", required=True, help="Input PDB file")
    p.add_argument(
        "--chain", required=True,
        help="Chain(s) to extract interface for. "
             "One chain (e.g. A): extract everything within cutoff of chain A. "
             "Two chains (e.g. A,B): extract interface between A and B only."
    )
    p.add_argument(
        "--cutoff", type=float, default=5.0,
        help="Distance cutoff in Angstroms (default: 5.0)"
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output PDB file (default: <input>_interface.pdb)"
    )
    args = p.parse_args()

    chains = [c.strip() for c in args.chain.split(',')]
    if len(chains) > 2:
        p.error("--chain accepts at most 2 chains (e.g. A or A,B)")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('s', args.input)

    interface_residues = get_interface_residues(structure, chains, args.cutoff)

    if not interface_residues:
        print("WARNING: no interface residues found", file=sys.stderr)
        sys.exit(0)

    output = args.output
    if output is None:
        base = args.input.rsplit('.pdb', 1)[0]
        output = f"{base}_interface.pdb"

    io = PDBIO()
    io.set_structure(structure)
    io.save(output, InterfaceSelect(interface_residues))

    # Summary
    from collections import Counter
    chain_counts = Counter(r.get_parent().id for r in interface_residues)
    summary = ", ".join(f"chain {c}: {n} residues" for c, n in sorted(chain_counts.items()))
    print(f"Extracted {len(interface_residues)} interface residues ({summary}) → {output}", file=sys.stderr)


if __name__ == "__main__":
    main()
