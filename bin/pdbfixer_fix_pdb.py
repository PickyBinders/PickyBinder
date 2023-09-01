#!/usr/bin/env python

from pdbfixer import PDBFixer
from openmm.app import PDBFile
import sys

input_pdb = sys.argv[1]
protein_name = sys.argv[2]


def fix_pdb(pdb, pdb_name):
    fixer = PDBFixer(pdb)

    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    for key in list(keys):
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    fixer.removeHeterogens(keepWater=True)

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(7.0)

    #outfile = pdb.split('.pdb')[0] + '_Hs.pdb'
    outfile = pdb_name + '_Hs.pdb'

    PDBFile.writeFile(fixer.topology, fixer.positions, open(outfile, 'w+'))


def main():
    fix_pdb(input_pdb, protein_name)


if __name__ == '__main__':
    main()
