#!/usr/bin/env python

from ost import io
import sys

pdb_files = sys.argv[1:]

for pdb_file in pdb_files:
    name_sdf = pdb_file.split('.pdb')[0] + '.sdf'
    try:
        pdb_ent = io.LoadPDB(pdb_file, read_conect=True)
        io.SaveEntity(pdb_ent, name_sdf, 'sdf')
    except Exception as e:
        print(name_sdf + ': ' + str(e))
        