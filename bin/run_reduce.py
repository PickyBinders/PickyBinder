#!/usr/bin/env python

import os
import sys
import subprocess

pdb_input = sys.argv[1]
receptor_name = sys.argv[2]
sing_img = sys.argv[3]
het_dict = sys.argv[4]


def run_reduce(singularity_img, het_dict_file, input_pdb, output_path):
    """Reduce (add Hydrogens) with ADFR reduce.

    The reduce command often returns an exit code 1. This function catches the
    CalledProcessError and if the output file looks OK (large enough)
    it will ignore the error. All other errors (other exit codes,
    error opening the file etc.) or in case the output file is empty will be
    raised.
    """
    # Add Hydrogens
    #hydrogenate_cmd = ["singularity", "exec", "-B", "/scicore/home/schwede", "-B", "/scratch", "-B", "/home", singularity_img, "reduce",
    #sing_dir = os.path.dirname(singularity_img)
    hydrogenate_cmd = ["singularity", "exec", "-B", "/scicore/home/schwede", "-B", "/scratch", "-B", "/home", singularity_img, "reduce",
                       "-FLIP",
                       "-DB", het_dict_file,
                       input_pdb
                       ]

    try:
        with open(output_path, "w+") as out_fd:
            subprocess.run(hydrogenate_cmd, stdout=out_fd, check=True)
        #subprocess.run(hydrogenate_cmd)
    except subprocess.CalledProcessError as err:
        if err.returncode == 1:
            print("ADFR reduce (hydrogenate) returned exit code 1.")
            # If the output file ends with END, and there are more than 1000
            # bytes before that, ignore the failure
            with open(output_path, "rb") as out_fd:
                pos = out_fd.seek(-4, os.SEEK_END)
                if pos > 1000:
                    print("Output file looks OK, continuing.")
                else:
                    print("ADFR reduce returned exit code 1 and file looks too short/truncated.")
                    raise
        else:
            print("ADFR reduce (hydrogenate) returned an error code %s", err.returncode)
            print("Command run: %s", " ".join(hydrogenate_cmd))
            raise


def main():
    out_pdb = receptor_name + '_Hs.pdb'
    run_reduce(sing_img, het_dict, pdb_input, out_pdb)


if __name__ == '__main__':
    main()
