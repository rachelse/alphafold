# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Protein data type."""
import sys
import os
sys.path.append('/home/seamustard52/repository/alphafold-rachelse')
import argparse
import concurrent.futures
from pathlib import Path

from alphafold.common import protein

def pdb2cif(pdb_file, cif_file, file_id=None):
    with open(pdb_file) as f:
        pdb_string = f.read()

    if file_id is None:
        file_id = os.path.basename(pdb_file).split(".")[0]

    prot = protein.from_pdb_string(pdb_string)
    print(f"Converting {pdb_file} to {cif_file}")
    cif = protein.to_mmcif(prot, file_id,  "Monomer")
    with open(cif_file, "w") as f:
        print(f"Writing cif to {cif_file}")
        f.write(cif)

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()

    argparser.add_argument("--pdb", "-i", type=str, 
                           help="Input pdb file or directory",
                           default = "bfvd/test_pdb",
                        #    default="/home/seamustard52/bfvd-analysis/pdbs_bfvd_logan"
                           )
    argparser.add_argument("--cif", "-o", type=str,
                            help="Output cif file or directory",
                            default = "bfvd/test_cif")
    argparser.add_argument("--cpus", "-c", type=int, default=os.cpu_count())
    
    args = argparser.parse_args()

    if not os.path.exists(args.pdb):
        raise ValueError("Input path does not exist")

    input_path = Path(args.pdb)

    if input_path.is_file():
        if args.cif.endswith(".cif"):
            output_path = Path(args.cif)
        else:
            raise ValueError("Output path must be a .cif file")
        output_path = Path(args.cif)
        pdb2cif(input_path, output_path)
    elif args.cpus ==1:
        if not os.path.exists(args.cif):
            os.makedirs(args.cif)
        
        pdb_files = list(input_path.glob("*.pdb"))
        cif_files = [os.path.join(args.cif, os.path.basename(pdb_file).split(".")[0] + ".cif") for pdb_file in pdb_files]
        for pdb_file, cif_file in zip(pdb_files, cif_files):
            pdb2cif(pdb_file, cif_file)
    else:
        if not os.path.exists(args.cif):
            os.makedirs(args.cif)
        
        pdb_files = list(input_path.glob("*.pdb"))
        cif_files = [os.path.join(args.cif, os.path.basename(pdb_file).split(".")[0] + ".cif") for pdb_file in pdb_files]


        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.cpus) as executor:

            futures = [executor.submit(pdb2cif, pdb_file, cif_file) for pdb_file, cif_file in zip(pdb_files, cif_files)]
            concurrent.futures.wait(futures)