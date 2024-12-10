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
sys.path.append('/home/seamustard52/repository/alphafold')

import collections
import dataclasses
import functools
import io
from typing import Any, Dict, List, Mapping, Optional, Tuple
from alphafold.common import mmcif_metadata
from alphafold.common import residue_constants
from alphafold.common import protein
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.Structure import Structure
import numpy as np

if __name__ == "__main__":
    pdb_dir = "/home/seamustard52/bfvd-analysis/pdbs_bfvd_logan"
    pdb_files = os.listdir(pdb_dir)
    # out_dir = "/home/seamustard52/bfvd-analysis/test"
    out_dir = "bfvd/test"


    for ex_in in pdb_files[:1]:
        ex_in = os.path.join(pdb_dir, ex_in)
        file_id = os.path.basename(ex_in).split(".")[0]
        ex_out = os.path.join(out_dir, file_id + ".cif")

        with open(ex_in) as f:
            pdb_string = f.read()

        prot = protein.from_pdb_string(pdb_string)
        cif = protein.to_mmcif(prot, file_id, "Monomer")
        print(cif)
        # with open(ex_out, "w") as f:
        #     f.write(cif)