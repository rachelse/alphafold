"""extra metadata for BFVD mmCIF files."""

import sys
sys.path.append('/home/seamustard52/repository/alphafold-rachel')

from typing import Mapping, Sequence
# from alphafold.common.protein import _get_entity_poly_seq
from alphafold.common import residue_constants
import numpy as np

def add_ma_qa_metric_local(old_cif: Mapping[str, Sequence[str]]) -> Mapping[str, Sequence[str]]:

  """Adds local pLDDT QA metric to the CIF file. (_ma_qa_metric_local)"""
  cif = {}

  # Extract Residue level pLDDT scores
  index = []
  last = 0
  for i in range(len(old_cif['_atom_site.label_seq_id'])):
    if old_cif['_atom_site.label_seq_id'][i] != last:
       index.append(i)
       last = old_cif['_atom_site.label_seq_id'][i]
  
  cif['_ma_qa_metric_local.label_asym_id'] = [old_cif['_atom_site.label_asym_id'][i] for i in index]
  cif['_ma_qa_metric_local.labe_comp_id'] = [old_cif['_atom_site.label_comp_id'][i] for i in index]
  cif['_ma_qa_metric_local.label_seq_id'] = [old_cif['_atom_site.label_seq_id'][i] for i in index]
  cif['_ma_qa_metric_local.metric_id'] = ['2'] * len(index)
  cif['_ma_qa_metric_local.metric_value'] = [old_cif['_atom_site.B_iso_or_equiv'][i] for i in index]
  cif['_ma_qa_metric_local.model_id'] = [old_cif['_atom_site.pdbx_PDB_model_num'][i] for i in index]
  cif['_ma_qa_metric_local.ordinal_id'] = cif['_ma_qa_metric_local.label_seq_id']

  return cif

def get_entity_poly_seq_one(chain_index : np.ndarray, residue_index : np.ndarray, aatype : np.ndarray
                            ) -> Mapping[str, Sequence[str]]:

  entity_seq = {}
  last_chain = ""
  last_res = ""
  for i in range(len(chain_index)):
    if chain_index[i] not in entity_seq:
      entity_seq[chain_index[i]] = ""
    if last_chain == chain_index[i] and (residue_index[i] - last_res > 1):
      entity_seq[chain_index[i]] += "X" * (residue_index[i] - last_res -1)
    
    aa_code = residue_constants.restypes[int(aatype[i])]
    entity_seq[chain_index[i]] += aa_code
    last_chain = chain_index[i]
    last_res = residue_index[i]

  # insert new lines every 80 characters
  for chain in entity_seq.keys():
    entity_seq[chain] = "\n".join([entity_seq[chain][i:i+80] for i in range(0, len(entity_seq[chain]), 80)])    
  return entity_seq