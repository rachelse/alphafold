"""extra metadata for BFVD mmCIF files."""

import sys
sys.path.append('/home/seamustard52/repository/alphafold-rachel')

from typing import Mapping, Sequence
# from alphafold.common.protein import _get_entity_poly_seq
from alphafold.common import residue_constants
import numpy as np
import pandas as pd
import collections

def get_ma_qa_metric_local(old_cif: Mapping[str, Sequence[str]]) -> Mapping[str, Sequence[str]]:

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

def get_pdbx_audit_revision():
  """Returns the _pdbx_audit_revision_history category. 
  history: (ordinal, data_content_type, major_revision, minor_revision, revision_date)
  details: (ordinal, revision_ordinal, data_content_type, provider, type, description)"""
  versions = [
    {
      'data_content_type': 'Structure model', 'provider': 'BFVD', 
      'major_revision': '1', 'minor_revision': '0', 'revision_date': '2024-11-01', 
      'type': 'Initial release', 'description': 'Initial release'
     },
  ]
  cif = collections.defaultdict(list)

  for i, version in enumerate(versions, start=1):
    cif['_pdbx_audit_revision_history.ordinal'].append(str(i))
    cif['_pdbx_audit_revision_history.data_content_type'].append(version['data_content_type'])
    cif['_pdbx_audit_revision_history.major_revision'].append(version['major_revision'])
    cif['_pdbx_audit_revision_history.minor_revision'].append(version['minor_revision'])
    cif['_pdbx_audit_revision_history.revision_date'].append(version['revision_date'])

    cif['_pdbx_audit_revision_details.ordinal'].append(str(i))
    cif['_pdbx_audit_revision_details.revision_ordinal'].append(str(i))
    cif['_pdbx_audit_revision_details.data_content_type'].append(version['data_content_type'])
    cif['_pdbx_audit_revision_details.provider'].append(version['provider'])
    cif['_pdbx_audit_revision_details.type'].append(version['type'])
    cif['_pdbx_audit_revision_details.description'].append(version['description'])
  return cif

def get_ma_target_ref_db_details(old_cif: Mapping[str, Sequence[str]]) -> Mapping[str, Sequence[str]]:
  """Returns the _ma_target_ref_db_details category. 
  details: (target_entity_id, 
            db_name, db_name_other_details, db_code, db_accession,
            seq_db_isoform, seq_db_align_begin, seq_db_align_end, 
            gene_name, ncbi_taxonomy_id, organism__scientific,
            seq_db_sequence_version_date, seq_db_sequence_checksum)"""
  #TODO : Provide metadata with arguments
  uniprot_path = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/uniprot-acc_length_taxid_organism_src_id_description_gene_sampled.tsv"
  bfvd_path = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/bfvd_logan-entry_acc_start_end_len_plddt_taxid_organism_src.tsv"

  uniprot_data = pd.read_csv(uniprot_path, sep="\t", 
                           names=["acc", "length", "taxid", "organism", "src", "id", "description", "gene"],
                           dtype={"acc": str, "length": int, "taxid": str, "organism": str, "src": str, "id": str, "description": str, "gene": str},
                           index_col=0
                           ).fillna("?").T.to_dict()
  
  bfvd_data = pd.read_csv(bfvd_path, sep="\t",
                        names=["entry", "acc", "start", "end", "len", "plddt", "taxid", "organism", "src"], 
                        dtype={"entry": str, "acc": str, "start": int, "end": int, "len": int, "plddt": float, "taxid": str, "organism": str, "src": str},
                        index_col=0).T.to_dict()

  cif = collections.defaultdict(list)
  filename = old_cif['_entry.id']
  idx = filename.find("_UNRELAXED_RANK_001")
  entry = filename[:idx]
  acc = entry.split("_")[0]
  print(filename)
  for id in range(len(old_cif['_entity.id'])):

    cif['_ma_target_ref_db_details.target_entity_id'].append(old_cif['_entity.id'][id])
    cif['_ma_target_ref_db_details.db_name'].append('UNP')

    if bfvd_data[entry]["src"] == "UNIPARC":
      db_details = "UNIPARC"
      db_code = "?"
      gene = "?"
    else:
      db_details = "?"
      db_code = uniprot_data[acc]['id']
      gene = uniprot_data[acc]['gene']

    cif['_ma_target_ref_db_details.db_name_other_details'].append(db_details)
    cif['_ma_target_ref_db_details.db_code'].append(db_code)
    cif['_ma_target_ref_db_details.db_accession'].append(acc)
    cif['_ma_target_ref_db_details.seq_db_isoform'].append('?')
    cif['_ma_target_ref_db_details.seq_db_align_begin'].append(str(bfvd_data[entry]['start']))
    cif['_ma_target_ref_db_details.seq_db_align_end'].append(str(bfvd_data[entry]['end']))
    cif['_ma_target_ref_db_details.gene_name'].append(gene)
    cif['_ma_target_ref_db_details.ncbi_taxonomy_id'].append(str(bfvd_data[entry]['taxid']))
    cif['_ma_target_ref_db_details.organism__scientific'].append(bfvd_data[entry]['organism'])
    cif['_ma_target_ref_db_details.seq_db_sequence_version_date'].append('2023-02') #TODO
    cif['_ma_target_ref_db_details.seq_db_sequence_checksum'].append('?') #TODO

  return cif