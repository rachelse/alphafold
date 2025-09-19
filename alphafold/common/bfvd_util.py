"""extra metadata for BFVD mmCIF files."""

import sys
sys.path.append('/home/seamustard52/repository/alphafold-rachel')

from typing import Mapping, Sequence
# from alphafold.common.protein import _get_entity_poly_seq
from alphafold.common import residue_constants
import numpy as np
import pandas as pd
import collections
from typing import List
import alphafold.common.bfvd_constants as bfvd_constants

def add_atom_site(old_cif: Mapping[str, Sequence[str]]):
  """Adds the _atom_site category to the CIF file."""
  acc = old_cif['_ma_target_ref_db_details.db_accession'][0] # NOTE:This only works for monomer
  db_name = old_cif['_ma_target_ref_db_details.db_name'][0]
  seq_start = int(old_cif['_ma_target_ref_db_details.seq_db_align_begin'][0])

  old_cif['_atom_site.pdbx_sifts_xref_db_acc'] = [acc] * len(old_cif['_atom_site.label_seq_id'])
  old_cif['_atom_site.pdbx_sifts_xref_db_name'] = [db_name] * len(old_cif['_atom_site.label_seq_id'])
  old_cif['_atom_site.pdbx_sifts_xref_db_num'] = [str(seq_start + int(i) - 1) for i in old_cif['_atom_site.label_seq_id']]
  old_cif['_atom_site.pdbx_sifts_xref_db_res'] = [residue_constants.restype_3to1[aa] for aa in old_cif['_atom_site.label_comp_id']]

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
  cif['_ma_qa_metric_local.label_comp_id'] = [old_cif['_atom_site.label_comp_id'][i] for i in index]
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

def get_pdbx_audit_revision(versions: List) -> Mapping[str, Sequence[str]]:
  """Returns the _pdbx_audit_revision_history category. 
  history: (ordinal, data_content_type, major_revision, minor_revision, revision_date)
  details: (ordinal, revision_ordinal, data_content_type, provider, type, description)"""

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

def get_ma_target_ref_db_details(old_cif: Mapping[str, Sequence[str]],
                                ) -> Mapping[str, Sequence[str]]:
  """Returns the _ma_target_ref_db_details category. 
  details: (target_entity_id, 
            db_name, db_name_other_details, db_code, db_accession,
            seq_db_isoform, seq_db_align_begin, seq_db_align_end, 
            gene_name, ncbi_taxonomy_id, organism__scientific,
            seq_db_sequence_version_date, seq_db_sequence_checksum)"""

  uniprot_data = pd.read_csv(bfvd_constants._UNIPROT_DATA, sep="\t", 
                           names=["acc", "length", "taxid", "organism", "src", "id", "description", "gene"],
                           dtype={"acc": str, "length": int, "taxid": str, "organism": str, "src": str, "id": str, "description": str, "gene": str},
                           index_col=0
                           ).fillna("?").T.to_dict()
  
  bfvd_data = pd.read_csv(bfvd_constants._BFVD_DATA, sep="\t",
                        names=["entry", "acc", "start", "end", "len", "plddt", "taxid", "organism", "src"], 
                        dtype={"entry": str, "acc": str, "start": int, "end": int, "len": int, "plddt": float, "taxid": str, "organism": str, "src": str},
                        index_col=0).T.to_dict()

  cif = collections.defaultdict(list)
  filename = old_cif['_entry.id']
  idx = filename.find("_UNRELAXED_RANK_001", 0, len(filename))
  entry = filename[:idx]
  entry = filename
  acc = entry.split("_")[0]
  
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

    taxid = bfvd_data[entry]['taxid']
    organism = bfvd_data[entry]['organism']
    try:
      valid = len(taxid)
    except:
      taxid = "?"
      organism = "?"
      
    cif['_ma_target_ref_db_details.db_name_other_details'].append(db_details)
    cif['_ma_target_ref_db_details.db_code'].append(db_code)
    cif['_ma_target_ref_db_details.db_accession'].append(acc)
    cif['_ma_target_ref_db_details.seq_db_isoform'].append('?')
    cif['_ma_target_ref_db_details.seq_db_align_begin'].append(str(bfvd_data[entry]['start']))
    cif['_ma_target_ref_db_details.seq_db_align_end'].append(str(bfvd_data[entry]['end']))
    cif['_ma_target_ref_db_details.gene_name'].append(gene)
    cif['_ma_target_ref_db_details.ncbi_taxonomy_id'].append(taxid)
    cif['_ma_target_ref_db_details.organism__scientific'].append(organism)
    cif['_ma_target_ref_db_details.seq_db_sequence_version_date'].append('2023-02') #TODO
    cif['_ma_target_ref_db_details.seq_db_sequence_checksum'].append('?') #TODO

  return cif

def get_citation(citation: List) -> Mapping[str, Sequence[str]]:
  """Returns the _citation category. 
  details: (id, title, journal_abbrev, journal_volume, page_first, page_last, year, pdbx_database_id_PubMed, pdbx_database_id_DOI, journal_id_ISSN, journal_id_ASTM, journal_id_CSD, country, book_publisher, authors)"""
  cif = collections.defaultdict(list)
  for i, cite in enumerate(citation, start=1):
    cif['_citation.id'].append(str(i))
    cif['_citation.title'].append(cite['title'])
    cif['_citation.journal_full'].append(cite['journal_full'])
    cif['_citation.journal_volume'].append(cite['journal_volume'])
    cif['_citation.page_first'].append(cite['page_first'])
    cif['_citation.page_last'].append(cite['page_last'])
    cif['_citation.year'].append(cite['year'])
    cif['_citation.pdbx_database_id_PubMed'].append(cite['pdbx_database_id_PubMed'])
    cif['_citation.pdbx_database_id_DOI'].append(cite['pdbx_database_id_DOI'])
    cif['_citation.journal_id_ISSN'].append(cite['journal_id_ISSN'])
    cif['_citation.journal_id_ASTM'].append(cite['journal_id_ASTM'])
    cif['_citation.journal_id_CSD'].append(cite['journal_id_CSD'])
    cif['_citation.country'].append(cite['country'])
    cif['_citation.book_publisher'].append(cite['book_publisher'])
  return cif

def get_citation_author(citation: List) -> Mapping[str, Sequence[str]]:
  """Returns the _citation_author category. 
  details: (citation_id, name, ordinal)"""
  cif = collections.defaultdict(list)
  ordinal = 1
  for id, cite in enumerate(citation, start=1):
    for author in cite['authors']:
      cif['_citation_author.citation_id'].append(str(id))
      cif['_citation_author.name'].append(author)
      cif['_citation_author.ordinal'].append(str(ordinal))
      ordinal += 1
  return cif

def get_ma_protocol_step() -> Mapping[str, Sequence[str]]:
  """Returns the _ma_protocol_step category. 
  details: (ordinal_id, protocol_id, step_id, 'method_type')"""
  cif = collections.defaultdict(list)
  cif['_ma_protocol_step.ordinal_id'] = ['1', '2']
  cif['_ma_protocol_step.protocol_id'] = ['1', '1']
  cif['_ma_protocol_step.step_id'] = ['1', '2']
  cif['_ma_protocol_step.method_type'] = ['coevolution MSA', 'modeling']
  return cif

def get_software(softwares: List) -> Mapping[str, Sequence[str]]:
  """Returns the _software category. 
  details: (ordinal, name, version, classification)"""
  cif = collections.defaultdict(list)
  for i, software in enumerate(softwares, start=1):
    cif['_software.pdbx_ordinal'].append(str(i))
    cif['_software.name'].append(software['name'])
    cif['_software.version'].append(software['version'])
    cif['_software.type'].append(software['type'])
    cif['_software.description'].append(software['description'])
    cif['_software.classification'].append(software['classification'])

    cif['_ma_software_group.ordinal_id'].append(str(i))
    cif['_ma_software_group.software_id'].append(str(i))
    cif['_ma_software_group.group_id'].append('1')

  return cif