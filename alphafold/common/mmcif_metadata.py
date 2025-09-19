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

"""mmCIF metadata."""

from typing import Mapping, Sequence
from alphafold import version
from alphafold.common import bfvd_util, bfvd_constants
import numpy as np


_DISCLAIMER = """?""" 

# _DISCLAIMER = """ALPHAFOLD DATA, COPYRIGHT (2021) DEEPMIND TECHNOLOGIES LIMITED.
# THE INFORMATION PROVIDED IS THEORETICAL MODELLING ONLY AND CAUTION SHOULD BE
# EXERCISED IN ITS USE. IT IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND,
# WHETHER EXPRESSED OR IMPLIED. NO WARRANTY IS GIVEN THAT USE OF THE INFORMATION
# SHALL NOT INFRINGE THE RIGHTS OF ANY THIRD PARTY. DISCLAIMER: THE INFORMATION IS
# NOT INTENDED TO BE A SUBSTITUTE FOR PROFESSIONAL MEDICAL ADVICE, DIAGNOSIS, OR
# TREATMENT, AND DOES NOT CONSTITUTE MEDICAL OR OTHER PROFESSIONAL ADVICE. IT IS
# AVAILABLE FOR ACADEMIC AND COMMERCIAL PURPOSES, UNDER CC-BY 4.0 LICENCE."""

# Authors of the Nature methods paper we reference in the mmCIF.
# _MMCIF_PAPER_AUTHORS = (
#     'Kim, Rachel Seongeun',
#     'Levy Karin, Eli',
#     'Mirdita, Milot',
#     'Chikhi, Rayan',
#     'Steinegger, Martin',
# )

# Authors of the mmCIF - we set them to be equal to the authors of the paper.
_MMCIF_AUTHORS = bfvd_constants._BFVD_CITATION["authors"]
_CITATIONS = [bfvd_constants._BFVD_CITATION, bfvd_constants._COLABFOLD_CITATION, bfvd_constants._ALPHAFOLD_CITATION]


def add_metadata_to_mmcif(
    old_cif: Mapping[str, Sequence[str]], model_type: str,
) -> Mapping[str, Sequence[str]]:
  """Adds AlphaFold metadata in the given mmCIF."""
  cif = {}

  # ModelCIF conformation dictionary.
  cif['_audit_conform.dict_name'] = ['mmcif_ma.dic']
  cif['_audit_conform.dict_version'] = ['1.3.9']
  cif['_audit_conform.dict_location'] = [
      'https://raw.githubusercontent.com/ihmwg/ModelCIF/master/dist/'
      'mmcif_ma.dic'
  ]

  # License and disclaimer.
  cif['_pdbx_data_usage.id'] = ['1', '2']
  cif['_pdbx_data_usage.type'] = ['license', 'disclaimer']
  cif['_pdbx_data_usage.details'] = [
      'Data in this file is available under a CC-BY-4.0 license.',
      _DISCLAIMER,
  ]
  cif['_pdbx_data_usage.url'] = [
      'https://creativecommons.org/licenses/by/4.0/',
      '?',
  ]
  cif['_pdbx_data_usage.name'] = ['CC-BY-4.0', '?']

  # Structure author details.
  cif['_audit_author.name'] = []
  cif['_audit_author.pdbx_ordinal'] = []
  for author_index, author_name in enumerate(_MMCIF_AUTHORS, start=1):
    cif['_audit_author.name'].append(author_name)
    cif['_audit_author.pdbx_ordinal'].append(str(author_index))

  # Paper author details.
  cif.update(bfvd_util.get_citation_author(_CITATIONS))

#   # Paper citation details.
  cif.update(bfvd_util.get_citation(_CITATIONS))

  # Type of data in the dataset including data used in the model generation.
  cif['_ma_data.id'] = ['1']
  cif['_ma_data.name'] = ['Model']
  cif['_ma_data.content_type'] = ['model coordinates']

  # Description of number of instances for each entity.
  cif['_ma_target_entity_instance.asym_id'] = old_cif['_struct_asym.id']
  cif['_ma_target_entity_instance.entity_id'] = old_cif[
      '_struct_asym.entity_id'
  ]
  cif['_ma_target_entity_instance.details'] = ['.'] * len(
      cif['_ma_target_entity_instance.entity_id']
  )

  # Details about the target entities.
  cif['_ma_target_entity.entity_id'] = cif[
      '_ma_target_entity_instance.entity_id'
  ]
  cif['_ma_target_entity.data_id'] = ['1'] * len(
      cif['_ma_target_entity.entity_id']
  )
  cif['_ma_target_entity.origin'] = ['.'] * len(
      cif['_ma_target_entity.entity_id']
  )

  # DOING: _ma_target_ref_db_details
  cif.update(bfvd_util.get_ma_target_ref_db_details(old_cif))

#   cif.update(bfvd_util.get_pdbx_audit_revision(bfvd_constants._VERSION))
  
#   # Details of the models being deposited.
#   cif['_ma_model_list.ordinal_id'] = ['1']
#   cif['_ma_model_list.model_id'] = ['1']
#   cif['_ma_model_list.model_group_id'] = ['1']
#   cif['_ma_model_list.model_name'] = ['Top ranked model']

#   cif['_ma_model_list.model_group_name'] = [
#       f'ColabFold {model_type} v{bfvd_constants._SOFTWARE[0]["version"]} model'
#   ]
#   cif['_ma_model_list.data_id'] = ['1']
#   cif['_ma_model_list.model_type'] = ['Ab initio model']

#   # Software used.
#   cif.update(bfvd_util.get_software(bfvd_constants._SOFTWARE))

#   # Method description to conform with ModelCIF.
#   cif.update(bfvd_util.get_ma_protocol_step())

#   # Details of the metrics use to assess model confidence.
#   cif['_ma_qa_metric.id'] = ['1', '2']
#   cif['_ma_qa_metric.name'] = ['pLDDT', 'pLDDT']
#   # Accepted values are distance, energy, normalised score, other, zscore.
#   cif['_ma_qa_metric.type'] = ['pLDDT', 'pLDDT']
#   cif['_ma_qa_metric.mode'] = ['global', 'local']
#   cif['_ma_qa_metric.software_group_id'] = ['1', '1']

#   # Global model confidence metric value.
#   cif['_ma_qa_metric_global.ordinal_id'] = ['1']
#   cif['_ma_qa_metric_global.model_id'] = ['1']
#   cif['_ma_qa_metric_global.metric_id'] = ['1']
#   global_plddt = np.mean(
#       [float(v) for v in old_cif['_atom_site.B_iso_or_equiv']]
#   )
#   cif['_ma_qa_metric_global.metric_value'] = [f'{global_plddt:.2f}']

#   ma_qa_metric_local = bfvd_util.get_ma_qa_metric_local(old_cif)
#   cif.update(ma_qa_metric_local)

#   cif['_atom_type.symbol'] = sorted(set(old_cif['_atom_site.type_symbol']))

  return cif
