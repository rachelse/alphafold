CODEN = {
  # JOURNAL_ID_ASTM
    "Nucleic Acids Research": "NARHAD",
    "Nature": "NATUAS",
    "Nature Methods": "NAMNBI",
}

PUBLISHER = {
    "Nucleic Acids Research": "Oxford University Press",
    "Nature": "Nature Publishing Group",
    "Nature Methods": "Nature Publishing Group",
}

COUNTRY = {
    "Nucleic Acids Research": "UK",
    "Nature": "UK",
    "Nature Methods": "US",
}

CSD = {
    "Nature": "0006",

}

_BFVD_CITATION = {
    "title": "BFVD—a large repository of predicted viral protein structures",
    "journal_full": "Nucleic Acids Research",
    "journal_volume": "?",
    "page_first": "gkae1119",
    "page_last": "gkae1119",
    "year": "2024",
    "pdbx_database_id_PubMed": "39574394",
    "pdbx_database_id_DOI": "10.1093/nar/gkae1119",
    "journal_id_ISSN": "0305-1048",
    "journal_id_ASTM": CODEN["Nucleic Acids Research"],
    "journal_id_CSD": CSD.get("Nucleic Acids Research", "?"),
    "country": COUNTRY["Nucleic Acids Research"],
    "book_publisher": PUBLISHER["Nucleic Acids Research"],
    "authors": (
        'Kim, Rachel Seongeun',
        'Levy Karin, Eli',
        'Mirdita, Milot',
        'Chikhi, Rayan',
        'Steinegger, Martin',
    )
}

_COLABFOLD_CITATION = {
    "title": "ColabFold: Making Protein folding accessible to all",
    "journal_full": "Nature Methods",
    "journal_volume": "19",
    "page_first": "679",
    "page_last": "682",
    "year": "2022",
    "pdbx_database_id_PubMed": "35637307",
    "pdbx_database_id_DOI": "10.1038/s41586-021-03819-2",
    "journal_id_ISSN": "1548-7105",
    "journal_id_ASTM": CODEN["Nature Methods"],
    "journal_id_CSD": CSD.get("Nature Methods", "?"),
    "country": COUNTRY["Nature Methods"],
    "book_publisher": PUBLISHER["Nature Methods"],
    "authors": (
        'Mirdita, Milot',
        'Schütze, Konstantin',
        'Moriwaki, Yoshitaka',
        'Heo, Lim',
        'Ovchinnikov, Sergey',
        'Steinegger, Martin',
    )
}

_ALPHAFOLD_CITATION = {
    "title": "Highly accurate protein structure prediction with AlphaFold",
    "journal_full": "Nature",
    "journal_volume": "596",
    "page_first": "583",
    "page_last": "589",
    "year": "2021",
    "pdbx_database_id_PubMed": "34265844",
    "pdbx_database_id_DOI": "10.1038/s41586-021-03819-2",
    "journal_id_ISSN": "0028-0836",
    "journal_id_ASTM": CODEN["Nature"],
    "journal_id_CSD": CSD.get("Nature", "?"),
    "country": COUNTRY["Nature"],
    "book_publisher": PUBLISHER["Nature"],
    "authors": (
        'Jumper, John',
        'Evans, Richard',
        'Pritzel, Alexander',
        'Green, Tim',
        'Figurnov, Michael',
        'Ronneberger, Olaf',
        'Tunyasuvunakool, Kathryn',
        'Bates, Russ',
        'Žídek, Augustin',
        'Potapenko, Anna',
        'Bridgland, Alex',
        'Meyer, Clemens',
        'Kohl, Simon A. A.',
        'Ballard, Andrew J.',
        'Cowie, Andrew',
        'Romera-Paredes, Bernardino',
        'Nikolov, Stanislav',
        'Jain, Rishub',
        'Adler, Jonas',
        'Back, Trevor',
        'Petersen, Stig',
        'Reiman, David',
        'Clancy, Ellen',
        'Zielinski, Michal',
        'Steinegger, Martin',
        'Pacholska, Michalina',
        'Berghammer, Tamas',
        'Silver, David',
        'Vinyals, Oriol',
        'Senior, Andrew W.',
        'Kavukcuoglu, Koray',
        'Kohli, Pushmeet',
        'Hassabis, Demis',
    )
}

VERSION = [
  {
    'data_content_type': 'Structure model', 'provider': 'BFVD', 
    'major_revision': '1', 'minor_revision': '0', 'revision_date': '2024-11-01', 
    'type': 'Initial release', 'description': 'Initial release'
    },
]