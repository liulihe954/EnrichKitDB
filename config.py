UPDATE = False

RAW_DATA_INFO = {

    # genome annotation (ensembl)
    'annotation': {
        # Cow -  Bos taurus
        'Bos_taurus.ARS-UCD1.2.107.gtf.gz': "http://ftp.ensembl.org/pub/release-107/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.107.gtf.gz",
        # Goat - Capra hircus
        'Capra_hircus.ARS1.107.gtf.gz': "http://ftp.ensembl.org/pub/release-107/gtf/capra_hircus/Capra_hircus.ARS1.107.gtf.gz",
        # Sheep - Ovis aries
        'Ovis_aries.Oar_v3.1.107.gtf.gz': "http://ftp.ensembl.org/pub/release-107/gtf/ovis_aries/Ovis_aries.Oar_v3.1.107.gtf.gz",
        # Chicken - Gallus gallus
        'Gallus_gallus_gca000002315v5.GRCg6a.107.gtf.gz': "http://ftp.ensembl.org/pub/release-107/gtf/gallus_gallus_gca000002315v5/Gallus_gallus_gca000002315v5.GRCg6a.107.gtf.gz",
        # Pig - Duroc - Sus scrofa
        'Sus_scrofa.Sscrofa11.1.107.gtf.gz': "http://ftp.ensembl.org/pub/release-107/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.107.gtf.gz",
        # Horse Equus caballus
        'Equus_caballus.EquCab3.0.107.gtf.gz': "http://ftp.ensembl.org/pub/release-107/gtf/equus_caballus/Equus_caballus.EquCab3.0.107.gtf.gz",

    },

    # id mapper information
    'id_mapper': {

        # match ensembl gene id to entrez id
        'gene2ensembl.gz': 'http://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz',

        # match entrez id to gene symbol name
        'gene_info.gz': 'http://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz',

        # Vertebrate Gene Nomenclature Committee: match vertebrate ensembl gene id to official vgnc gene symbol
        'vgnc_protein_coding_gene.txt': 'http://ftp.ebi.ac.uk/pub/databases/genenames/vgnc/tsv/all/locus_groups/all_protein-coding_gene_All.txt',

        # HUGO Gene Nomenclature Committee: match vertebrate ensembl gene id to human orthologous (ensembl + entrez + symbol)
        'hgnc_protein_coding_gene.txt': 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt',

        # Cow -  Bos taurus
        'bta_ensembl.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "btaurus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>\'',
        # Goat - Capra hircus
        'cap_ensembl.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "chircus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>\'',
        # Sheep - Ovis aries
        'ovi_ensembl.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "oaries_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>\'',
        # Chicken - Gallus gallus
        'gal_ensembl.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "ggallus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>\'',
        # Pig - Duroc - Sus scrofa
        'sus_ensembl.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "sscrofa_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>\'',
        # Horse Equus caballus
        'equ_ensembl.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "ecaballus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>\'',
    },

    # go - using biomart api
    'go': {
        # Cow -  Bos taurus
        'go_bta.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "btaurus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "go_id" /><Attribute name = "name_1006" /></Dataset></Query>\'',
        # Goat - Capra hircus
        'go_cap.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "chircus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "go_id" /><Attribute name = "name_1006" /></Dataset></Query>\'',
        # Sheep - Ovis aries
        'go_ovi.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "oaries_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "go_id" /><Attribute name = "name_1006" /></Dataset></Query>\'',
        # Chicken - Gallus gallus
        'go_gal.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "ggallus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "go_id" /><Attribute name = "name_1006" /></Dataset></Query>\'',
        # Pig - Duroc - Sus scrofa
        'go_sus.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "sscrofa_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "go_id" /><Attribute name = "name_1006" /></Dataset></Query>\'',
        # Horse Equus caballus
        'go_equ.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "ecaballus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "go_id" /><Attribute name = "name_1006" /></Dataset></Query>\'',

    },

    # interpro - using biomart api
    'interpro': {
        # Cow -  Bos taurus
        'interpro_bta.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "btaurus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "interpro" /><Attribute name = "interpro_description" /></Dataset></Query>\'',
        # Goat - Capra hircus
        'interpro_cap.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "chircus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "interpro" /><Attribute name = "interpro_description" /></Dataset></Query>\'',
        # Sheep - Ovis aries
        'interpro_ovi.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "oaries_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "interpro" /><Attribute name = "interpro_description" /></Dataset></Query>\'',
        # Chicken - Gallus gallus
        'interpro_gal.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "ggallus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "interpro" /><Attribute name = "interpro_description" /></Dataset></Query>\'',
        # Pig - Duroc - Sus scrofa
        'interpro_sus.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "sscrofa_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "interpro" /><Attribute name = "interpro_description" /></Dataset></Query>\'',
        # Horse Equus caballus
        'interpro_equ.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "ecaballus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "interpro" /><Attribute name = "interpro_description" /></Dataset></Query>\'',

    },

    # kegg
    'kegg': {
        # Cow -  Bos taurus
        'bta_gene_sets.txt': 'https://rest.kegg.jp/link/pathway/bta',  # gene sets: link gene to pathway
        'bta_pathway_description.txt': 'https://rest.kegg.jp/list/pathway/bta',  # kegg gene sets: match pathway id to readable term descriptions
        # Goat - Capra hircus
        'cap_gene_sets.txt': 'https://rest.kegg.jp/link/pathway/chx',  # gene sets: link gene to pathway
        'cap_pathway_description.txt': 'https://rest.kegg.jp/list/pathway/chx',  # kegg gene sets: match pathway id to readable term descriptions
        # Sheep - Ovis aries
        'ovi_gene_sets.txt': 'https://rest.kegg.jp/link/pathway/oas',  # gene sets: link gene to pathway
        'ovi_pathway_description.txt': 'https://rest.kegg.jp/list/pathway/oas',  # kegg gene sets: match pathway id to readable term descriptions
        # Chicken - Gallus gallus
        'gal_gene_sets.txt': 'https://rest.kegg.jp/link/pathway/gga',  # gene sets: link gene to pathway
        'gal_pathway_description.txt': 'https://rest.kegg.jp/list/pathway/gga',  # kegg gene sets: match pathway id to readable term descriptions
        # Pig - Duroc - Sus scrofa
        'sus_gene_sets.txt': 'https://rest.kegg.jp/link/pathway/ssc',  # gene sets: link gene to pathway
        'sus_pathway_description.txt': 'https://rest.kegg.jp/list/pathway/ssc',  # kegg gene sets: match pathway id to readable term descriptions
        # Horse Equus caballus
        'equ_gene_sets.txt': 'https://rest.kegg.jp/link/pathway/ecb',  # gene sets: link gene to pathway
        'equ_pathway_description.txt': 'https://rest.kegg.jp/list/pathway/ecb',  # kegg gene sets: match pathway id to readable term descriptions

    },

    # reactome gene sets: gene + pathway + descriptions
    'reactome': {
        'NCBI2Reactome_All_Levels.txt': 'https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt',
    },

    # mesh
    'mesh': {
        'MeSH.db.sqlite': 'https://bioconductorhubs.blob.core.windows.net/annotationhub/AHMeSHDbs/v001/MeSH.db.sqlite',  # gene sets: gene id link to mesh id in bta

        # Cow -  Bos taurus
        'MeSH.Bta.eg.db.sqlite': 'https://bioconductorhubs.blob.core.windows.net/annotationhub/AHMeSHDbs/v001/MeSH.Bta.eg.db.sqlite',  # mesh gene sets: mesh id + term description
        # Goat - Capra hircus
        # NA

        # Sheep - Ovis aries
        # NA

        # Chicken - Gallus gallus
        'MeSH.Gga.eg.db.sqlite': 'https://bioconductorhubs.blob.core.windows.net/annotationhub/AHMeSHDbs/v001/MeSH.Gga.eg.db.sqlite',  # mesh gene sets: mesh id + term description

        # Pig - Duroc - Sus scrofa
        'MeSH.Ssc.eg.db.sqlite': 'https://bioconductorhubs.blob.core.windows.net/annotationhub/AHMeSHDbs/v001/MeSH.Ssc.eg.db.sqlite',  # mesh gene sets: mesh id + term description

        # Horse Equus caballus
        'MeSH.Eca.eg.db.sqlite': 'https://bioconductorhubs.blob.core.windows.net/annotationhub/AHMeSHDbs/v001/MeSH.Eca.eg.db.sqlite',  # mesh gene sets: mesh id + term description

    },

    # msigdb gene sets
    'msigdb': {
        'msigdb_v7.5.1.xml': 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/msigdb_v7.5.1.xml',
    },

    # tftargets
    'tftargets':{
        'tftargets.rda':'https://github.com/slowkow/tftargets/raw/master/data/tftargets.rda',
    }
}
