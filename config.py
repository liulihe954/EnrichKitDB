UPDATE = False

RAW_DATA_INFO = {

# genome annotation (ensembl)
'annotation': {
    'Bos_taurus.ARS-UCD1.2.107.gtf.gz':"http://ftp.ensembl.org/pub/release-107/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.107.gtf.gz",
},

# id mapper information
'id_mapper': {

    # match ensembl gene id to entrez id
    'gene2ensembl.gz': 'http://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz',

    # Vertebrate Gene Nomenclature Committee: match vertebrate ensembl gene id to official vgnc gene symbol
   'vgnc_protein_coding_gene.txt':'http://ftp.ebi.ac.uk/pub/databases/genenames/vgnc/tsv/all/locus_groups/all_protein-coding_gene_All.txt',

    # HUGO Gene Nomenclature Committee: match vertebrate ensembl gene id to human orthologous (ensembl + entrez + symbol)
   'hgnc_protein_coding_gene.txt':'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt',
},

# go - using biomart api
'go': {
    'go_bta.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "btaurus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "go_id" /><Attribute name = "name_1006" /></Dataset></Query>\'',

},

# interpro - using biomart api
'interpro': {
    'interpro_bta.txt': '\'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "btaurus_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "interpro" /><Attribute name = "interpro_description" /></Dataset></Query>\'',
},

# kegg 
'kegg':{
    # gene sets: link gene to pathway
    'bta_gene_sets.txt':'https://rest.kegg.jp/link/pathway/bta',

    # kegg gene sets: match pathay id to readable term descriptions
    'bta_pathway_description.txt': 'https://rest.kegg.jp/list/pathway/bta'
},

# reactome gene sets: gene + pathway + descriptions
'reactome': {
    'NCBI2Reactome_All_Levels.txt':'https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt',
},

# mesh 
'mesh': {

    # gene sets: gene id link to mesh id in bta
    'MeSH.db.sqlite':'https://bioconductorhubs.blob.core.windows.net/annotationhub/AHMeSHDbs/v001/MeSH.db.sqlite',

    # mesh gene sets: mesh id + term description
    'MeSH.Bta.eg.db.sqlite':'https://bioconductorhubs.blob.core.windows.net/annotationhub/AHMeSHDbs/v001/MeSH.Bta.eg.db.sqlite',
},

# msigdb gene sets
'msigdb': {
    'msigdb_v7.5.1.xml':'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/msigdb_v7.5.1.xml',
}
}
