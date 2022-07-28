from urllib.request import urlopen
from collections import defaultdict
import gzip
from os.path import exists
import pandas as pd
import certifi
import ssl

class id_parser():
    
    '''A class that can parse gene id from public databases.
    
    1. ensembl id to entrez id - ftp://ftp.ncbi.nih.gov/gene/DATA//gene2ensembl.gz
    2. ensembl to vgnc symbol: ftp.ebi.ac.uk/pub/databases/genenames/vgnc/tsv/all/locus_groups/all_protein-coding_gene_All.txt
    3. human orthologs - ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt
    
    return a single pandas dataframe - ID Mapper
    '''
    
    def __init__(self, item_list):
        # public database api
        self.gene2ensembl_local = item_list[0]
        self.vgnc_local = item_list[1]
        self.hgnc_local = item_list[2]
        
        #self.gene2ensembl = 'http://ftp.ncbi.nih.gov/gene/DATA//gene2ensembl.gz'
        #self.vgnc = 'http://ftp.ebi.ac.uk/pub/databases/genenames/vgnc/tsv/all/locus_groups/all_protein-coding_gene_All.txt'
        #self.hgnc = 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt' 
    
    
    def parseEntrezID(self):
        #
        entrez_mapper = defaultdict(list)


        print('Now working on ensembl id to entrez id using local copy...')
        print('About to read from ', self.gene2ensembl_local)
        
        #ssl._create_default_https_context = ssl._create_unverified_context
        #streamed_file = urlopen(self.gene2ensembl)
        #print('request status: ', streamed_file.status)
        
        # source
        opener = gzip.open if self.gene2ensembl_local.endswith('.gz') else open
        
        count_tmp = 0
        
        print('Started reading lines from', self.gene2ensembl_local)
        
        # process data
        with opener(self.gene2ensembl_local) as f:
            
            for line in f:     
                if line.decode('utf8').startswith('#'): continue
                else:
                    count_tmp += 1
                    if count_tmp % 500000 == 0: print(f'parsed {count_tmp} records...')
                    
                    fields = line.decode('utf8').rstrip().split('\t')
                    
                    tax_id = fields[0];
                    tmp_entrez = fields[1];
                    tmp_ensembl = fields[2]

                    if (tmp_ensembl.startswith('')):
                        
                        if len(entrez_mapper[(tax_id,tmp_ensembl)]) == 0:    
                            
                            entrez_mapper[(tax_id,tmp_ensembl)].append(tmp_entrez)
                            
                        else:
                            
                            if tmp_entrez in entrez_mapper[(tax_id,tmp_ensembl)]: continue
                                
                            else:
                                entrez_mapper[(tax_id,tmp_ensembl)].append(tmp_entrez)
        ##
        output = {
            'tax_id': [],
            'gene_id': [],
            'entrez_id': []
        }

        for k, v in entrez_mapper.items():
            for i in range(len(v)):
                output['tax_id'].append(k[0])
                output['gene_id'].append(k[1])
            output['entrez_id'].extend(v)

        print('Done!')
        print()
        
        return pd.DataFrame.from_dict(output)
    
    #
    def parseHumanOrtholog(self):
        
        # container
        vgnc_mapper = {
        'gene_id': [],
        'vgnc_id': [],
        'vgnc_symbol': [],
        'hgnc_orthologs': []
        }
        
        #
        print('Now working on ensembl id to vgnc symbol and human ortholog...')
        
        #print('About to start streaming from ', self.vgnc)
        #ssl._create_default_https_context = ssl._create_unverified_context
        #streamed_file = urlopen(self.vgnc)
        #print('request status: ', streamed_file.status)
        
        print('About to start read from ', self.vgnc_local)
        
        opener = gzip.open if self.vgnc_local.endswith('.gz') else open
        
        count_tmp = 0
        
        print('Started reading lines from', self.vgnc_local)
        
        # process data
        with opener(self.vgnc_local) as f:
            
            for line in f:

                count_tmp += 1
               
                if count_tmp % 10000 == 0: print(f'parsed {count_tmp} records...')

                tmp_line = line.split('\t')

                if tmp_line[0] == 'taxon_id': continue

                vgnc_mapper['gene_id'].append(tmp_line[20])
                vgnc_mapper['vgnc_id'].append(tmp_line[1])
                vgnc_mapper['vgnc_symbol'].append(tmp_line[2])
                vgnc_mapper['hgnc_orthologs'].append(tmp_line[25].replace('"','').rstrip())

        print('Done!')
        print()
        
        return pd.DataFrame.from_dict(vgnc_mapper)
    

    def parseHumanIDs(self):
        
        hgnc_ref = {
            'hgnc_orthologs': [],
            'human_gene_id': [],
            'human_entrez_id': [],
            'hgnc_symbol': []}
        
        
        print('Now working on retrieving extra human ortholog info...')
        #print('About to start streaming from ', self.hgnc)
        #ssl._create_default_https_context = ssl._create_unverified_context
        #streamed_file = urlopen(self.hgnc)
        #print('request status: ', streamed_file.status)
        
        print('About to start read from ', self.hgnc_local)
        
        opener = gzip.open if self.hgnc_local.endswith('.gz') else open
        
        count_tmp = 0
        
        print('Started reading lines from', self.hgnc_local)
        
        # process data
        with opener(self.hgnc_local) as f:
            
            for line in f:
                count_tmp += 1

                if count_tmp % 10000 == 0: print(f'parsed {count_tmp} records...')

                tmp_line = line.split('\t')

                if tmp_line[0] == 'hgnc_id': continue

                hgnc_ref['hgnc_orthologs'].append(tmp_line[0])
                hgnc_ref['human_gene_id'].append(tmp_line[19])
                hgnc_ref['human_entrez_id'].append(tmp_line[18])
                hgnc_ref['hgnc_symbol'].append(tmp_line[1])
        
        print('Done!')
        print()
        
        return pd.DataFrame.from_dict(hgnc_ref)
    
#     def attach_extra_info(self):
        
#         IDMapper = pd.DataFrame().assign(
#             ek_id = self.base_gene_list['ek_id'],
#             gene_id = self.base_gene_list['gene_id']
#         )
        
#         print('Now try to retrieve information.')
#         print('')
        
#         #
#         entrezid_out = self.parseEntrezID()
        
#         human_ortholog_out = self.parseHumanOrtholog()
        
#         hgnc_ref_out = self.parseHumanIDs()
        
#         #
#         IDMapper = IDMapper.merge(entrezid_out, on = 'gene_id', how = 'left')
#         IDMapper = IDMapper.merge(human_ortholog_out, on = 'gene_id', how = 'left')
#         IDMapper = IDMapper.merge(hgnc_ref_out, on = 'hgnc_orthologs', how = 'left')
        
#         print('Done!')
        
#         return IDMapper