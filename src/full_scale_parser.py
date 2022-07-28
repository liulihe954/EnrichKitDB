import sys
import os
from src.gtf_parser import gtf_parser
from src.id_parser import id_parser
import pandas as pd
import sqlite3
import pandas as pd
import numpy as np
import shutil
from os.path import splitext, basename
import xml.etree.ElementTree as ET


class full_scale_parser():
    
    def __init__(self, BASE):
        self.BASE = BASE
        
    def do_parse(self):
        self.parse_gtf()
        print('\n')
        self.parse_id_mapper()
        print('\n')
        self.parse_go()
        print('\n')
        self.parse_interpro()
        print('\n')
        self.parse_kegg()
        print('\n')
        self.parse_reactome()
        print('\n')
        self.parse_mesh()
        print('\n')
        self.parse_msigdb()
        print('\n')
        
        
    def parse_gtf(self):
        source_path = os.path.join(self.BASE,'data','raw','annotation')
        source_list = os.listdir(source_path)

        for item in source_list:
            # form the import source
            input_path = os.path.join(source_path,item)
            # form the output source
            output_path = input_path.replace('raw','tmp').replace(item,item.split('.')[0].lower()) + '_gtf.txt'
            # parse
            tmp = gtf_parser(input_path).dataframe()
            # output to tmp
            tmp.to_csv(output_path, index = False, encoding='utf-8')
    
    def parse_id_mapper(self):
        #
        source_path = os.path.join(self.BASE,'data','raw','id_mapper')
        source_list = os.listdir(source_path)

        input_list = [None] * 3

        for source in source_list:
            tmp_path = os.path.join(source_path,source)
            if source.startswith('vgnc'): input_list[1] = tmp_path
            elif source.startswith('hgnc'): input_list[2] = tmp_path
            else: input_list[0] = tmp_path

        # instantiate
        parser = id_parser(input_list)

        # redirect output

        # gene2ensembl          
        parser.parseEntrezID().to_csv(input_list[0].replace('raw','tmp').replace('gz','txt'), index = False, encoding='utf-8')

        # vgnc
        parser.parseHumanOrtholog().to_csv(input_list[1].replace('raw','tmp').replace('gz','txt'), index = False, encoding='utf-8')

        # hgnc
        parser.parseHumanIDs().to_csv(input_list[2].replace('raw','tmp').replace('gz','txt'), index = False, encoding='utf-8')

        print('Done')
        print('')

    def parse_go(self):
        print('Now parse go.')
        source_path = os.path.join(self.BASE,'data','raw','go')
        output_path = source_path.replace('raw','tmp')

        # simply copy
        cmd = 'cp ' + os.path.join(source_path,'* ') + os.path.join(output_path,'')
        
        os.system(cmd)
        print('Done')
        print('')
    
    
    def parse_interpro(self):
        print('Now parse interpro.')
        source_path = os.path.join(self.BASE,'data','raw','interpro')
        output_path = source_path.replace('raw','tmp')

        # simply copy
        cmd = 'cp ' + os.path.join(source_path,'* ') + os.path.join(output_path,'')
        os.system(cmd)
        print('Done')
    
    def parse_kegg(self):
        print('Now parse kegg.')
        source_path = os.path.join(self.BASE,'data','raw','kegg')
        source_list = os.listdir(source_path)

        for item in source_list:
            output_path_tmp = os.path.join(source_path.replace('raw','tmp'),item)
            if 'gene_sets' in item:
                kegg_involve = []
                with open(os.path.join(source_path,item)) as f:
                    for line in f:
                        tmp = line.rstrip().split('\t')
                        kegg_involve.append([
                            tmp[0].split(':')[1],
                            tmp[1].split(':')[1]
                        ])
                #
                kegg_involve_df = pd.DataFrame(kegg_involve, columns = ['gene_id', 'pathway_id'])
                kegg_involve_df.to_csv(output_path_tmp, index = False, encoding='utf-8')
                print(f'Done parsing {item}')

            elif 'pathway_description' in item:
                kegg_pathway_list = []
                with open(os.path.join(source_path,item)) as f:
                    for line in f:
                        tmp = line.rstrip().split('\t')    
                        kegg_pathway_list.append([
                            tmp[0].split(':')[1],
                            tmp[1]
                        ])
                    #
                kegg_pathway_list_df = pd.DataFrame(kegg_pathway_list, columns = ['pathway_id', 'pathway_description'])
                kegg_pathway_list_df.to_csv(output_path_tmp, index = False, encoding='utf-8')

                print(f'Done parsing {item}')

    def parse_reactome(self):
        print('Now parse reactome.')
        source_path = os.path.join(self.BASE,'data','raw','reactome')
        source_list = os.listdir(source_path)

        for item in source_list:
            output_path_tmp = os.path.join(source_path.replace('raw','tmp'),item)
            reactome_involve = []
            with open(os.path.join(source_path,item)) as f:
                for line in f:
                    tmp = line.rstrip().split('\t')
                    reactome_involve.append([tmp[0],tmp[1],tmp[3],tmp[5]])
            print(f'Done parsing {item}')
            #
            reactome_involve_df = pd.DataFrame(reactome_involve, columns = ['gene_id', 'pathway_id','pathway_description','species'])
            reactome_involve_df.to_csv(output_path_tmp, index = False, encoding='utf-8')
    
    def parse_mesh(self):
        source_path = os.path.join(self.BASE,'data','raw','mesh')
        source_list = os.listdir(source_path)

        for item in source_list:
            db_path = os.path.join(source_path,item)
            output_path_tmp = os.path.join(source_path.replace('raw','tmp'),item.replace('sqlite','txt'))

            if item.split(".")[1] == 'db':
                print(f'Query data from {item}')
                #
                conn = sqlite3.connect(db_path)
                cursor = conn.cursor()
                cursor.execute("""SELECT MESHID, MESHTERM FROM DATA WHERE CATEGORY IN ('D','G');""")
                mesh_term_match = cursor.fetchall()
                mesh_term_match_df = pd.DataFrame(mesh_term_match, columns = ['pathway_id', 'pathway_description']).drop_duplicates()
                mesh_term_match_df.to_csv(output_path_tmp, index = False, encoding='utf-8')
                print(f'Obtained {len(mesh_term_match_df)} records from {item}')
            else:
                print(f'Query data from {item}')
                conn = sqlite3.connect(db_path)
                cursor = conn.cursor()
                cursor.execute("""SELECT GENEID, MESHID FROM DATA WHERE MESHCATEGORY IN ('D','G');""")
                involve_list = cursor.fetchall()
                involve_list_df = pd.DataFrame(involve_list, columns = ['gene_id','pathway_id'])
                involve_list_df.to_csv(output_path_tmp, index = False, encoding='utf-8')
                print(f'Obtained {len(involve_list_df)} records from {item}')

        print('Done!')
    
            
    def parse_msigdb(self):
        print('Now parse msigdb.')
        source_path = os.path.join(self.BASE,'data','raw','msigdb')
        source_list = os.listdir(source_path)

        for item in source_list:
            #
            input_path = os.path.join(source_path, item)
            output_path = input_path.replace('raw','tmp').replace('xml','txt')

            #
            tree = ET.parse(input_path)
            root = tree.getroot()
            cols = ['STANDARD_NAME','SYSTEMATIC_NAME','ORGANISM','EXTERNAL_DETAILS_URL','CHIP',
                    'CATEGORY_CODE','SUB_CATEGORY_CODE','CONTRIBUTOR','CONTRIBUTOR_ORG',
                    'DESCRIPTION_BRIEF','MEMBERS_EZID','DESCRIPTION_FULL']
            convert = [{attrib:xe.attrib[attrib] for attrib in cols} for xe in root]
            df = pd.DataFrame(columns=cols)
            df = df.append(convert,ignore_index = True)[['SYSTEMATIC_NAME','DESCRIPTION_BRIEF','MEMBERS_EZID']]

            # post process
            # 1. pathway id + pathway decription
            df_pathway_list = df[['SYSTEMATIC_NAME','DESCRIPTION_BRIEF']].rename(columns={"SYSTEMATIC_NAME": "pathway_id", "DESCRIPTION_BRIEF": "pathway_description"})

            print('pathway list done.')
            df_pathway_list.to_csv(output_path.replace('txt','pathway_list') + '.txt', index = False, encoding='utf-8')

            # 2. involvement list

            involve_list = []
            for index, row in df.iterrows():
                tmp_pathway_id = row[0]
                tmp_gene_set = row[2].split(',')
                for item in tmp_gene_set:
                    involve_list.append([item, tmp_pathway_id])
            print('involvement list done.')
            df_involve_list = pd.DataFrame(involve_list, columns = ['gene_id','pathway_id'])          

            df_involve_list.to_csv(output_path.replace('txt','involve_list') + '.txt', index = False, encoding='utf-8')
            print('Done.')
            print()        