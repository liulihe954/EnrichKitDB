import sys
import pandas as pd
import numpy as np
import shutil
from os.path import splitext, basename
import os
from src.gtf_parser import gtf_parser
from src.id_parser import id_parser
import pandas as pd
import sqlite3
import requests
import xml.etree.ElementTree as ET
from src.gtf_splitter import gtf_splitter


class full_scale_parser():

    def __init__(self, BASE):
        self.BASE = BASE
        self.ABBR_DICT = {
            'Bos_taurus': 'bta',
            'Capra_hircus': 'cap',
            'Equus_caballus': 'equ',
            'Gallus_gallus_gca000002315v5': 'gal',
            'Ovis_aries': 'ovi',
            'Sus_scrofa': 'sus'
        }
        self.SPECIES_DICT = {
            'bta': ['9913', 'ENSBTAG', 'Bta', 'Bos taurus', 0],
            'cap': ['9925', 'ENSCHIG', None, None, 1],
            'equ': ['9796', 'ENSECAG', 'Eca', None, 2],
            'gal': ['9031', 'ENSGALG', 'Gga', 'Gallus gallus', 3],
            'ovi': ['9940', 'ENSOARG', None, None, 4],
            'sus': ['9823', 'ENSSSCG', 'Ssc', 'Sus scrofa', 5]
        }

    def do_parse(self):
        self.parse_id_mapper()
        print('\n')

        self.parse_gtf()
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
        source_path = os.path.join(self.BASE, 'data', 'raw', 'annotation')
        source_list = os.listdir(source_path)

        for item in source_list:
            #
            species_tmp = self.ABBR_DICT[item.split('.')[0]]
            # form the import source
            input_path = os.path.join(source_path, item)

            # parse
            tmp = gtf_parser(input_path).dataframe()

            # split and format
            spl = gtf_splitter(self.BASE, species_tmp, tmp)

            tmp_skinny = spl.migrate_tables()

#             tmp_skinny = pd.read_csv(os.path.join(self.BASE, 'data', 'tmp', 'annotation', ''.join(['genes_', species_tmp, '.csv'])))

            # output id mapper
            print('About to compose id mapper for ' + species_tmp)

            self.compose_id_mapper(species_tmp, tmp_skinny)

            print('Compose id mapper for ' + species_tmp, + ' done!')

    def compose_id_mapper(self, species, gtf_skinny):

        # all source paths
        gene2ensembl_path_tmp = os.path.join(self.BASE, 'data', 'tmp', 'id_mapper', 'gene2ensembl.txt')
        vgnc_path_tmp = os.path.join(self.BASE, 'data', 'tmp', 'id_mapper', 'vgnc_protein_coding_gene.txt')
        hgnc_path_tmp = os.path.join(self.BASE, 'data', 'tmp', 'id_mapper', 'hgnc_protein_coding_gene.txt')

        if species in self.SPECIES_DICT:

            # entrez id: differentiate the given SPECIES by tax id
            gene2ensembl_iter_csv = pd.read_csv(gene2ensembl_path_tmp, iterator=True, chunksize=500000)
            df_gene2ensembl_tmp = pd.concat([chunk[chunk['tax_id'] == int(self.SPECIES_DICT[species][0])] for chunk in gene2ensembl_iter_csv])

            #
            vgnc_iter_csv = pd.read_csv(vgnc_path_tmp, iterator=True, chunksize=10000)
            df_vgnc_tmp = pd.concat([chunk[chunk['gene_id'].str.contains(self.SPECIES_DICT[species][1], na=False)] for chunk in vgnc_iter_csv])

            #
        df_hgnc_tmp = pd.read_csv(hgnc_path_tmp)

        # take union as the index
        tmp1 = pd.DataFrame(gtf_skinny['gene_id'].drop_duplicates())
        tmp1['source'] = 4

        tmp2 = pd.DataFrame(df_gene2ensembl_tmp['gene_id'].drop_duplicates())
        tmp2['source'] = 2

        tmp3 = pd.DataFrame(df_vgnc_tmp['gene_id'].drop_duplicates())
        tmp3['source'] = 1

        dfs = [tmp1, tmp2, tmp3]
        gene_id_idx = pd.concat(dfs)
        gene_id_out = pd.DataFrame(gene_id_idx.groupby(by=['gene_id'])['source'].sum())

        # merge
        id_mapper = pd.merge(gene_id_out, df_gene2ensembl_tmp[['gene_id', 'entrez_id']], on='gene_id', how='left')
        id_mapper = pd.merge(id_mapper, df_vgnc_tmp, on='gene_id', how='left')
        id_mapper = pd.merge(id_mapper, df_hgnc_tmp, on='hgnc_orthologs', how='left')

        # # post processing
        # id_mapper['entrez_id'] = id_mapper['entrez_id'].astype('Int64')
        # id_mapper['human_entrez_id'] = id_mapper['human_entrez_id'].astype('Int64')

        # append two extra columns - gene symbol from ensembl and gene symbol from ncbi
        id_mapper_old = id_mapper
        ncbi_symbol_df = pd.read_csv(os.path.join(self.BASE, 'data', 'tmp', 'id_mapper', ''.join([species, '_ncbi.txt'])))
        ensembl_symbol_df = pd.read_csv(os.path.join(self.BASE, 'data', 'tmp', 'id_mapper', ''.join([species, '_ensembl.txt'])))

        id_mapper_new = pd.merge(id_mapper_old, ensembl_symbol_df, on='gene_id', how='left')
        id_mapper_final = pd.merge(id_mapper_new, ncbi_symbol_df, on='entrez_id', how='outer')
        id_mapper_final.pop('source')

        ensembl_symbol_col = id_mapper_final.pop('ensembl_symbol')
        ncbi_symbol_col = id_mapper_final.pop('ncbi_symbol')

        id_mapper_final.insert(1, 'ensembl_symbol', ensembl_symbol_col)
        id_mapper_final.insert(3, 'ncbi_symbol', ncbi_symbol_col)
        id_mapper_final.insert(0, "species", self.SPECIES_DICT[species][4])

        # id_mapper_final['source'] = id_mapper_final['source'].astype('Int64')
        id_mapper_final['species'] = id_mapper_final['species'].astype('Int64')
        id_mapper_final['entrez_id'] = id_mapper_final['entrez_id'].astype('Int64')
        id_mapper_final['entrez_id'] = id_mapper_final['entrez_id'].astype('Int64')
        id_mapper_final['human_entrez_id'] = id_mapper_final['human_entrez_id'].astype('Int64')

        id_mapper_final.to_csv(os.path.join(self.BASE, 'data', 'tmp', 'id_mapper', ''.join(['id_mapper_', species, '.txt'])), index=False, encoding='utf-8')

    def check_ensembl_symbol_helper(self, gene_id):
        print('checking ', gene_id)

        try:
            server = "https://rest.ensembl.org"
            ext = ''.join(['/xrefs/id/', gene_id, '?'])

            r = requests.get(server + ext, headers={"Content-Type": "application/json"})

            if not r.ok:
                r.raise_for_status()
                sys.exit()

            decoded = r.json()

            symbol_tmp = decoded[0]['display_id']

            if symbol_tmp is not None and symbol_tmp != gene_id:
                return symbol_tmp
            else:
                return None
        except Exception:
            return None

    def parse_id_mapper(self):
        #
        source_path = os.path.join(self.BASE, 'data', 'raw', 'id_mapper')
        source_list = os.listdir(source_path)

        input_list = [None] * 4

        for source in source_list:
            tmp_path = os.path.join(source_path, source)
            if source.startswith('vgnc'):
                input_list[1] = tmp_path
            elif source.startswith('hgnc'):
                input_list[2] = tmp_path
            elif source.startswith('gene_info'):
                input_list[3] = tmp_path
            elif source.startswith('gene2ensembl'):
                input_list[0] = tmp_path
            elif source.endswith('ensembl.txt'):
                print('now formatting ensembl symbols for ' + source)
                tmp = pd.read_csv(tmp_path, sep='\t', header=None)
                tmp.columns = ['gene_id', 'ensembl_symbol']
                tmp.to_csv(tmp_path.replace('raw', 'tmp'), index=False, encoding='utf-8')

        # instantiate
        parser = id_parser(input_list, self.BASE)

        # redirect output

        # gene2ensembl
        parser.parseEntrezID().to_csv(input_list[0].replace('raw', 'tmp').replace('gz', 'txt'), index=False, encoding='utf-8')

        # vgnc
        parser.parseHumanOrtholog().to_csv(input_list[1].replace('raw', 'tmp').replace('gz', 'txt'), index=False, encoding='utf-8')

        # hgnc
        parser.parseHumanIDs().to_csv(input_list[2].replace('raw', 'tmp').replace('gz', 'txt'), index=False, encoding='utf-8')

        # gene_info: ncbi to gene symbol
        parser.parseGeneInfo()

        print('Done')
        print('')

    def parse_go(self):
        print('Now parse go.')
        source_path = os.path.join(self.BASE, 'data', 'raw', 'go')
        output_path = source_path.replace('raw', 'tmp')

        for item in os.listdir(source_path):
            prefix_tmp = item.split('_')[1][: 3]

            print('parse go: currently parsing ' + item)

            tmp = pd.read_csv(os.path.join(source_path, item), sep='\t', header=None).dropna().drop_duplicates()
            tmp.columns = ['gene_id', 'pathway_id', 'pathway_description']

            tmp_involve = tmp[['gene_id', 'pathway_id']].dropna().drop_duplicates()
            tmp_pathway_description = tmp[['pathway_id', 'pathway_description']].dropna().drop_duplicates()

            tmp_involve.to_csv(os.path.join(output_path, ''.join([prefix_tmp, '_involve.txt'])), index=False, encoding='utf-8')
            tmp_pathway_description.to_csv(os.path.join(output_path, ''.join([prefix_tmp, '_pathway_description.txt'])), index=False, encoding='utf-8')

        # # simply copy
        # cmd = 'cp ' + os.path.join(source_path, '* ') + os.path.join(output_path, '')

        # os.system(cmd)
        print('Done')
        print('')

    def parse_interpro(self):
        print('Now parse interpro.')
        source_path = os.path.join(self.BASE, 'data', 'raw', 'interpro')
        output_path = source_path.replace('raw', 'tmp')

        for item in os.listdir(source_path):
            prefix_tmp = item.split('_')[1][: 3]

            print('parse interpro: currently parsing ' + item)
            tmp = pd.read_csv(os.path.join(source_path, item), sep='\t', header=None).dropna().drop_duplicates()

            tmp.columns = ['gene_id', 'pathway_id', 'pathway_description']

            tmp_involve = tmp[['gene_id', 'pathway_id']].dropna().drop_duplicates()
            tmp_pathway_description = tmp[['pathway_id', 'pathway_description']].dropna().drop_duplicates()

            tmp_involve.to_csv(os.path.join(output_path, ''.join([prefix_tmp, '_involve.txt'])), index=False, encoding='utf-8')

            tmp_pathway_description.to_csv(os.path.join(output_path, ''.join([prefix_tmp, '_pathway_description.txt'])), index=False, encoding='utf-8')

        # # simply copy
        # cmd = 'cp ' + os.path.join(source_path, '* ') + os.path.join(output_path, '')
        # os.system(cmd)

        print('Done')

    def parse_kegg(self):
        print('Now parse kegg.')
        source_path = os.path.join(self.BASE, 'data', 'raw', 'kegg')
        source_list = os.listdir(source_path)

        for item in source_list:
            output_path_tmp = os.path.join(source_path.replace('raw', 'tmp'), item)
            if 'gene_sets' in item:
                kegg_involve = []
                with open(os.path.join(source_path, item)) as f:
                    for line in f:
                        tmp = line.rstrip().split('\t')
                        kegg_involve.append([
                            tmp[0].split(':')[1],
                            tmp[1].split(':')[1]
                        ])
                #
                kegg_involve_df = pd.DataFrame(kegg_involve, columns=['gene_id', 'pathway_id'])
                kegg_involve_df.to_csv(output_path_tmp.replace('sets', 'involve'), index=False, encoding='utf-8')
                print(f'Done parsing {item}')

            elif 'pathway_description' in item:
                kegg_pathway_list = []
                with open(os.path.join(source_path, item)) as f:
                    for line in f:
                        tmp = line.rstrip().split('\t')
                        kegg_pathway_list.append([
                            tmp[0],
                            tmp[1]
                        ])
                    #
                kegg_pathway_list_df = pd.DataFrame(kegg_pathway_list, columns=['pathway_id', 'pathway_description'])
                kegg_pathway_list_df.to_csv(output_path_tmp, index=False, encoding='utf-8')

                print(f'Done parsing {item}')

    def parse_reactome(self):
        print('Now parse reactome.')
        source_path = os.path.join(self.BASE, 'data', 'raw', 'reactome')
        source_list = os.listdir(source_path)

        for item in source_list:
            reactome_involve = []
            with open(os.path.join(source_path, item)) as f:
                for line in f:
                    tmp = line.rstrip().split('\t')
                    reactome_involve.append([tmp[0], tmp[1], tmp[3], tmp[5]])

            #
            reactome_involve_df = pd.DataFrame(reactome_involve, columns=['gene_id', 'pathway_id', 'pathway_description', 'species'])

            species_category = [('Bos taurus', 'bta'), ('Gallus gallus', 'gal'), ('Sus scrofa', 'sus')]

            for item in species_category:
                print(f'Now parsing {item}')
                species_tmp = item[0]
                reactome_df_tmp = reactome_involve_df[reactome_involve_df['species'] == species_tmp]
                #
                reactome_df_tmp_involve = reactome_df_tmp[['gene_id', 'pathway_id']].dropna().drop_duplicates()
                reactome_df_tmp_involve.to_csv(os.path.join(source_path.replace('raw', 'tmp'), ''.join([item[1], '_involve.txt'])), index=False, encoding='utf-8')
                #
                reactome_df_tmp_pathway_description = reactome_df_tmp[['pathway_id', 'pathway_description']].dropna().drop_duplicates()
                reactome_df_tmp_pathway_description.to_csv(os.path.join(source_path.replace('raw', 'tmp'), ''.join([item[1], '_pathway_description.txt'])), index=False, encoding='utf-8')
            print(f'Done parsing {item}')

    def parse_mesh(self):
        source_path = os.path.join(self.BASE, 'data', 'raw', 'mesh')
        source_list = os.listdir(source_path)

        for item in source_list:

            if item.split('.')[1] == 'Ssc':
                prefix_tmp = 'ssc'
            elif item.split('.')[1] == 'Bta':
                prefix_tmp = 'bta'
            elif item.split('.')[1] == 'Eca':
                prefix_tmp = 'equ'
            elif item.split('.')[1] == 'Gga':
                prefix_tmp = 'gal'
            else:
                prefix_tmp = 'all'

            db_path = os.path.join(source_path, item)

            if item.split(".")[1] == 'db':
                print(f'Query data from {item}')
                #
                conn = sqlite3.connect(db_path)
                cursor = conn.cursor()
                cursor.execute("""SELECT MESHID, MESHTERM FROM DATA WHERE CATEGORY IN ('D','G');""")
                mesh_term_match = cursor.fetchall()
                mesh_term_match_df = pd.DataFrame(mesh_term_match, columns=['pathway_id', 'pathway_description']).drop_duplicates()

                output_path_tmp = os.path.join(source_path.replace('raw', 'tmp'), ''.join([prefix_tmp, '_pathway_description.txt']))

                mesh_term_match_df.to_csv(output_path_tmp, index=False, encoding='utf-8')
                print(f'Obtained {len(mesh_term_match_df)} records from {item}')
                conn.close()
            else:
                print(f'Query data from {item}')
                conn = sqlite3.connect(db_path)
                cursor = conn.cursor()
                cursor.execute("""SELECT GENEID, MESHID FROM DATA WHERE MESHCATEGORY IN ('D','G');""")
                involve_list = cursor.fetchall()
                involve_list_df = pd.DataFrame(involve_list, columns=['gene_id', 'pathway_id'])

                output_path_tmp = os.path.join(source_path.replace('raw', 'tmp'), ''.join([prefix_tmp, '_involve.txt']))

                involve_list_df.to_csv(output_path_tmp, index=False, encoding='utf-8')
                print(f'Obtained {len(involve_list_df)} records from {item}')
                conn.close()

        print('Done!')

    def parse_msigdb(self):
        print('Now parse msigdb.')
        source_path = os.path.join(self.BASE, 'data', 'raw', 'msigdb')
        source_list = os.listdir(source_path)

        for item in source_list:
            #
            input_path = os.path.join(source_path, item)
            output_path = input_path.replace('raw', 'tmp').replace('xml', 'txt')

            #
            tree = ET.parse(input_path)
            root = tree.getroot()
            cols = ['STANDARD_NAME', 'SYSTEMATIC_NAME', 'ORGANISM', 'EXTERNAL_DETAILS_URL', 'CHIP',
                    'CATEGORY_CODE', 'SUB_CATEGORY_CODE', 'CONTRIBUTOR', 'CONTRIBUTOR_ORG',
                    'DESCRIPTION_BRIEF', 'MEMBERS_EZID', 'DESCRIPTION_FULL']
            convert = [{attrib: xe.attrib[attrib] for attrib in cols} for xe in root]
        
            df = pd.DataFrame(columns=cols)
            df = pd.concat([df, pd.DataFrame(convert)], ignore_index=True)[['SYSTEMATIC_NAME', 'DESCRIPTION_BRIEF', 'MEMBERS_EZID']]
            # df = df.append(convert, ignore_index=True)[['SYSTEMATIC_NAME', 'DESCRIPTION_BRIEF', 'MEMBERS_EZID']]

            # post process
            # 1. pathway id + pathway decription
            df_pathway_list = df[['SYSTEMATIC_NAME', 'DESCRIPTION_BRIEF']].rename(columns={"SYSTEMATIC_NAME": "pathway_id", "DESCRIPTION_BRIEF": "pathway_description"})

            print('pathway list done.')
            df_pathway_list.to_csv(output_path.replace('txt', 'pathway_list') + '.txt', index=False, encoding='utf-8')

            # 2. involvement list

            involve_list = []
            for index, row in df.iterrows():
                tmp_pathway_id = row[0]
                tmp_gene_set = row[2].split(',')
                for item in tmp_gene_set:
                    involve_list.append([item, tmp_pathway_id])
            print('involvement list done.')
            df_involve_list = pd.DataFrame(involve_list, columns=['gene_id', 'pathway_id'])

            df_involve_list.to_csv(output_path.replace('txt', 'involve_list') + '.txt', index=False, encoding='utf-8')
            print('Done.')
            print()
