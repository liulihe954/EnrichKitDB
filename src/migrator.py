import os
import pandas as pd
from datetime import datetime


class migrator:

    def __init__(self, BASE, SPECIES):
        self.BASE = BASE
        self.DB_UNITS = ['go', 'interpro', 'kegg', 'mesh', 'reactome', 'msigdb']
        self.SPECIES = SPECIES

    # read in everything

    def read_gtf(self):

        # condition on species to figure out the right table
        if self.SPECIES == 'bta':
            path_tmp = os.path.join(self.BASE, 'data', 'tmp', 'annotation', 'bos_taurus_gtf.txt')

        # read raw (keep a skinny version for id mapper composition)
        gtf_raw = pd.read_csv(path_tmp, low_memory=False)

        gtf_skinny = gtf_raw[gtf_raw.feature.isin(['gene'])][['gene_id', 'seqname', 'start', 'end', 'strand', 'gene_biotype']]

        # clean version: keep those one chr 1-29 and X
        gtf_clean = gtf_raw[gtf_raw.seqname.isin([
            '1', '2', '3', '4', '5',
            '6', '7', '8', '9', '10',
            '11', '12', '13', '14', '15',
            '16', '17', '18', '19', '20',
            '21', '22', '23', '24', '25',
            '26', '27', '28', '29', 'X', ])]

    #     # extract gene info
        gtf_clean_genes = gtf_clean[gtf_clean.feature.isin(['gene'])][['gene_id', 'seqname', 'start', 'end', 'strand', 'gene_biotype']]

    #     # extract exon info
        gtf_clean_exons = gtf_clean[gtf_clean.feature.isin(['exon'])][['gene_id', 'exon_id', 'seqname', 'start', 'end', 'exon_number', 'gene_biotype', 'strand', 'transcript_id']]

    #     # extract other feature info
        gtf_clean_features = gtf_clean[~gtf_clean.feature.isin(['gene', 'exon'])][['gene_id', 'feature', 'seqname', 'start', 'end', 'gene_biotype', 'strand']]

        ##
        self.gtf_skinny = gtf_skinny
        self.gtf_clean_genes = gtf_clean_genes
        self.gtf_clean_exons = gtf_clean_exons
        self.gtf_clean_features = gtf_clean_features

    def compose_id_mapper(self):

        # all source paths
        gene2ensembl_path_tmp = os.path.join(self.BASE, 'data', 'tmp', 'id_mapper', 'gene2ensembl.txt')
        vgnc_path_tmp = os.path.join(self.BASE, 'data', 'tmp', 'id_mapper', 'vgnc_protein_coding_gene.txt')
        hgnc_path_tmp = os.path.join(self.BASE, 'data', 'tmp', 'id_mapper', 'hgnc_protein_coding_gene.txt')

        if self.SPECIES == 'bta':

            # entrez id: differentiate the given SPECIES by tax id
            gene2ensembl_iter_csv = pd.read_csv(gene2ensembl_path_tmp, iterator=True, chunksize=500000)
            df_gene2ensembl_bta = pd.concat([chunk[chunk['tax_id'] == 9913] for chunk in gene2ensembl_iter_csv])

            #
            vgnc_iter_csv = pd.read_csv(vgnc_path_tmp, iterator=True, chunksize=10000)
            df_vgnc_bta = pd.concat([chunk[chunk['gene_id'].str.startswith('ENSBTAG', na=False)] for chunk in vgnc_iter_csv])

            #
        df_hgnc_bta = pd.read_csv(hgnc_path_tmp)

        # take union as the index
        tmp1 = pd.DataFrame(self.gtf_skinny['gene_id'].drop_duplicates())
        tmp1['source'] = 4

        tmp2 = pd.DataFrame(df_gene2ensembl_bta['gene_id'].drop_duplicates())
        tmp2['source'] = 2

        tmp3 = pd.DataFrame(df_vgnc_bta['gene_id'].drop_duplicates())
        tmp3['source'] = 1

        dfs = [tmp1, tmp2, tmp3]
        gene_id_idx = pd.concat(dfs)
        gene_id_out = pd.DataFrame(gene_id_idx.groupby(by=['gene_id'])['source'].sum())

        # merge
        id_mapper = pd.merge(gene_id_out, df_gene2ensembl_bta[['gene_id', 'entrez_id']], on='gene_id', how='left')
        id_mapper = pd.merge(id_mapper, df_vgnc_bta, on='gene_id', how='left')
        id_mapper = pd.merge(id_mapper, df_hgnc_bta, on='hgnc_orthologs', how='left')

        # post processing
        id_mapper['entrez_id'] = id_mapper['entrez_id'].astype('Int64')
        id_mapper['human_entrez_id'] = id_mapper['human_entrez_id'].astype('Int64')

        id_mapper.reset_index(drop=True, inplace=True)
        id_mapper['ek_gene_id'] = id_mapper.index
        first_column = id_mapper.pop('ek_gene_id')
        id_mapper.insert(0, 'ek_gene_id', first_column)

        self.id_mapper = id_mapper

    def replace_id(self, df,id_mapper_tiny_no_dup):
        output = pd.merge(df, id_mapper_tiny_no_dup, how='left', on='gene_id')
        first_column = output.pop('ek_gene_id')
        output.insert(0, 'ek_gene_id', first_column)
        df_final = output.drop(columns=['gene_id'])

        return df_final

    def process_db_unit(self):

        output = []
        tmp_path = os.path.join(self.BASE, 'data', 'tmp')
        for db_unit in self.DB_UNITS:
            if self.SPECIES == 'bta':
                if db_unit in ['go', 'interpro']:
                    path_tmp = os.path.join(tmp_path, db_unit, db_unit + '_bta.txt')
                    df_tmp = pd.read_csv(path_tmp, sep='\t', header=None, names=['gene_id', 'pathway_id', 'pathway_description'])
                    #
                    pathway_list_df = df_tmp[['pathway_id', 'pathway_description']].drop_duplicates().dropna()
                    involve_list_df = df_tmp[['gene_id', 'pathway_id']].dropna()

                    output.append([pathway_list_df, involve_list_df])

                if db_unit in ['kegg', 'mesh', 'msigdb']:

                    if db_unit == 'kegg':
                        pathway_list_path_tmp = os.path.join(tmp_path, db_unit, 'bta_pathway_description.txt')
                        involve_list_df_path_tmp = os.path.join(tmp_path, db_unit, 'bta_gene_sets.txt')

                    elif db_unit == 'mesh':
                        pathway_list_path_tmp = os.path.join(tmp_path, db_unit, 'MeSH.db.txt')
                        involve_list_df_path_tmp = os.path.join(tmp_path, db_unit, 'MeSH.Bta.eg.db.txt')

                    elif db_unit == 'msigdb':
                        pathway_list_path_tmp = os.path.join(tmp_path, db_unit, 'msigdb_v7.5.1.pathway_list.txt')
                        involve_list_df_path_tmp = os.path.join(tmp_path, db_unit, 'msigdb_v7.5.1.involve_list.txt')

                    pathway_list_df = pd.read_csv(pathway_list_path_tmp)
                    involve_list_df = pd.read_csv(involve_list_df_path_tmp)
                    involve_list_df['gene_id'] = involve_list_df['gene_id'].astype('Int64')

                    output.append([pathway_list_df, involve_list_df])

                if db_unit in ['reactome']:
                    path_tmp = os.path.join(tmp_path, db_unit, 'NCBI2Reactome_All_Levels.txt')
                    iter_csv = pd.read_csv(path_tmp, iterator=True, chunksize=50000)
                    df_tmp = pd.concat([chunk[chunk['species'] == 'Bos taurus'] for chunk in iter_csv])

                    pathway_list_df = df_tmp[['pathway_id', 'pathway_description']].drop_duplicates().dropna()
                    involve_list_df = df_tmp[['gene_id', 'pathway_id']].dropna()

                    output.append([pathway_list_df, involve_list_df])

        self.twolist_before_convert = output

    def create_tables(self):
        # meta table
        current_time = datetime.now().strftime("%D")
        pathway_meta = [
            [1, 'bta', 'biomart_go', 'gene_id', current_time],
            [2, 'bta', 'biomart_interpro', 'gene_id', current_time],
            [3, 'bta', 'kegg', 'entrez_id', current_time],
            [4, 'bta', 'mesh', 'entrez_id', current_time],
            [5, 'bta', 'reactome', 'entrez_id', current_time],
            [6, 'bta', 'msigdb', 'human_entrez_id', current_time],
        ]
        self.pathway_meta_df = pd.DataFrame(pathway_meta, columns=['pathway_meta_id', 'species', 'name', 'id_type', 'update_time'])

        # all pathways
        pathway_collection = []

        for index, pathway_section in enumerate(self.twolist_before_convert):

            pathway_tmp = pathway_section[0]
            pathway_tmp['pathway_meta_id'] = index + 1
            pathway_collection.append(pathway_tmp)

        pathway_collection_df = pd.concat(pathway_collection)
        pathway_collection_df.reset_index(drop=True, inplace=True)
        pathway_collection_df['ek_pathway_id'] = pathway_collection_df.index
        first_column = pathway_collection_df.pop('ek_pathway_id')
        pathway_collection_df.insert(0, 'ek_pathway_id', first_column)

        self.all_pathway_after_convert = pathway_collection_df

        ####
        converted_involve_collection = []

        for index, pathway_section in enumerate(self.twolist_before_convert):

            tmp = pathway_section[1].drop_duplicates()

            if index in [0, 1]:  # 0 1 gene id
                tmp_out = pd.merge(tmp, self.id_mapper[['ek_gene_id', 'gene_id']].drop_duplicates('gene_id'), how='left', on='gene_id')
                tmp_out2 = pd.merge(tmp_out, pathway_collection_df[['pathway_id', 'ek_pathway_id']].drop_duplicates('pathway_id'), how='left', on='pathway_id')
                tmp_out2['ek_pathway_id'] = tmp_out2['ek_pathway_id'].astype('Int64')
                tmp_final = tmp_out2[['ek_gene_id', 'ek_pathway_id']]
                converted_involve_collection.append(tmp_final)
                pass
            elif index in [2, 3, 4]:  # 2 3 4 entrez id
                tmp_mapper = self.id_mapper[['ek_gene_id', 'entrez_id']].drop_duplicates('entrez_id').dropna()
                tmp_out = pd.merge(tmp, tmp_mapper, how='left', left_on='gene_id', right_on='entrez_id')
                tmp_out2 = pd.merge(tmp_out, pathway_collection_df[['pathway_id', 'ek_pathway_id']].drop_duplicates('pathway_id'), how='left', on='pathway_id')
                tmp_out2['ek_gene_id'] = tmp_out2['ek_gene_id'].astype('Int64')
                tmp_out2['ek_pathway_id'] = tmp_out2['ek_pathway_id'].astype('Int64')
                tmp_out2.dropna()
                tmp_final = tmp_out2[['ek_gene_id', 'ek_pathway_id']]
                converted_involve_collection.append(tmp_final)
            elif index in [5]:  # 5 human_entrez_id
                tmp_mapper = self.id_mapper[['ek_gene_id', 'human_entrez_id']].drop_duplicates('human_entrez_id').dropna()
                tmp_out = pd.merge(tmp, tmp_mapper, how='left', left_on='gene_id', right_on='human_entrez_id')
                tmp_out2 = pd.merge(tmp_out, pathway_collection_df[['pathway_id', 'ek_pathway_id']].drop_duplicates('pathway_id'), how='left', on='pathway_id')
                tmp_out2['ek_gene_id'] = tmp_out2['ek_gene_id'].astype('Int64')
                tmp_out2['ek_pathway_id'] = tmp_out2['ek_pathway_id'].astype('Int64')
                tmp_out2.dropna()
                tmp_final = tmp_out2[['ek_gene_id', 'ek_pathway_id']]
                converted_involve_collection.append(tmp_final)

        involve = pd.concat(converted_involve_collection)
        involve = involve.dropna()
        involve.reset_index(drop=True, inplace=True)

        #
        self.all_involve_after_convert = involve

    def migrate_tables(self):
        print('Now to format gtf files.....')
        self.read_gtf()
        print('ok...')

        print('Now to formatid mapper files.....')
        self.compose_id_mapper()
        print('ok...')

        print('Now to format db units files.....')
        self.process_db_unit()
        print('ok...')

        print('Now to create files.....')
        self.create_tables()
        print('ok...')

        #
        id_mapper_tiny = self.id_mapper[['ek_gene_id', 'gene_id']]
        id_mapper_tiny_no_dup = id_mapper_tiny.drop_duplicates('gene_id')
        self.gtf_clean_genes = self.replace_id(self.gtf_clean_genes, id_mapper_tiny_no_dup)
        self.gtf_clean_exons = self.replace_id(self.gtf_clean_exons, id_mapper_tiny_no_dup)
        self.gtf_clean_features = self.replace_id(self.gtf_clean_features, id_mapper_tiny_no_dup)
        #
        print('Now to output all tables.....')
        if self.SPECIES == 'bta':
            #
            path_tmp = os.path.join(self.BASE, 'data', 'output')
            
            self.gtf_clean_exons['exon_number'] = self.gtf_clean_exons['exon_number'].astype('Int64')
            self.gtf_clean_exons = self.gtf_clean_exons.drop(columns = ['seqname'])
            self.gtf_clean_features = self.gtf_clean_features.drop(columns = ['seqname'])
            self.gtf_clean_genes['seqname'] = self.gtf_clean_genes['seqname'].replace(['X'],30)
            
            self.gtf_clean_genes['seqname'] = pd.to_numeric(self.gtf_clean_genes['seqname'], errors='coerce')
            #self.gtf_clean_genes['seqname'].astype('Int64')
            
            #
            print(f'About to output gene list with dimension {self.gtf_clean_genes.shape}')
            self.gtf_clean_genes.to_csv(os.path.join(path_tmp, 'genes_bta.csv'), index=False, encoding='utf-8')
            
            #
            print(f'About to output exon list with dimension {self.gtf_clean_exons.shape}')
            self.gtf_clean_exons.to_csv(os.path.join(path_tmp, 'exons_bta.csv'), index=False, encoding='utf-8')
            
            #
            print(f'About to output feature list with dimension {self.gtf_clean_features.shape}')
            self.gtf_clean_features.to_csv(os.path.join(path_tmp, 'features_bta.csv'), index=False, encoding='utf-8')
            
            
            #
            print(f'About to output pathway db meta list with dimension {self.pathway_meta_df.shape}')
            self.pathway_meta_df.to_csv(os.path.join(path_tmp, 'pathway_meta_bta.csv'), index=False, encoding='utf-8')
            
            print(f'About to output id mapper list with dimension {self.id_mapper.shape}')
            self.id_mapper['entrez_id'] = self.id_mapper['entrez_id'].astype('Int64')
            self.id_mapper['human_entrez_id'] = self.id_mapper['human_entrez_id'].astype('Int64')
            self.id_mapper.to_csv(os.path.join(path_tmp, 'id_mapper_bta.csv'), index=False, encoding='utf-8')
            
            
            print(f'About to output all pathway list with dimension {self.all_pathway_after_convert.shape}')
            self.all_pathway_after_convert.to_csv(os.path.join(path_tmp, 'pathways_bta.csv'), index=False, encoding='utf-8')
            
            print(f'About to output all gene-pathway involve list with dimension {self.all_involve_after_convert.shape}')
            self.all_involve_after_convert.to_csv(os.path.join(path_tmp, 'involve_bta.csv'), index=False, encoding='utf-8')

        print('Done!')
