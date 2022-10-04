import os
import pandas as pd
from datetime import datetime
from collections import defaultdict


class migrator:

    def __init__(self, BASE, SPECIES):
        self.BASE = BASE
        self.DB_UNITS = ['go', 'interpro', 'kegg', 'mesh', 'reactome', 'msigdb']
        self.SPECIES_DICT = {
            'bta': ['9913', 'ENSBTAG', 'Bta', 'Bos taurus'],
            'cap': ['9925', 'ENSCHIG', None, None],
            'equ': ['9796', 'ENSECAG', 'Eca', None],
            'gal': ['9031', 'ENSGALG', 'Gga', 'Gallus gallus'],
            'ovi': ['9940', 'ENSOARG', None, None],
            'sus': ['9823', 'ENSSSCG', 'Ssc', 'Sus scrofa']
        }
        self.SPECIES = SPECIES

    # caller
    def migrate_tables(self):
        # call manipulate_tables
        self.manipulate_tables()

    # handler
    def manipulate_tables(self):
        print('Now to format gtf files.....')
        self.read_gtf()
        print('ok...')

#         print('Now to format computated features.....')
#         self.compute_feature()
#         print('ok...')

        print('Now to format id mapper files.....')
        self.compose_id_mapper()
        print('ok...')

        print('Now to format db units files.....')
        self.process_db_unit()
        print('ok...')

        print('Now to format gene-pathway involvement files.....')
        self.create_involve_tables()
        print('ok...')

        #
        id_mapper_tiny = self.id_mapper[['ek_gene_id', 'gene_id']]
        id_mapper_tiny_no_dup = id_mapper_tiny.drop_duplicates('gene_id')
        self.gtf_clean_genes = self.replace_id(self.gtf_clean_genes, id_mapper_tiny_no_dup)
        self.gtf_clean_exons = self.replace_id(self.gtf_clean_exons, id_mapper_tiny_no_dup)
        self.gtf_clean_features = self.replace_id(self.gtf_clean_features, id_mapper_tiny_no_dup)

#       #
        self.flatten_exon()

        print('Now to format computated features.....')
        self.compute_feature()
        print('ok...')

        #
        print('Now to output all tables.....')
        if self.SPECIES in self.SPECIES_DICT:
            #
            path_tmp = os.path.join(self.BASE, 'data', 'output')

            self.gtf_clean_exons['exon_number'] = self.gtf_clean_exons['exon_number'].astype('Int64')
            self.gtf_clean_exons = self.gtf_clean_exons.drop(columns=['seqname'])
            self.gtf_clean_features = self.gtf_clean_features.drop(columns=['seqname'])

            self.gtf_clean_genes['seqname'] = self.gtf_clean_genes['seqname'].replace(['X'], 30)

            self.gtf_clean_genes['seqname'] = pd.to_numeric(self.gtf_clean_genes['seqname'], errors='coerce')
            # self.gtf_clean_genes['seqname'].astype('Int64')

            #
            print(f'About to output gene list with dimension {self.gtf_clean_genes.shape}')
            self.gtf_clean_genes.to_csv(os.path.join(path_tmp, ''.join(['genes_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            #
            print(f'About to output exon list with dimension {self.gtf_clean_exons.shape}')
            self.gtf_clean_exons_falt.to_csv(os.path.join(path_tmp, ''.join(['exons_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            #
            print(f'About to output feature list with dimension {self.gtf_clean_features.shape}')
            self.gtf_clean_features.to_csv(os.path.join(path_tmp, ''.join(['features_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            #
            print(f'About to output pathway db meta list with dimension {self.pathway_meta_df.shape}')
            self.pathway_meta_df.to_csv(os.path.join(path_tmp, ''.join(['pathway_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            print(f'About to output id mapper list with dimension {self.id_mapper.shape}')
            self.id_mapper['entrez_id'] = self.id_mapper['entrez_id'].astype('Int64')
            self.id_mapper['human_entrez_id'] = self.id_mapper['human_entrez_id'].astype('Int64')
            self.id_mapper.to_csv(os.path.join(path_tmp, ''.join(['id_mapper_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            print(f'About to output all pathway list with dimension {self.all_pathway_after_convert.shape}')
            self.all_pathway_after_convert.to_csv(os.path.join(path_tmp, ''.join(['pathways_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            print(f'About to output all gene-pathway involve list with dimension {self.all_involve_after_convert.shape}')
            self.all_involve_after_convert.to_csv(os.path.join(path_tmp, ''.join(['involve_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            print(f'About to output all gene-pathway involve list with dimension {self.all_involve_after_convert.shape}')
            self.computatd_features.to_csv(os.path.join(path_tmp, ''.join(['computed_feature_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

        print('Done!')
        print('')

    ## read in everything ##

    def read_gtf(self):
        # read in gtfs
        if self.SPECIES in self.SPECIES_DICT:
            path_tmp = os.path.join(self.BASE, 'data', 'tmp', 'annotation', ''.join([self.SPECIES, '_gtf.txt']))

        # read raw (keep a skinny version for id mapper composition)
        gtf_raw = pd.read_csv(path_tmp, low_memory=False)

        gtf_skinny = gtf_raw[gtf_raw.feature.isin(['gene'])][['gene_id', 'seqname', 'start', 'end', 'strand', 'gene_biotype']]

        # clean version: remove all non-chr items
        gtf_clean = gtf_raw[gtf_raw.seqname.isin([
            '1', '2', '3', '4', '5',
            '6', '7', '8', '9', '10',
            '11', '12', '13', '14', '15',
            '16', '17', '18', '19', '20',
            '21', '22', '23', '24', '25',
            '26', '27', '28', '29', '30',
            '31', '32', '33', '34', '35',
            '36', '37', '38', '39', '40',
            'X', 'Y', 'W', 'Z'])]

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

        if self.SPECIES in self.SPECIES_DICT:

            # entrez id: differentiate the given SPECIES by tax id
            gene2ensembl_iter_csv = pd.read_csv(gene2ensembl_path_tmp, iterator=True, chunksize=500000)
            df_gene2ensembl_tmp = pd.concat([chunk[chunk['tax_id'] == self.SPECIES_DICT[self.SPECIES][0]] for chunk in gene2ensembl_iter_csv])

            #
            vgnc_iter_csv = pd.read_csv(vgnc_path_tmp, iterator=True, chunksize=10000)
            df_vgnc_tmp = pd.concat([chunk[chunk['gene_id'].str.startswith(self.SPECIES_DICT[self.SPECIES][1], na=False)] for chunk in vgnc_iter_csv])

            #
        df_hgnc_tmp = pd.read_csv(hgnc_path_tmp)

        # take union as the index
        tmp1 = pd.DataFrame(self.gtf_skinny['gene_id'].drop_duplicates())
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

        # post processing
        id_mapper['entrez_id'] = id_mapper['entrez_id'].astype('Int64')
        id_mapper['human_entrez_id'] = id_mapper['human_entrez_id'].astype('Int64')

        # id_mapper.reset_index(drop=True, inplace=True)
        # id_mapper['ek_gene_id'] = id_mapper.index
        # first_column = id_mapper.pop('ek_gene_id')
        # id_mapper.insert(0, 'ek_gene_id', first_column)

        self.id_mapper = id_mapper

    def replace_id(self, df, id_mapper_tiny_no_dup):
        output = pd.merge(df, id_mapper_tiny_no_dup, how='left', on='gene_id')
        first_column = output.pop('ek_gene_id')
        output.insert(0, 'ek_gene_id', first_column)
        df_final = output.drop(columns=['gene_id'])

        return df_final

    # helper to dedup exons
    def flatten_exon(self):
        exon_df = self.gtf_clean_exons.copy(deep=True)
        exon_df = exon_df.groupby('ek_gene_id')
        # logic here
        print('About to flatten exons list...')
        container = []
        for name, group in exon_df:
            if name % 1000 == 0:
                print(f'flattend {name} genes...')
            tmp_df = group
            tmp_merger = defaultdict(list)
            for row in tmp_df.iterrows():
                tmp_transcript_exon = row[1][1]
                tmp_transcript_collection = tmp_merger[tmp_transcript_exon]
                tmp_transcript_item = row[1][8]
                if tmp_transcript_item not in tmp_transcript_collection:
                    tmp_merger[tmp_transcript_exon].append(tmp_transcript_item)
            merger_df = pd.DataFrame(list([k, '+'.join(v)] for k, v in tmp_merger.items()), columns=['exon_id', 'transcript_id'])
            tmp_df = tmp_df.drop(['seqname', 'exon_number', 'transcript_id'], axis=1)
            tmp_output = pd.merge(tmp_df, merger_df, how='left', on='exon_id').drop_duplicates().reset_index(drop=True)  # reset_index(drop=True,)

            container.append(tmp_output)

        print('Done.')

        self.gtf_clean_exons_falt = pd.concat(container)

    def process_db_unit(self):

        output = []
        tmp_path = os.path.join(self.BASE, 'data', 'tmp')
        for db_unit in self.DB_UNITS:
            if self.SPECIES in self.SPECIES_DICT.keys():
                if db_unit in ['go', 'interpro']:
                    path_tmp = os.path.join(tmp_path, db_unit, db_unit + ''.join(['_', self.SPECIES, '.txt']))
                    df_tmp = pd.read_csv(path_tmp, sep='\t', header=None, names=['gene_id', 'pathway_id', 'pathway_description']).drop_duplicates().dropna()
                    #
                    pathway_list_df = df_tmp[['pathway_id', 'pathway_description']].drop_duplicates().dropna()
                    involve_list_df = df_tmp[['gene_id', 'pathway_id']].drop_duplicates().dropna()

                    output.append([pathway_list_df, involve_list_df])

                if db_unit in ['kegg', 'mesh', 'msigdb']:

                    if db_unit == 'kegg':
                        pathway_list_path_tmp = os.path.join(tmp_path, db_unit, ''.join([self.SPECIES, '_pathway_description.txt']))
                        involve_list_df_path_tmp = os.path.join(tmp_path, db_unit, ''.join([self.SPECIES, '_gene_sets.txt']))

                    elif db_unit == 'mesh':
                        if self.SPECIES not in ['ovi', 'cap']:
                            pathway_list_path_tmp = os.path.join(tmp_path, db_unit, 'MeSH.db.txt')
                            involve_list_df_path_tmp = os.path.join(tmp_path, db_unit, ''.join(['MeSH.', self.SPECIES_DICT[self.SPECIES][2], '.eg.db.txt']))
                        else:
                            pathway_list_path_tmp = None
                            involve_list_df_path_tmp = None

                    elif db_unit == 'msigdb':
                        pathway_list_path_tmp = os.path.join(tmp_path, db_unit, 'msigdb_v7.5.1.pathway_list.txt')
                        involve_list_df_path_tmp = os.path.join(tmp_path, db_unit, 'msigdb_v7.5.1.involve_list.txt')

                    if pathway_list_path_tmp is not None:
                        pathway_list_df = pd.read_csv(pathway_list_path_tmp).drop_duplicates().dropna()
                        involve_list_df = pd.read_csv(involve_list_df_path_tmp).drop_duplicates().dropna()

                        involve_list_df['keep'] = map(lambda x: x.isdigit(), involve_list_df['gene_id'])
                        involve_list_df = involve_list_df[involve_list_df['keep'] == True]
                        involve_list_df.gene_id = involve_list_df.gene_id.astype('float').astype('Int64')

                        output.append([pathway_list_df, involve_list_df])

                if db_unit in ['reactome']:
                    path_tmp = os.path.join(tmp_path, db_unit, 'NCBI2Reactome_All_Levels.txt')
                    iter_csv = pd.read_csv(path_tmp, iterator=True, chunksize=50000)

                    df_tmp = pd.concat([chunk[chunk['species'] == self.SPECIES_DICT[self.SPECIES][3]] for chunk in iter_csv])

                    pathway_list_df = df_tmp[['pathway_id', 'pathway_description']].drop_duplicates().dropna()
                    involve_list_df = df_tmp[['gene_id', 'pathway_id']].drop_duplicates().dropna()

                    output.append([pathway_list_df, involve_list_df])

        self.twolist_before_convert = output

    def create_involve_tables(self):

        # meta table
        current_time = datetime.now().strftime("%D")
        pathway_meta = [
            # Bos taurus
            [1, 'bta', 'biomart_go', 'gene_id', current_time],
            [2, 'bta', 'biomart_interpro', 'gene_id', current_time],
            [3, 'bta', 'kegg', 'entrez_id', current_time],
            [4, 'bta', 'mesh', 'entrez_id', current_time],
            [5, 'bta', 'reactome', 'entrez_id', current_time],
            [6, 'bta', 'msigdb', 'human_entrez_id', current_time],
            # Ovis aries
            [7, 'ovi', 'biomart_go', 'gene_id', current_time],
            [8, 'ovi', 'biomart_interpro', 'gene_id', current_time],
            [9, 'ovi', 'kegg', 'entrez_id', current_time],
            [10, 'ovi', 'msigdb', 'human_entrez_id', current_time],
            # Capra hircus
            [11, 'cap', 'biomart_go', 'gene_id', current_time],
            [12, 'cap', 'biomart_interpro', 'gene_id', current_time],
            [13, 'cap', 'kegg', 'entrez_id', current_time],
            [14, 'cap', 'msigdb', 'human_entrez_id', current_time],
            # Equus caballus
            [15, 'equ', 'biomart_go', 'gene_id', current_time],
            [16, 'equ', 'biomart_interpro', 'gene_id', current_time],
            [17, 'equ', 'kegg', 'entrez_id', current_time],
            [18, 'equ', 'mesh', 'entrez_id', current_time],
            [19, 'equ', 'msigdb', 'human_entrez_id', current_time],

            # Sus scrofa
            [20, 'sus', 'biomart_go', 'gene_id', current_time],
            [21, 'sus', 'biomart_interpro', 'gene_id', current_time],
            [22, 'sus', 'kegg', 'entrez_id', current_time],
            [23, 'sus', 'mesh', 'entrez_id', current_time],
            [24, 'sus', 'reactome', 'entrez_id', current_time],
            [25, 'sus', 'msigdb', 'human_entrez_id', current_time],
            # Gallus gallus
            [26, 'gal', 'biomart_go', 'gene_id', current_time],
            [27, 'gal', 'biomart_interpro', 'gene_id', current_time],
            [28, 'gal', 'kegg', 'entrez_id', current_time],
            [29, 'gal', 'mesh', 'entrez_id', current_time],
            [30, 'gal', 'reactome', 'entrez_id', current_time],
            [31, 'gal', 'msigdb', 'human_entrez_id', current_time]
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

    def compute_feature(self):
        #
        exons_raw = self.gtf_clean_exons_falt.copy(deep=True)
        exons_raw = exons_raw.groupby('ek_gene_id')
        #
        upstream_len = downstream_len = 10000
        donor_len = acceptor_len = 50
        all_computed_feature = []
        #
        for name, group in exons_raw:
            #
            tmp_df = group.sort_values(by=['start', 'end'])
            #
            if name % 1000 == 1:
                print(f'annotated {name} genes...')

            # for each gene
            length = len(tmp_df.index)
            row = 0
            next_row = 0
            tmp_ek_gene_id = tmp_df.iloc[row][0]
            tmp_gene_type = tmp_df.iloc[row][4]
            tmp_gene_strand = tmp_df.iloc[row][5]
            #
            tmp_compute_feature = []

            tmp_compute_feature.append(['upstream',
                                        tmp_df.iloc[row][2] - 1 - upstream_len,
                                        tmp_df.iloc[row][2] - 1,
                                        tmp_gene_type,
                                        tmp_gene_strand])
            #
            while next_row < length - 1:
                tmp_intron_start = tmp_df.iloc[row][3] + 1
                while tmp_df.iloc[next_row][2] - 1 <= tmp_intron_start and next_row < length - 1:
                    next_row += 1
                tmp_intron_end = tmp_df.iloc[next_row][2] - 1
                row = next_row

                # logic for donor and splice
                if tmp_intron_end - tmp_intron_start <= 2 * donor_len:
                    tmp_compute_feature.append([
                        'short intron',
                        tmp_intron_start,
                        tmp_intron_end,
                        tmp_gene_type,
                        tmp_gene_strand])
                else:
                    # donor
                    tmp_compute_feature.append(
                        ['splice donor',
                         tmp_intron_start,
                         tmp_intron_start + donor_len,
                         tmp_gene_type,
                         tmp_gene_strand])

                    # intron
                    tmp_compute_feature.append(
                        ['intron',
                         tmp_intron_start + donor_len + 1,
                         tmp_intron_end - donor_len - 1,
                         tmp_gene_type,
                         tmp_gene_strand])

                    # acceptor
                    tmp_compute_feature.append(
                        ['splice acceptor',
                         tmp_intron_end - donor_len,
                         tmp_intron_end,
                         tmp_gene_type,
                         tmp_gene_strand])
            tmp_compute_feature.append(
                ['downstream',
                 tmp_df.iloc[length - 1][3] + 1,
                 tmp_df.iloc[length - 1][3] + 1 + upstream_len,
                 tmp_gene_type,
                 tmp_gene_strand])

            tmp_out = pd.DataFrame(tmp_compute_feature, columns=['feature', 'start', 'end', 'gene_biotype', 'strand'])
            tmp_out['ek_gene_id'] = tmp_ek_gene_id
            first_col = tmp_out.pop('ek_gene_id')
            tmp_out.insert(0, "ek_gene_id", first_col)
            all_computed_feature.append(tmp_out)

        # output
        self.computatd_features = pd.concat(all_computed_feature)
