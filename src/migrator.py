import os
import pandas as pd
import numpy as np
from datetime import datetime
import glob2
from src.update_id_mapper_row import update_row
import concurrent.futures

class migrator:
    def __init__(self, BASE, max_threads):
        self.BASE = BASE
        self.max_threads = max_threads - 1
        self.DB_UNITS = ['msigdb', 'go', 'interpro', 'kegg', 'mesh', 'reactome']
        self.SPECIES_DICT = {
            'bta': [0, 'Bos_taurus', 'ARS-UCD1.2.107'],
            'cap': [1, 'Capra_hircus', 'ARS1.107'],
            'equ': [2, 'Equus_caballus', 'EquCab3.0.107'],
            'gal': [3, 'Gallus_gallus', 'gca000002315v5.GRCg6a.107'],
            'ovi': [4, 'Ovis_aries', 'Oar_v3.1.107'],
            'sus': [5, 'Sus_scrofa', 'Sscrofa11.1.107'], }

    def generate_tables(self):
        self.do_migrate()

        output_path = os.path.join(self.BASE, 'data', 'output')
        print('Now writing species.txt')
        self.all_species_table.to_csv(os.path.join(output_path, 'species.txt'), index=False)

        print('Now writing id_mapper.txt')
        self.id_mapper_all_updated.to_csv(os.path.join(output_path, 'id_mapper.txt'), index=False)

        print('Now writing genes.txt')
        self.genes_all.to_csv(os.path.join(output_path, 'genes.txt'), index=False)

        print('Now writing exons.txt')
        self.exons_all.to_csv(os.path.join(output_path, 'exons.txt'), index=False)

        print('Now writing features.txt')
        self.features_all.to_csv(os.path.join(output_path, 'features.txt'), index=False)

        print('Now writing computed_features.txt')
        self.computed_feature_all.to_csv(os.path.join(output_path, 'computed_features.txt'), index=False)

        print('Now writing pathway_meta.txt')
        self.pathway_meta.to_csv(os.path.join(output_path, 'pathway_meta.txt'), index=False)

        print('Now writing pathway.txt')
        self.all_pathway.to_csv(os.path.join(output_path, 'pathway.txt'), index=False)

        print('Now writing involve.txt')
        self.all_involve.to_csv(os.path.join(output_path, 'involve.txt'), index=False)
        self.all_involve_msigdb.to_csv(os.path.join(output_path, 'involve_msigdb.txt'), index=False)

    def do_migrate(self):
        print('About to compose species table.')
        self.compose_species_table()
        print('Done.\n')

        print('About to compose id_mapper table.')
        self.compose_id_mapper()
        print('Done.\n')

        print('About to compose genes table.')
        self.compose_genes()
        print('Done.\n')

        print('About to compose exons table.')
        self.compose_exons()
        print('Done.\n')

        print('About to compose features table.')
        self.compose_features()
        print('Done.\n')

        print('About to compose features table.')
        self.compose_compute_features()
        print('Done.\n')

        print('About to compose pathway meta table.')
        self.compose_pathway_meta()
        print('Done.\n')

        print('About to compose pathway table.')
        self.compose_pathway()
        print('Done.\n')

        print('About to compose involve table.')
        self.compose_involve()
        print('Done.\n')

    def compose_species_table(self):
        all_species_table = pd.DataFrame([
            [0, 'bta', 'Bos_taurus', 'ARS-UCD1.2.107'],
            [1, 'cap', 'Capra_hircus', 'ARS1.107'],
            [2, 'equ', 'Equus_caballus', 'EquCab3.0.107'],
            [3, 'gal', 'Gallus_gallus', 'gca000002315v5.GRCg6a.107'],
            [4, 'ovi', 'Ovis_aries', 'Oar_v3.1.107'],
            [5, 'sus', 'Sus_scrofa', 'Sscrofa11.1.107'],
        ])
        all_species_table.columns = ['ek_species', 'name_short', 'name_long', 'current_gtf']

        # mark
        self.all_species_table = all_species_table

    def process_row(self, row):
        # Wrapper for processing a single row
        # This can include print statements or additional logic
        print(f'Processing row with gene_id: {row.get("gene_id", "N/A")}')
        return update_row(row)
    
    def manual_check_id(self):

        # define paths
        input_path = os.path.join(self.BASE, 'data', 'tmp', 'id_mapper')
        files_to_concatenate = glob2.glob(os.path.join(input_path, 'id_mapper_*.txt'))
        
        # compose a full id_mapper_pre_check
        all_id_mapper_df = []
        for file in files_to_concatenate:
            df = pd.read_csv(file, engine='python')
            all_id_mapper_df.append(df)
        if not all_id_mapper_df:
            print("No files to process.")
        
        #
        self.id_mapper_all_raw = pd.concat(all_id_mapper_df, ignore_index=True)
        self.id_mapper_all_raw.astype(str)

        # do manual check then output a "id_mapper_all_updated.txt" (In-place Check)
        
        # Convert DataFrame to list of dictionaries for parallel processing
        rows = self.id_mapper_all_raw.to_dict(orient='records')
        # threads pool
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_threads) as executor:
            # Submit all rows for processing
            results = executor.map(self.process_row, rows)

        # Updating the DataFrame with the results
        for i, updated_row in enumerate(results):
            self.id_mapper_all_raw.iloc[i] = updated_row
        
        # write to csv in case other issue happen
        self.id_mapper_all_updated.to_csv(os.path.join(self.BASE,'data', 'tmp', 'id_mapper', 'id_mapper_all_updated.txt'), index=False)
        # for index, row in self.id_mapper_all_raw.iterrows():
        #     print('manual_check_id for ' + str(index))
        #     # Perform the checks and updates as per your original logic here
        #     updated_row = update_row(row)

        #     # If 'gene_id' is not missing, do nothing
        #     if not pd.isna(row['gene_id']):
        #         pass

        #     # Update the row in the DataFrame
        #     self.id_mapper_all_raw.loc[index] = updated_row
        # Process rows in parallel using ThreadPoolExecutor

    def compose_id_mapper(self):
        # input_path = os.path.join(self.BASE, 'data', 'tmp', 'id_mapper')

        ## check id to generate id_mapper_all_updated.txt
        self.manual_check_id()

        ##
        id_mapper_all_updated = pd.read_csv(os.path.join(self.BASE,'data', 'tmp', 'id_mapper', 'id_mapper_all_updated.txt'), low_memory=False)
        id_mapper_all_updated.fillna(-1, inplace=True)
        id_mapper_all_updated['entrez_id'] = id_mapper_all_updated['entrez_id'].astype(int)
        id_mapper_all_updated['human_entrez_id'].replace('None', -1, inplace=True)
        id_mapper_all_updated['human_entrez_id'] = id_mapper_all_updated['human_entrez_id'].astype(float)
        id_mapper_all_updated['human_entrez_id'] = id_mapper_all_updated['human_entrez_id'].astype(int)
        # id_mapper_all_updated.pop('Unnamed: 0')
        id_mapper_all_updated = id_mapper_all_updated.sort_values(by=['species', 'gene_id'])
        id_mapper_all_updated.reset_index(drop=True, inplace=True)
        id_mapper_all_updated['ek_gene_id'] = id_mapper_all_updated.index

        # id_mapper_all_updated = id_mapper_all_updated.astype({'entrez_id': float})
        # id_mapper_all_updated = id_mapper_all_updated.astype({'entrez_id': int})
        # id_mapper_all_updated = id_mapper_all_updated.astype({'entrez_id': str})

        # id_mapper_all_updated = id_mapper_all_updated.astype({'human_entrez_id': float})
        # id_mapper_all_updated = id_mapper_all_updated.astype({'human_entrez_id': int})
        # id_mapper_all_updated = id_mapper_all_updated.astype({'human_entrez_id': str})
        id_mapper_all_updated.astype(str)
        id_mapper_all_updated.replace(-1, '', inplace=True)

        id_mapper_all_updated = id_mapper_all_updated[['ek_gene_id', 'gene_id', 'ensembl_symbol', 'entrez_id', 'ncbi_symbol', 'vgnc_id', 'vgnc_symbol', 'hgnc_orthologs', 'human_gene_id', 'human_entrez_id', 'hgnc_symbol', 'species']]

        print('Done reading id_mapper_all_updated.txt.')

        id_mapper_tiny_no_dup = id_mapper_all_updated[['ek_gene_id', 'gene_id']].drop_duplicates().dropna()

        # mark
        self.id_mapper_all_updated = id_mapper_all_updated
        self.id_mapper_tiny_no_dup = id_mapper_tiny_no_dup

    def compose_genes(self):
        input_path = os.path.join(self.BASE, 'data', 'tmp', 'annotation')
        all_out = []
        for item in os.listdir(input_path):
            if item.startswith('genes_'):
                print('Reading genes information for species: {}...'.format(item))
                cur_name_short = item.split('_')[-1].split('.')[0]
                tmp_df = pd.read_csv(os.path.join(input_path, item))
                tmp_df.insert(0, "species", self.SPECIES_DICT[cur_name_short][0])
                all_out.append(tmp_df)
        all_out_df = pd.concat(all_out)
        all_out_df = all_out_df.sort_values(by=['species', 'gene_id'])

        all_out_df_tagged = pd.merge(all_out_df, self.id_mapper_tiny_no_dup, how='left', on='gene_id')
        first_column = all_out_df_tagged.pop('ek_gene_id')
        all_out_df_tagged.insert(0, 'ek_gene_id', first_column)
        genes_all = all_out_df_tagged.drop(columns=['gene_id'])

        genes_all = genes_all.astype({'ek_gene_id': float})
        genes_all = genes_all.dropna().astype({'ek_gene_id': int})
        genes_all = genes_all.astype({'ek_gene_id': str})

        genes_all = genes_all[['seqname', 'start', 'end', 'strand', 'gene_biotype', 'ek_gene_id', 'species']]

        # mark
        self.genes_all = genes_all

    def compose_exons(self):

        input_path = os.path.join(self.BASE, 'data', 'tmp', 'annotation')
        all_out = []
        for item in os.listdir(input_path):
            if item.startswith('exons_'):
                print('Reading exon information for species: {}...'.format(item))
                tmp_df = pd.read_csv(os.path.join(input_path, item))
                all_out.append(tmp_df)
        all_out_df = pd.concat(all_out)
        all_out_df = all_out_df.sort_values(by=['gene_id'])

        all_out_df_tagged = pd.merge(all_out_df, self.id_mapper_tiny_no_dup, how='left', on='gene_id')
        first_column = all_out_df_tagged.pop('ek_gene_id')
        all_out_df_tagged.insert(0, 'ek_gene_id', first_column)
        exons_all = all_out_df_tagged.drop(columns=['gene_id'])
        exons_all = exons_all.sort_values(by=['ek_gene_id'])

        exons_all = exons_all.astype({'ek_gene_id': float})
        exons_all = exons_all.dropna().astype({'ek_gene_id': int})
        exons_all = exons_all.astype({'ek_gene_id': str})

        exons_all = exons_all[['exon_id', 'start', 'end', 'gene_biotype', 'strand', 'transcript_id', 'ek_gene_id']]

        # mark
        self.exons_all = exons_all

    def compose_features(self):
        input_path = os.path.join(self.BASE, 'data', 'tmp', 'annotation')
        all_out = []
        for item in os.listdir(input_path):
            if item.startswith('features_'):
                print('Reading feature information for species: {}...'.format(item))
                tmp_df = pd.read_csv(os.path.join(input_path, item))
                all_out.append(tmp_df)
        all_out_df = pd.concat(all_out)
        all_out_df = all_out_df.sort_values(by=['gene_id'])

        all_out_df_tagged = pd.merge(all_out_df, self.id_mapper_tiny_no_dup, how='left', on='gene_id')
        first_column = all_out_df_tagged.pop('ek_gene_id')
        all_out_df_tagged.insert(0, 'ek_gene_id', first_column)
        features_all = all_out_df_tagged.drop(columns=['gene_id'])
        features_all = features_all.sort_values(by=['ek_gene_id'])

        features_all = features_all.astype({'ek_gene_id': float})
        features_all = features_all.dropna().astype({'ek_gene_id': int})
        features_all = features_all.astype({'ek_gene_id': str})

        features_all = features_all[['feature', 'start', 'end', 'gene_biotype', 'strand', 'ek_gene_id']]

        # mark
        self.features_all = features_all

    def compose_compute_features(self):
        input_path = os.path.join(self.BASE, 'data', 'tmp', 'annotation')
        all_out = []
        for item in os.listdir(input_path):
            if item.startswith('computed_feature_'):
                print('Reading computed feature information for species: {}...'.format(item))
                tmp_df = pd.read_csv(os.path.join(input_path, item))
                all_out.append(tmp_df)
        all_out_df = pd.concat(all_out)
        all_out_df = all_out_df.sort_values(by=['gene_id'])

        all_out_df_tagged = pd.merge(all_out_df, self.id_mapper_tiny_no_dup, how='left', on='gene_id')
        first_column = all_out_df_tagged.pop('ek_gene_id')
        all_out_df_tagged.insert(0, 'ek_gene_id', first_column)
        computed_feature_all = all_out_df_tagged.drop(columns=['gene_id'])

        computed_feature_all = computed_feature_all.astype({'ek_gene_id': float})
        computed_feature_all = computed_feature_all.dropna().astype({'ek_gene_id': int})
        computed_feature_all = computed_feature_all.astype({'ek_gene_id': str})

        computed_feature_all = computed_feature_all[['feature', 'start', 'end', 'gene_biotype', 'strand', 'ek_gene_id']]

        # mark
        self.computed_feature_all = computed_feature_all

    def compose_pathway_meta(self):
        current_time = datetime.now().strftime("%D")
        all_db_id_type = {
            'msigdb': 'human_entrez_id',
            'go': 'gene_id',
            'interpro': 'gene_id',
            'kegg': 'entrez_id',
            'mesh': 'entrez_id',
            'reactome': 'entrez_id',
        }
        out = []
        for db in self.DB_UNITS:
            print('Formatting pathway meta information for: {}...'.format(db))
            pathways_dir_path = os.path.join(self.BASE, 'data', 'tmp', db)
            involve_list = os.listdir(pathways_dir_path)
            if db == 'msigdb':
                for species in self.SPECIES_DICT.keys():
                    tmp_tuple = [species, db, all_db_id_type[db], current_time]
                    out.append(tmp_tuple)
            else:
                for species_involve in involve_list:
                    if species_involve.endswith('_involve.txt'):
                        tmp_tuple = [species_involve.split('_')[0], db, all_db_id_type[db], current_time]
                        out.append(tmp_tuple)
        pathway_meta = pd.DataFrame(out)
        pathway_meta.columns = ['species', 'name', 'id_type', 'update_time']
        pathway_meta.reset_index(drop=True, inplace=True)
        pathway_meta.insert(0, 'pathway_meta_id', pathway_meta.index + 1)

        # mark
        self.pathway_meta = pathway_meta

    def compose_pathway(self):
        out = []
        for db in self.DB_UNITS:
            print('Formatting pathway information for: {}...'.format(db))
            pathways_dir_path = os.path.join(self.BASE, 'data', 'tmp', db)
            pathway_description_list = os.listdir(pathways_dir_path)
            if db == 'msigdb':
                for description in pathway_description_list:
                    if description.endswith('pathway_list.txt'):
                        tmp_df = pd.read_csv(os.path.join(pathways_dir_path, description))
                        tmp_df.insert(2, 'pathway_meta', 1)
                        out.append(tmp_df)

                next
            else:
                for description in pathway_description_list:
                    if description.endswith('_description.txt'):
                        tmp_db = db
                        tmp_species = description.split('_')[0]
                        if tmp_db == 'mesh' and tmp_species == 'all':
                            tmp_id_list = self.pathway_meta.loc[self.pathway_meta.name == db]['pathway_meta_id']
                            for id_mesh in tmp_id_list:
                                tmp_df = pd.read_csv(os.path.join(pathways_dir_path, description))
                                tmp_df.insert(2, 'pathway_meta', id_mesh)
                                out.append(tmp_df)
                        else:
                            tmp_id = self.pathway_meta.loc[(self.pathway_meta.species == tmp_species) & (self.pathway_meta.name == db)]['pathway_meta_id']
                            tmp_df = pd.read_csv(os.path.join(pathways_dir_path, description))
                            tmp_df.insert(2, 'pathway_meta', tmp_id.values[0])
                            out.append(tmp_df)

        all_pathway = pd.concat(out)
        all_pathway.reset_index(drop=True, inplace=True)
        all_pathway.insert(0, 'ek_pathway_id', all_pathway.index)

        # mark
        self.all_pathway = all_pathway

    def compose_involve(self):
        out = []
        self.id_mapper_all_updated = self.id_mapper_all_updated.astype(str)
        for db in self.DB_UNITS:
            print('Formatting involve information for: {}...'.format(db))
            pathways_dir_path = os.path.join(self.BASE, 'data', 'tmp', db)
            pathway_involve_list = os.listdir(pathways_dir_path)
            if db == 'msigdb':
                for involve in pathway_involve_list:
                    if involve.endswith('involve_list.txt'):
                        print('currently compose involve list for msigdb')
                        tmp_df = pd.read_csv(os.path.join(pathways_dir_path, involve))
                        tmp_df = tmp_df.astype({'gene_id': float})
                        tmp_df = tmp_df.dropna().astype({'gene_id': int})
                        tmp_df = tmp_df.astype({'gene_id': str})
                        tmp_df = tmp_df.rename(columns={"gene_id": "ek_gene_id"})
                        # print('tmp_df - ', tmp_df.head)

                        target_pathway_meta = 1

                        pathway_map_tmp = self.all_pathway[self.all_pathway['pathway_meta'] == target_pathway_meta]

                        tmp_df_2 = pd.merge(tmp_df, pathway_map_tmp[['ek_pathway_id', 'pathway_id']],
                                            how='left',
                                            on='pathway_id')
                        out_tmp = tmp_df_2[['ek_gene_id', 'ek_pathway_id']]

                        box = out_tmp
                        box['ek_gene_id'] = box['ek_gene_id'].astype(int)
                        for species_tmp in self.SPECIES_DICT.keys():
                            print('Getting gene id for msigdb for species - ', species_tmp)
                            ek_species_tmp = self.all_species_table.loc[self.all_species_table.name_short == species_tmp]['ek_species'].values[0]
                            species_mapper_tmp = self.id_mapper_all_updated.loc[self.id_mapper_all_updated.species == str(ek_species_tmp)]
                            t = species_mapper_tmp[['gene_id', 'human_entrez_id']].drop_duplicates()
                            t.replace('', np.nan, inplace=True)
                            tmp_mapper = t.dropna()
                            tmp_mapper['human_entrez_id'] = tmp_mapper['human_entrez_id'].astype(int)
                            mapping = {tmp_mapper.columns[0]: species_tmp}
                            tmp_mapper = tmp_mapper.rename(columns=mapping)
                            involve_list_raw_merge = pd.merge(box, tmp_mapper, how='left', left_on='ek_gene_id', right_on='human_entrez_id')
                            involve_list_raw_merge.pop('human_entrez_id')
                            box = involve_list_raw_merge
                        #
                        all_involve_msigdb = box
                        all_involve_msigdb = all_involve_msigdb.rename(columns={all_involve_msigdb.columns[0]: 'human_entrez_id'})
                        all_involve_msigdb.reset_index(drop=True, inplace=True)

                next
            else:
                for involve in pathway_involve_list:
                    if involve.endswith('involve.txt'):
                        tmp_db = db
                        tmp_species = involve.split('_')[0]
                        print('currently compose involve list of ', tmp_db, ' database for species ', tmp_species)
                        tmp_df = pd.read_csv(os.path.join(pathways_dir_path, involve)).drop_duplicates()
                        tmp_df = tmp_df.astype(str)

                        id_type = self.pathway_meta.loc[(self.pathway_meta.name == tmp_db) & (self.pathway_meta.species == tmp_species)]['id_type'].values[0]

                        mapper_tmp = self.id_mapper_all_updated[['ek_gene_id', id_type]].dropna().drop_duplicates()

                        tmp_df_1 = pd.merge(tmp_df, mapper_tmp,
                                            how='left',
                                            left_on='gene_id', right_on=id_type).dropna()

                        target_pathway_meta = self.pathway_meta.loc[(self.pathway_meta.name == tmp_db) & (self.pathway_meta.species == tmp_species)]['pathway_meta_id'].values[0]
                        pathway_map_tmp = self.all_pathway[self.all_pathway['pathway_meta'] == target_pathway_meta]

                        tmp_df_2 = pd.merge(tmp_df_1, pathway_map_tmp[['ek_pathway_id', 'pathway_id']],
                                            how='left',
                                            on='pathway_id')  # .dropna()

                        out_tmp = tmp_df_2[['ek_gene_id', 'ek_pathway_id']]
                        out.append(out_tmp)
        all_involve = pd.concat(out)
        all_involve.reset_index(drop=True, inplace=True)

        # mark
        self.all_involve = all_involve
        self.all_involve_msigdb = all_involve_msigdb
