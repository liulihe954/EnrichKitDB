import os
import pandas as pd
from datetime import datetime


class migrator:

    def __init__(self, BASE):
        self.BASE = BASE
        self.DB_UNITS = ['msigdb', 'go', 'interpro', 'kegg', 'mesh', 'reactome']
        self.SPECIES_DICT = {
            'bta': [0, 'Bos_taurus', 'ARS-UCD1.2.107'],
            'cap': [1, 'Capra_hircus', 'ARS1.107'],
            'ovi': [2, 'Ovis_aries', 'Oar_v3.1.107'],
            'gal': [3, 'Gallus_gallus', 'gca000002315v5.GRCg6a.107'],
            'sus': [4, 'Sus_scrofa', 'Sscrofa11.1.107'],
            'equ': [5, 'Equus_caballus', 'EquCab3.0.107'], }

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
            [2, 'ovi', 'Ovis_aries', 'Oar_v3.1.107'],
            [3, 'gal', 'Gallus_gallus', 'gca000002315v5.GRCg6a.107'],
            [4, 'sus', 'Sus_scrofa', 'Sscrofa11.1.107'],
            [5, 'equ', 'Equus_caballus', 'EquCab3.0.107'],
        ])
        all_species_table.columns = ['name_short', 'ek_species', 'name_long', 'current_gtf']

        # mark
        self.all_species_table = all_species_table

    def compose_id_mapper(self):
        input_path = os.path.join(self.BASE, 'data', 'tmp', 'id_mapper')
        id_mapper_all_updated = pd.read_csv(os.path.join(input_path, 'id_mapper_all_updated.txt'), low_memory=False)

        id_mapper_all_updated.pop('Unnamed: 0')
        id_mapper_all_updated.fillna('')
        id_mapper_all_updated = id_mapper_all_updated.sort_values(by=['species', 'gene_id'])
        id_mapper_all_updated.reset_index(drop=True, inplace=True)
        id_mapper_all_updated['ek_gene_id'] = id_mapper_all_updated.index

        # id_mapper_all_updated = id_mapper_all_updated.astype({'entrez_id': float})
        # id_mapper_all_updated = id_mapper_all_updated.astype({'entrez_id': int})
        # id_mapper_all_updated = id_mapper_all_updated.astype({'entrez_id': str})

        # id_mapper_all_updated = id_mapper_all_updated.astype({'human_entrez_id': float})
        # id_mapper_all_updated = id_mapper_all_updated.astype({'human_entrez_id': int})
        # id_mapper_all_updated = id_mapper_all_updated.astype({'human_entrez_id': str})

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
                tmp_tuple = ['all', db, all_db_id_type[db], current_time]
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
                        # print('out_tmp - ', out_tmp.head)
                        out.append(out_tmp)
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
                        if id_type == 'entrez_id':
                            mapper_tmp = mapper_tmp.astype({id_type: float})
                            mapper_tmp = mapper_tmp.dropna().astype({id_type: int})
                            mapper_tmp = mapper_tmp.astype({id_type: str})

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
