import os
import pandas as pd
from collections import defaultdict


class gtf_splitter():

    def __init__(self, BASE, SPECIES, gtf):
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
        self.gtf = gtf

    # caller
    def migrate_tables(self):
        # call manipulate_tables
        self.manipulate_tables()
        return self.gtf_skinny

    # handler
    def manipulate_tables(self):
        # init
        self.read_gtf()

        print('Now to format flatten_exon.....')
        self.flatten_exon()

        print('Now to format computated features.....')
        self.compute_feature()
        print('ok...')

        print('Now to output all tables.....')
        if self.SPECIES in self.SPECIES_DICT:
            #
            path_tmp = os.path.join(self.BASE, 'data', 'tmp', 'annotation')

            self.gtf_clean_exons['exon_number'] = self.gtf_clean_exons['exon_number'].astype('float').astype('Int64')

            self.gtf_clean_exons = self.gtf_clean_exons.drop(columns=['seqname'])
            self.gtf_clean_features = self.gtf_clean_features.drop(columns=['seqname'])

            if self.SPECIES == 'bta':
                self.gtf_clean_genes['seqname'] = self.gtf_clean_genes['seqname'].replace(['X'], 30)
            elif self.SPECIES == 'gal':
                self.gtf_clean_genes['seqname'] = self.gtf_clean_genes['seqname'].replace(['W'], 34).replace(['Z'], 35)
            elif self.SPECIES == 'sus':
                self.gtf_clean_genes['seqname'] = self.gtf_clean_genes['seqname'].replace(['X'], 19).replace(['Y'], 20)
            elif self.SPECIES == 'equ':
                self.gtf_clean_genes['seqname'] = self.gtf_clean_genes['seqname'].replace(['X'], 32)
            elif self.SPECIES == 'cap':
                self.gtf_clean_genes['seqname'] = self.gtf_clean_genes['seqname']
            elif self.SPECIES == 'ovi':
                self.gtf_clean_genes['seqname'] = self.gtf_clean_genes['seqname'].replace(['X'], 27)

            self.gtf_clean_genes['seqname'] = pd.to_numeric(self.gtf_clean_genes['seqname'], errors='coerce')

            #
            print(f'About to output gene list with dimension {self.gtf_clean_genes.shape}')
            self.gtf_clean_genes.to_csv(os.path.join(path_tmp, ''.join(['genes_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            #
            print(f'About to output exon list with dimension {self.gtf_clean_exons.shape}')
            self.gtf_clean_exons_falt.to_csv(os.path.join(path_tmp, ''.join(['exons_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            print(f'About to output feature list with dimension {self.gtf_clean_features.shape}')
            self.gtf_clean_features.to_csv(os.path.join(path_tmp, ''.join(['features_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

            print(f'About to output feature list with dimension {self.computatd_features.shape}')
            self.computatd_features.to_csv(os.path.join(path_tmp, ''.join(['computed_feature_', self.SPECIES, '.csv'])), index=False, encoding='utf-8')

        print('Done!')
        print('')

    ## read in everything ##

    def read_gtf(self):

        gtf_raw = self.gtf
        gtf_raw['start'] = gtf_raw['start'].astype('float').astype('Int64')
        gtf_raw['end'] = gtf_raw['end'].astype('float').astype('Int64')

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

        self.gtf_skinny = gtf_raw[gtf_raw.feature.isin(['gene'])][['gene_id', 'seqname', 'start', 'end', 'strand', 'gene_biotype']]

    #     # extract gene info
        gtf_clean_genes = gtf_clean[gtf_clean.feature.isin(['gene'])][['gene_id', 'seqname', 'start', 'end', 'strand', 'gene_biotype']]

    #     # extract exon info
        gtf_clean_exons = gtf_clean[gtf_clean.feature.isin(['exon'])][['gene_id', 'exon_id', 'seqname', 'start', 'end', 'exon_number', 'gene_biotype', 'strand', 'transcript_id']]

    #     # extract other feature info
        gtf_clean_features = gtf_clean[~gtf_clean.feature.isin(['gene', 'exon'])][['gene_id', 'feature', 'seqname', 'start', 'end', 'gene_biotype', 'strand']]

        ##
        self.gtf_clean_genes = gtf_clean_genes
        self.gtf_clean_exons = gtf_clean_exons
        self.gtf_clean_features = gtf_clean_features

    # helper to dedup exons
    def flatten_exon(self):
        exon_df = self.gtf_clean_exons.copy(deep=True)
        exon_df = exon_df.groupby('gene_id')
        # logic here
        print('About to flatten exons list...')
        container = []

        count = 0
        for name, group in exon_df:
            count += 1

            if count % 2000 == 0:
                print(f'flattend {count} genes...')

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

    def compute_feature(self):
        #
        exons_raw = self.gtf_clean_exons_falt.copy(deep=True)
        exons_raw = exons_raw.groupby('gene_id')

        #
        upstream_len = downstream_len = 10000
        donor_len = acceptor_len = 50
        all_computed_feature = []
        #
        count = 0
        for name, group in exons_raw:
            count += 1
            #
            tmp_df = group.sort_values(by=['start', 'end'])
            #
            if count % 2000 == 1:
                print(f'annotated {count} genes...')

            # for each gene
            length = len(tmp_df.index)
            row = 0
            next_row = 0
            tmp_gene_id = tmp_df.iloc[row][0]
            tmp_gene_type = tmp_df.iloc[row][4]
            tmp_gene_strand = tmp_df.iloc[row][5]
            #
            tmp_compute_feature = []

            tmp_compute_feature.append([tmp_gene_id,
                                        'upstream',
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
                        tmp_gene_id,
                        'short intron',
                        tmp_intron_start,
                        tmp_intron_end,
                        tmp_gene_type,
                        tmp_gene_strand])
                else:
                    # donor
                    tmp_compute_feature.append(
                        [tmp_gene_id,
                         'splice donor',
                         tmp_intron_start,
                         tmp_intron_start + donor_len,
                         tmp_gene_type,
                         tmp_gene_strand])

                    # intron
                    tmp_compute_feature.append(
                        [tmp_gene_id,
                         'intron',
                         tmp_intron_start + donor_len + 1,
                         tmp_intron_end - donor_len - 1,
                         tmp_gene_type,
                         tmp_gene_strand])

                    # acceptor
                    tmp_compute_feature.append(
                        [tmp_gene_id,
                         'splice acceptor',
                         tmp_intron_end - acceptor_len,
                         tmp_intron_end,
                         tmp_gene_type,
                         tmp_gene_strand])
            tmp_compute_feature.append(
                [tmp_gene_id,
                 'downstream',
                 tmp_df.iloc[length - 1][3] + 1,
                 tmp_df.iloc[length - 1][3] + 1 + downstream_len,
                 tmp_gene_type,
                 tmp_gene_strand])

            tmp_out = pd.DataFrame(tmp_compute_feature, columns=['gene_id', 'feature', 'start', 'end', 'gene_biotype', 'strand'])

            all_computed_feature.append(tmp_out)

        # output
        self.computatd_features = pd.concat(all_computed_feature)
