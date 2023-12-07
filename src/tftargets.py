from collections import defaultdict
import pandas as pd
import os

class tftarget:
    def __init__(self, BASE):
        self.BASE = BASE
        self.preprocess_r_path = os.path.join(BASE, 'src', 'process_tf_enrichkit.R')
        self.input_path = os.path.join(BASE, 'data', 'tmp','tftargets','tftargets.csv')
        self.output_path = os.path.join(BASE, 'data', 'output','tf_targets.txt')
        self.hgnc_path = os.path.join(BASE, 'data', 'raw','id_mapper','hgnc_protein_coding_gene.txt')

    def process_data_R(self):
        r_script = 'Rscript '+ self.preprocess_r_path + ' ' + self.BASE
        # print(self.preprocess_r_path)
        os.system(r_script)

    def process_data_py(self):
        #
        hgnc_ref = {
            'hgnc_orthologs': [],
            'human_gene_id': [],
            'human_entrez_id': [],
            'hgnc_symbol': []
        }

        count_tmp = 0
        with open(self.hgnc_path) as f:
            for line in f:
                count_tmp += 1
        
                if count_tmp % 10000 == 0:
                    print(f'parsed {count_tmp} records...')
        
                tmp_line = line.split('\t')
        
                if tmp_line[0] == 'hgnc_id':
                    continue
        
                hgnc_ref['hgnc_orthologs'].append(tmp_line[0])
                hgnc_ref['human_gene_id'].append(tmp_line[19])
                hgnc_ref['human_entrez_id'].append(tmp_line[18])
                hgnc_ref['hgnc_symbol'].append(tmp_line[1])
        hgnc_ref_df = pd.DataFrame(hgnc_ref)

        ##
        tf_df = pd.read_csv(self.input_path)
        tf_df_new = tf_df.copy(deep=True)
        #
        converted_rows = []
        miss_count = 0
        all_count = 0
        all_genes = []
        for index, row in tf_df_new.iterrows():
            source = row['Source']
            if source in ['ENCODE','TRED']:
                target_str = row['Gene']
                temp_list = target_str.split(',')
                temp_list_convert = []
                for item in temp_list:
                    all_count += 1
        #             hgnc_temp = hgnc_ref.loc[hgnc_ref['human_entrez_id'] == str(item)]['hgnc_symbol'].values[0]
                    try:
                        hgnc_temp = hgnc_ref_df.loc[hgnc_ref_df['human_entrez_id'] == str(item)]['hgnc_symbol'].values[0]
                        # print(hgnc_temp)
                        temp_list_convert.append(hgnc_temp)
                        if not hgnc_temp in all_genes:
                            all_genes.append(hgnc_temp)
                    except:
                        miss_count += 1
                convert_str = ','.join(temp_list_convert)
                tf_df_new.at[index,'Gene'] = convert_str
                
        tf_df_new_order = tf_df_new[['Source','TF','Gene']]
        
        #
        last_row_dict = {
            'Source':'ALL',
            'TF':'ALL',
            'Gene': ','.join(all_genes)
            
        }
        last_row_df = pd.DataFrame(last_row_dict, index=[0])
        output_df = pd.concat([tf_df_new_order, last_row_df], ignore_index=True)
        output_df.to_csv(self.output_path, index=False)