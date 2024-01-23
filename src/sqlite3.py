import sqlite3
import csv
import os

class create_sqlite:
    def __init__(self,BASE, db_name):
        self.BASE = BASE
        self.df_path = os.path.join(self.BASE, 'sqlite', db_name) #'EnrichKitDB.sqlite'
        self.file_names = [("Species","species.txt"),
                           ("ID_Mapper","id_mapper.txt"),
                           ("Gene","genes.txt"),
                           ("Exon","exons.txt"),
                           ("ComputedFeatures","computed_features.txt"),
                           ("Feature","features.txt"),
                           ("Pathway_Meta","pathway_meta.txt"),
                           ("Pathway","pathway.txt"),
                           ("Involve","involve.txt"),         
                           ("InvolveM","involve_msigdb.txt"),
                           ("TF_Gene","tf_targets.txt")]
        self.create_table_statements = [
            "CREATE TABLE Species (ek_species INTEGER PRIMARY KEY, name_short TEXT, name_long TEXT, current_gtf TEXT);",
            "CREATE TABLE ID_Mapper (ek_gene_id INTEGER PRIMARY KEY, species INTEGER, gene_id TEXT, ensembl_symbol TEXT, entrez_id TEXT, ncbi_symbol TEXT, vgnc_id TEXT, vgnc_symbol TEXT, hgnc_orthologs TEXT, human_gene_id TEXT, human_entrez_id TEXT, hgnc_symbol TEXT, FOREIGN KEY(species) REFERENCES Species(ek_species));",
            "CREATE TABLE Gene (seqname TEXT, start INTEGER, end INTEGER, strand TEXT, gene_biotype TEXT, ek_gene_id INTEGER, species INTEGER, FOREIGN KEY(ek_gene_id) REFERENCES ID_Mapper(ek_gene_id), FOREIGN KEY(species) REFERENCES Species(ek_species));",
            "CREATE TABLE Exon (ek_gene_id INTEGER, exon_id TEXT, start INTEGER, end INTEGER, gene_biotype TEXT, strand TEXT, transcript_id TEXT, FOREIGN KEY(ek_gene_id) REFERENCES ID_Mapper(ek_gene_id));",
            "CREATE TABLE ComputedFeatures (ek_gene_id INTEGER, feature TEXT, start INTEGER, end INTEGER, gene_biotype TEXT, strand TEXT, FOREIGN KEY(ek_gene_id) REFERENCES ID_Mapper(ek_gene_id));",
            "CREATE TABLE Feature (ek_gene_id INTEGER, feature TEXT, start INTEGER, end INTEGER, gene_biotype TEXT, strand TEXT, FOREIGN KEY(ek_gene_id) REFERENCES ID_Mapper(ek_gene_id));",
            "CREATE TABLE Pathway_Meta (pathway_meta_id INTEGER PRIMARY KEY, species TEXT, name TEXT, id_type TEXT, update_time TEXT);",
            "CREATE TABLE Pathway (ek_pathway_id INTEGER PRIMARY KEY, pathway_id TEXT, pathway_description TEXT, pathway_meta INTEGER, FOREIGN KEY(pathway_meta) REFERENCES Pathway_Meta(pathway_meta_id));",
            "CREATE TABLE Involve (ek_gene_id INTEGER, ek_pathway_id INTEGER, FOREIGN KEY(ek_gene_id) REFERENCES ID_Mapper(ek_gene_id), FOREIGN KEY(ek_pathway_id) REFERENCES Pathway(ek_pathway_id));",
            "CREATE TABLE InvolveM (human_entrez_id INTEGER, ek_pathway_id INTEGER, bta TEXT, cap TEXT, equ TEXT, gal TEXT, ovi TEXT, sus TEXT, FOREIGN KEY(ek_pathway_id) REFERENCES Pathway(ek_pathway_id));",
            "CREATE TABLE TF_Gene (source TEXT, TF TEXT, Gene TEXT);",
            "CREATE INDEX idx_gene_start ON Gene(start);",
            "CREATE INDEX idx_gene_end ON Gene(end);",
            "CREATE INDEX idx_exon_start ON Exon(start);",
            "CREATE INDEX idx_exon_end ON Exon(end);",
            "CREATE INDEX idx_computedfeatures_start ON ComputedFeatures(start);",
            "CREATE INDEX idx_computedfeatures_end ON ComputedFeatures(end);",
            "CREATE INDEX idx_feature_start ON Feature(start);",
            "CREATE INDEX idx_feature_end ON Feature(end);"
            ]
    def create_tables(self):
        if not os.path.exists(self.df_path):
            db_conn = sqlite3.connect(self.df_path)
            cursor = db_conn.cursor()
            for statement in self.create_table_statements:
                cursor.execute(statement)
            db_conn.commit()
            cursor.close()
            db_conn.close()

    def load_all_tables(self):
        # list table names and paths
        for file_tuple in self.file_names:
            file_path = os.path.join(self.BASE, 'data', 'output', file_tuple[1])
            #
            print('Output to sqliate: about to working ', file_tuple[1])
            #
            db_conn = sqlite3.connect(self.df_path)
            self.load_data_into_table(db_conn, file_tuple[0], file_path)
            

    def load_data_into_table(self, db_conn, table_name, file_path):
        #
        cursor = db_conn.cursor()
        with open(file_path, 'r') as file:
            reader = csv.reader(file)
            headers = next(reader)
            placeholders = ','.join('?' * len(headers))
            insert_sql = f'INSERT INTO {table_name} ({','.join(headers)}) VALUES ({placeholders})'

            for row in reader:
                cursor.execute(insert_sql, row)

        db_conn.commit()
        cursor.close()
