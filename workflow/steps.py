import sys
from src.download import downloader
from src.full_scale_parser import full_scale_parser
from src.migrator import migrator
from src.tftargets import tftarget


# parse
#from parse import parser

# migrate
#from migrate import migrater

# populate
# from populate import injector


class work_flow():

    def __init__(self, RAW_DATA_INFO, UPDATE, BASE, SPECIES):
        #
        self.RAW_DATA_INFO = RAW_DATA_INFO
        self.UPDATE = UPDATE
        self.BASE = BASE
        self.SPECIES = SPECIES

    def download(self):
        #
        d = downloader(self.RAW_DATA_INFO, self.UPDATE, self.BASE)
        #
        d.download_raw()

    def parse(self):
        #
        p = full_scale_parser(self.BASE)
        
        p.do_parse()
        # can also do single work
        #   .parse_gtf()
        #   .parse_id_mapper()
        #   .parse_go() 
        #   .parse_interpro()
        #   .parse_kegg()
        #   .parse_reactome()
        #   .parse_mesh()
        #   .parse_msigdb()
        #   .parse_reactome()
        #   .parse_id_mapper()
        #   .parse_kegg()

    def migrate(self):
        m = migrator(self.BASE)
        m.generate_tables()

    def tftarget(self):  
        t = tftarget(self.BASE)
        t.process_data_R()
        t.process_data_py()
        
    def populate():
        pass
