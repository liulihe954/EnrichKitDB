import os
from config import RAW_DATA_INFO, UPDATE, MAX_THREADS, CHECK_WEB
from workflow.steps import work_flow
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#
BASE = os.getcwd()

SPECIES = ['bta', 'cap', 'equ', 'gal', 'ovi', 'sus']
units = ['annotation', 'id_mapper', 'go', 'kegg', 'interpro', 'mesh', 'reactome', 'msigdb','tftargets']
SQLITE_NAME = 'EnrichKitDB_test.sqlite'


path_list = [
    os.path.join(BASE, 'sqlite'),
    os.path.join(BASE, 'data')
]

path_list.extend([os.path.join(path_list[1], 'raw', x) for x in units])
path_list.extend([os.path.join(path_list[1], 'tmp', x) for x in units])
path_list.extend([os.path.join(path_list[1], 'output')])
print(f'Now try to create {len(path_list)} paths for the DB...')
print()
for item in path_list:
    if not os.path.exists(item):
        os.makedirs(item)
        print(f'created path {item}')
    else:
        print(f'{item} already exists.')
print()

# internal mapping

SPECIES_DICT = {
    'bta': [0, 'Bos_taurus', 'ARS-UCD1.3.113'],
    'cap': [1, 'Capra_hircus', 'ARS1.113'],
    'equ': [2, 'Equus_caballus', 'EquCab3.0.113'],
    'gal': [3, 'Gallus_gallus', 'bGalGal1.mat.broiler.GRCg7b.113'],
    'ovi': [4, 'Ovis_aries', 'ARS-UI_Ramb_v2.0.113'],
    'sus': [5, 'Sus_scrofa', 'Sscrofa11.1.113'],
    }

print("--------------------------------")
print("  INITIALIZING A NEW WORK FLOW  ")
print("--------------------------------")
current = work_flow(RAW_DATA_INFO, UPDATE, BASE, SPECIES, MAX_THREADS, SQLITE_NAME, SPECIES_DICT, CHECK_WEB)

print('DONE!')
print("------------------------------\n\n")

if UPDATE:
    # 1. download raw data from source if the purpose is to UPDATE
    print("--------------------------------")
    print("  DOWNLOADING FILES FROM SOURCE ")
    print("--------------------------------")
    # current.download()
    print('DONE!\n')
    print("------------------------------\n\n")

# 2. parse raw data and generate priliminary tables
print("--------------------------------")
print("PARSING CURRENT VERSION RAW DATA")
print("--------------------------------")
# current.parse()
print('DONE!\n')
print("------------------------------\n\n")

# 3. assign UID(ek_id) and generate schema
print("--------------------------------")
print("    GENERATING DATABASE SCHEMA  ")
print("--------------------------------")
current.migrate()
print('DONE!\n')
print("------------------------------\n\n")

# 4. process tf targets
print("--------------------------------")
print("    PROCESSING TFTARGETS        ")
print("--------------------------------")
current.tftarget()
print('DONE!\n')
print("------------------------------\n\n")

# 5. sqlite database
print("----------------------------------")
print(" Create a portable sqlite database")
print("----------------------------------")
# current.populate()
print('DONE!\n')
print("------------------------------\n\n")

