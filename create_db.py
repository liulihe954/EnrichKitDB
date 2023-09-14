import os
from config import RAW_DATA_INFO, UPDATE
from workflow.steps import work_flow

# Create structure
BASE = os.getcwd()

SPECIES = ['sus']
units = ['annotation', 'id_mapper', 'go', 'kegg', 'interpro', 'mesh', 'reactome', 'msigdb']

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


print("--------------------------------")
print("  INITIALIZING A NEW WORK FLOW  ")
print("--------------------------------")
current = work_flow(RAW_DATA_INFO, UPDATE, BASE, SPECIES)

print('DONE!')
print("------------------------------\n\n")

if UPDATE:
    # 1. download raw data from source if the purpose is to UPDATE
    print("--------------------------------")
    print("  DOWNLOADING FILES FROM SOURCE ")
    print("--------------------------------")
    current.download()
    print('DONE!\n')
    print("------------------------------\n\n")

# 2. parse raw data and generate priliminary tables
print("--------------------------------")
print("PARSING CURRENT VERSION RAW DATA")
print("--------------------------------")
current.parse()
print('DONE!\n')
print("------------------------------\n\n")

# 3. assign UID(ek_id) and generate schema
print("--------------------------------")
print("    GENERATING DATABASE SCHEMA  ")
print("--------------------------------")
current.migrate()
print('DONE!\n')
print("------------------------------\n\n")

# 4. insert into sqlite and generate portable instance
# current.populate()
