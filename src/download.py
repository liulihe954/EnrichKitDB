import os

class downloader():
    #
    def __init__(self, RAW_DATA_INFO, UPDATE, BASE):
        #
        self.RAW_DATA_INFO = RAW_DATA_INFO
        self.UPDATE = UPDATE
        self.BASE = BASE
        
        
    def download_raw(self):
        #
        if self.UPDATE:
            for path, collection in self.RAW_DATA_INFO.items():
                #
                tmp_file_destination = os.path.join(self.BASE,'data','raw', path)
                #
                for file_name, url in collection.items():
                    #
                    file_destination = os.path.join(tmp_file_destination, file_name)
                    #
                    print(f'Now try to download {file_name} from source...')
                    print('\n In case of server down, use control + c to skip the current item...\n')
                    print('\n Remember which file(s) were missed...\n')
#                     print((f'wget -O {file_destination} {url}'))
                    os.system(f'wget -O {file_destination} {url}')
                    print('Done!\n\n')
        else:
            print('Please toggle the UPDATE parameter in config.py to proceed.')