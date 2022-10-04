from collections import defaultdict
import gzip
import pandas as pd
import re
from urllib.request import urlopen
import certifi
import ssl

# Code reference:
# https://gist.github.com/slowkow/8101481


class gtf_parser:

    """An class to prase the genome annotation file (.gtf)
    """

    def __init__(self, input_file):
        self.input_file = input_file
        self.GTF_HEADER = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
        self.R_SEMICOLON = re.compile(r'\s*;\s*')
        self.R_COMMA = re.compile(r'\s*,\s*')
        self.R_KEYVALUE = re.compile(r'(\s+|\s*=\s*)')

    def lines(self):
        """reading the pgzipped GTF file from source (e.g., ensembl.org) and generate a dict for each line.
        """

#         Additional feature: stream in from web source

        #print('Started streaming file from ', self.stream_url)
        #ssl._create_default_https_context = ssl._create_unverified_context
        #streamed_file = urlopen(self.stream_url)
        #print('request status: ', streamed_file.status)

        opener = gzip.open if self.input_file.endswith('.gz') else open

        print('Started reading lines from', self.input_file)

        with opener(self.input_file) as f:
            for line in f:
                if line.decode('utf8').startswith('#'):
                    continue
                else:
                    yield self.parse(line)

    def parse(self, line):
        """Parse a single GTF line and return a dict.
        """
        result = {}

        fields = line.decode('utf8').rstrip().split('\t')

        for i, col in enumerate(self.GTF_HEADER):
            result[col] = self._get_value(fields[i])

        # INFO field consists of kv paris like " key1=value; key2=value; ...".
        infos = [x for x in re.split(self.R_SEMICOLON, fields[8]) if x.strip()]

        for i, info in enumerate(infos, 1):
            try:  # It should be key="value".
                key, _, value = re.split(self.R_KEYVALUE, info, 1)
            except ValueError:  # But sometimes it is just "value".
                key = 'INFO{}'.format(i)
                value = info
            if value:  # Ignore the field if there is no value.
                result[key] = self._get_value(value)

        return result

    ##
    def _get_value(self, value):

        if not value:
            return None
        value = value.strip('"\'')  # Strip double and single quotes.

        # Return a list if the value has a comma.
        if ',' in value:
            value = re.split(self.R_COMMA, value)

        # These values are equivalent to None.
        elif value in ['', '.', 'NA']:
            return None

        return value

    def dataframe(self):
        """Convert to a return a pandas.DataFrame.
        """
        # Each column is a list stored as a value in this dict.
        result = defaultdict(list)
        print('About to start parsing...')
        ct = 0
        for i, line in enumerate(self.lines()):

            if i % 200000 == 0:
                ct += 1
                print('parsed', 200 * ct, 'k records...')

            for key in line.keys():
                # This key has not been seen yet, so set it to None for all previous lines.
                if key not in result:
                    result[key] = [None] * i

            # Ensure this row has some value for each column.
            for key in result.keys():
                result[key].append(line.get(key, None))

        print('Done!')

        return pd.DataFrame(result)
