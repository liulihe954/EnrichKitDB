import pandas as pd
import requests
import sys
import xmltodict
from lxml import etree
from io import StringIO
import argparse


def check_ensembl_symbol_helper(gene_id):
    print('checking ens symbol for ', gene_id)
    server = "https://rest.ensembl.org"
    ext = ''.join(['/xrefs/id/', gene_id, '?'])

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        return None

    decoded = r.json()

    symbol_tmp = decoded[0]['display_id']

    if symbol_tmp is not None and symbol_tmp != gene_id:
        return symbol_tmp
    else:
        return None


def from_ens_2_entrez(gene_id):
    print('checking ncbi entrez id and symbol for ', gene_id)
    target_url = 'https://www.ncbi.nlm.nih.gov/gene/?term=' + gene_id
    response = requests.get(target_url)
    dict_data = xmltodict.parse(response.content)
    package = dict_data['html']['body']['div']['div']['form']['div'][0]['div'][3]['div']['div'][5]
    entrez_id = package['div'][0]['div'][0]['span']['#text'].split(' ')[2][:-1]
    name = package['div'][1]['div'][0]['div']['div']['dl']['dd'][0]['#text']
    return entrez_id, name


def get_human_orth_ncbi(gene_id):
    print('checking ncbi human orthlogus for ', gene_id)
    target_url = 'https://www.ncbi.nlm.nih.gov/gene/?term=' + gene_id
    response = requests.get(target_url)
    dict_data = xmltodict.parse(response.content)
    package = dict_data['html']['body']['div']['div']['form']['div'][0]['div'][3]['div']['div'][5]
    for item in package['div'][1]['div'][0]['div']['div']['dl']['dd'][7]['a']:
        if item['#text'] == 'human':
            human_entrez_id = item['@href'].split('/')[-1]
    return human_entrez_id


def parse_human_url(human_entrez_id):
    print('parsing ncbi human orthlogus for', human_entrez_id)
    target_url = 'https://www.ncbi.nlm.nih.gov/gene/' + human_entrez_id
    response = requests.get(target_url)
    dict_data = xmltodict.parse(response.content)
    package = dict_data['html']['body']['div']['div']['form']['div'][0]['div'][3]['div']['div'][5]['div'][1]['div'][0]

    human_symbol = package['div']['div']['dl']['dd'][0]['#text']
    hgnc = package['div']['div']['dl']['dd'][2]['a']['#text']
    ens_id = package['div']['div']['dl']['dd'][3]['a'][0]['#text']
    if len(ens_id.split(':')) == 2:
        ens_id = ens_id.split(':')[1]

    return human_symbol, hgnc, ens_id


def it_child(node):
    for element in node.iterchildren():
        if element.tag == 'a' and 'href' in element.attrib and element.attrib['href'].startswith('http://www.ncbi.nlm.nih.gov/entrez'):
            return element.text
        it_child(element)


def from_ens_2_entrez_use_enblorg(tmp_gene_id):
    target_url = 'https://uswest.ensembl.org/Gene/Summary?g=' + tmp_gene_id
    page = requests.get(target_url)
    html = page.content.decode("utf-8")
    tree = etree.parse(StringIO(html), parser=parser)
    root = tree.getroot()
    return it_child(root)


def check_human_gene_id(gene_id):
    server = "https://rest.ensembl.org"
    ext = "/homology/id/" + gene_id + "?type=orthologues;format=condensed"
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()

    decoded = r.json()
    return decoded


def check_ensembl_symbol_helper(gene_id):
    print('checking ', gene_id)
    try:
        server = "https://rest.ensembl.org"
        ext = ''.join(['/xrefs/id/', gene_id, '?'])

        r = requests.get(server + ext, headers={"Content-Type": "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()

        symbol_tmp = decoded[0]['display_id']

        if symbol_tmp is not None and symbol_tmp != gene_id:
            return symbol_tmp
        else:
            return None
    except Exception:
        return None


def main(args):
    data_file_path = args.dataset
    # Where you run your code.
    id_mapper_final = pd.read_csv(data_file_path)

    # for index, row in id_mapper_final.iterrows():
    #     tmp_gene_id = id_mapper_final.at[index, 'gene_id']
    #     if tmp_gene_id is not None and tmp_gene_id != '':
    #         returned_symbol = check_ensembl_symbol_helper(tmp_gene_id)
    #         if id_mapper_final.at[index, 'ensembl_symbol'] != returned_symbol:
    #             id_mapper_final.at[index, 'ensembl_symbol'] = returned_symbol

    for index, row in id_mapper_final.iterrows():
        tmp_gene_id = id_mapper_final.at[index, 'gene_id']
        if tmp_gene_id is not None and tmp_gene_id != '':

            # check ensembl_symbol
            try:
                returned_ens_symbol = check_ensembl_symbol_helper(tmp_gene_id)
            except Exception:
                pass

            if id_mapper_final.at[index, 'ensembl_symbol'] != returned_ens_symbol:
                id_mapper_final.at[index, 'ensembl_symbol'] = returned_ens_symbol

    #         print(id_mapper_final.at[index, 'entrez_id'] == 'nan')
    #         print(len(id_mapper_final.at[index, 'entrez_id']))
    #         print(pd.isan(id_mapper_final.at[index, 'entrez_id']), '?')
            # check entrez id
            if id_mapper_final.at[index, 'entrez_id'] == 'nan':
                returned_entrez_id = None
                returned_entrez_symbol = None
                # ens 2 entrez: trust ncbi
                try:
                    returned_entrez_id, returned_entrez_symbol = from_ens_2_entrez(tmp_gene_id)
    #                 id_mapper_final.at[index, 'entrez_id'] = returned_entrez_id
    #                 id_mapper_final.at[index, 'ncbi_symbol'] = returned_entrez_symbol
                except Exception:
                    pass
                if returned_entrez_id == None or returned_entrez_symbol == None:
                    # ens 2 entrez: trust ensembl
                    try:
                        returned_entrez_id = from_ens_2_entrez_use_enblorg(tmp_gene_id)
                        _, returned_entrez_symbol = from_ens_2_entrez(returned_entrez_id)
                    except Exception:
                        pass

                if not returned_entrez_id is None and not returned_entrez_symbol is None:

                    id_mapper_final.at[index, 'entrez_id'] = returned_entrez_id
                    id_mapper_final.at[index, 'ncbi_symbol'] = returned_entrez_symbol

            # check human orthologus
            if not id_mapper_final.at[index, 'entrez_id'] == 'nan' and id_mapper_final.at[index, 'human_gene_id'] == 'nan':
                #             try:
                print(str(id_mapper_final.at[index, 'entrez_id']), 'ddd')
                human_entrez_id = get_human_orth_ncbi(str(int(id_mapper_final.at[index, 'entrez_id'])))

                id_mapper_final.at[index, 'human_entrez_id'] = human_entrez_id

                human_symbol, hgnc, ens_id = parse_human_url(str(human_entrez_id))

                id_mapper_final.at[index, 'hgnc_orthologs'] = hgnc

                id_mapper_final.at[index, 'human_gene_id'] = ens_id

                id_mapper_final.at[index, 'hgnc_symbol'] = human_symbol
    #             except Exception:
    #                 pass

    id_mapper_final.to_csv(''.join(['new_', data_file_path]), index=False, encoding='utf-8')


# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module
if __name__ == "__main__":
    # Note: this example shows named command line arguments.  See the argparse
    # documentation for positional arguments and other examples.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dataset',
                        help='input gene expression data file path',
                        type=str)
    args = parser.parse_args()
    # Note: this simply calls the main function above, which we could have
    # given any name
    main(args)
