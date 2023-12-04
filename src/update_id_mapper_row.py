import pandas as pd
import requests
import sys
import xmltodict
from lxml import etree
from io import StringIO

pd.set_option('mode.chained_assignment',None)

def check_ensembl_symbol_helper(gene_id):
    # print('checking ens symbol for ', gene_id)
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
    # print('checking ncbi entrez id and symbol for ', gene_id)
    target_url = 'https://www.ncbi.nlm.nih.gov/gene/?term=' + gene_id
    response = requests.get(target_url)
    dict_data = xmltodict.parse(response.content)
    package = dict_data['html']['body']['div']['div']['form']['div'][0]['div'][3]['div']['div'][5]
    entrez_id = package['div'][0]['div'][0]['span']['#text'].split(' ')[2][:-1]
    name = package['div'][1]['div'][0]['div']['div']['dl']['dd'][0]['#text']
    return entrez_id, name


def get_human_orth_ncbi(gene_id):
    # print('checking ncbi human orthlogus for ', gene_id)
    target_url = 'https://www.ncbi.nlm.nih.gov/gene/?term=' + gene_id
    response = requests.get(target_url)
    dict_data = xmltodict.parse(response.content)
    package = dict_data['html']['body']['div']['div']['form']['div'][0]['div'][3]['div']['div'][5]
    for item in package['div'][1]['div'][0]['div']['div']['dl']['dd'][7]['a']:
        if item['#text'] == 'human':
            human_entrez_id = item['@href'].split('/')[-1]
    return human_entrez_id


def parse_human_url(human_entrez_id):
    # print('parsing ncbi human orthlogus for', human_entrez_id)
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
    parser = etree.HTMLParser(recover=True)
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
    # print('checking ensembl symbol for', gene_id)
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
def from_hgncid_to_hgncsymb(hgnc_id):
    server = 'https://rest.genenames.org/'
    ext = 'search/hgnc_id/' + str(hgnc_id.split(":")[-1])

    r = requests.get(server + ext)

    if not r.ok:
        r.raise_for_status()

#     decoded = r.content.decode()
    dict_data = xmltodict.parse(r.content)
    return dict_data['response']['result']['doc']['str'][1]['#text']

def human_from_entrez_to_ens(returned_entrez_id):
    server = 'https://www.ncbi.nlm.nih.gov/'
    ext = 'gene/' + str(returned_entrez_id)

#     parser = etree.HTMLParser(recover=True)
    page = requests.get(server + ext)

    if not page.ok:
        page.raise_for_status()
    dict_data = xmltodict.parse(page.content)
    package = dict_data['html']['body']['div']['div']['form']['div'][0]['div'][3]['div']['div'][5]
    returned_ensid = package['div'][1]['div'][0]['div']['div']['dl']['dd'][3]['a'][0]['#text']

    return returned_ensid.split(':')[-1]

def it_child2(node, found):
    if node == None:
        return
    if len(found) == 0:
        for element in node.iterchildren():
            if 'data-section' in element.attrib and element.attrib['data-section'] == 'Featured':
                found.append(element.attrib['data-item-id'])
            it_child2(element, found)
    else:
        return found
    return found

def human_from_symbol_2entrez(human_symbol):
    server = 'https://www.ncbi.nlm.nih.gov/'
    ext = 'search/all/?term=' + 'TDH'
    parser = etree.HTMLParser(recover=True)
    page = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not page.ok:
        page.raise_for_status()

    html = page.content.decode("utf-8")
    tree = etree.parse(StringIO(html), parser=parser)
    root = tree.getroot()

    returned_entrez_id = None
    returned_ens_id = None

    try:
        found = []
        returned_entrez_id = it_child2(root,found)[0].split(':')[-1]
    except Exception:
        pass

    try:
        returned_ens_id = human_from_entrez_to_ens(returned_entrez_id)
    except Exception:
        pass

    return returned_entrez_id, returned_ens_id


def update_row(input_row):
    # convert to series
    row = pd.Series(input_row)

    if row.isna().sum() == 0:
        print(row['gene_id'], ' is complete.')
    else:
        if not pd.isna(row['gene_id']):
            if pd.isna(row['ensembl_symbol']):
                try:
                    row['ensembl_symbol'] = check_ensembl_symbol_helper(row['gene_id'])
                except Exception:
                    pass

            if pd.isna(row['entrez_id']):
                try:
                    returned_entrez_id, returned_entrez_symbol = from_ens_2_entrez(row['gene_id'])
                    if returned_entrez_id is not None and returned_entrez_symbol is not None:
                        row['entrez_id'] = returned_entrez_id
                        row['ncbi_symbol'] = returned_entrez_symbol
                except Exception:
                    pass

                if pd.isna(row['entrez_id']):
                    try:
                        returned_entrez_id = from_ens_2_entrez_use_enblorg(row['gene_id'])
                        if returned_entrez_id is not None:
                            _, returned_entrez_symbol = from_ens_2_entrez(returned_entrez_id)
                            row['entrez_id'] = returned_entrez_id
                            row['ncbi_symbol'] = returned_entrez_symbol
                    except Exception:
                        pass

        if not pd.isna(row['entrez_id']) and pd.isna(row['human_gene_id']):
            try:
                human_entrez_id = get_human_orth_ncbi(str(int(float(row['entrez_id']))))
                if human_entrez_id:
                    row['human_entrez_id'] = human_entrez_id
                    human_symbol, hgnc, ens_id = parse_human_url(str(human_entrez_id))
                    row['hgnc_orthologs'] = hgnc
                    row['human_gene_id'] = ens_id
                    row['hgnc_symbol'] = human_symbol
            except Exception:
                pass

        if not pd.isna(row['hgnc_orthologs']) and pd.isna(row['hgnc_symbol']):
            try:
                row['hgnc_symbol'] = from_hgncid_to_hgncsymb(row['hgnc_orthologs'])
            except Exception:
                pass

        if not pd.isna(row['hgnc_symbol']) and (pd.isna(row['human_gene_id']) or pd.isna(row['human_entrez_id'])):
            try:
                returned_entrez_id_human, returned_ens_id_human = human_from_symbol_2entrez(row['hgnc_symbol'])
                if returned_entrez_id_human is not None:
                    row['human_entrez_id'] = returned_entrez_id_human
                if returned_ens_id_human is not None:
                    row['human_gene_id'] = returned_ens_id_human
            except Exception:
                pass

    return row.to_dict()