import json
from ete3 import PhyloTree, TreeStyle

from skbio.diversity.alpha import fisher_alpha as alpha_function

def normalise_row(row):
    '''The older json files don't have the correct types'''
    int_fields = ('taxon', 'parent', 'count_norm')
    for k in int_fields:
        v = row.get(k)
        if v is not None:
            row[k] = int(v)

def clean_json(fp):
    '''Some of the files appear to have HTML junk at the start and end of the file.'''
    data = ''.join(line.rstrip() for line in fp)
    start = data.find('{')
    end = data.rfind('}')
    return data[start:end + 1]

def load_json(fp):
    data = json.loads(clean_json(fp))
    taxonomy = {}
    count_total = 0
    counts = []

    for row in data['ubiome_bacteriacounts']:
        normalise_row(row)
        counts.append(row['count_norm'])
        t = PhyloTree()
        t.name = row['tax_name']
        t.add_features(**row)
        taxonomy[row['taxon']] = t

    root = taxonomy[min(taxonomy.keys())]
    count_total = root.count_norm
    root.alpha = alpha_function(counts)

    for t in taxonomy.values():
        t.add_feature('count_pct', float(t.count_norm) / count_total * 100)
        parent = t.parent
        tp = taxonomy.get(parent)
        if tp is not None:
            tp.add_child(t)
    print('loaded {} into tree depth {} diversity {:.2f}'.format(len(taxonomy), len(root), root.alpha))
    return root

def get_bacteria_count(root, name):
    nodes = root.search_nodes(name=name)
    if len(nodes) == 0:
        return 0
    return nodes[0].count_pct

def check_bacteria_(root, name, lower=None, upper=None):
    nodes = root.search_nodes(name=name)
    if len(nodes) == 0:
        raise KeyError(name)
    count_pct = nodes[0].count_pct
    print('{:5.2f} '.format(count_pct, c='.'), end='')
    if lower is not None and count_pct < lower:
        return -1
    if upper is not None and count_pct > upper:
        return 1
    return 0

def check_bacteria(root, name, lower=None, upper=None):
    print('Checking for {:.<40}'.format(name, c='.'), end='')
    try:
        retval = check_bacteria_(root, name, lower=lower, upper=upper)
    except KeyError:
        print('[NOT FOUND]')
        return
    if retval < 0:
        s = '[LOW]'
    elif retval > 0:
        s = '[HIGH]'
    else: s = '[OK]'
    print(s)

__all__ = ['normalise_row', 'clean_json', 'load_json', 'check_bacteria',
           'get_bacteria_count']
