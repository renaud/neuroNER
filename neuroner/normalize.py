import glob
import re
from itertools import groupby
from operator import itemgetter

from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

import oboparser
import os
import sys

'''
Cleans up raw Sherlok annotations
in: raw Sherlok output
out: only keep ontologyIds and if none is present concatenate the type and text
'''

def clean_annotations(n, orig_neuron_str = None):
    clean = []
    filt_attrib_list = ['Neuron', 'PreNeuron', 'PostNeuron', 'Electrophysiology', 'ProteinTrigger']
    for (begin, end, txt, type_, props) in n:
        if 'ontologyId' in props:
            clean.append((begin, end, props[u'ontologyId']))
        elif type_ not in filt_attrib_list:
            clean.append((begin, end, '{}:{}'.format(type_, txt)))
    clean = sorted(clean, key = lambda tup: tup[0])

    # remove unwanted UNKN_REGION ontology terms if a better ontology term available
    # also defer to NCBI_GENE over HBP_NEUROTRANSMITTER
    pos_list = [c[0] for c in clean]
    new_clean = clean
    for p in pos_list:
        indices = [i for i, x in enumerate(pos_list) if x == p]
        if len(indices) > 1:
            try:
                conflicting_terms = [clean[j] for j in indices]
                for c in conflicting_terms:
                    if 'UNKN' in c[2]:
                        new_clean.remove(c)
                    if 'HBP_NEUROTRANSMITTER' in c[2]:
                        new_clean.remove(c)
            except Exception:
                continue
    clean = new_clean

    # required because sherlok doesn't always return all missing terms
    if orig_neuron_str:
        # do some stuff to return missing terms - really crappy code
        nl = []
        for a in clean:
            for i in range(a[0],a[1]):
                nl.append(i)
            nl.append(i+1)
        unfound_inds = []
        for i in range(1,len(orig_neuron_str)):
            if i in nl:
                continue
            else:
                unfound_inds.append(i)
        data = unfound_inds
        ranges = []
        for k, g in groupby(enumerate(data), lambda (i,x):i-x):
            group = map(itemgetter(1), g)
            ranges.append((group[0], group[-1]))
        for r in ranges:
            l = [r[0], r[1], 'Missing:' + orig_neuron_str[(r[0]-1):(r[1]+1)].strip()]
            clean.append(l)

        clean = sorted(clean, key= lambda tup: tup[0])

    clean = [c[2] for c in clean]

    # convert projection annotations to new ontology
    proj_annots = [22,113,7322,7323,7324,7325]
    proj_annot_strs = ['HBP_PROJECTION:%s' % p for p in proj_annots]
    region_annot_strs = ['UNKN_REGION:%s' % p for p in proj_annots]
    new_clean = clean
    for i,c in enumerate(clean):
        for j,p in enumerate(region_annot_strs):
            if c == p:
                new_clean[i] = proj_annot_strs[j]
    clean = new_clean

    # filter out neuron triggers because not meaningful
    clean = [c for c in clean if 'NeuronTrigger' not in c]

    return clean

onto_root = os.path.dirname(sys.modules['neuroner'].__file__) + '/resources/bluima/neuroner/'
assert os.path.exists(onto_root), 'could not find onto_root at ' + onto_root

# load all ontologies into one large dictionary
def load_ontologies():
    """Loads all of the ontologies into a nice dictionary data structure"""

    # a massive dictionary containing key : dictionary mappings between HBP ontology id's and .obo ontology terms
    big_onto = {}
    mcc = MouseConnectivityCache()
    aba_onto = mcc.get_ontology()

    file_name_list = [f for f in glob.glob(onto_root + "*.robo")]
    file_name_list.extend([f for f in glob.glob(onto_root + "*.obo")])
    for fn in file_name_list:
        for o in oboparser.parse(fn):
            if 'id' in o:
                big_onto[o['id']] = o
    for k in big_onto.keys():
        if 'ABA' in k:
            new_o = big_onto[k]
            aba_id = int(k[11:])
            new_o['acronym'] = [aba_onto[aba_id]['acronym'].item()]
            big_onto[k] = new_o
    return big_onto


big_onto = load_ontologies()

def normalize_annots(onto_list, shorten = False):
    """Convenience function for turning a neuroNER ontology list back into ontology names or acronyms
    """
    normalized_term = []
    for l in onto_list:
        if l in big_onto:
            d = big_onto[l]
            if shorten:
                if 'acronym' in d:
                    normalized_term.append(d['acronym'][0])
                else:
                    normalized_term.append(d['name'])
            else:
                normalized_term.append(d['name'])
        else:
            new_term = re.sub('^\w+:', '', l)
            normalized_term.append(new_term)
    return ' '.join(normalized_term)