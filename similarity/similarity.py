'''
Computes the semantic similarity between two neuron mentions.
The main function is similarity()
Computation is delegated to similarity_inter and similarity_intra
'''

import similarity_inter, similarity_intra
from operator import itemgetter
from itertools import groupby
import re

from sherlok import Sherlok # pip install --upgrade sherlok
s = Sherlok('neuroner')

WEIGHTS = { #TODO: implement weights
    'Layer': 1.0,
    'ProteinProp': 1.0
}

'''
Computes the intra and inter semantic similarity between two neurons
in: n1@str, n2@str: the two neurons to measure similarity
out: (score:float, [([matching_properties], explanation@str)])
'''
def similarity(n1, n2, weights=WEIGHTS):
    assert type(n1) is str and len(n1) > 0, 'n1 cannot be empty'
    assert type(n2) is str and len(n2) > 0, 'n2 cannot be empty'
    n1c = _cleanup(s.annotate(n1).annotations)
    n2c = _cleanup(s.annotate(n2).annotations)
    return similarity2(n1c, n2c)

'''
Cleans up raw Sherlok annotations
in: raw Sherlok output
out: only keep ontologyIds and if none is present concatenate the type and text
'''
def _cleanup(n, orig_neuron_str = None):
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
    # filter out neuron triggers because not meaningful
    clean = [c for c in clean if 'NeuronTrigger' not in c]

    return clean

# load all ontologies into one large dictionary
big_onto = similarity_intra.load_ontologies()

def _normalize(onto_list, shorten = False):
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



'''
Computes the intra and inter semantic similarity between two neurons
in:  n1@[], n2@[]: the two neurons to measure similarity
out: (score:float, [([matching_properties, explanation@str])])
'''
def similarity2(n1, n2, weights=WEIGHTS):
    s_intra = similarity_intra._similarity_intra(n1, n2, weights)
    s_inter = similarity_inter._similarity_inter(n1, n2)
    print('s_intra', s_intra, 's_inter', s_inter)
    return (s_intra[0] + s_inter[0], s_intra[1] + s_inter[1])

