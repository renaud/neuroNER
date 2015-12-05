'''
Computes the semantic similarity between two neuron mentions.
The main function is similarity()
Computation is delegated to similarity_inter and similarity_intra
'''

import similarity_inter, similarity_intra

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
def _cleanup(n):
    clean = []
    filt_attrib_list = ['Neuron', 'PreNeuron', 'PostNeuron', 'Electrophysiology', 'ProteinTrigger', 'NeuronTrigger']
    for (begin, end, txt, type_, props) in n:
        if 'ontologyId' in props:
            clean.append(props[u'ontologyId'])
        elif type_ not in filt_attrib_list:
            clean.append('{}:{}'.format(type_, txt))
    return clean


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

