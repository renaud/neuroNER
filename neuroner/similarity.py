'''
Computes the semantic similarity between two neuron mentions.
The main function is similarity()
Computation is delegated to similarity_inter and similarity_intra
'''

from sherlok import Sherlok # pip install --upgrade sherlok

import similarity_inter
import similarity_intra
from normalize import clean_annotations

s = Sherlok('neuroner')

WEIGHTS = { #TODO: implement weights
            #TODO: merge with BASE_MULTIPLIER implemented in similarity_intra
    'Layer': 1.0,
    'ProteinProp': 1.0
}

'''
Computes the intra and inter semantic similarity between two neurons
in: n1@str, n2@str: the two neurons to measure similarity
out: (score:float, [([matching_properties], explanation@str)])
'''
def similarity(n1, n2, weights=WEIGHTS, symmetric=True, use_inter_similarity=True):
    assert type(n1) is str and len(n1) > 0, 'n1 cannot be empty'
    assert type(n2) is str and len(n2) > 0, 'n2 cannot be empty'
    n1c = clean_annotations(s.annotate(n1).annotations)
    n2c = clean_annotations(s.annotate(n2).annotations)
    return similarity2(n1c, n2c, weights, symmetric, use_inter_similarity)





'''
Computes the intra and inter semantic similarity between two neurons
in:  n1@[], n2@[]: the two neurons to measure similarity
out: (score:float, [([matching_properties, explanation@str])])
'''
def similarity2(n1, n2, weights=WEIGHTS, symmetric=True, use_inter_similarity=True):
    # TODO assert type(n1) is array
    s_intra = similarity_intra._similarity_intra(n1, n2, weights, symmetric)
    if not use_inter_similarity:
        return s_intra
    else:
        s_inter = similarity_inter._similarity_inter(n1, n2)
        #print('s_intra', s_intra, 's_inter', s_inter)
        return (s_intra[0] + s_inter[0], s_intra[1] + s_inter[1])

