'''
Computes the similarity between neuron properties within a property class.
E.g. compares 'layer 1a' with 'layer 1-2' and 'Hypothalamus' with 'Midbrain'
Implemented in XXXXSimilarity classes, e.g. LayerSimilarity, BrainR
'''
import os.path
import oboparser
from config import cfg
from itertools import chain

onto_root = cfg['onto_root']
assert os.path.exists(onto_root), 'could not find onto_root at ' + onto_root


class LayerSimilarity(object):
    """LayerSimilarity: similarity """
    PREFIX = 'HBP_LAYER'
    def __init__(self):
        self.onto_id2layer_numbers = {}
        # 'HBP_LAYER:0000123' -> 123
        def onto_id2layer_number(o):
            return int(o[10:])
        for o in oboparser.parse(onto_root + 'hbp_layer_ontology.robo'):
            if 'id' in o:
                onto_id = o['id'] # HBP_LAYER:0000001
                assert len(onto_id) == 17, 'invalid onto_id: {} in {}'.format(onto_id, o)
                if onto_id2layer_number(onto_id) < 8: # HBP_LAYER:0000001 to 7
                    # just reslove to layer id
                    self.onto_id2layer_numbers[onto_id] = [onto_id2layer_number(onto_id)]
                elif 'union_of' in o:
                    # put the layer numbers of all 'union_of' ids
                    self.onto_id2layer_numbers[onto_id] = [onto_id2layer_number(u) for u in o['union_of']]
                elif 'is_a' in o: # e.g. L5a is_a L5
                    # to simplify, resolve to L5
                    self.onto_id2layer_numbers[onto_id] = [onto_id2layer_number(o['is_a'][0])]
                else:
                    raise Exception('invalid layer entry: {}'.format(o))

    '''
    returns 1 if both neuron share a similar layer, e.g. 'layer2/3' and 'layer 2A', else 0
    '''
    def similarity(self, n1, n2):
        # HBP_LAYER:0000101 (layer 1-2) --> [1,2]
        def neuron2layer_numbers(neuron):
            layer_numbers = []
            for n in neuron:
                if n.startswith('HBP_LAYER:'):
                    for layer_number in self.onto_id2layer_numbers[n]:
                        layer_numbers.append(layer_number)
            return layer_numbers
        n1_layer_numbers = neuron2layer_numbers(n1)
        n2_layer_numbers = neuron2layer_numbers(n2)
        # on same layer?
        common_layer = set(n1_layer_numbers).intersection(set(n2_layer_numbers))
        if len(common_layer) > 0:
            return (1.0, (['HBP_LAYER:000000{}'.format(common_layer.pop())], 'located on same Layer') )
        else:
            return (0, []) # no layers in both neurons




class BrainRegionSimilarity(object):
    """BrainRegionSimilarity: similarity """
    def __init__(self):
        # forSJ: the place to load ABA ontology
        pass

    # input: n1 and n2, each a list of ontology ids, e.g. ['HBP_EPHYS:0000080', 'HBP_LAYER:0000001', 'NCBI_GENE:19293']
    # returns a tuple with:
    #        1) a score (float) reflecting how similar these neuron are regarding brain regions
    #        2) an array of explanations, each a tuple of the form (['matching_ontology_ids'], 'human readable explanation')
    #        forSJ: just return [] for 2) if it's confusing...
    def similarity(self, n1, n2):
        return (0, []) # forSJ: placeholder: no similarity


similarities = [LayerSimilarity(), BrainRegionSimilarity()]

def _similarity_intra(n1, n2, weights):
    # perfect similarity/equality
    if n1 == n2:
        return (1,[])
    else:
        # dispatch to each similarities; aggregate score & explanations
        (sim_intra, explanations) = (0, [])
        for s in similarities:
            (s_sim, s_expl) = s.similarity(n1, n2) #TODO weights
            sim_intra += s_sim; explanations.append(s_expl)
        return (sim_intra, explanations)

