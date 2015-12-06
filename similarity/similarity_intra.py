'''
Computes the similarity between neuron properties within a property class.
E.g. compares 'layer 1a' with 'layer 1-2' and 'Hypothalamus' with 'Midbrain'
Implemented in XXXXSimilarity classes, e.g. LayerSimilarity, BrainR
'''
import os.path
import oboparser
from config import cfg
from itertools import chain
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
# requires pip'ing: pip install https://github.com/AllenInstitute/AllenSDK/archive/v0.10.1.tar.gz

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
        self.onto_id2parent_regions = {}

        # only find parents up to these major brain regions
        big_reg_acros = ['Isocortex', 'OLF', 'STR', 'PAL', 'TH', 'HY', 'MB', 'PAL', 'MY', 'CB', 'HPF', 'CTXsp']

        mcc = MouseConnectivityCache()
        onto = mcc.get_ontology()
        df = onto.df
        struct_ids = df.id
        for struct_id in struct_ids:
            structure_path_str = onto[struct_id].structure_id_path.item()
            region_path = structure_path_str.split('/')
            parent_list = []
            for r in reversed(region_path):
                if r:
                    parent_list.append(r)
                    tdf = df[df.id == int(r)]
                    acronym = tdf.acronym.item()
                    if acronym in big_reg_acros:
                        break
            self.onto_id2parent_regions[struct_id] = parent_list

    def onto_id2region_id(self, o):
        return int(o[11:])

    # input: n1 and n2, each a list of ontology ids, e.g. ['HBP_EPHYS:0000080', 'HBP_LAYER:0000001', 'NCBI_GENE:19293']
    # returns a tuple with:
    #        1) a score (float) reflecting how similar these neuron are regarding brain regions
    #        2) an array of explanations, each a tuple of the form (['matching_ontology_ids'], 'human readable explanation')
    def similarity(self, n1, n2):
        def neuron2region_ids(neuron):
            region_ids = []
            for n in neuron:
                if n.startswith('ABA_REGION:'):
                    region_id = self.onto_id2region_id(n)
                    region_ids.append(self.onto_id2parent_regions[region_id])
            return list(chain(*region_ids))
        n1_parent_regions = neuron2region_ids(n1)
        n2_parent_regions = neuron2region_ids(n2)
        # print n1_parent_regions
        # print n2_parent_regions

        common_regions = [item for item in n1_parent_regions if item in n2_parent_regions]

        if common_regions:
            if n1_parent_regions[0] == n2_parent_regions[0]:
                return 1.0, (['ABA_REGION:{}'.format(n1_parent_regions[0])], 'exact same brain region')
            elif n1_parent_regions[0] in n2_parent_regions:
                return 1.0, (['ABA_REGION:{}'.format(n1_parent_regions[0])], 'sharing a common brain region')
            elif n2_parent_regions[0] in n1_parent_regions:
                return 1.0, (['ABA_REGION:{}'.format(n2_parent_regions[0])], 'sharing a common brain region')
            elif len(common_regions) > 0:
                return .5, (['ABA_REGION:{}'.format(common_regions.pop())], 'sibling regions')
        else:
            return (0, []) # no regions in both neurons


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

