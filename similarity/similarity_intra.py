'''
Computes the similarity between neuron properties within a property class.
E.g. compares 'layer 1a' with 'layer 1-2' and 'Hypothalamus' with 'Midbrain'
Implemented in XXXXSimilarity classes, e.g. LayerSimilarity, BrainR
'''
import os.path
import oboparser
from config import cfg
from itertools import chain
import glob
import numpy as np

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
                    # not happy with this - SJT
                    self.onto_id2layer_numbers[onto_id] = [onto_id2layer_number(onto_id), onto_id2layer_number(o['is_a'][0])]
                else:
                    raise Exception('invalid layer entry: {}'.format(o))

    '''
    returns .5 if both neurons share a similar layer, e.g. 'layer2/3' and 'layer 2A', else 0
    '''
    def similarity(self, n1, n2):
        BASE_MULTIPLIER = 2
        # HBP_LAYER:0000101 (layer 1-2) --> [1,2]
        def neuron2layer_numbers(neuron):
            layer_numbers = []
            for n in neuron:
                if n.startswith(self.PREFIX):
                    for layer_number in self.onto_id2layer_numbers[n]:
                        layer_numbers.append(layer_number)
            return layer_numbers
        n1_layer_numbers = neuron2layer_numbers(n1)
        n2_layer_numbers = neuron2layer_numbers(n2)
        sum_layer_count = len(n1_layer_numbers) + len(n2_layer_numbers) + 1

        # on same layer?
        common_layers = list(set(n1_layer_numbers).intersection(set(n2_layer_numbers)))
        match_percent = float(len(common_layers) + 1)/sum_layer_count * 2
        for i,c in enumerate(common_layers):
            common_layers[i] = 'HBP_LAYER:000000%d' % c
        if len(common_layers) > 0:
            #return (1.0, (['HBP_LAYER:000000{}'.format(common_layer.pop())], 'located on same layer') )
            return (match_percent * BASE_MULTIPLIER, (common_layers, 'shares layers') )
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
        BASE_MULTIPLIER = 1.5
        def neuron2region_ids(neuron):
            region_ids = []
            parent_region_ids = []
            for n in neuron:
                if n.startswith('ABA_REGION:'):
                    region_id = self.onto_id2region_id(n)
                    region_ids.append(region_id)
            for region_id in region_ids:
                parent_region_ids.append(self.onto_id2parent_regions[region_id])
            return list(chain(*parent_region_ids))

        #def assess_sub_region(r1, r2):
        n1_parent_regions = neuron2region_ids(n1)
        n2_parent_regions = neuron2region_ids(n2)
        # print n1_parent_regions
        # print n2_parent_regions

        common_regions = [item for item in n1_parent_regions if item in n2_parent_regions]

        if common_regions:
            if n1_parent_regions[0] == n2_parent_regions[0]:
                return 1.0 * BASE_MULTIPLIER, (['ABA_REGION:{}'.format(n1_parent_regions[0])], 'exact same brain region')
            elif n1_parent_regions[0] in n2_parent_regions:
                return .875 * BASE_MULTIPLIER, (['ABA_REGION:{}'.format(n1_parent_regions[0])], 'sharing a common brain region')
            elif n2_parent_regions[0] in n1_parent_regions:
                return .875 * BASE_MULTIPLIER, (['ABA_REGION:{}'.format(n2_parent_regions[0])], 'sharing a common brain region')
            elif len(common_regions) > 0:
                return .75 * BASE_MULTIPLIER, (['ABA_REGION:{}'.format(common_regions.pop())], 'sibling regions')
        else:
            return (0, []) # no regions in both neurons



class MorphologySimilarity(object):
    """LayerSimilarity: similarity """
    PREFIX = 'HBP_MORPHOLOGY'

    '''
    returns number of shared terms if both neurons share same morphology, e.g. 'Pyr' and 'pyramidal', else 0
    '''
    def similarity(self, n1, n2):
        def neuron2morphology(neuron):
            morphologies = []
            for n in neuron:
                if n.startswith(self.PREFIX):
                    morphologies.append(n)
            return morphologies

        n1_morphologies = neuron2morphology(n1)
        n2_morphologies = neuron2morphology(n2)

        common_morphologies = list(set(n1_morphologies).intersection(set(n2_morphologies)))
        if len(common_morphologies) > 0:
            return (len(list(common_morphologies)), (common_morphologies, 'shares morphology') )
        else:
            return (0, []) # no morphologies common to both neurons

class MouseLineSimilarity(object):
    """LayerSimilarity: similarity """
    PREFIX = 'MOUSE_LINE'

    '''
    returns number of shared terms if both neurons share same mouse line, else 0
    '''
    def similarity(self, n1, n2):
        BASE_MULTIPLIER = 2
        def neuron2line(neuron):
            morphologies = []
            for n in neuron:
                if n.startswith(self.PREFIX):
                    morphologies.append(n)
            return morphologies

        n1_lines = neuron2line(n1)
        n2_lines = neuron2line(n2)

        common_lines = list(set(n1_lines).intersection(set(n2_lines)))
        if len(common_lines) > 0:
            return (BASE_MULTIPLIER * len(list(common_lines)), (common_lines, 'shares mouse_line') )
        else:
            return (0, []) # no morphologies common to both neurons

class ProteinSimilarity(object):
    """LayerSimilarity: similarity """
    PREFIX = 'NCBI_GENE'

    '''
    returns number of shared terms, else 0
    '''
    def similarity(self, n1, n2):
        def neuron2terms(neuron):
            matching_terms = []
            for n in neuron:
                if n.startswith(self.PREFIX):
                    matching_terms.append(n)
            return matching_terms

        n1_terms = neuron2terms(n1)
        n2_terms = neuron2terms(n2)

        common_terms = list(set(n1_terms).intersection(set(n2_terms)))
        if len(common_terms) > 0:
            return (len(common_terms), (common_terms, 'shares proteins') )
        else:
            return (0, []) # no morphologies common to both neurons

class NeurotransmitterSimilarity(object):
    """LayerSimilarity: similarity """
    PREFIX = 'HBP_NEUROTRANSMITTER'

    '''
    returns number of shared terms, else 0
    '''
    def similarity(self, n1, n2):
        def neuron2terms(neuron):
            matching_terms = []
            for n in neuron:
                if n.startswith(self.PREFIX):
                    matching_terms.append(n)
            return matching_terms

        n1_terms = neuron2terms(n1)
        n2_terms = neuron2terms(n2)

        common_terms = list(set(n1_terms).intersection(set(n2_terms)))
        if len(common_terms) > 0:
            return (len(common_terms), (common_terms, 'shares neurotransmitters') )
        else:
            return (0, []) # no morphologies common to both neurons

class ProjectionSimilarity(object):
    """LayerSimilarity: similarity """
    PREFIX = 'HBP_PROJECTION'

    '''
    returns number of shared terms, else 0
    '''
    def similarity(self, n1, n2):
        BASE_MULTIPLIER = 4
        def neuron2terms(neuron):
            matching_terms = []
            for n in neuron:
                if n.startswith(self.PREFIX):
                    matching_terms.append(n)
            return matching_terms

        n1_terms = neuron2terms(n1)
        n2_terms = neuron2terms(n2)

        common_terms = list(set(n1_terms).intersection(set(n2_terms)))
        if len(common_terms) > 0:
            return (len(common_terms) * BASE_MULTIPLIER, (common_terms, 'shares projection patterns') )
        else:
            return (0, []) # no regions common to both neurons

class UnknRegionSimilarity(object):
    """LayerSimilarity: similarity """
    PREFIX = 'UNKN_REGION'

    '''
    returns number of shared terms, else 0
    '''
    def similarity(self, n1, n2):
        BASE_MULTIPLIER = .5
        def neuron2terms(neuron):
            matching_terms = []
            for n in neuron:
                if n.startswith(self.PREFIX):
                    matching_terms.append(n)
            return matching_terms

        n1_terms = neuron2terms(n1)
        n2_terms = neuron2terms(n2)

        common_terms = list(set(n1_terms).intersection(set(n2_terms)))
        if len(common_terms) > 0:
            return (len(common_terms) * BASE_MULTIPLIER, (common_terms, 'shares general regions') )
        else:
            return (0, []) # no regions common to both neurons

similarities = [LayerSimilarity(), BrainRegionSimilarity(), ProjectionSimilarity(), UnknRegionSimilarity(),
                MorphologySimilarity(), NeurotransmitterSimilarity(), ProteinSimilarity(), MouseLineSimilarity()]


def _similarity_intra(n1, n2, weights, symmetric):
    # assumes that the first entered neuron is the target neuron, and normalizes score based on highest possible
    #   for first entered neuron

    # calculate base similarity
    (sim_intra_base_first, exp) = _calc_base_similarity_intra(n1, n1, weights)
    (sim_intra_base_second, exp) = _calc_base_similarity_intra(n2, n2, weights)
    if sim_intra_base_first == 0.0:
        sim_intra_base_first = 1.0

    # perfect similarity/equality
    if n1 == n2:
        return (1, [])
    else:
        # dispatch to each similarities; aggregate score & explanations
        (sim_intra, explanations) = _calc_base_similarity_intra(n1, n2, weights)

        # normalize score based highest possible score
        sim_intra_norm = sim_intra / sim_intra_base_first

        if symmetric:
            sim_intra_norm = sim_intra / ((sim_intra_base_first + sim_intra_base_second) / 2)

        return (sim_intra_norm, explanations)

def _calc_base_similarity_intra(n1, n2, weights):
    (sim_intra, explanations) = (0.0, [])
    for s in similarities:
        (s_sim, s_expl) = s.similarity(n1, n2) #TODO weights
        sim_intra += s_sim; explanations.append(s_expl)
    return (sim_intra, explanations)


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


