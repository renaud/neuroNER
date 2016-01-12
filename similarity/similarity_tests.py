import unittest;

from similarity import similarity, _cleanup
from sherlok import Sherlok # pip install --upgrade sherlok


class TestSimilarity(unittest.TestCase):

    def test_cleanup(self):
        s = Sherlok('neuroner')
        an = s.annotate('layer 4 pyramidal long large neuron').annotations
        clean = _cleanup(an)
        self.assertEqual(clean, ['HBP_LAYER:0000004', u'HBP_MORPHOLOGY:0000001', 'Missing:long', 'Size:large'])

    def test_exact_similiarity(self):
        s = similarity('layer 4 neuron', 'layer 4 neuron', use_inter_similarity=False)
        self.assertEqual(s[0], 1.0)

    def test_inter_similiarity_PV(self):
        s = similarity('PV neuron', 'fast-spiking neuron', use_inter_similarity=True)
        self.assertEqual(s[0], 0.9, 'inter similarity works for PV and fast-spiking')
        s_reverse = similarity('fast-spiking neuron', 'PV neuron')
        self.assertEqual(s_reverse[0], 0.9, 'inter similarity works in both directions')

from similarity_intra import *

class TestLayerSimilarity(unittest.TestCase):
    ls = LayerSimilarity()

    def test_exact_similiarity(self):
        sim = self.ls.similarity(['HBP_LAYER:0000001'], ['HBP_LAYER:0000001'])
        self.assertEqual(sim[1], (['HBP_LAYER:0000001'], 'shares layers') )

    def test_no_similiarity(self):
        sim = self.ls.similarity(['HBP_LAYER:0000001'], ['HBP_LAYER:0000002'])
        self.assertEqual(sim, (0, []))

    def test_semantic_similiarity(self):
        # L1 and L4 <--> L1/2
        sim = self.ls.similarity(['HBP_LAYER:0000001', 'HBP_LAYER:0000041'], ['HBP_LAYER:0000101'])
        self.assertEqual(sim[1], (['HBP_LAYER:0000001'], 'shares layers') )

    def test_semantic_similiarity2(self):
        # layer 1-3 <--> layer 3b
        sim = self.ls.similarity(['HBP_LAYER:0000102'], ['HBP_LAYER:0000031'])
        self.assertEqual(sim[1], (['HBP_LAYER:0000003'], 'shares layers') )



class TestBrainRegionSimilarity(unittest.TestCase):
    brs = BrainRegionSimilarity()

    # pvz and mez share hyp as a common parent
    hyp = 'ABA_REGION:1097' # Hypothalamus
    pvz = 'ABA_REGION:157'  # Periventricular zone
    mez = 'ABA_REGION:467'  # Hypothalamic medial zone

    def test_no_similiarity(self):
        sim = self.brs.similarity([self.hyp], ['ABA_REGION:1024']) # grooves
        self.assertEqual(sim, (0, [])), 'these two brain regions are not similar at all, therefore a score of 0'

    def test_exact_similiarity(self):
        sim = self.brs.similarity([self.hyp], [self.hyp])
        self.assertEqual(sim[1], ([self.hyp], 'exact same brain region') )

    def test_parent_child_relationship(self):
        sim = self.brs.similarity([self.hyp], [self.mez])
        self.assertEqual(sim[1], ([self.hyp],  'sharing a common brain region') )

    def test_sharing_direct_parents(self):
        sim = self.brs.similarity([self.pvz], [self.mez])
        self.assertEqual(sim[1], ([self.hyp], 'sibling regions') )



if __name__ == '__main__':
    unittest.main()
