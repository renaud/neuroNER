import unittest

import oboparser
from config import cfg

class TestParser(unittest.TestCase):

    def test_parse(self):
        obo_file = cfg['onto_root'] + 'hbp_layer_ontology.robo'
        obo = list(oboparser.parse(obo_file))
        self.assertTrue(len(obo) > 10)
        #for o in obo:
        #    print o
        print(obo[0])
        #import pdb; pdb.set_trace()


if __name__ == '__main__':
    unittest.main()
