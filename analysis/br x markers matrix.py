
# coding: utf-8

# In[10]:

'''
Computes matrix brain_region x marker_genes
'''
import ijson
from collections import defaultdict

matrix_dict = defaultdict(dict)

#f = open('neuroner_20160122s_index_sample2.json')
f = open('/Users/richarde/data/neuroner/neuroner_20160122s_index.json') # data in https://dl.dropboxusercontent.com/u/975350/tmp/neuroner_20160122s_index.json.tgz
docs = ijson.items(f, 'item')
cnt = 0
for d in docs:
    pmid = d['_id']
    for neuron in d['_source']['neuron']:
        #print(neuron)
        bregions = [nprop['onto_id'] for nprop in neuron['neuron_properties'] if nprop['neuron_type'] == 'brainregion']
        genes =    [nprop['onto_id'] for nprop in neuron['neuron_properties'] if nprop['neuron_type'] == 'protein']

        #print 'bregions', bregions, 'genes', genes
        for bn in bregions:
            for g in genes:
                cnt += 1
                #print bn, g
                try:
                    matrix_dict[bn][g] += 1
                except Exception, e:
                    matrix_dict[bn][g] = 1
                    
print cnt, 'co-occurences'


# In[13]:

import pandas as pd

mat = pd.DataFrame(matrix_dict)

writer = pd.ExcelWriter('matrix_dict.xlsx')
mat.to_excel(writer, 'sheet1')


# In[ ]:




# In[ ]:



