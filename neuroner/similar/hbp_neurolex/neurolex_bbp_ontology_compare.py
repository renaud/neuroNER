
# coding: utf-8

# In[1]:

get_ipython().magic(u'load_ext autoreload')
get_ipython().magic(u'autoreload 2')


# # Comparing NeuroLex and BBP ontologies
# 

# In[2]:

import sys
sys.path.append('../')
from neuroner.neuroner import clean_annotations, similarity2

from sherlok import Sherlok
neuroner = Sherlok('neuroner')


# In[3]:

# PARSE OBO
import oboparser, re

hbp_obo_file  = 'hbp_cell_ontology.obo'
nlex_obo_file = 'neurolex.obo'

# a simple function to pull out the cell names
SYNONOYM_NAME = re.compile(r'"(.*?)"').search
def get_cell_names(obo_onto_file):
    cell_names = {}    
    for stanza in oboparser.parse(obo_onto_file):
        if stanza['@type'] != "Term": continue
        id = stanza["id"]
        if id == '': continue
        name = stanza["name"]
        synonyms = []
        for s in stanza["synonym"]:
            synonyms.append(SYNONOYM_NAME(s).group(1))
        cell_names[id] = (name, synonyms)
    return cell_names

hbp_cell_names  = get_cell_names(hbp_obo_file)
nlex_cell_names = get_cell_names(nlex_obo_file)
print 'hbp_cell_names', len(hbp_cell_names)
print 'nlex_cell_names', len(nlex_cell_names)
assert hbp_cell_names['HBP_CELL:0000016'][0] == 'neuron expressing Parvalbumin'
assert hbp_cell_names['HBP_CELL:0000016'][1][1] == 'Pvalb+ cell'


# In[ ]:




# In[4]:

# PREPROCESS NEUROLEX and HBP
def preprocess(cell_names, prefix=''):
    cell_names_processed = {}
    for id, entity in cell_names.items():
        name, synonyms = entity
        print('name::', name, "id:", id, "prefix:", prefix)
        variants = []
        variants.append(clean_annotations(neuroner.annotate(prefix + name).annotations)) # name
        for s in synonyms:
            print('         syn::', s)
            variants.append(clean_annotations(neuroner.annotate(prefix + s).annotations)) # synonyms
        cell_names_processed[id] = variants
    return cell_names_processed

hbp_cell_names_processed  = preprocess(hbp_cell_names, prefix='neocortex ')
nlex_cell_names_processed = preprocess(nlex_cell_names)
print('done :-)')


# In[5]:

import pickle
pickle.dump(hbp_cell_names_processed,  open('hbp_cell_names_processed.pckl',  'wb'))
pickle.dump(nlex_cell_names_processed, open('nlex_cell_names_processed.pckl', 'wb'))


# In[6]:

print hbp_cell_names['HBP_CELL:0000064']
print hbp_cell_names_processed['HBP_CELL:0000064']


# In[17]:

print nlex_cell_names['nifext_56']
print nlex_cell_names_processed['nifext_56']


# In[18]:

similarity2(hbp_cell_names_processed['HBP_CELL:0000064'][0], nlex_cell_names_processed['nifext_56'][0])


# In[20]:

#%debug
hbp_hits={}
for hbp_id, hbp_cell_variants in hbp_cell_names_processed.items():
    hbp_variant = hbp_cell_variants[0] # just the 1st one TODO
    if len(hbp_variant) == 0:
        print('NO VARIANT FOR ', hbp_id, hbp_cell_names[hbp_id])
    else:
        hits = []
        for nlex_id, nlex_cell_variants in nlex_cell_names_processed.items():
            for variant in nlex_cell_variants:
                if len(variant) == 0: continue
                #print("sim: HBP::",hbp_variant, "NLEX::",variant, "nlex_id",nlex_id)
                sim = similarity2(hbp_variant, variant, use_inter_similarity=True)
                if sim[0] > 0:
                    #print('hit for',hbp_variant, 'WITH',variant[0])
                    hits.append((nlex_id,) + sim)
        hbp_hits[hbp_id] = hits
print('hits', len(hbp_hits))
#hit: (nlex_id, sim_score, explanations...)


# In[21]:

hbp_hits['HBP_CELL:0000064']


# In[60]:

#%debug

import collections
def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def get_key(item):
    return item[1]

with open('hbp-neurolex.tsv', 'w') as outf:
    outf.write('BBP label	BBP rdf	Neurolex name	%	Neurolex ID\n')
    for hbp_id, hits in hbp_hits.items():
        if len(hits) > 0:
            hits_sorted = sorted(hits, key=get_key, reverse=True)

            #print('{} ({})'.format(hbp_cell_names[hbp_id][0], hbp_id))
            already_printed = []
            for nlex_id, score, explain in hits_sorted[:5]:
                
                # only if more than Neocortex
                explain_str = ' '.join(flatten(explain))
                
                if nlex_id not in already_printed and explain_str != 'ABA_REGION:315 exact same brain region':
                    already_printed.append(nlex_id)
                    #print explain_str
                    #print '* {} ({}, {})'.format(nlex_cell_names[nlex_id][0], round(score, 2), nlex_id)
                    outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(hbp_cell_names[hbp_id][0], hbp_id, nlex_cell_names[nlex_id][0], round(score, 2), nlex_id, explain_str))


# In[45]:

cell_= 'Nest Basket Cell' #'nest basket cell'
clean_annotations(neuroner.annotate(cell_).annotations)


# In[47]:

sorted(hbp_hits['HBP_CELL:0000061'], key=get_key, reverse=True)[:5]


# In[ ]:



