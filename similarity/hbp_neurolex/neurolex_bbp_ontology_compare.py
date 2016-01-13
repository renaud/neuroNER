
# coding: utf-8

# In[1]:

get_ipython().magic(u'load_ext autoreload')
get_ipython().magic(u'autoreload 2')


# # Comparing NeuroLex and BBP ontologies
# 

# In[2]:

from similarity import _cleanup, _normalize,   similarity2

from sherlok import Sherlok
neuroner = Sherlok('neuroner')


# In[3]:

# PARSE OBO
import oboparser, re

hbp_obo_file = 'hbp_neurolex/hbp_cell_ontology.obo'
nlex_obo_file = 'hbp_neurolex/neurolex.obo'

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
def preprocess(cell_names):
    cell_names_processed = {}
    for id, entity in cell_names.items():
        name, synonyms = entity
        print('name::', name, "id:", id)
        variants = []
        variants.append(_cleanup(neuroner.annotate(name).annotations)) # name
        for s in synonyms:
            print('         syn::',s)
            variants.append(_cleanup(neuroner.annotate(name).annotations)) # synonyms
        cell_names_processed[id] = variants
    return cell_names_processed

hbp_cell_names_processed  = preprocess(hbp_cell_names)
nlex_cell_names_processed = preprocess(nlex_cell_names)
print('done :-)')


# In[5]:

import pickle
pickle.dump(hbp_cell_names_processed,  open('hbp_cell_names_processed.pckl',  'wb'))
pickle.dump(nlex_cell_names_processed, open('nlex_cell_names_processed.pckl', 'wb'))


# In[13]:

print hbp_cell_names_processed['HBP_CELL:0000064']


# In[21]:

#print nlex_cell_names['nlx_cell_1006021']
#print nlex_cell_names_processed['nlx_cell_1006021']


# In[43]:

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
                sim = similarity2(hbp_variant, variant[0], use_inter_similarity=True)
                if sim[0] > 0:
                    hits.append((nlex_id,) + sim)
        hbp_hits[hbp_id] = hits
print('hits', len(hbp_hits))
#hit: (nlex_id, sim_score, explanations...)


# In[56]:

#hbp_hits['HBP_CELL:0000028']


# In[55]:

def get_key(item):
    return item[0]

for hbp_id, hits in hbp_hits.items():
    if len(hits) > 0:
        hits_sorted = sorted(hits, key=get_key)
        
        print('{} ({})'.format(hbp_cell_names[hbp_id][0], hbp_id))
        already_printed = []
        for nlex_id, score, explain in hits_sorted[:5]:
            if nlex_id not in already_printed:
                already_printed.append(nlex_id)
                print '* {} ({}, {})'.format(nlex_cell_names[nlex_id][0], round(score, 2), nlex_id)
        print ''


# In[ ]:



