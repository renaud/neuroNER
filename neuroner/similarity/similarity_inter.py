'''
Computes the similarity between neuron properties of different property classes.
E.g. compares neurotransmitters (GABA) with proteins (GAD67, and ephys (FS) with proteins (PV)

Similarity rules inter (between) neuron property classes are explicited in the following format:

(['onto_id1', 'onto_id2', ...],     #0: expectations for first neuron
 ['onto_id3'],                      #1: expectations for second neuron
 is_bidirectional,                  #2: can the rule be applied in both directions?
 probability,                       #3: [0 - 1.0]
 'citation',                        #4: a paper supporting this rule
 'explanation'                      #5: explanation or citation snippet
)
'''

_rules = [

# FS <--> PV
(['HBP_EPHYS:0000080', 'HBP_EPHYS_TRIGGER:0000003'], ['NCBI_GENE:19293'], True, 0.9, 'TODO', 'implicit correspondance'),

# Chandelier are FS
(['HBP_MORPHOLOGY:0000007'], ['HBP_EPHYS:0000080', 'HBP_EPHYS_TRIGGER:0000003'], False, 0.9, 'TODO', 'implicit correspondance'),

]


def _similarity_inter(n1, n2):
    (sim_inter, explanations) = (0, [])
    for rule in _rules:
        #print 'rule[0]', rule[0]
        #print 'rule[1]', rule [1]
        sim_this_rule = _sim_rule(rule, n1, n2)
        if rule[2]: # is_bidirectional -> apply rule the other direction
            sim_this_rule += _sim_rule(rule, n2, n1)
        if sim_this_rule > 0:
            sim_inter += sim_this_rule
            explain_text = '{} ({})'.format(rule[5], rule[4])
            explanations.append( (rule[0] + rule[1], explain_text))

    # cap similarity at max 1
    return (min(1, sim_inter), explanations)


# check whether n1 and n2 satisfy that rule
def _sim_rule(rule, n1, n2):
    #print '_sim_rule', rule, n1, n2
    for expectation in rule[0]:
        if expectation not in n1:
            #print('expectation', expectation, 'not in', n1)
            return 0 # n1 does not satisfy expectations
    for expectation in rule[1]:
        if expectation not in n2:
            #print('expectation', expectation, 'not in', n2)
            return 0 # n1 does not satisfy expectations
    return rule[3] # rule probability



