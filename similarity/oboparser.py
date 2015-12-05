#!/usr/bin/env python

"""oboparse.py: Some functions to read the OBO format
(http://www.geneontology.org/GO.format.obo-1_2.shtml) and load it
(lazily) in a sequence of Python dicts"""

__author__ = "Stelios Sfakianakis"

from collections import defaultdict
import fileinput
import json


def parse_lines(lines):
    for line in lines:
        i = line.find('!')
        if i >= 0:
            line = line[0:i]
        line = line.strip()
        if line:
            yield line


def read_stanzas(lines):
    def unescape(str):
        return str.strip()
    single_valued = set(['id', 'name', 'is_obsolete', 'obsolete_since',
                        'ontology', 'format-version', 'data-version', 'date',
                        'comment', 'def', 'created_in'])
    stanza = defaultdict(list)
    stanza['@type'] = 'Header'
    for line in lines:
        if line.startswith('['):
            yield stanza
            stanza = defaultdict(list)
            stanza['@type'] = line[1:-1]
        else:
            k, v = line.split(':', 1)
            key = unescape(k)
            val = unescape(v)
            if key in single_valued:
                stanza[key] = val
            else:
                stanza[key].append(val)
    yield stanza

def parse(f):
    return read_stanzas(parse_lines(open(f, 'r')))

if __name__ == '__main__':
    sz = read_stanzas(parse_lines(fileinput.input()))
    for s in sz:
        print json.dumps(s)
