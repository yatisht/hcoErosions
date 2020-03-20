import pickle
import argparse
import numpy as np
import sys


CHAIN_ID_DIR = 'chain_ids/'
MAPPINGS_DIR = 'mappings/'
MD_DIR = 'md/'

def get_species_md_dict(reference, query):
    codingSeqFile = MAPPINGS_DIR + reference + '.' + query + \
            '.codingMutations.out'
    
    geneLengths = {}
    geneDelRate = {}
    geneMutRate = {}
    geneChainAssigned = {}
    
    fH = open(MAPPINGS_DIR + reference + '.referenceGenes.out', 'r')
    lines = fH.readlines()
    for line in lines:
        words = line.split()
        geneLengths[words[0]] = len(words[1])
        geneChainAssigned[words[0]] = 0 
    fH.close()
    
    chain_ids = pickle.load(open(CHAIN_ID_DIR + query + '.p', 'rb'))
    for k in chain_ids.keys():
        geneChainAssigned[k] = 1
    
    fH = open(codingSeqFile, 'r')
    lines = fH.readlines()
    fH.close()
    
    for line in lines:
        words = line.split()
        geneName = words[1]
        geneLen = geneLengths[geneName]
        numNonSyn = sum(1 for c in words[3] if c.islower())
        numDel = 0.0
        numDel10 = 0.0
        exons = words[2].split('|')
        for exon in exons:
            entries = exon.split(',')
            num_entries = len(entries)
            for i in range(1, num_entries, 2):
                entry = entries[i]
                indels = entry.split('-')
                if (indels[1].isdigit()):
                    numDel += int(indels[1])/3.0
                if (indels[0].isdigit() and indels[1].isdigit()):
                    if (int(indels[1]) > 30):
                        numDel10 += int(indels[1])/3.0
        geneDelRate[geneName] = (1.0*numNonSyn)/(geneLen - numDel + 0.1)
        geneMutRate[geneName] = (1.0*numDel10)/geneLen
    
    
    a = []
    a_chain = []
    geneNames = geneLengths.keys()
    
    
    for geneName in geneNames:
        a.append ([geneMutRate[geneName], geneDelRate[geneName]])
        if (geneChainAssigned[geneName] == 1):
            a_chain.append ([geneMutRate[geneName], geneDelRate[geneName]])
    
    
    x = np.array(a)
    x_chain = np.array(a_chain)
    y = np.mean(x_chain, axis=0)
    C = np.cov(x_chain, rowvar=False)
    C_inv = np.linalg.inv(C)
    Md = np.array([(z - y).dot(C_inv).dot(z - y) for z in x])
    
    
    species_dict = {}
    for k, geneName in enumerate(geneNames):
        species_dict[geneName]=[[geneDelRate[geneName], geneMutRate[geneName], \
                                   Md[k], geneChainAssigned[geneName]]]
    
    return species_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate Mahalanobis Distance'
                                     ' for a species list.')
    parser.add_argument("-reference", type=str,
                        help="reference assembly name")

    parser.add_argument("-speciesList", type=str,
                        help="file containing list of query species")

    if len(sys.argv) <= 4:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())
    r = args['reference']
    qF = args['speciesList']
    
    md_dict = {}
    with open(qF, 'r') as f:
        for line in f:
            query = line.rstrip()
            d = get_species_md_dict (r, query)
            md_dict[query] = d
    pickle.dump(md_dict, open(MD_DIR+'md_dict.p', 'w'))



