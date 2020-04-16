import pickle
import argparse
import numpy as np
import sys


CHAIN_ID_DIR = 'chain_ids/'
MAPPINGS_DIR = 'mappings/'
MD_DIR = 'md/'

def get_species_md_dict(reference, query, transcript_gene_dict):
    codingSeqFile = MAPPINGS_DIR + reference + '.' + query + \
            '.codingMutations.out'
    
    transcriptLengths = {}
    transcriptDelRate = {}
    transcriptMutRate = {}
    transcriptChainAssigned = {}
    
    fH = open(MAPPINGS_DIR + reference + '.referenceGenes.out', 'r')
    lines = fH.readlines()
    for line in lines:
        words = line.split()
        transcriptLengths[words[0]] = len(words[1])
        transcriptChainAssigned[words[0]] = 0 
    fH.close()

    if 'ENST00000618881' in transcriptLengths.keys():
        print len(transcriptLengths.keys())
    
    chain_ids = pickle.load(open(CHAIN_ID_DIR + query + '.p', 'rb'))
    for k in chain_ids.keys():
        transcriptChainAssigned[k] = 1
    
    fH = open(codingSeqFile, 'r')
    lines = fH.readlines()
    fH.close()
    
    for line in lines:
        words = line.split()
        transcriptName = words[1]
        transcriptLen = transcriptLengths[transcriptName]
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
        transcriptMutRate[transcriptName] = (1.0*numNonSyn)/(transcriptLen \
                                                             - numDel + 0.1)
        transcriptDelRate[transcriptName] = (1.0*numDel10) / transcriptLen
    
    
    a = []
    a_chain = []
    transcriptNames = transcriptLengths.keys()
    
    
    for transcriptName in transcriptNames:
        a.append ([transcriptMutRate[transcriptName], \
                   transcriptDelRate[transcriptName]])
        if (transcriptChainAssigned[transcriptName] == 1):
            a_chain.append ([transcriptMutRate[transcriptName], \
                             transcriptDelRate[transcriptName]])
    
    
    x = np.array(a)
    x_chain = np.array(a_chain)
    y = np.mean(x_chain, axis=0)
    C = np.cov(x_chain, rowvar=False)
    C_inv = np.linalg.inv(C)
    Md = np.array([(z - y).dot(C_inv).dot(z - y) for z in x])
    
    
    species_dict = {}
    for k, transcriptName in enumerate(transcriptNames):
        if (transcriptName in transcript_gene_dict.keys()):
            geneName = transcript_gene_dict[transcriptName]
            l = species_dict.get(geneName, [])
            l.append ([transcriptDelRate[transcriptName], \
                       transcriptMutRate[transcriptName], \
                       Md[k], transcriptChainAssigned[transcriptName]])
            species_dict[geneName] = l
    
    return species_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate Mahalanobis Distance'
                                     ' for a species list.')
    parser.add_argument("-reference", type=str,
                        help="reference assembly name")
    parser.add_argument("-speciesList", type=str,
                        help="file containing list of query species")
    parser.add_argument("-geneTranscriptIds", type=str,
            help="reference gene set complete trancript ids")    

    if len(sys.argv) <= 6:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())
    r = args['reference']
    qF = args['speciesList']
    
    gene_transcript_filename = args['geneTranscriptIds']
    
    transcript_gene_dict = {}
    with open(gene_transcript_filename, 'r') as f:
        for line in f:
            words = line.split()
            transcript_gene_dict[words[1]] = words[0]
    print len(transcript_gene_dict.keys())
    
    
    md_dict = {}
    with open(qF, 'r') as f:
        for line in f:
            query = line.rstrip()
            d = get_species_md_dict (r, query, transcript_gene_dict)
            md_dict[query] = d
    pickle.dump(md_dict, open(MD_DIR+'md_dict.p', 'w'))



