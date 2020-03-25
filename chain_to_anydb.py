import anydbm
import argparse
import os
import sys
import pickle
import gzip

CHAIN_DIR = 'chains/' 
CHAIN_ID_DIR = 'chain_ids/'

def chainHeaderToDb(chain_ids, chainIn, dbOut):
    db = anydbm.open(dbOut, 'n')
    with gzip.open(chainIn, 'r') as f:
        for line in f:
            if (line[0:5] == 'chain'):
                words = line.split()
                chain_id = int(words[12])
                if words[12] in chain_ids:
                    #print chain_id, line
                    db[str(chain_id)] = line
    db.close()

def chainLinkToDb(chain_ids, chainIn, dbOut):
    db = anydbm.open(dbOut, 'n')
    chain_id = 0
    link_id = 0
    dq = 0
    dr = 0
    match = 0
    rpos = 0
    qpos = 0
    frame = 1
    strand = '+'
    in_chain = False
    print_chain = False
    with gzip.open(chainIn, 'r') as f:
        for line in f:
            words = line.split()
            if (line[0:5] == 'chain'):
                chain_id = int(words[12])
                rpos = int(words[5])
                strand = words[9]
                if (strand == '+'):
                    qpos = int(words[10])
                else:
                    qpos = int(words[8]) - int(words[10])
                frame = 0
                link_id = 0
                in_chain =True
                if words[12] in chain_ids:
                    print_chain = True
                else:
                    print_chain = False

            elif (in_chain):
                if (len(words) > 0):
                    if (len(words) == 3):
                        match, dr, dq = [int(c) for c in line.split()]
                    else:
                        match, dr, dq = [int(words[0]), 0, 0]
                    key = str(chain_id)+'_'+str(link_id)
                    if (strand == '+'):
                        link = ' '.join([str(rpos), str(qpos), str(match),
                                         str(frame), str(chain_id),
                                         str(link_id)]) 
                        qpos += match + dq
                    else:
                        qpos -= match
                        link = ' '.join([str(rpos), str(qpos), str(match),
                                         str(frame), str(chain_id),
                                         str(link_id)]) 
                        qpos -= dq
                    
                    rpos += match + dr
                    frame = (frame + dq + 3 - (dr % 3)) % 3
                    if print_chain:
                        #print key, link
                        db[key] = link 
                    link_id += 1
    db.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts each chain to db '
            'via anydbm.')
    parser.add_argument("-reference", type=str,
            help="reference assembly name")
    parser.add_argument("-query", type=str,
            help="query assembly name")
    
    if len(sys.argv) <= 4:
        parser.print_help()
        sys.exit(1)
    
    args = vars(parser.parse_args())
    reference = args['reference']
    query = args['query']

    chainIn = CHAIN_DIR + reference + '.' + query + '.all.chain.gz'
    dbOut = CHAIN_DIR + reference + '.' + query + '.all.chain.db'
    linkDbOut = CHAIN_DIR + reference + '.' + query + '.all.chain.link.db'
    
    chain_ids = set([])
    chain_id_dict = pickle.load(open(CHAIN_ID_DIR + query + '.p')) 
    for k,v in chain_id_dict.items():
        chain_ids.add(v[0])
    
    chainLinkToDb(chain_ids, chainIn, linkDbOut)
    chainHeaderToDb(chain_ids, chainIn, dbOut)

