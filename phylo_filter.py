
from treelib import Node, Tree
import pickle
import sys
import argparse
import csv

MD_DIR = 'md/' 
REF_GENES_DIR = 'ref_genes/' 

def create_tree(tree_filename):
    tree = Tree()
    f = open(tree_filename)
    lines = f.readlines()
    first_line = lines[0]
    l = first_line.rstrip()
    s1 = l.split(',')
    s2 = [s.replace('(', '').replace(')', '') for s in s1]
    stack = [(w.count('('), w.count(')')) for w in s1]
    num_open = sum([s[0] for s in stack])
    num_close = sum([s[1] for s in stack])

    if ((num_open != num_close) or (num_open != len(s2)-1)):
        print 'ERROR: PhyloTree in incorrect format!'
        sys.exit()

    curr_node = '0'
    parent_stack = []

    for (k, species) in enumerate(s2):
        no = stack[k][0]
        nc = stack[k][1]
        for i in range(no):
            curr_node = str(int(curr_node)+1)
            if len(parent_stack) == 0:
                tree.create_node(curr_node, curr_node)
            else:
                tree.create_node(curr_node, curr_node, parent=parent_stack[-1])
            parent_stack.append(curr_node)
        tree.create_node(species, species, parent=parent_stack[-1])
        for i in range(nc):
            parent_stack.pop()

    return tree

def indicator(feature_vectors):
    intact_thresh = 5
    outlier_thresh = 15
    l = []
    all_del = True
    for feature_vec in feature_vectors:
        md = feature_vec[2]
        chain_assigned = feature_vec[3]
        non_gaps = feature_vec[1]
        if feature_vec[0] < 1:
            all_del = False
        if (chain_assigned < 1): # or (feature_vec[2] == 1):
            l.append(0)
        elif md > outlier_thresh and feature_vec[0] > 0.1:
            l.append(-1)
        elif (md > intact_thresh) or (chain_assigned == 0) or (non_gaps > 0.2): 
            l.append(0)
        else:
            l.append(1)
    # If any transcript is intact
    if sum([i > 0 for i in l]) > 0:
        return 1
    # If all transcripts are outliers
    elif (sum(l) == -len(l)) and (all_del == False):# and (chain_assigned == 1):
        return -1
    else:
        return 0

def LCA(tree, node1, node2):
    parent_list = [node1]
    p = node1
    while (p != tree.root):
        p_node = tree.parent(p)
        p = p_node.identifier
        parent_list.append(p)
    p = node2
    while (p not in parent_list) and (p != tree.root):
        p_node = tree.parent(p)
        p = p_node.identifier
    return p

def dfs(gene, md_dict, ref, tree, node, leaf_1, leaf_2, first, all_ret):
    children = tree.children(node)
    n = tree.get_node(node)
    l1 = tree.get_node(leaf_1)
    l2 = tree.get_node(leaf_2)
    intact = 0
    lost = 0
    undetermined = 0
    if n.is_leaf():
        if node == ref:
            ind = 1
        else:
            ind = get_indicator(md_dict, node, gene)
        if ind == 1:
            return (1, 0, 0, [], [], all_ret)
        elif ind == -1:
            return (0, 1, 0, [n.tag], [], all_ret)
        else:
            return (0, 0, 1, [], [n.tag], all_ret)
    lost_children = []
    lost_leaves = []
    undetermined_leaves = []
    for child in children:
        (i, l, u, ll, uu, ret) = dfs(gene, md_dict, ref, tree, \
                                     child.identifier, leaf_1, leaf_2, 0, [])
        for r in ret:
            all_ret.append(r)
        intact += i
        lost += l
        undetermined += u
        for x in ll: 
            lost_leaves.append(x)
        for x in uu: 
            undetermined_leaves.append(x)
        if (i == 0) and (l > 1) and (l >= u):
            lost_children.append([child.identifier, (i, l, u), ll + uu])
    for lost_child in lost_children:
        all_ret.append(lost_child)
    if (first):
        return all_ret 
    else:
        return (intact, lost, undetermined, lost_leaves, undetermined_leaves, \
                all_ret)

def get_indicator (md_dict, species, gene_id):
    if (gene_id in md_dict[species].keys()):
        return indicator(md_dict[species][gene_id])
    else:
        return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Apply phylogenetic filter'
                                     ' to discover high-confidence orthologous'
                                     ' gene erosions.')
    parser.add_argument("-reference", type=str,
                        help="reference assembly name")

    parser.add_argument("-speciesList", type=str,
                        help="file containing list of query species")

    parser.add_argument("-phyloTree", type=str,
                        help="file containing phylogenetic tree of species "\
                        "list in newick format")
    
    parser.add_argument("-geneTranscriptIds", type=str,
            help="reference gene set complete trancript ids")    

    parser.add_argument("-outFile", type=str,
                        help="output file containing hiconfErosion predictions"\
                        " (csv format)")

    if len(sys.argv) <= 10:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())
    ref = args['reference']
    qF = args['speciesList']
    tree_filename = args['phyloTree']
    gene_transcript_filename = args['geneTranscriptIds']
    out_filename = args['outFile']

    tree = create_tree(tree_filename)
    print 'Input tree:'
    tree.show()
    
    csvfile = open(out_filename, 'w')
    writer = csv.writer(csvfile)

    leaves = [l.identifier for l in  tree.leaves(tree.root)]
    md_dict = pickle.load(open(MD_DIR+'md_dict.p'))

    #    transcript_gene_dict = {}
    gene_names = set([])

    with open(gene_transcript_filename, 'r') \
            as f:
        for line in f:
            words = line.split()
            gene_names.add(words[0])
#            transcript_gene_dict[words[1]] = words[0]
    #print transcript_gene_dict

    print 'Finding hiconfErosion predictions.'
    print '#Gene, Intact-1, Intact-2, hiconfErosion species'
    writer.writerow(['Gene', 'Intact-1', 'Intact-2', 'hiconfErosion species (-'\
                     ' separated)'])
#    for transcript in transcript_gene_dict.keys():
    for gene_name in gene_names:
#        gene_name = transcript_gene_dict[transcript]
        intact_leaves = [k for k in leaves if (k != ref) and \
                         (get_indicator(md_dict, k, gene_name) == 1)]
        intact_leaves.append(ref)
        #print gene_name, intact_leaves
        intact_ancestors = []
        intact_ancestor_leaves = {}
        for leaf_1 in [ref]:
            for leaf_2 in intact_leaves:
                if (leaf_1 != leaf_2):
                    ancestor = LCA(tree, leaf_1, leaf_2)
                    intact_ancestors.append(ancestor)
                    intact_ancestor_leaves[ancestor] = [leaf_1, leaf_2]
        s = set(intact_ancestors)
        #print transcript, s
        
        for node in s:
            l = intact_ancestor_leaves[node]
            leaf_1 = l[0]
            leaf_2 = l[1]
            l1 = tree.get_node(leaf_1)
            l2 = tree.get_node(leaf_2)
            all_ret = dfs(gene_name, md_dict, ref, tree, node, leaf_1, \
                          leaf_2, 1, [])
            all_ids = [ret[0] for ret in all_ret]
            sub_tree_ids = []
            for curr_id in all_ids:
                sub_t = tree.subtree(curr_id)
                node_ids = [tree[node].tag for node in \
                            sub_t.expand_tree(mode=Tree.DEPTH) if \
                            tree[node].tag != curr_id]
                for node_id in node_ids:
                    sub_tree_ids.append(node_id)
            unique_ids = list(set(all_ids) - set(sub_tree_ids))
            for ret in all_ret:
                if ret[0] in unique_ids:
                    print gene_name, l1.tag, l2.tag, ','.join(ret[2])
                    writer.writerow([gene_name, l1.tag, l2.tag,
                                     '-'.join(ret[2])])

