import math
import re
import copy
import random

def get_tree_from_file(treefile="tree.txt"):
    '''
    :param treefile: file of a tree
    :return: array of nodes, map of nodes to their children, map of edges between nodes,
            q_a

    '''

    nodes = ["Root"]
    q_a = [1/4, 1/4, 1/4, 1/4]
    child_dict = dict()
    edges_dict = dict()

    with open(treefile) as f:
        for line in f:
            line = line.rstrip('\n')
            child, parent, weight = line.split(' ')
            weight = float(weight)

            if parent not in child_dict:
                child_dict[parent] = [child]
            else:
                child_dict[parent]+=[child]
            edges_dict[(parent,child)] = weight
            nodes.append(child)

    return [nodes, child_dict, edges_dict, q_a]
a = get_tree_from_file()
#print(a[1])

def subtree_prob(a,b,t,alpha, bases ="ACGT"):
    '''
    calculates probabilty of one subtree (one child)
    probabilities of children are multiplied later
    (apparently its faster
    http://compbio.fmph.uniba.sk/vyuka/mbi/index.php/CI08
    Zlozitost, zlepsenie
    )

    :param a: parent node
    :param b: child node

    '''
    
    new_value = []
    for i, base_val in enumerate(bases):
        b_sum = 0
        for bi, b_vaL in enumerate(bases):
            b_sum+=(b[bi] * jukes_matrix_prob(base_to_int(base_val),bi,t,alpha))
        new_value.append(b_sum)
    return new_value

def jukes_matrix_prob(b1,b2,t,alpha):
    #base 1, base2, time, alpha
    if b1==b2:
        return (1+3*(math.e**(-4.0/3.0 * alpha * t)))/4
    else:
        return (1 - (math.e**(-4.0/3.0 * alpha * t)))/4


def get_leaves(child_dict, nodes):
    leaves = []
    for i in nodes:
        if i not in child_dict:
            leaves.append(i)
    return leaves

def get_leaf_parents(nodes, leaves, child_dict):
    leaf_parents = []
    for node in nodes:
        if node in child_dict:
            for child in child_dict[node]:
                if child not in leaves:
                    break
                if node not in leaf_parents:
                    leaf_parents.append(node)
    return leaf_parents


def base_to_int(base, bases = "ACGT"):
    return bases.index(base)

def felsenstein(tree, alpha, col,bases = "ACGT"):
    '''
    :param tree: return of get_tree_from_file
    :param alpha: rate of change
    :param col: alignment column

    :return: probability, float (0,1)
    '''
    t=copy.deepcopy(tree)
    nodes, child_dict, edges_dict, q_a = t
    leaves = get_leaves(child_dict, nodes)
    node_values = dict()

    #create node values
    for i, node in enumerate(nodes):
        start_val = 0
        #if info is missing use [1,1,1,1...,1]
        if "-" in node:
            start_val = 1
        val_arr = [start_val]*len(bases)
        if node not in leaves:
            node_values[node] = val_arr
        #if node is a leaf 
        else:
            #if info is missing use [1,1,1,1...,1]
            if col[leaves.index(node)] in "N-":
                node_values[node]=[1]*len(bases)
            #else [0,0,....1,0,...0,0]
            else:
                val_arr[base_to_int(col[leaves.index(node)])] = 1
                node_values[node] = val_arr
    #until we disconnect everything
    while len(leaves)<len(nodes):
        #get parents of leaves
        leaf_parents = get_leaf_parents(nodes,leaves,child_dict)
        for parent in leaf_parents:
            parent_vals = node_values[parent]
            leaf1 = child_dict[parent][0]
            leaf2 = child_dict[parent][1]
            #(sum A[y,b]P[B|a,t])
            b_val = subtree_prob(parent_vals,node_values[leaf1],edges_dict[(parent,leaf1)], alpha)
            c_val = subtree_prob(parent_vals,node_values[leaf2],edges_dict[(parent,leaf2)],alpha)
            #A[v,a] = ....
            for i in range(len(b_val)):
                parent_vals[i]=b_val[i]*c_val[i]

            node_values[parent] = parent_vals
            child_dict.pop(parent)
            # if parent=="Root":
            #     final_prob = 0
            #     for i, prob in enumerate(node_values["Root"]):
            #         final_prob+=prob*q_a[i]
            #     return final_prob

        leaves = get_leaves(child_dict,nodes)

    final_prob = 0
    for i, prob in enumerate(node_values["Root"]):
        final_prob+=prob*q_a[i]
    return final_prob


def make_cols(file="cftr.txt"):
    with open(file) as f:
        a  = f.read()
    #split podla >...\n
    rows = re.split(">[A-Za-z]*\n",a)
    rows = [i.replace("\n","") for i in rows[1:]]
    #array aby sedel ordering
    species = re.findall(">[A-Za-z]*\n",a)
    species = [i[1:-1] for i in species]
    return (species,rows)


def find_alpha(tree, species_info,start=0,stop=None):
    '''
    :param tree: tree
    :param species info: list of names of species, list of rows that represent whole genetic sequence
                        (\n removed etc.)
    :param start, stop: task d) best alpha for window
    :return: a* for window [start,stop)
    '''

    species, rows = species_info
    t = [i.copy() for i in tree]
    alpha_star = 0.1
    best_prob = float('-inf')
    alpha = 0.1
    probs = []
    haf = dict()
    if stop==None:
        stop = len(rows[0])
    while alpha <=2:
        alpha = round(alpha,2)
        probs = []
        for i in range(start,stop):
            col = "".join([x[i] for x in rows])
            prob = felsenstein(t,alpha,col)
            if prob!=0:
                probs.append(math.log(prob))
            else:
                probs.append(0)
        if sum(probs)>=best_prob:
            best_prob=sum(probs)
            alpha_star=alpha
        haf[alpha]=sum(probs)
        alpha+=0.1
    return alpha_star



def animal_genome_string(animal_string,filename="ctfr.txt"):
    #turns cftr data into 1 long line string
    #ran it through regex validator so it should be fine
    with open("cftr.txt") as file:
        f = file.read()
        animal = re.search(f">{animal_string}\n([ACGTN]|-|\n)*",f)
        animal = animal.group()
        animal = re.sub(f"(\n|>{animal_string})","",animal)        
    return animal
def remap_exon_ind(ind,string):
    #bazy císlujeme od 1
    out = ind-1
    i=0
    out = 0
    #nepocítame pomlcky
    while i!=ind:
        if (string[i]=="-"):
            ind+=1
        i+=1
    return ind

def make_exons_list(filename="exons.txt"):
    #turn exons.txt into array of binary tupples
    exons = []
    with open(filename) as file:
        for line in file:
            tmp = line.rstrip("\n").split()
            exons.append([int(i) for i in tmp])
    return exons



def make_alphas_arr(filename="out.txt"):
    #out.txt from task d has the format of float;float;float;....float;float
    with open(filename) as file:
        f = file.read()
        alpha_arr = f.split(";")
        alpha_arr = [float(i) for i in alpha_arr]
    return alpha_arr
    
def alphas_for_exons(t):
    human_string = animal_genome_string("Human")
    alpha_arr = make_alphas_arr()
    #print(len(alpha_arr))
    exons = make_exons_list()
    new_exons = []
    exon_windows = set()
    non_exon_windows = set()
    all_windows = set()
    for exon in exons:
        new_start = remap_exon_ind(exon[0],human_string)
        new_end = remap_exon_ind(exon[1],human_string)
        new_exons.append([new_start,new_end])
    for i in range(0,len(human_string),100):
        all_windows.add(i)
        for e_start, e_end in new_exons:
            if (i in range(e_start,e_end) or i+100 in range(e_start,e_end)):
                exon_windows.add(i)
                break
    all_windows = list(all_windows)
    # print(exon_windows)
    exonX = [i for i in list(exon_windows) if i//100 < len(alpha_arr)]
    exonY = [alpha_arr[i//100] for i in all_windows[:] if i in exon_windows and i//100 < len(alpha_arr)]
    non_exonX = [i//1 for i in all_windows[:] if i not in exon_windows and i//100 < len(alpha_arr)]
    # print(len(alpha_arr))
    # print(len(all_windows))
    non_exonY = [alpha_arr[(i)//100] for i in all_windows[:] if i not in exon_windows and i//100 < len(alpha_arr)]

    return exonX, exonY, non_exonX, non_exonY


def compute_alpha_for_windows(tree=a):
    window_dict = []
    cols = make_cols()
    print(len(cols[1][1]))
    out=[]
    for i in range(0,len(cols[1][0])-100,100):
        try:
            out.append(find_alpha(a,cols,i,i+100))
            with open("out.txt","w") as file:
                g = [str(i) for i in out]
                file.write(";".join(g)) 

        except:
            break
    with open("out.txt","w") as file:
        g = [str(i) for i in out]
        file.write(";".join(g)) 



haf = alphas_for_exons(a)
for i in haf:
    print(i[:10])