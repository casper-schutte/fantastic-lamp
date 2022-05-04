import odgi
import pandas as pd

# In order to use this script, you might need to run these two commands in the terminal:
# env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=lib python3 -c 'import odgi'
# export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2

# Put the path to the .og and .gfa files here:
# og_path = "Test/test_1.og"
og_path = "yeast+edits.og"
gaf_path = "GE00001631-DOT_H11_S191_R2_001.subset.gaf"
# gaf_path = "Test/test.gaf"
gr = odgi.graph()
gr.load(og_path)
path_names = []
node_ids = []
gr.for_each_path_handle(lambda p: path_names.append(gr.get_path_name(p)))
gr.for_each_path_handle(lambda p: node_ids.append(p))


def step_str(step):
    path_str = gr.get_path_name(gr.get_path_handle_of_step(step))
    dir_str = "+" if not gr.get_is_reverse(gr.get_handle_of_step(step)) else "-"
    return path_str + dir_str


all_path = []


def show_steps(handle):
    steps = []
    gr.for_each_step_on_handle(handle, lambda step: steps.append(step))
    all_path.append([
        gr.get_id(handle), " ".join([step_str(s) for s in steps])
    ])


# Everything above this line was borrowed and slightly adapted from the odgi tutorial on GitHub
gr.for_each_handle(show_steps)
# The list "all_path" contains a list in the form of [node_id, path_name]


def read_gaf(gaf_file_name):
    """
    This function takes a gaf file and returns a list of the nodes that reads in the
    file were mapped to.
    :param gaf_file_name:
    :return:
    """
    headers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    gaf_file = pd.read_csv(gaf_file_name, sep="\t", header=None, names=headers)
    # gaf_file = gaf_file.dropna()
    # Using dropna() broke something? Will just have to deal with the missing values manually
    # further down the line.
    nodes = []
    for j in range(0, len(gaf_file)):
        # Collect the node IDs in the appropriate column
        nodes.append(gaf_file[6][j])
    # print(f"nodes: {nodes[:10]}")
    return nodes


gaf_nodes = read_gaf(gaf_path)

# The code below names the columns in the gaf file, in case you want to use more than just the nodes that
# had reads mapped to them.
# specify headers: headers = ["q_name", "q_len", "q_start", "q_end", "rel_orient", "path", "path_len", "path_start",
# "path_end", "res_matches", "aln_block_len", "mapq"]


def make_paths(path):
    """
    This function separates the names of the paths.
    In the input, the names are in one string, separated by spaces.
    :param path:
    :return:
    """
    my_list = []
    for i in path:
        my_list.append([i[0], i[1].split(" ")])
    return my_list


pan_path = make_paths(all_path)


def separate_paths(paths):
    """
    This function returns all the path names in a uniform format. This is necessary because the
    some of the names are in a list with more than 1 element.
    :param paths:
    :return:
    """
    sep_path_names = []
    for x in paths:
        if len(x[1]) > 1:
            for y in x[1]:
                if y not in sep_path_names:
                    sep_path_names.append(y)
        else:
            if x[1][0] not in sep_path_names:
                sep_path_names.append(x[1][0])
    return sep_path_names


def get_path_names(paths):
    """
    This function returns 2 lists, one with the path names for the reference
    and one for the homology arms.
    :param paths:
    :return:
    """
    # Need to write this function such that it would work for any example, not only
    # this specific one. Could take the names from the original files, need to think
    # about this.
    h_arms = []
    ref_paths = []
    for path_name in paths:
        if path_name.startswith("h"):
            h_arms.append(path_name)
        else:
            ref_paths.append(path_name)
    return h_arms, ref_paths


hom_path, ref_path = get_path_names(separate_paths(pan_path))


def create_node_dict(path):
    """
    creates a dictionary with the node IDs as keys and the path names as values
    :param path:
    :return:
    """
    path_dict = {}
    for i in path:
        path_dict[i[0]] = i[1]

    return path_dict


node_dict = create_node_dict(pan_path)


def get_path(node):
    """
    Careful, this function returns a list, even if there is only one element.
    Use a length or type check when this function is called.
    :param node:
    :return:
    """
    return node_dict.get(node)


# I now have a dictionary of node_id: path_name that I can use to query whether a node (from the reads)
# is in a reference path, homology arm, or both.

# I now need to figure out how to make edges from the paths.


def find_legit_edges(reads):
    """
    This function takes a list of nodes to which reads mapped (from the GAF file)
    and picks out the relevant edges (only reads that map to more than one node, for now).
    The dictionary with the edges still needs to be created in a subsequent step.
    :param reads:
    :return:
    """
    con_nodes = []
    for x in reads:
        if x.count('<') > 1 or x.count('>') > 1:
            x = x.replace('<', '>')
            con_nodes.append(x[1:].split('>'))

    # print(len(con_nodes))
    return con_nodes


my_edges = find_legit_edges(gaf_nodes)


def create_edge_dict(edges):
    """
    This function takes a list of edges and creates a dictionary with the tuple of 2 concurrent nodes
    as the key and the number of reads mapping to that edge as the value. Tuples are sorted such that
    they are in ascending order. This means that the first node in the tuple is always the smaller node,
    regardless of original orientation of the read.
    :param edges:
    :return:
    """
    edge_list = []
    edge_dict = {}
    for nodes in edges:
        if len(nodes) == 2:
            temp = (int(nodes[0]), int(nodes[1]))
            temp = sorted(temp)
            edge_list.append(tuple(temp))
        elif len(nodes) > 2:
            for j in range(len(nodes)-1):
                temp = (int(nodes[j]), int(nodes[j+1]))
                temp = sorted(temp)
                edge_list.append(tuple(temp))

    for i in edge_list:
        if i in edge_dict:
            edge_dict[i] += 1
        else:
            edge_dict[i] = 1
    return edge_dict


read_edges = create_edge_dict(my_edges)


def reverse_dict_search(dictionary, value):
    """
    This function takes a dictionary and a value and returns a list of keys that have that value.
    :param dictionary:
    :param value:
    :return:
    """
    # Old method:
    # return [k for k, v in dictionary.items() if v == value[n] for n in range(len(value))]
    key = []
    for k, v in dictionary.items():
        if value in v:
            key.append(k)
    return key


# testing
# print(f'{hom_path[619]}: {reverse_dict_search(node_dict, hom_path[619])}')
# print(f'{ref_path[1]}: {reverse_dict_search(node_dict, [ref_path[1]])}')

# Next step: write a function that iterates over each edge in each path (ref and homology) and tallies the number
# of reads that map in each path.


def create_edges(path):
    """
    This function takes a list of nodes and creates a list of edges. This function will be called to
    create the edges of the paths so that they can be compared to the edges that had reads mapped to them.
    :param path:
    :return:
    """
    edges = []
    for j in range(len(path) - 1):
        temp = (int(path[j]), int(path[j + 1]))
        temp = sorted(temp)
        edges.append(tuple(temp))

    return edges


# Reminder of where what data is. The names of the paths are in ref_path and hom_path. The actual paths are called
# with reverse_dict_search(node_dict, PATH[i]). They then need to be made into edges with create_edges(). I can
# then compare the edges of the paths to the edges that had reads mapped to them.

def tally_reads_in_path(path, edge_dict):
    """
    This function takes a path and a dictionary of edges and tallies the number of reads that map to each edge in
    the path.
    :param path:
    :param edge_dict:
    :return:
    """
    read_tally = 0
    graph_edges = create_edges(reverse_dict_search(node_dict, path))
    for edge in graph_edges:
        if edge in edge_dict:
            read_tally += edge_dict[edge]
    return read_tally


def get_tallies(paths, edge_dict):
    """
    This function takes a list of paths and a dictionary of edges and tallies the number of reads that mapped
    to each path.
    :param paths:
    :param edge_dict:
    :return:
    """
    tally_count = 0
    for i in paths:
        print(f'{i}: {tally_reads_in_path(i, edge_dict)}')
        # print(f'{i}')

        tally_count += tally_reads_in_path(i, edge_dict)
    return tally_count


total_ref = get_tallies(ref_path, read_edges)
total_hom = get_tallies(hom_path, read_edges)
print(f"Homology arm: {total_hom}")
print(f"Reference: {total_ref}")

print(reverse_dict_search(node_dict, "homology_arm_35426:0-138+"))

# Next, make a table with the homology arm coverage and ref coverage for each homology arm.


def make_coverage_table(path_ids):
    """
    This function will take a list of path names and return a table of the coverage of each path.
    Format of the table is to be decided.
    :param path_ids:
    :return:
    """
