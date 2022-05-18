import odgi
import pandas as pd
import csv
import sys

# In order to use this script, you might need to run these two commands in the terminal:
# env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=lib python3 -c 'import odgi'
# export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2

# Put the path to the .og and .gfa files here:
# og_path = "Test/test_1.og"
og_path = "yeast+edits.og"
# gaf_path = "GE00001631-DOT_H11_S191_R2_001.subset.gaf"
# gaf_path = "Test/test.gaf"
gaf_path = "GE00001631-DOT_A07_S103_R1_001.test.gaf"
# name_for_output = "GE00001631-DOT_A07_S103_R1_001.subset"
name_for_output = "test1"
# For automation:
if sys.argv[1:]:
    gaf_path = sys.argv[1]
    name_for_output = sys.argv[2]


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
# The list "all_path" contains a list in the form of [node_id, path_name]


def show_steps(handle):
    steps = []
    gr.for_each_step_on_handle(handle, lambda step: steps.append(step))
    all_path.append([
        gr.get_id(handle), " ".join([step_str(s) for s in steps])
    ])


# Everything above this line was borrowed and slightly adapted from the odgi tutorial on GitHub
gr.for_each_handle(show_steps)


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
    print(f"nodes: {nodes[:10]}")
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
# The list "pan_path" contains lists in the form of [node_id, [path_name, path_name, ...]]


def separate_paths(paths):
    """
    This function returns all the path names in a uniform format. This is necessary because
    some names are in a list with more than 1 element.
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
# These are just lists of the path names for the homology arms and the reference.


def create_node_dict(path):
    """
    creates a dictionary with the node IDs as keys and the path names as values
    :param path:
    :return:
    """
    temp_dict = {}
    for i in path:
        temp_dict[i[0]] = i[1]

    return temp_dict


node_dict = create_node_dict(pan_path)
# This dictionary can be used to look up the path names for a given node.


def get_path(node):
    """
    Careful, this function returns a list, even if there is only one element.
    Use a length or type check when this function is called.
    :param node:
    :return:
    """
    return node_dict.get(node)


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
            x = x.replace('<', '>')  # This makes it easier to split the string.
            # This is something that could change if we are interested in the orientation of the reads.
            con_nodes.append(x[1:].split('>'))
            # The [1:] above takes away the leading empty string that is left behind by the split function.
    return con_nodes


my_edges = find_legit_edges(gaf_nodes)


def create_edge_dict(edges):
    """
    This function takes a list of edges and creates a dictionary with the tuple of 2 concurrent nodes
    as the key and the number of reads mapping to that edge as the value. Tuples are sorted such that
    they are in ascending order. This means that the first node in the tuple is always the smaller node,
    regardless of original orientation of the read. This is for uniformity of the data as the list of
    edges in the reference and homology arm paths are also created in the same fashion.
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
            for j in range(len(nodes) - 1):
                temp = (int(nodes[j]), int(nodes[j + 1]))
                temp = sorted(temp)
                edge_list.append(tuple(temp))

    for i in edge_list:
        if i in edge_dict:
            edge_dict[i] += 1
        else:
            edge_dict[i] = 1
    return edge_dict


read_edges = create_edge_dict(my_edges)


def create_reverse_dict(dictionary):
    """
    This function takes as input a dictionary where there the keys are nodes and the values are a list of paths.
    The output is a dictionary where the keys are paths and the values are a list of nodes. This dictionary has far
    fewer entries than the input dictionary, making it much more efficient to look up the nodes for a given path.
    :param dictionary:
    :return:
    """
    graph_dict = {}
    for k, v in dictionary.items():
        for x in v:
            if x in graph_dict:
                graph_dict[x].append(k)
            else:
                graph_dict[x] = [k]
    return graph_dict


path_dict = create_reverse_dict(node_dict)
print(f"path_dict: {len(path_dict)}")
print(f"pan_path: {len(pan_path)}")
# Notice how much smaller the search space for the keys is in the path_dict.


def create_edges(path):
    """
    This function takes a list of nodes and creates a list of edges. This function will be called to
    create the edges of the paths so that they can be compared to the edges that had reads mapped to them.
    :param path:
    :return:
    """
    edges = []
    path = sorted(path)
    for j in range(len(path) - 1):
        temp = (int(path[j]), int(path[j + 1]))
        temp = sorted(temp)
        edges.append(tuple(temp))

    return edges


def create_edge_tally_dict(dictionary, paths):
    """
    This function takes a dictionary of edges {(node1, node2): number of reads} and a list of path names.
    It then creates a dictionary that maps the path names to the number of edges in the path that had reads mapped
    to them.
    :param dictionary:
    :param paths:
    :return:
    """
    tally_dictionary = {}
    for i in paths:
        for j in i:
            tally = 0
            # print(f"{j}: {path_dict.get(j)}")
            # print(f"path: {j}, tally: {create_edges(reverse_dict_search(node_dict, j))}")
            # tally_dictionary[j] = len(create_edges(reverse_dict_search(node_dict, j)))
            temp = create_edges(path_dict.get(j))
            for k in temp:
                if dictionary.get(k) is not None:
                    tally += dictionary.get(k)
            tally_dictionary[j] = tally

    print(len(tally_dictionary))
    return tally_dictionary


edge_tally = create_edge_tally_dict(read_edges, get_path_names(separate_paths(pan_path)))


def tally_reads_in_path(path):
    """
    This function takes a path and returns the number of edges in the path that had reads mapped to them.
    :param path:
    :return:
    """
    return edge_tally[path]


def get_tallies(paths):
    """
    This function takes a list of paths and a dictionary of edges and tallies the number of reads that mapped
    to each path.
    :param paths:
    :return:
    """
    tally_count = 0
    for i in paths:
        # print(f'{i}: {tally_reads_in_path(i, edge_dict)}')
        # print(f'{i}')

        tally_count += tally_reads_in_path(i)
    return tally_count


# Next, make a table with the homology arm coverage and ref coverage for each homology arm.
# I need to group the corresponding homology arm and ref-homology arm together.


def num_edges(path):
    """
    This function takes a path and returns the number of edges in the path.
    :param path:
    :return:
    """

    return len(create_edges(path_dict.get(path)))


print(f"{'homology_arm_34728-'}: {num_edges('homology_arm_34728-')}")
# Verified this works correctly.


def group_paths(paths):
    """
    This function takes a list of homology arm paths and groups the corresponding reference paths (over the
    homology arm path range) together. If either the homology arm or the reference path is not found to have zero
    edges, it is NOT grouped.
    :param paths:
    :return:
    """
    grouped_paths = []
    # This list will be in the form: [homology_arm_path, ref_path]. Remember that the ref_path is a sub-path over the
    # homology arm path's range.
    for hpath in paths:
        if num_edges(hpath) != 0 and num_edges(f"ref_{hpath}") != 0:
            grouped_paths.append([hpath, num_edges(hpath), f"ref_{hpath}", num_edges(f"ref_{hpath}")])
    return grouped_paths


paths_for_coverage = group_paths(hom_path)


def make_coverage_table(path_ids):
    """
    This function will take a list of path names and return a table of the coverage of each path.
    Format of the table is to be decided.
    :param path_ids:
    :return:
    """
    coverage_list = []
    # this list will be in the form: [[hom_arm_name, hom_arm_coverage, ref_arm_coverage], [etc]]
    # Expanded form: [[hom_arm_name, hom_arm_coverage, ref_subpath_coverage, #_of_hom_arm_edges,
    # count_for_hom_arm, #_of_ref_subpath_edges, count_for_ref_subpath], [etc]]
    for path in path_ids[1:]:
        if path[1] != 0 and path[3] != 0:
            coverage_list.append(
                [path[0], tally_reads_in_path(path[0]) / path[1],
                 tally_reads_in_path(path[2]) / path[3], path[1], tally_reads_in_path(path[0]),
                 path[3], tally_reads_in_path(path[2])])
    return coverage_list


cov_list = make_coverage_table(paths_for_coverage)
# The line below sorts the list by the homology arm coverage, making it easy to see which hom_arm has the highest
# coverage.
sorted_cov_list = sorted(cov_list, key=lambda row: row[1], reverse=True)


def write_to_tsv(coverage_list):
    """
    This function will take a list of coverage data and write it to a tsv file.
    """
    with open(f"{name_for_output}.tsv", "wt") as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        tsv_writer.writerow(["Homology arm", "Homology arm coverage", "Reference coverage", "Homology arm edges",
                             "Count for homology arm", "Reference edges", "Count for reference"])
        for i in coverage_list:
            tsv_writer.writerow(i)


write_to_tsv(sorted_cov_list)
print("Done!")
# Need to clean up this script and make it more readable. Also, make it more efficient.
