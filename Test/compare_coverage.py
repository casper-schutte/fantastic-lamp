import odgi
import csv
import argparse


# In order to use this script, you might need to run these two commands in the terminal:
# env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=lib python3 -c 'import odgi'
# export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2
# export LD_PRELOAD=/home/ec2-user/.conda/pkgs/libjemalloc-5.2.1-h9c3ff4c_6/lib/l$


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--og-path", required=True)
    parser.add_argument("--gaf-path", required=True)
    parser.add_argument("--out-path", required=True)
    parser.add_argument("--info-path", required=False)
    parser.add_argument("--test-example", required=False)
    return parser.parse_args()


args = parse_args()
gr = odgi.graph()
gr.load(args.og_path)
path_names = []
node_ids = []
gr.for_each_path_handle(lambda p: path_names.append(gr.get_path_name(p)))
gr.for_each_path_handle(lambda p: node_ids.append(p))


def step_str(step):
    path_str = gr.get_path_name(gr.get_path_handle_of_step(step))
    dir_str = "+" if not gr.get_is_reverse(gr.get_handle_of_step(step)) else "-"
    return path_str + dir_str


all_path = []


# The list "all_path" contains a list in the form of [node_id, path_names]


def show_steps(handle):
    steps = []
    gr.for_each_step_on_handle(handle, lambda step: steps.append(step))
    all_path.append([
        gr.get_id(handle), " ".join([step_str(s) for s in steps])
    ])


# Most everything above this line was borrowed and slightly adapted from the odgi tutorial on GitHub
gr.for_each_handle(show_steps)


def read_gaf(gaf_file_name):
    """
    This function takes a gaf file and returns a list of the nodes that reads in the
    file were mapped to.
    :param: gaf_file_name
    :return: list of nodes
    """
    nodes = []
    for row in open(gaf_file_name):
        row = row.split()
        mapq = int(row[11])
        if mapq < 30:
            continue
        nodes.append(row[5])
    return nodes


gaf_nodes = read_gaf(args.gaf_path)


def make_paths(path):
    """
    This function separates the names of the paths.
    In the input, the names are in one string, separated by spaces.
    :param: path
    :return: list of paths
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
    # The bash script ensures that the correct names are given (the homology arm paths all start with "h")
    h_arms = []
    ref_paths = []
    for path_name in paths:
        if path_name.startswith("hom"):
            h_arms.append(path_name)
        elif path_name.startswith("ref_h"):
            ref_paths.append(path_name)
    return h_arms, ref_paths


hom_path, ref_hom_path = get_path_names(separate_paths(pan_path))


# These are just lists of the path names for the homology arms and the reference homology arms.


def create_node_dict(path):
    """
    creates a dictionary with the node IDs as keys and the path names as values. If it detects that a node is shared
    between a homology-arm path and a corresponding reference-homology-arm path, it will add the path name "shared" to
    the node in this dictionary.
    :param path:
    :return:
    """
    temp_dict = {}
    for i in path:
        shared_score = 0
        for j in i[1]:
            if j.startswith("hom"):
                shared_score += 1
            if j.startswith("ref_h"):
                shared_score += 1
        if shared_score < 2:
            temp_dict[i[0]] = i[1]
        else:
            # i[1].append("shared")
            temp_dict[i[0]] = i[1]
    return temp_dict


# Need to investigate the efficiency of this function/method of separating the shared and unshared paths.


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
    and picks out the relevant edges (only reads that map to more than one node).
    The dictionary with the edges still needs to be created in a subsequent step.
    :param reads:
    :return:
    """
    con_nodes = []
    for x in reads:
        if x.count('<') > 1 or x.count('>') > 1:
            x = x.replace('<', '>')
            # This makes it easier to split the string.
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
        if v is not None:
            for x in v:
                if x in graph_dict:
                    graph_dict[x].append(k)
                else:
                    graph_dict[x] = [k]
    return graph_dict


path_dict = create_reverse_dict(node_dict)


def create_edges(path):
    """
    This function takes a list of nodes and creates a list of edges. This function will be called to
    create the edges of the paths so that they can be compared to the edges that had reads mapped to them.
    :param path:
    :return:
    """
    edges = []
    if path is not None:
        path = sorted(path)
        for j in range(len(path) - 1):
            temp = (int(path[j]), int(path[j + 1]))
            temp = sorted(temp)
            edges.append(tuple(temp))

    return edges


shared_edges = []


def create_shared_edges(hpaths):
    """
    This function takes each homology arm, and it's corresponding reference homology arm, and creates a list of edges
    that are shared between the two paths. This new list can be checked against to exclude shared edges further down
    the line.
    """
    counter = 0
    ref_edges = []
    h_edges = []
    for i in hpaths[1:]:
        ref_edges.append(create_edges(path_dict.get(f"ref_{i}")))
        h_edges.append(create_edges(path_dict.get(i)))
        for j in h_edges[counter]:
            if j in ref_edges[counter]:
                shared_edges.append(j)
        counter += 1


create_shared_edges(hom_path)
# print(shared_edges)


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
            temp = create_edges(path_dict.get(j))
            for k in temp:
                if k not in shared_edges and dictionary.get(k) is not None:
                    tally += dictionary.get(k)
            tally_dictionary[j] = tally
    return tally_dictionary


edge_tally = create_edge_tally_dict(read_edges, get_path_names(separate_paths(pan_path)))


def tally_reads_in_path(path):
    """
    This function takes a path and returns the number of edges in the path that had reads mapped to them.
    :param path:
    :return:
    """
    return edge_tally[path]


def num_edges(path):
    """
    This function takes a path and returns the number of edges in the path. It excludes edges that are shared between
    the reference and the homology arm.
    :param path:
    :return:
    """
    temp = []
    for i in create_edges(path_dict.get(path)):
        if i not in shared_edges:
            temp.append(i)
    return len(temp)


def group_paths(paths):
    """
    This function takes a list of homology arm paths and groups the corresponding reference paths (over the
    homology arm path range) together. If either the homology arm or the reference path is not found to have zero
    edges, it is NOT grouped. The list is in the form [homology arm path, reference path].
    :param paths:
    :return:
    """
    grouped_paths = []
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
    # this list will be in the form: [[hom_arm_name, hom_arm_coverage, ref_subpath_coverage, #_of_hom_arm_edges,
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
#print(sorted_cov_list[:10])


def write_to_tsv(coverage_list):
    """
    This function will take a list of coverage data and write it to a tsv file.
    """
    with open(args.out_path, "wt") as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        tsv_writer.writerow(["Homology arm", "Homology arm coverage", "Reference coverage", "Homology arm edges",
                             "Count for homology arm", "Reference edges", "Count for reference"])
        for i in coverage_list:
            tsv_writer.writerow(i)


write_to_tsv(sorted_cov_list)
print("Done!")
