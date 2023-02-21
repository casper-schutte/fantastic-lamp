import odgi
import csv
import argparse


# In order to use this script, you might need to run these two commands in the terminal:
# env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=lib python3 -c 'import odgi'
# export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2 (or wherever your jemalloc is located)
# If you are running the bash script, these commands are automatically initiated.


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--og-path", required=True)
    parser.add_argument("--gaf-path", required=True)
    parser.add_argument("--out-path", required=True)
    parser.add_argument("--info-path", required=False)
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
    file were mapped to. Only reads with MAPQ higher than 30 are considered.
    :param: gaf_file_name
    :return: list of nodes
    """
    nodes = []
    with open(gaf_file_name) as f:
        for row in f:
            row = row.split()
            qual = int(row[11])
            path = row[5]
            if qual < 30:
                continue
            nodes.append(path)

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
    and one with the path names for homology arms.
    :param paths:
    :return:
    """
    # The bash script ensures that the correct names are given (the homology arm paths all start with "homology_arm_*" and
    # the reference homology arm paths all start with "ref_homology_arm_*")
    h_arms = []
    ref_paths = []
    for path_name in paths:
        if path_name.startswith("hom"):
            h_arms.append(path_name)
        elif path_name.startswith("ref_h"):
            ref_paths.append(path_name)
    return h_arms, ref_paths


hom_path, ref_hom_path = get_path_names(separate_paths(pan_path))
# This creates the lists of the path names for the homology arms and the reference homology arms.


def create_node_dict(path):
    """
    creates a dictionary with the node IDs as keys and the path names as values. If it detects that a node is shared
    between a homology-arm path and a corresponding reference-homology-arm path, it will add the path name "shared" to
    the node in this dictionary. These shared nodes are discarded when calculating coverage.
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
            # Double check this function? Is this correct?
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
    for i in hpaths[1:]:
        ref_edges = create_edges(path_dict.get(f"ref_{i}"))
        h_edges = create_edges(path_dict.get(i))
        for j in h_edges:
            if j in ref_edges:
                shared_edges.append(j)


create_shared_edges(hom_path)


def read_name_dict():
    """
    This function mapes edges (keys) to the read names (values) in a dictionary.
    :return:
    """
    read_info = {}
    with open(args.gaf_path) as f:
        for row in f:
            row = row.split()
            qual = int(row[11])
            name = row[0]
            path = row[5]
            if qual < 30:
                continue
            if path.count(">") > 1 or path.count("<") > 1:
                if f"{name}A" in read_info.keys():
                    read_info[f"{name}B"] = create_edges(find_legit_edges([path])[0])
                else:
                    read_info[f"{name}A"] = create_edges(find_legit_edges([path])[0])

    new_dict = {}
    for k, v in read_info.items():
        if v is not None:
            for x in v:
                if x not in shared_edges:
                    if x in new_dict:
                        new_dict[x].append(k)
                    else:
                        new_dict[x] = [k]
                else:
                    if x not in new_dict:
                        new_dict[x] = None
    return new_dict


read_dict = read_name_dict()


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
                    if dictionary.get(k) > tally:
                        tally = dictionary.get(k)
            tally_dictionary[j] = tally
    return tally_dictionary


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


read_table = {}
count_table = {}


def group_paths(paths):
    """
    This function takes a list of homology arm paths and groups the corresponding reference paths (over the
    homology arm path range) together. If either the homology arm or the reference path is found to have zero
    edges, it is NOT grouped. The list is in the form [homology arm path, reference path].
    :param paths:
    :return:
    """
    grouped_paths = []
    for hpath in paths:
        if num_edges(hpath) != 0 and num_edges(f"ref_{hpath}") != 0:
            grouped_paths.append([hpath, num_edges(hpath), f"ref_{hpath}", num_edges(f"ref_{hpath}")])
            temp = create_edges(path_dict.get(hpath))
            temp2 = [hpath, []]
            for i in temp:
                if read_dict.get(i) is not None:
                    for j in read_dict.get(i):
                        if j not in temp2[1] and j not in shared_edges:
                            temp2[1].append(j)
            # read_table and count_table map the name of a path to the names of reads that map to them and to the
            # number of reads that map to them, respectively.
            read_table[temp2[0]] = temp2[1]
            count_table[temp2[0]] = len(temp2[1])
            ref_edges = create_edges(path_dict.get(f"ref_{hpath}"))
            temp3 = [f"ref_{hpath}", []]
            for i in ref_edges:
                if read_dict.get(i) is not None:
                    for j in read_dict.get(i):
                        if j not in temp3[1] and j not in shared_edges:
                            temp3[1].append(j)
            read_table[temp3[0]] = temp3[1]
            count_table[temp3[0]] = len(temp3[1])
    return grouped_paths


paths_for_coverage = group_paths(hom_path)


def tally_reads_in_path(path):
    """
    This function takes a path and returns the number of edges in the path that had reads mapped to them.
    :param path:
    :return:
    """
    return count_table.get(path)


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


def write_to_tsv(coverage_list):
    """
    This function will take a list of coverage data and write it to a tsv file. It also writes a file containing information 
    # mapping each read to the homology arms to which it maps, which is used for deubugging and for additional information.
    """
    with open(args.out_path, "wt") as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        tsv_writer.writerow(["Homology arm", "Homology arm coverage", "Reference coverage", "Homology arm edges",
                             "Count for homology arm", "Reference edges", "Count for reference"])
        for i in coverage_list:
            tsv_writer.writerow(i)

    with open("my_read_info.tsv", "wt") as tsv_file:
        fieldnames = ["Reads", "Homology arm"]
        tsv_writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        tsv_writer.writeheader()
        for i in read_table:
            for j in read_table[i]:
                tsv_writer.writerow({"Reads": j, "Homology arm": i})


write_to_tsv(sorted_cov_list)
print("Done!")
