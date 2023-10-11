import os
import shutil
import statistics
import subprocess
import sys
import pandas as pd
import networkx as nx
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
import numpy as np

from collections import Counter
from Bio import AlignIO, SeqRecord, Seq, SeqIO
from Bio.Align import MultipleSeqAlignment
from ete3 import Tree
from scipy.stats import entropy, skew
from EBG.Features.parsimony_mutation_counter import count_subst_freqs
from EBG.utils.utils import setup_logger, check_raxml_availability


class FeatureComputer:
    """
        Performs all feature compuatations necessary for EBG.

        Attributes
        ----------
        logger : logger
            for communcation with command line
        msa_filepath : str
            absolute path to the .fasta file
        model_filepath : str
            absolute path to the RAxML-NG model file
        tree_filepath : str
            absolute path to the .newick tree file
        current_directory : os.path
            current working directory
        output_prefix : str
             name of the ouput EBG produces
        raxml_ng_path : str
            path to the RAxML-NG executable

        Methods
        -------
        traverse_and_add_edges(self, node_, graph):
            Creates a networkx graph out of ete3 tree
        compute_substitution_features(self):
            Computes features from the parsimony substitution counts per site
        compute_split_features_tree(self):
            Computes features out of possible bipartitions of the tree
        compute_parsimony_support_features(self):
            Computes support features from parsimony starting trees
        compute_parsimony_support(self):
            Computes parsimony support itself
        compute_parsimony_bootstrap_support_features(self):
            Computes features from parsimony bootstrap support
       compute_parsimony_bootstrap_support(self):
            Computes parsimony bootstrap support itself
        """
    def __init__(self, msa_filepath, tree_filepath, model_filepath, output_prefix, raxml_ng_path, redo):
        self.logger = setup_logger("FeatureComputer")
        self.msa_filepath = msa_filepath
        self.model_filepath = model_filepath
        self.tree_filepath = tree_filepath
        self.current_directory = os.path.abspath(os.curdir)
        self.output_prefix = output_prefix
        if raxml_ng_path == "raxml-ng":
            self.raxml_path = check_raxml_availability()
            if self.raxml_path == None:
                self.logger.error("RAxML-NG not found in path variable. Please make sure it exists or pass the full path as input parameter (-raxmlng PATH)")
        else:
            self.raxml_path = raxml_ng_path
            if not os.path.isfile(self.raxml_path):
                self.logger.error(
                    f"RAxML-NG not found in the provided path {raxml_ng_path}. Please make sure it exists and provide the absolute path. Exiting EBG.")
                sys.exit()

        tmp_folder_path = os.path.abspath(os.path.join(os.curdir, output_prefix))

        if os.path.exists(tmp_folder_path):
            if redo:
                self.logger.warning(
                    f"Found exisiting result folder: {tmp_folder_path}, it will be overwritten because of redo mode")
                shutil.rmtree(tmp_folder_path)
            else:
                self.logger.error(
                    f"Found exisiting result folder: {tmp_folder_path}, please rename the folder or run in -redo mode. Exiting EBG.")
                sys.exit()

        os.makedirs(tmp_folder_path)
        tmp_folder_path_work = os.path.abspath(os.path.join(os.curdir, output_prefix, "tmp"))
        os.makedirs(tmp_folder_path_work)
        os.chdir(tmp_folder_path_work)

    def traverse_and_add_edges(self, node_, graph) -> nx.Graph:
        """
        Recursively constructs a networkx graph out of a ete3 tree

                Parameters:
                        :param node_: ete3 tree
                        :param graph: instance of networkx graph

                Returns:
                        :return graph: ete3 tree as graph

        """
        for child in node_.children:
            edge_weight = node_.get_distance(child)
            graph.add_edge(node_.name, child.name, weight=edge_weight)
            self.traverse_and_add_edges(child, graph)
        return graph

    def compute_substitution_features(self) -> dict:
        """
        Computes coefficient of variation, skewness, max and mean of the parsimony substitution count per site

                Parameters:

                Returns:
                        :return dict: dictionary with the four summary statistics as values

        """
        with open(self.tree_filepath, "r") as tree_file:
            tree_data = tree_file.read()
            tree = Tree(tree_data)
            leaf_count = len(tree.get_leaves())

            result = count_subst_freqs(self.tree_filepath, self.msa_filepath)

            max_subst_freq = max(result) / leaf_count
            avg_subst_freq = (sum(result) / len(result)) / leaf_count

            if statistics.mean(result) != 0:
                cv_subst_freq = statistics.stdev(result) / statistics.mean(result)
            else:
                cv_subst_freq = 0
            if max(result) == min(result):
                normalized_result = result
            else:
                normalized_result = [(x - min(result)) / (max(result) - min(result)) for x in
                                     result]
            sk_subst_freq = skew(normalized_result)

        self.logger.info("Finished computing parsimony substitution features!")
        return {"max_substitution_frequency": max_subst_freq,
                "mean_substitution_frequency": avg_subst_freq,
                "cv_substitution_frequency": cv_subst_freq,
                "skw_substitution_frequency": sk_subst_freq}

    def compute_split_features_tree(self) -> pd.DataFrame:
        """
        Splits the tree at each inner branch. Computes the mean closeness centrality ratio and the branch length ratio between the produced subtrees.

                Parameters:

                Returns:
                        :return pd.DataFrame: dataframe with the branchId and both features as columns

        """
        with open(self.tree_filepath, "r") as tree_file:
            tree_str = tree_file.read()
            phylo_tree = Tree(tree_str)
            branch_id_counter = 0

            for node in phylo_tree.traverse():
                branch_id_counter += 1
                if not node.is_leaf():
                    node.__setattr__("name", branch_id_counter)

            results = []
            phylo_tree_original = phylo_tree.copy()  # get reference copy
            for node in phylo_tree.traverse():
                if not node.is_leaf():
                    phylo_tree = phylo_tree_original.copy()  # copy reference
                    left_subtree = node.detach()
                    right_subtree = phylo_tree

                    branch_lengths_left = [node.dist for node in left_subtree.traverse() if not node.is_root()]
                    branch_lengths_right = [node.dist for node in right_subtree.traverse() if not node.is_root()]
                    sum_bl_left = sum(branch_lengths_left)
                    sum_bl_right = sum(branch_lengths_right)
                    bl_ratio = min(sum_bl_left, sum_bl_right) / max(sum_bl_right, sum_bl_left)

                    G_left = nx.DiGraph()
                    G_left_ = self.traverse_and_add_edges(left_subtree, graph=G_left)
                    closeness_centrality_left = list(nx.closeness_centrality(G_left_, distance='weight').values())
                    mean_clo_sim_left = sum(closeness_centrality_left) / len(closeness_centrality_left)

                    G_right = nx.DiGraph()
                    G_right_ = self.traverse_and_add_edges(right_subtree, graph=G_right)
                    closeness_centrality_right = list(nx.closeness_centrality(G_right_, distance='weight').values())
                    mean_clo_sim_right = sum(closeness_centrality_right) / len(closeness_centrality_right)

                    try:
                        mean_clo_sim_ratio = min(mean_clo_sim_left, mean_clo_sim_right) / max(mean_clo_sim_left,
                                                                                              mean_clo_sim_right)
                    except ZeroDivisionError:
                        mean_clo_sim_ratio = -1

                    result = (node.name, bl_ratio, mean_clo_sim_ratio)
                    results.append(result)
            self.logger.info("Finished computing tree split features!")

            return pd.DataFrame(results,
                                columns=["branchId", "branch_length_ratio_split", "mean_closeness_centrality_ratio"])

    def compute_split_features_msa(self):
        with open(self.tree_filepath, "r") as tree_file:
            tree_str = tree_file.read()
            phylo_tree = Tree(tree_str)
            branch_id_counter = 0

            for node in phylo_tree.traverse():
                branch_id_counter += 1
                if not node.is_leaf():
                    node.__setattr__("name", branch_id_counter)

            results = []

            for node in phylo_tree.traverse("postorder"):
                if (not node.is_root()) and (not node.is_leaf()):

                    list_a = []
                    list_a_dist_branch = []
                    list_a_dist_topo = []
                    list_b = []
                    list_b_dist_branch = []
                    list_b_dist_topo = []
                    for leaf in phylo_tree.get_leaves():
                        if leaf in node.get_leaves():
                            list_a.append(leaf.name)
                            list_a_dist_branch.append(leaf.get_distance(target=phylo_tree.get_tree_root()))
                            list_a_dist_topo.append(
                                leaf.get_distance(topology_only=True, target=phylo_tree.get_tree_root()))
                        else:
                            list_b.append(leaf.name)
                            list_b_dist_branch.append(leaf.get_distance(target=phylo_tree.get_tree_root()))
                            list_b_dist_topo.append(
                                leaf.get_distance(topology_only=True, target=phylo_tree.get_tree_root()))

                    split_mean_dist_branch_a = statistics.mean(list_a_dist_branch)
                    split_std_dist_branch_a = np.std(list_a_dist_branch)

                    split_mean_dist_branch_b = statistics.mean(list_b_dist_branch)
                    split_std_dist_branch_b = np.std(list_b_dist_branch)

                    split_std_ratio_branch = min(split_std_dist_branch_a, split_std_dist_branch_b) / max(
                        split_std_dist_branch_a, split_std_dist_branch_b)
                    split_mean_ratio_branch = min(split_mean_dist_branch_a, split_mean_dist_branch_b) / max(
                        split_mean_dist_branch_a, split_mean_dist_branch_b)

                    split_std_dist_topo_a = np.std(list_a_dist_topo)

                    split_std_dist_topo_b = np.std(list_b_dist_topo)

                    split_std_ratio_topo = min(split_std_dist_topo_a, split_std_dist_topo_b) / max(
                        split_std_dist_topo_a, split_std_dist_topo_b)

                    alignment = AlignIO.read(self.msa_filepath, 'fasta')
                    alignment_a = MultipleSeqAlignment([])
                    alignment_b = MultipleSeqAlignment([])
                    for record in alignment:
                        if record.id in list_a:
                            alignment_a.append(record)
                        elif record.id in list_b:
                            alignment_b.append(record)

                    freqs_b = []
                    freqs_a = []

                    entropy_differences = []

                    for i in range(len(alignment_a[0])):
                        column_a = alignment_a[:, i]

                        column_b = alignment_b[:, i]

                        combined_values = column_a + column_b
                        all_keys = set(combined_values)

                        counter_a = Counter({key: 0 for key in all_keys})
                        counter_b = Counter({key: 0 for key in all_keys})

                        counter_a.update(column_a)
                        counter_b.update(column_b)

                        sorted_keys = sorted(all_keys)

                        counter_a = Counter({key: counter_a[key] for key in sorted_keys})
                        counter_b = Counter({key: counter_b[key] for key in sorted_keys})

                        freqs_a.append(counter_a)
                        freqs_b.append(counter_b)

                    for site_freq_a, site_freq_b in zip(freqs_a, freqs_b):
                        total_count_a = sum(site_freq_a.values())
                        total_count_b = sum(site_freq_b.values())
                        try:
                            normalized_freq_a = {k: v / total_count_a for
                                                 k, v in site_freq_a.items()}
                        except ZeroDivisionError:
                            normalized_freq_a = {k: 0 for
                                                 k, v in site_freq_a.items()}
                        try:
                            normalized_freq_b = {k: v / total_count_b for
                                                 k, v in site_freq_b.items()}
                        except:
                            normalized_freq_b = {k: 0 for
                                                 k, v in site_freq_b.items()}

                        site_freq_a_array = np.array(list(normalized_freq_a.values()))

                        site_freq_b_array = np.array(list(normalized_freq_b.values()))

                        entropy_a = entropy(site_freq_a_array)
                        entropy_b = entropy(site_freq_b_array)

                        entropy_difference = abs(entropy_a - entropy_b)
                        entropy_differences.append(entropy_difference)

                    split_std_entropy_diff = np.std(entropy_differences)

                    result = (node.name, split_std_entropy_diff, split_std_ratio_topo, split_mean_ratio_branch,
                              split_std_ratio_branch)
                    results.append(result)
            self.logger.info("Finished computing msa split features!")
            return pd.DataFrame(results, columns=["branchId", "std_split_site_entropy_difference",
                                                  "std_split_hop_distance_ratio", "mean_split_branch_length_ratio",
                                                  "std_split_branch_length_ratio"])

    def compute_parsimony_support_features(self, support_path):
        """
        Computes multiple summary statistics over the parsimony support of inner branches.

                Parameters:
                        :param support_path: path to the support file generated by RAxML-NG

                Returns:
                        :return pd.DataFrame: dataframe with the branchId and all features as columns

        """
        with open(support_path, "r") as support_file:
            tree_str = support_file.read()
            tree = Tree(tree_str)
            branch_id_counter = 0
            results = []
            farthest_branch = tree.get_farthest_leaf(topology_only=False)[1]
            for node in tree.traverse():
                branch_id_counter += 1
                node.__setattr__("name", branch_id_counter)
                length = node.dist
                length_relative = length / farthest_branch
                if node.support is not None and not node.is_leaf():
                    num_children = sum([1 for child in node.traverse()])
                    number_nodes = sum([1 for node in tree.traverse()])
                    number_children_relative = num_children / number_nodes
                    childs_inner = [node_child for node_child in node.traverse() if not node_child.is_leaf()]
                    parents_inner = node.get_ancestors()
                    supports_childs = []
                    weighted_supports_childs = []
                    for child in childs_inner:
                        supports_childs.append(child.support)
                        weighted_supports_childs.append(child.support * child.dist)

                    supports_parents = []
                    weighted_supports_parents = []
                    for parent in parents_inner:
                        supports_parents.append(parent.support)
                        weighted_supports_parents.append(parent.support * parent.dist)

                    if len(weighted_supports_childs) >= 1:
                        min_pars_supp_child_w = min(weighted_supports_childs)
                        max_pars_supp_child_w = max(weighted_supports_childs)
                    else:
                        min_pars_supp_child_w = -1
                        max_pars_supp_child_w = -1

                    if len(supports_childs) >= 2:
                        std_pars_supp_child = np.std(supports_childs)
                    else:
                        std_pars_supp_child = -1

                    if len(weighted_supports_childs) >= 1:
                        mean_pars_supp_children_w = statistics.mean(weighted_supports_childs)
                    else:
                        mean_pars_supp_children_w = -1

                    if len(weighted_supports_parents) >= 1:
                        mean_pars_supp_parents_w = statistics.mean(weighted_supports_parents)
                    else:
                        mean_pars_supp_parents_w = -1

                    results.append((node.name, node.support / 100,
                                    min_pars_supp_child_w,
                                    max_pars_supp_child_w, mean_pars_supp_parents_w, length, length_relative,
                                    min(supports_childs), std_pars_supp_child, number_children_relative,
                                    mean_pars_supp_children_w
                                    ))
        self.logger.info("Finished computing Parsimony features from 1000 trees!")
        return pd.DataFrame(results, columns=["branchId", "parsimony_support", "min_pars_support_children_weighted",
                                              "max_pars_support_children_weighted",
                                              "mean_pars_support_parents_weighted", "branch_length",
                                              "norm_branch_length", "min_pars_support_children",
                                              "std_pars_support_children", "number_children_relative",
                                              "mean_pars_support_children_weighted"])

    def compute_parsimony_support(self):
        """
        Calls RAxML-NG for creating 1000 parsimony starting trees, then it computes the support induced by them on the tree at self.tree_filepath

                Parameters:

                Returns:
                        :return string: absolute path to the created support file

        """
        output_prefix = "parsimony_tmp_1000"

        raxml_command = [
            self.raxml_path,
            "--start",
            "--model", self.model_filepath,
            "--tree", "pars{1000}",
            "--msa", self.msa_filepath,
            "--redo",
            "--prefix", output_prefix,
            "--threads", "auto{60}",
            "--log", "ERROR",
        ]

        subprocess.run(raxml_command, shell=False)

        parsimonies_filepath = os.path.join(os.curdir, output_prefix + ".raxml.startTree")
        raxml_command = [self.raxml_path,
                         "--support",
                         "--tree", self.tree_filepath,
                         "--bs-trees", parsimonies_filepath,
                         "--redo",
                         f"--prefix", output_prefix,
                         "--threads", "auto{60}",
                         "--log", "ERROR",
                         ]

        subprocess.run(raxml_command, shell=False)

        return self.compute_parsimony_support_features(
            os.path.abspath(parsimonies_filepath.replace("startTree", "support")))

    def compute_parsimony_bootstrap_support_features(self, support_path):
        """
        Computes multiple features from the parsimony bootstrap support file

                Parameters:
                        :param support_path: absolut path to the support file
                Returns:
                        :return pd.Dataframe: dataframe with branchId and all features as columns

        """
        results = []

        with open(support_path, "r") as support_file:
            tree_str = support_file.read()
            tree = Tree(tree_str)

            all_supports = []
            for node in tree.traverse():
                if not node.is_leaf():
                    all_supports.append(node.support)

            skewness_bootstrap_pars_support_tree = skew(all_supports)

            branch_id_counter = 0
            for node in tree.traverse():
                branch_id_counter += 1
                node.__setattr__("name", branch_id_counter)
                if node.support is not None and not node.is_leaf():

                    parents_inner = node.get_ancestors()
                    childs_inner = [node_child for node_child in node.traverse() if not node_child.is_leaf()]

                    supports_parents = []
                    supports_children_weighted = []
                    supports_children = []
                    for parent in parents_inner:
                        supports_parents.append(parent.support)

                    for child in childs_inner:
                        supports_children_weighted.append(child.support * child.dist)
                        supports_children.append(child.support)

                    if len(supports_children) >= 2:
                        std_pars_bootsupp_child = np.std(supports_children)
                    else:
                        std_pars_bootsupp_child = -1

                    if len(supports_parents) >= 1:
                        mean_pars_bootsupp_parents = statistics.mean(supports_parents)
                        std_pars_bootsupp_parents = np.std(supports_parents)
                    else:
                        mean_pars_bootsupp_parents = -1
                        std_pars_bootsupp_parents = -1

                    results.append(
                        (node.name, node.support / 100, mean_pars_bootsupp_parents, std_pars_bootsupp_parents,
                         skewness_bootstrap_pars_support_tree,
                         min(supports_children_weighted), max(supports_children_weighted), std_pars_bootsupp_child))
        self.logger.info("Finished computing Parsimony Bootstrap features!")

        return pd.DataFrame(results,
                            columns=["branchId", "parsimony_bootstrap_support", "mean_pars_bootstrap_support_parents",
                                     "std_pars_bootstrap_support_parents", "skewness_bootstrap_pars_support_tree",
                                     "min_pars_bootstrap_support_children_w", "max_pars_bootstrap_support_children_w",
                                     "std_pars_bootstrap_support_children"])

    def compute_parsimony_bootstrap_support(self):
        """
        Performs the parsimony bootstrap. Samples alignment sites with replacement 200 times. For each samples MSA
        RAxML-NG is called to infer the corresponding parsimony starting tree. Stores those trees in a final file.
        Calculates the mean Robinson-Foulds distance on the fly.
        Finally computes the support those trees induce on the tree at self.tree_path and the corresponding features.

                Parameters:
                Returns:
                        :return pd.DataFrame: dataframe with branchId and the parsimony bootstrap support features

        """
        alignment = AlignIO.read(self.msa_filepath, "fasta")
        sequence_data = [list(record.seq) for record in alignment]
        alignment_array = np.array(sequence_data)
        original_ids = [record.id for record in alignment]

        trees_path = os.path.join(os.curdir, "parsimony_bootstraps_tmp.txt")

        for x in range(200):
            if (x % 20 == 0) and (x != 0):
                self.logger.info(f"Finished computing {x} from 200 parsimony bootstraps ... ")

            sampled_columns = np.random.choice(alignment_array.shape[1], size=alignment_array.shape[1],
                                               replace=True)

            replicate_alignment = alignment_array[:, sampled_columns]

            seq_records = [SeqRecord.SeqRecord(Seq.Seq(''.join(seq)), id=original_ids[i], description="") for i, seq
                           in
                           enumerate(replicate_alignment)]

            msa_new = AlignIO.MultipleSeqAlignment(seq_records)

            new_msa_path = os.path.join(os.curdir, "parsimony_bootstrap_tmp_" + str(x) + ".fasta")

            output_prefix = "parsimony_bootstrap_tmp_" + str(
                x)

            SeqIO.write(msa_new, new_msa_path, "fasta")

            raxml_command = [
                self.raxml_path,
                "--start",
                "--model", self.model_filepath,
                "--tree", "pars{1}",
                "--msa", new_msa_path,
                "--redo",
                "--prefix", output_prefix,
                "--threads", "auto{60}",
                "--log", "ERROR"
            ]
            try:
                subprocess.run(raxml_command, shell=False)
            except:
                self.logger.error()

            result_tree_path = os.path.join(os.curdir, output_prefix + ".raxml.startTree")

            try:
                with open(os.path.abspath(result_tree_path), 'r') as tree_file:
                    newick_tree = tree_file.read()
            except FileNotFoundError:
                print("boot tree not found")
                continue

            with open(trees_path, 'a') as trees_file:
                if not os.path.exists(os.path.abspath(trees_path)):
                    trees_file.write(newick_tree)
                else:
                    trees_file.write(newick_tree)

        self.logger.info("Finished computing 200 Parsimony Bootstraps")
        raxml_command = [self.raxml_path,
                         "--support",
                         "--tree", self.tree_filepath,
                         f"--bs-trees", trees_path,
                         "--redo",
                         "--prefix", output_prefix,
                         "--log", "ERROR",
                         "--threads", "auto{60}",
                         ]
        subprocess.run(raxml_command, shell=False)

        raxml_command = [self.raxml_path,
                         "--rfdist",
                         "--tree", trees_path,
                         "--redo",
                         "--prefix", output_prefix,
                         "--threads", "auto{60}",
                         "--log", "ERROR",
                         ]
        result = subprocess.run(raxml_command, shell=False)
        rf_dist_path = trees_path.replace(".txt", "_199.raxml.rfDistances").replace("bootstraps", "bootstrap")
        rf_distances = []

        with open(rf_dist_path, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if len(columns) >= 4:
                    try:
                        float_value = float(columns[-1])
                        rf_distances.append(float_value)
                    except ValueError:
                        print(f"Skipping invalid value in line: {line}")

        result = self.compute_parsimony_bootstrap_support_features(
            os.path.join(os.path.abspath(os.curdir), "parsimony_bootstrap_tmp_199.raxml.support"))
        result["mean_norm_rf_distance"] = statistics.mean(rf_distances)

        return result
