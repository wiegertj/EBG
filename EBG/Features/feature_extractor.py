import time
import os

from EBG.Features.feature_computer import FeatureComputer
from EBG.utils.utils import setup_logger, check_file_exists


class FeatureExtractor:
    """
           Coordinates the feature computation.

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
           feature_computer : FeatureComputer
               object that encapsulates all feature computations

           Methods
           -------
           extract_features(self):
               triggers feature computation in the right order, merges result and tracks execution time
           """
    def __init__(self, msa_file_path, tree_file_path, model_file_path, output_prefix, raxml_ng_path, redo):
        self.msa_file_path = msa_file_path
        self.tree_file_path = tree_file_path
        self.model_file_path = model_file_path
        self.feature_computer = FeatureComputer(msa_file_path, tree_file_path, model_file_path, output_prefix, raxml_ng_path, redo)
        self.logger = setup_logger("FeatureExtractor")
        check_file_exists(msa_file_path, "Fasta MSA file", self.logger)
        check_file_exists(tree_file_path, "Newick tree file", self.logger)
        check_file_exists(model_file_path, "RAxML-NG model file", self.logger)
        self.current_directory = os.path.abspath(os.curdir)

    def extract_features(self):
        """
        triggers feature computation in the right order, merges result and tracks execution time

                Parameters:

                Returns:
                        :return pd.DataFrame: dataframe with all features

        """
        self.logger.info("Starting feature extraction ... ")
        start_time = time.time()
        parsimony_features_bootstrap = self.feature_computer.compute_parsimony_bootstrap_support()
        elapsed_time = time.time() - start_time
        self.logger.info(f"Elpased time: {round(elapsed_time, 2)} seconds")
        start_time = time.time()
        split_features_tree = self.feature_computer.compute_split_features_tree()
        elapsed_time = time.time() - start_time
        self.logger.info(f"Elpased time: {round(elapsed_time, 2)} seconds")
        start_time = time.time()
        substitution_features = self.feature_computer.compute_substitution_features()
        elapsed_time = time.time() - start_time
        self.logger.info(f"Elpased time: {round(elapsed_time, 2)} seconds")
        start_time = time.time()
        parsimony_features = self.feature_computer.compute_parsimony_support()
        elapsed_time = time.time() - start_time
        self.logger.info(f"Elpased time: {round(elapsed_time, 2)} seconds")
        merged_df = split_features_tree.merge(parsimony_features, on="branchId", how="inner")
        merged_df = merged_df.merge(parsimony_features_bootstrap, on="branchId", how="inner")

        for key, value in substitution_features.items():
            merged_df[key] = value
        self.logger.info("Feature extraction finished!")

        return merged_df
