import logging
import os
import warnings
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
from ete3 import Tree
from scipy.stats import entropy
from EBG.Features.feature_extractor import FeatureExtractor
from EBG.utils.utils import setup_logger, format_predictions, get_model


class Predictor:
    """
        Class for making predictions and writing the results.

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
        prediction_model_median : LightGBM.Booster
            model predicting the median of the bootstrap
        prediction_model_lower5 : LightGBM.Booster
            model predicting the 5% lower bound with quantile regression
        prediction_model_lower10 : LightGBM.Booster
            model predicting the 10% lower bound with quantile regression
        classifier_70 : LightGBM.Booster
            model for predicting if the SBS value of a branch exceeds 70
        classifier_75 : LightGBM.Booster
            model for predicting if the SBS value of a branch exceeds 75
        classifier_80 : LightGBM.Booster
            model for predicting if the SBS value of a branch exceeds 80
        classifier_85 : LightGBM.Booster
            model for predicting if the SBS value of a branch exceeds 85
        prediction : pd.DataFrame
            dataframe with all prediction results
        t : str
            type of prediction to perform

    """
    def __init__(self, msa_filepath, tree_filepath, model_filepath, o="EBG_output", t="b", raxml_ng_path="raxml-ng", redo=False):
        self.current_directory = os.path.abspath(os.curdir)
        self.logger = setup_logger("Predictor")
        self.tree_filepath = tree_filepath
        self.prediction_model_median = get_model("median_model")
        self.prediction_model_lower5 = get_model("low_model_5")
        self.prediction_model_lower10 = get_model("low_model_10")
        self.classifier_70 = get_model("class_70")
        self.classifier_75 = get_model("class_75")
        self.classifier_80 = get_model("class_80")
        self.classifier_85 = get_model("class_85")
        self.feature_extractor = FeatureExtractor(msa_filepath, tree_filepath, model_filepath, o, raxml_ng_path, redo)
        self.output_prefix = o
        self.prediction = None
        self.type = t

    def predict(self):
        features = self.feature_extractor.extract_features()

        features_pred = features[
            ["branchId", "parsimony_bootstrap_support", "parsimony_support", "mean_substitution_frequency",
             "norm_branch_length",
             "branch_length", "mean_norm_rf_distance", "max_substitution_frequency",
             "skewness_bootstrap_pars_support_tree", "cv_substitution_frequency",
             "branch_length_ratio_split", "max_pars_bootstrap_support_children_w", "skw_substitution_frequency",
             "mean_pars_bootstrap_support_parents",
             "max_pars_support_children_weighted", "std_pars_bootstrap_support_parents", "min_pars_support_children",
             "min_pars_support_children_weighted",
             "number_children_relative", "mean_pars_support_children_weighted", "std_pars_bootstrap_support_children",
             "mean_closeness_centrality_ratio", "min_pars_bootstrap_support_children_w"]]

        if self.type == "b" or self.type == "r":
            prediction_lower5 = format_predictions(
                self.prediction_model_lower5.predict(features_pred, num_threads=1) * 100)
            prediction_lower10 = format_predictions(
                self.prediction_model_lower10.predict(features_pred, num_threads=1) * 100)
            median_prediction = format_predictions(
                self.prediction_model_median.predict(features_pred, num_threads=1) * 100)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                features["prediction_lower5"] = prediction_lower5
                features["prediction_lower10"] = prediction_lower10
                features["prediction_median"] = median_prediction
                features["prediction_median_lower5_distance"] = abs(
                    features["prediction_median"] - features["prediction_lower5"])
                features["prediction_median_lower10_distance"] = abs(
                    features["prediction_median"] - features["prediction_lower10"])
                features_pred["median_pred"] = [x / 100 for x in median_prediction]
                features_pred["lower_bound_10"] = [x / 100 for x in prediction_lower10]
                features_pred["lower_bound_5"] = [x / 100 for x in prediction_lower5]

        if self.type == "c" or self.type == "b":
            classification_70 = self.classifier_70.predict(features_pred, num_threads=1)
            classification_75 = self.classifier_75.predict(features_pred, num_threads=1)
            classification_80 = self.classifier_80.predict(features_pred, num_threads=1)
            classification_85 = self.classifier_85.predict(features_pred, num_threads=1)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                features["prediction_bs_over_70"] = classification_70
                features["prediction_bs_over_75"] = classification_75
                features["prediction_bs_over_80"] = classification_80
                features["prediction_bs_over_85"] = classification_85

                for index, row in features.iterrows():
                    entropy_row_80 = entropy([row["prediction_bs_over_80"], 1 - row["prediction_bs_over_80"]],
                                             base=2)
                    entropy_row_85 = entropy([row["prediction_bs_over_85"], 1 - row["prediction_bs_over_85"]],
                                             base=2)
                    entropy_row_75 = entropy([row["prediction_bs_over_75"], 1 - row["prediction_bs_over_75"]],
                                             base=2)
                    entropy_row_70 = entropy([row["prediction_bs_over_70"], 1 - row["prediction_bs_over_70"]],
                                             base=2)
                    features.loc[index, 'prediction_uncertainty_bs_over_80'] = entropy_row_80
                    features.loc[index, 'prediction_uncertainty_bs_over_85'] = entropy_row_85
                    features.loc[index, 'prediction_uncertainty_bs_over_75'] = entropy_row_75
                    features.loc[index, 'prediction_uncertainty_bs_over_70'] = entropy_row_70

        self.logger.info("Finished prediction!")

        features["dataset"] = self.output_prefix
        features.to_csv(os.path.abspath(os.path.join("..", self.output_prefix + "_features.csv")))

        self.prediction = features

    def map_results(self):
        classification_results = ["bs_over_70", "bs_over_75", "bs_over_80", "bs_over_85", "uncertainty_bs_over_70", "uncertainty_bs_over_75", "uncertainty_bs_over_80", "uncertainty_bs_over_85"]
        regression_results = ["lower5", "lower10", "median"]
        map_trees = []
        if self.type == "r" or self.type == "b":
            map_trees.extend(regression_results)
        if self.type == "c" or self.type == "b":
            map_trees.extend(classification_results)

        with open(self.tree_filepath, "r") as tree_file:
            tree_str = tree_file.read()
            phylo_tree_reference = Tree(tree_str)
            branch_id_counter = 0
            for node in phylo_tree_reference.traverse():
                branch_id_counter += 1
                if not node.is_leaf():
                    node.__setattr__("name", branch_id_counter)

        for value in map_trees:
            current_tree = phylo_tree_reference.copy()
            if value in regression_results:
                for node in current_tree.traverse():
                    if not node.is_leaf():
                        matching_row = self.prediction[self.prediction['branchId'] == node.name]
                        if not matching_row.empty:
                            node.__setattr__("support", matching_row.iloc[0]['prediction_' + value])
                current_tree.write(outfile=os.path.join(self.current_directory, self.output_prefix,
                                                        self.output_prefix + "_" + value + "_support_prediction.newick"),
                                   format=0)
            else:
                for node in current_tree.traverse():
                    if not node.is_leaf():
                        matching_row = self.prediction[self.prediction['branchId'] == node.name]
                        if not matching_row.empty:
                            target = matching_row.iloc[0]['prediction_' + value]
                            node.__setattr__("support", target)
                current_tree.write(outfile=os.path.join(self.current_directory, self.output_prefix,
                                                        self.output_prefix + "_" + value + "_support_prediction.newick"),
                                   format=0)

        self.logger.info(f"Stored trees with prediction result at: " + os.path.abspath(os.path.join(self.current_directory, self.output_prefix)))


