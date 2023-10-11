import math
import os
import unittest

from EBG.Features.feature_computer import FeatureComputer

test_msa_path = os.path.abspath(os.path.curdir + "/data/test.fasta")
test_model_path = os.path.abspath(os.path.curdir + "/data/test.bestModel")
test_tree_path = os.path.abspath(os.path.curdir + "/data/test.bestTree")


class FeatureComputerTests(unittest.TestCase):
    def setUp(self):
        os.chdir(os.path.dirname(__file__))
        self.feature_computer = FeatureComputer(msa_filepath=test_msa_path, model_filepath=test_model_path,
                                                tree_filepath=test_tree_path, output_prefix="test",
                                                raxml_ng_path="raxml-ng", redo=True)

    def tearDown(self):
        pass

    def test_substitution_features(self):

        substitution_features = self.feature_computer.compute_substitution_features()

        self.assertEqual(len(substitution_features), 4)

        for key, value in substitution_features.items():
            self.assertTrue(isinstance(value, float))

        self.assertTrue(math.isclose(substitution_features.get('max_substitution_frequency', 0), 0.2941, abs_tol=1e-3))
        self.assertTrue(math.isclose(substitution_features.get('mean_substitution_frequency', 0), 0.0363, abs_tol=1e-3))
        self.assertTrue(math.isclose(substitution_features.get('cv_substitution_frequency', 0), 1.6327, abs_tol=1e-3))
        self.assertTrue(math.isclose(substitution_features.get('skw_substitution_frequency', 0), 1.7165, abs_tol=1e-3))

    def test_split_features_tree(self):

        split_features_tree = self.feature_computer.compute_split_features_tree()
        self.assertEqual(split_features_tree.shape, (15, 3))

        self.assertTrue((0 <= split_features_tree['branch_length_ratio_split']).all() and (
                split_features_tree['branch_length_ratio_split'] <= 1).all())
        self.assertTrue((0 <= split_features_tree['mean_closeness_centrality_ratio']).all() and (
                split_features_tree['mean_closeness_centrality_ratio'] <= 1).all())

    def test_parsimony_support(self):

        self.feature_computer.compute_parsimony_support()
        self.assertTrue(os.path.exists(os.path.abspath(os.path.curdir) + "/parsimony_tmp_1000.raxml.support"))
        self.assertTrue(os.path.exists(os.path.abspath(os.path.curdir) + "/parsimony_tmp_1000.raxml.startTree"))

        try:
            with open(os.path.abspath(os.path.curdir) + "/parsimony_tmp_1000.raxml.startTree", 'r') as file:
                lines = file.readlines()
                self.assertEqual(len(lines), 1000)
        except AssertionError as e:
            print(f"Assertion failed: {e}")
        finally:
            os.remove(os.path.abspath(os.path.curdir) + "/parsimony_tmp_1000.raxml.support")
            os.remove(os.path.abspath(os.path.curdir) + "/parsimony_tmp_1000.raxml.startTree")
            os.remove(os.path.abspath(os.path.curdir) + "/parsimony_tmp_1000.raxml.log")
            os.remove(os.path.abspath(os.path.curdir) + "/parsimony_tmp_1000.raxml.rba")

    def test_parsimony_bootstrap(self):

        parsimony_bootstrap_features = self.feature_computer.compute_parsimony_bootstrap_support()
        self.assertEqual(parsimony_bootstrap_features.shape, (15, 9))
        self.assertTrue((0 <= parsimony_bootstrap_features['mean_norm_rf_distance']).all() and (
                parsimony_bootstrap_features['mean_norm_rf_distance'] <= 1).all())
        self.assertTrue((-1 <= parsimony_bootstrap_features['mean_pars_bootstrap_support_parents']).all() and (
                parsimony_bootstrap_features['mean_pars_bootstrap_support_parents'] <= 100).all())
        self.assertTrue((0 <= parsimony_bootstrap_features['min_pars_bootstrap_support_children_w']).all().all())
        self.assertTrue((0 <= parsimony_bootstrap_features['max_pars_bootstrap_support_children_w']).all().all())


if __name__ == '__main__':
    unittest.main()
