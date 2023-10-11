import unittest
import EBG.utils.utils


class FeatureComputerTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_check_raxml(self):
        self.assertEqual(EBG.utils.utils.check_raxml_availability(), "/usr/local/bin/raxml-ng")

    def test_format_predictions(self):
        raw_predictions = [55.3, 22.1, 1.00123, 103]
        formatted_predictions = EBG.utils.utils.format_predictions(raw_predictions)
        self.assertEqual(formatted_predictions, [55, 22, 1, 100])

    def test_file_check(self):
        with self.assertRaises(FileNotFoundError):
            EBG.utils.utils.check_file_exists("no_file.fasta", "msa file")


if __name__ == '__main__':
    unittest.main()