import argparse
import time
import os

os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
from EBG.Prediction.predictor import Predictor


def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Educated Bootstrap Guesser")
    parser.add_argument("-msa", required=True, type=str, help="PATH\t\tabsolute path to the MSA file in .fasta format ")
    parser.add_argument("-tree", required=True, type=str, help="PATH\t\tabsolute path to the tree file in newick format (.bestTree)")
    parser.add_argument("-model", required=True, type=str, help="PATH\t\tabsolute path to the raxml-ng model file (.bestModel)")
    parser.add_argument("-o", required=False, type=str, help="VALUE\t\toutput folder name, overwrites all files if already existing\tdefault: EBG_output")
    parser.add_argument("-t", required=False, type=str, default="b",
                        help="c | r | b\t\ttype of prediction: (c)lassification, (r)egression, (b)oth\tdefault: b")
    parser.add_argument("-raxmlng", required=False, type=str, default="raxml-ng",
                        help="PATH\t\tpath to RAxML-NG\tdefault: raxml-ng")
    parser.add_argument("-redo", action='store_true', required=False,
                        help="if set, all existing results will be overwritten")

    args = parser.parse_args()
    predictor = Predictor(args.msa, args.tree, args.model, args.o, args.t, args.raxmlng, args.redo)
    predictor.predict()
    predictor.map_results()
    total_time = time.time() - start_time
    print(f"Total elapsed time: {round(total_time, 2)} seconds")


if __name__ == "__main__":
    main()
