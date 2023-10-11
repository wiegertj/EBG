import pickle
import os


def get_model(name):
    with open(os.path.join(os.path.dirname(__file__), name + '.pkl'), 'rb') as model_file:
        model = pickle.load(model_file)
        return model
