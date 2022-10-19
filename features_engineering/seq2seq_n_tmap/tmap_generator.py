import tmap as tm
import numpy as np
from matplotlib import pyplot as plt
import json
import os

from PIL import Image
import base64

import pickle

from faerun import Faerun
from io import BytesIO

# input
file_name = "tmap_input_rct2_500n" #output file name
pickle_input_path = "/Users/jacky/Desktop/Current_projects/ML_CHEM/Tmap/tmap_input_rct2_500n"
image_input_path = "images"


images = []
image_labels = []

# open input_file
input_file = open(pickle_input_path,'rb')
input_dict = pickle.load(input_file)
n, weights, edge_list, label_list  = input_dict["n_limit"], input_dict["weights"], input_dict["edge_list"], input_dict["label_list"]


def generate_tmap(n, weights, edge_list, label_list = None, image_input_path = None):
    x, y, s, t, _ = tm.layout_from_edge_list(n, edge_list, create_mst=False)
    x_mst, y_mst, s_mst, t_mst, _ = tm.layout_from_edge_list(
        n, edge_list, create_mst=True
    )

    # show image
    if image_input_path:
        images = []
        image_labels = []
        for file in sorted(os.listdir(image_input_path)):
            encoded = base64.b64encode(open("images/" + file, "rb").read()).decode()
            image_labels.append( 'data:image/png;base64,{}'.format(encoded) )

    # all points are mono-colored
    c = np.zeros((len(x_mst),), dtype=int)

    # create faerun graph
    faerun = Faerun(clear_color="#111111", view="front", coords=False)

    if image_input_path:
        label_input = image_labels
    else:
        label_input = c

    # add data
    faerun.add_scatter(
        "CHEM_Tmap",
        {"x": x_mst, "y": y_mst,"c": c, "labels":label_input}, #, "c": labels, "labels": image_labels},
        colormap="tab20",
        shader="smoothCircle",
        point_scale=2.5,
        max_point_size=10,
        has_legend=True,
        categorical=True,
    )
    faerun.add_tree(
        file_name, {"from": s_mst, "to": t_mst}, point_helper="CHEM_Tmap", color="#666666"
    )

    faerun.plot(file_name, template="url_image")


if __name__ == "__main__":
    generate_tmap(n, weights, edge_list, label_list, image_input_path)
