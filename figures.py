""" Print figures by looking through all folders in output """

import matplotlib.pyplot as plt
import os

import mush


def find_folder():
    folder = input("Please give folder name for output: ")
    if not os.path.isdir(folder):
        print("Folder name not valid. Folder name changed to {}. Input was {}.".format("./output/",folder))
        folder = "./output/"
    else:
        folder += "/"
        print("Folder for output: {}.".format(folder))
    return folder

def all_figures(folder):
    # figure with statistics
    list_files = os.listdir(folder)
    timesteps = []
    for file in list_files:
        if file[-14:] == "statistics.txt":
            file_stat = folder + "/" + file
        elif file[-9:] == ".timestep":
            timesteps.append(folder + "/" + file)
        elif file[-5:] == ".yaml":
            parameter_file = folder + "/" + file
    mush.fig_stat(file_stat, save=True, output=folder)

if __name__ == "__main__":



    all_figures("output/exponent_1")