""" Print figures by looking through all folders in output """

import matplotlib.pyplot as plt
import numpy as np
import os
import re
from operator import itemgetter
import pandas as pd
import yaml

import mush

regex = re.compile(r'\d+')
def find_float(filename):
    numbers = regex.findall(filename) # gives 2 numbers as number with format "1.05" returns 1 and 05
    return float(numbers[0] + '.' + numbers[1])


def find_folder():
    folder = input("Please give folder name for output: ")
    if not os.path.isdir(folder):
        print("Folder name not valid. Folder name changed to {}. Input was {}.".format("./output/",folder))
        folder = "./output/"
    else:
        folder += "/"
        print("Folder for output: {}.".format(folder))
    return folder


def fig_stat(filename, save=False, output="./", print_all=True, print_list=[]):
    data = pd.read_csv(filename, sep=" ", index_col=False)
    if print_all:
        n_col = data.shape[1]
        fig, ax = plt.subplots(n_col-2, 2, figsize=[6, n_col*3]) #first column with iteration as x axis, 2nd column with time
        names = data.columns
        for i in range(n_col-2):
            ax[i, 0].plot(data['iteration_number'], data[names[i+2]])
            ax[i, 1].plot(data['time'], data[names[i+2]])
            ax[i, 0].set_ylabel(names[i+2])
        ax[n_col-3,0].set_xlabel("iteration_number")
        ax[n_col-3,1].set_xlabel("time")
    plt.savefig(output + filename[:-4]+".pdf")


def make_figure(filename, save=False, output="./", max_r=None):
    data = pd.read_csv(filename, sep=" ") # , index_col=False)
    fig, ax = plt.subplots(1,2, sharey=True)
    dr = data["radius"][1] - data["radius"][0]
    ax[0].plot(data["porosity"], data["radius"] + dr / 2.)
    ax[1].plot(data["velocity"], data["radius"] + dr)
    ax[0].set_xlabel("Porosity")
    ax[0].set_ylabel("Height (non-dim)")
    ax[1].set_xlabel("Solid velocity (non-dim)")
    if not max_r == None:
        ax[0].set_ylim([0., max_r])
        ax[1].set_ylim([0., max_r])
    if save:
        plt.savefig(output + filename[:-4] + '.pdf') # -4 to remove the .csv
        plt.close(fig)
    else: plt.show()


def all_lines(output_folder, save=False):
    list_files = os.listdir(output_folder)
    fig, ax = plt.subplots(1,2, sharey=True)
    ax[0].set_xlabel("Porosity")
    ax[0].set_ylabel("Height (non-dim)")
    ax[1].set_xlabel("Solid velocity (non-dim)")
    for file in list_files:
        #print(file[-9:])
        if file[-9:] == ".timestep":
            print(output_folder+file)
            data = pd.read_csv(output_folder +'/' + file, sep=" ")
            dr = data["radius"][1]-data["radius"][0]
            time = file[-14:-9]
            ax[0].plot(data["porosity"], data["radius"] + dr / 2., label=time)
            ax[1].plot(data["velocity"], data["radius"] + dr, label=time)
    #ax[1].legend(loc="upper left", bbox_to_anchor=(1,1))
    #plt.tight_layout()
    if save:
        plt.savefig(output_folder + "/all_fig"+ '.pdf') # -4 to remove the .csv
        plt.close(fig)
    else: plt.show()


def all_figures(folder):
    # figure with statistics
    list_files = os.listdir(folder)
    timesteps = {}
    for file in list_files:
        if file[-14:] == "statistics.txt":
            file_stat = folder + "/" + file
        elif file[-9:] == ".timestep":
            _name = folder + "/" + file
            _time = find_float(file)
            timesteps[_name] = _time
        elif file[-5:] == ".yaml":
            parameter_file = folder + "/" + file
    with open(parameter_file, 'r') as stream:
        try:
            options = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    print("Figures with global statistics: {}".format(file_stat))
    fig_stat(file_stat, save=True)
    # preparation for the figure
    fig, ax = plt.subplots(1,2, sharey=True)
    ax[0].set_xlabel("Porosity")
    ax[0].set_ylabel("Height (non-dim)")
    ax[1].set_xlabel("Solid velocity (non-dim)")
    i = 0
    n_i = 1 # if too many figures
    n = int(len(timesteps)/n_i)
    colors = plt.cm.viridis(np.linspace(0,1,n))
    for name, time in sorted(timesteps.items(), key = itemgetter(1)):
        print(name, time)
        # single figure 
        make_figure(name, save=True, output="./", max_r=options["Ric_adim"])
        # figure with all timesteps
        data = pd.read_csv(name, sep=" ")
        dr = data["radius"][1]-data["radius"][0]
        if i%n_i ==0: 
            ax[0].plot(data["porosity"], data["radius"] + dr / 2.,color=colors[i], label=time)
            ax[1].plot(data["velocity"], data["radius"] + dr, color=colors[i], label=time)
        i += 1
    plt.savefig(folder+"/all_figs.pdf")

if __name__ == "__main__":


    all_figures("output/exp_1_coeff__02")


    # snippet for ordering dictionnary and print values.
    # from operator import itemgetter  
    # for key, value in sorted(d.items(), key = itemgetter(1)):
    #     print(key, value)
