""" Print figures by looking through all folders in output """

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import re
from operator import itemgetter
import pandas as pd
import yaml
from scipy import optimize

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


def fig_stat(filename, save=False, output="", print_all=True, print_list=[]):
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
        if names[i+2] == "thickness_boundary":
            try:
                    maximum = data["thickness_boundary"].iloc[100]
                    maximum = data["thickness_boundary"].iloc[-1]*10
                    print(maximum)
            except Exception as e:
                    maximum = data["thickness_boundary"].iloc[-1]*10
            ax[i,0].set_ylim([0, maximum])
            ax[i,1].set_ylim([0, maximum])
                # ax[i,0].set_ylim([0, 10*data["radius"].iloc[-1]])
                # ax[i,1].set_ylim([0, 10*data["radius"].iloc[-1]])
        plt.savefig(output + filename[:-4]+".pdf")
    else:
        n_col = len(print_list)
        fig, ax = plt.subplots(n_col, 2, figsize=[6, n_col*4])
        for name in print_list:
            ax[i, 0].plot(data['iteration_number'], data[name])
            ax[i, 1].plot(data['time'], data[name])
            ax[i, 0].set_ylabel(name)
            if name == "thickness_boundary":
                try:
                    maximum = data["radius"].iloc[100]
                except Exception as e:
                    maximum = data["radius"].iloc[-1]
                ax[i,0].set_ylim([0, maximum])
                ax[i,1].set_ylim([0, maximum])


def make_figure(filename, save=False, output="", max_r=None):
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
    #fig, ax = plt.subplots(1,2, sharey=True)
    fig2, ax2 = plt.subplots()
    ax2.set_xlabel("timestep_print")
    ax2.set_ylabel("radius")
    #ax[0].set_xlabel("Porosity")
    #ax[0].set_ylabel("Height (non-dim)")
    #ax[1].set_xlabel("Solid velocity (non-dim)")
    for i, file in enumerate(list_files):
        #print(file[-9:])
        if file[-9:] == ".timestep":
            print(output_folder+file)
            data = pd.read_csv(output_folder +'/' + file, sep=" ")
            dr = data["radius"][1]-data["radius"][0]
            time = file[-14:-9]
            #ax[0].plot(data["porosity"], data["radius"] + dr / 2., label=time)
            #ax[1].plot(data["velocity"], data["radius"] + dr, label=time)
            ax2.scatter(i*np.ones(data.shape[0]), data["radius"] + dr / 2., c=data["porosity"], cmap=plt.cm.viridis)
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
    fig, ax = plt.subplots(2,1,figsize=[5,3])
    fig2, ax2 = plt.subplots()
    ax2.set_xlabel("timestep_print")
    ax2.set_ylabel("radius")
    i = 0
    n_i = 10 # if too many figures
    n = int(len(timesteps))
    colors = plt.cm.viridis(np.linspace(0,1,int(n)+1))
    # print(timesteps)


    name_max = max(timesteps.items(), key=itemgetter(1))[0]
    #print(name_max)
    data = pd.read_csv(name_max, sep=" ")
    dr = data["radius"][1]-data["radius"][0]
    Radius = np.array(data["radius"].values)+ dr / 2.
    #print(Radius)
    Time = np.array(sorted(timesteps.values()))
    #print((Time))
    X, Y = np.meshgrid(Time, Radius)
    Z = np.ones_like(X)
    print(Z.shape)

    for i, (name, time) in enumerate(sorted(timesteps.items(), key = itemgetter(1))):
        print(name, time)
        # single figure
        # make_figure(name, save=True, output="", max_r=options["Ric_adim"])
        # figure with all timesteps
        data = pd.read_csv(name, sep=" ")

        Porosity = np.array(data["porosity"].values)
        N_r = len(Porosity)
        print(N_r)
        Z[:N_r, i] = Porosity


        dr = data["radius"][1]-data["radius"][0]
        i += 1
        if i%n_i ==0:
            ax[0].plot(data["radius"] + dr / 2.,data["porosity"], color=colors[i], label=time)
            ax[1].plot(data["radius"] + dr,data["velocity"], color=colors[i], label=time)
            #ax2.scatter(i*np.ones(data.shape[0]), data["radius"] + dr / 2., c=data["porosity"], cmap=plt.cm.viridis)
        
    levels = np.linspace(0, 0.4, 100) # [0., 0.1, 0.2, 0.3, 0.4]
    sc = ax2.contourf(X, Y, Z, levels=levels, extend="max")
    sc2 = ax2.contour(X, Y, Z>=0.4, 1,colors='k')
    cbar = plt.colorbar(sc)
    cbar.set_label("Porosity")
    cbar.set_clim([0., 0.4])
    #plt.figure(fig.number)
    #plt.savefig(folder+"/all_figs.pdf")
    ax[0].set_xlim([0,10])
    ax[0].set_ylim([0,0.4])
    ax[0].set_xlabel("Radius")
    ax[0].set_ylabel("Porosity")
    plt.figure(fig2.number)
    plt.savefig(folder+"/all_figs_test.pdf")
    plt.figure(fig.number)
    plt.tight_layout()
    plt.savefig(folder+"/all_figs.pdf")



def fig_thickness(folder_main):

        fig, ax = plt.subplots(2, 1, sharex=True)

        size = np.array([1, 10, 20, 50, 100, 200])
        symbols = ["o", "v", "^", "s", "*", "d"]
        cm = plt.cm.magma(size/200)
        list_subfolder = os.listdir(folder_main)
        for subfolder_name in list_subfolder:
            list_files = os.listdir(folder_main+"/"+subfolder_name)
            for file in list_files:
                if file[-14:] == "statistics.txt":
                    file_stat = folder_main + "/" + subfolder_name + "/" + file
                if file[-5:] == ".yaml":
                    with open(folder_main + "/" + subfolder_name + "/" + file, 'r') as stream:
                        try:
                            param = yaml.safe_load(stream)
                            #print(param)
                        except yaml.YAMLError as exc:
                            print(exc)
            data = pd.read_csv(file_stat, sep=" ", index_col=False)
            index_size = np.argmin((size-param["Ric_adim"])**2)
            ax[0].scatter(param['coeff_velocity'], data["thickness_boundary"].iloc[-1]/param["Ric_adim"], c=param["Ric_adim"], linewidths=0, vmin=1, vmax=200, marker=symbols[index_size])            
            #print(param['coeff_velocity'], data["thickness_boundary"].iloc[-1]/param["Ric_adim"])
            ax[1].scatter(param['coeff_velocity'], data["thickness_boundary"].iloc[-1], c=param["Ric_adim"], linewidths=0, vmin=1, vmax=200, marker=symbols[index_size])
        #ax[0].set_ylim([0, 3])
        #ax[]
        ax[1].set_yscale("log")
        ax[1].set_xscale("log")
        ax[0].set_yscale("log")
        ax[0].set_xscale("log")
        ax[1].set_ylim([1e-4, 1e4])
        ax[0].set_ylim([1e-4, 1e4])
        ax[0].set_ylabel("$\delta$/$R_{ic}$")
        ax[1].set_ylabel("$\delta$")
        ax[1].set_xlabel("coeff growth velocity")
        plt.tight_layout()
        #ax[0].legend()
        #plt.show()


def fig_porosity(folder_main):

        fig, ax = plt.subplots(2,1, figsize=[4,6])

        size = np.array([1, 10, 20, 50, 100, 200])
        growth = np.array([1., 0., -1., -2.])
        cm = plt.cm.viridis((growth+2)/4)
        list_subfolder = os.listdir(folder_main)
        
        for subfolder_name in list_subfolder:
            list_files = os.listdir(folder_main+"/"+subfolder_name)
            for file in list_files:
                if file[-14:] == "statistics.txt":
                    file_stat = folder_main + "/" + subfolder_name + "/" + file
                if file[-5:] == ".yaml":
                    with open(folder_main + "/" + subfolder_name + "/" + file, 'r') as stream:
                        try:
                            param = yaml.safe_load(stream)
                            #print(param)
                        except yaml.YAMLError as exc:
                            print(exc)
            data = pd.read_csv(file_stat, sep=" ", index_col=False)
            index = np.argmin((10**growth-param["coeff_velocity"])**2)
            if param["Ric_adim"] == 100:
                marker = "--"
            else: marker = "-"
            if param["coeff_velocity"] == 1:
                label = r"$R_{ic}$ = "+"{}".format(param["Ric_adim"])
            else: label = ""
            ax[0].plot(data["radius"]/param["Ric_adim"], data["sum_phi"], marker,   color=cm[index], label=label)
            if param["Ric_adim"] == 10:
                label = r"$\dot{R}_{ic}$ = "+ "{}".format(param["coeff_velocity"])
            else: label = ""
            ax[1].plot(data["radius"]/param["Ric_adim"], data["thickness_boundary"], marker,  color=cm[index], label=label)
            
            #print(param['coeff_velocity'], data["thickness_boundary"].iloc[-1]/param["Ric_adim"])
            #ax[1].scatter(param['coeff_velocity'], data["sum_phi"].iloc[-1], c=param["Ric_adim"], linewidths=0, vmin=1, vmax=200, marker=symbols[index_size])
        #ax[0].set_ylim([0, 3])
        #ax[]
        ax[1].set_yscale("log")
        #ax[1].set_xscale("log")
        ax[1].set_ylim([1e-1, 1.e4])
        ax[1].set_xlim([0, 1])
        #ax[1].set_xlim([0.9e-2, 1])
        ax[0].legend()
        ax[1].legend()
        ax[1].set_xlabel("Radius/$R_{ic}$")
        ax[0].set_ylabel("<$\phi$>")
        ax[1].set_ylabel("$\delta$")
        plt.tight_layout()
        plt.savefig("phi_delta.pdf")
        #plt.show()



def fig_porosity_thickness(folder_main):
        columns = ["Ric_adim", "coeff_velocity", "exp", "sum_phi", "delta", "remarks"]
        df = pd.DataFrame(columns=columns)

        def add_value(df, ric, coeff, exp, phi, delta, remarks=""):
            df_add = pd.DataFrame({"Ric_adim":[ric], "coeff_velocity":[coeff], "exp":[exp], "sum_phi":[phi], "delta":[delta], "remarks":[remarks]})
            df = df.append(df_add)
            return df

        list_subfolder = os.listdir(folder_main)
        for subfolder_name in list_subfolder:
            list_files = os.listdir(folder_main+"/"+subfolder_name)
            for file in list_files:
                if file[-14:] == "statistics.txt":
                    file_stat = folder_main + "/" + subfolder_name + "/" + file
                if file[-5:] == ".yaml":
                    with open(folder_main + "/" + subfolder_name + "/" + file, 'r') as stream:
                        try:
                            param = yaml.safe_load(stream)
                            #print(param)
                        except yaml.YAMLError as exc:
                            print(exc)
            data = pd.read_csv(file_stat, sep=" ", index_col=False)
            if data["radius"].iloc[-1] < 0.99*param["Ric_adim"]:
                remarks = "run ended before completion. Radius {}/{}".format(data["radius"], param["Ric_adim"])
            else:
                remarks = ""
            df = add_value(df, param["Ric_adim"], param['coeff_velocity'], param['growth_rate_exponent'], 
                        data["sum_phi"].iloc[-1], data["thickness_boundary"].iloc[-1], remarks)

        df.to_csv("data.csv")

        fig, ax = plt.subplots(3, 1, sharex=True, figsize=[6, 6])
        cmap = plt.cm.viridis
        sc = ax[0].scatter(df["coeff_velocity"], df["delta"], c=df["Ric_adim"], marker='+', cmap=cmap)
        ax[1].scatter(df["coeff_velocity"], df["delta"]/df["Ric_adim"], c=df["Ric_adim"], marker='+', cmap=cmap)
        ax[2].scatter(df["coeff_velocity"], df["sum_phi"], c=df["Ric_adim"], marker='+', cmap=cmap)
        ax[0].set_yscale("log")
        ax[0].set_xscale("log")
        ax[1].set_yscale("log")
        ax[1].set_xscale("log")
        ax[2].set_yscale("log")
        ax[2].set_xscale("log")
        ax[0].set_ylim([1e-2, 1e4])
        ax[1].set_ylim([1e-4, 1e2])
        ax[2].set_ylim([1e-3, 0.5])
        ax[0].set_xlim([0.9e-3, 1.1e1])
        ax[0].set_ylabel("$\delta$")
        ax[1].set_ylabel("$\delta$/$R_{ic}$")
        ax[2].set_ylabel("<$\phi$>")
        ax[2].set_xlabel("Growth rate")
        #fig.subplots_adjust(right=0.8)
        #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cbar = plt.colorbar(sc, ax=ax.ravel().tolist())
        cbar.set_label("$R_{ic}$")
        #plt.tight_layout()
        plt.savefig("scaling_laws.pdf")

        fig, ax = plt.subplots(2, 1, sharex=True, figsize=[6, 4])
        cmap = plt.cm.viridis
        sc = ax[0].scatter(df["coeff_velocity"], df["delta"], c=df["Ric_adim"], marker='+', cmap=cmap)
        ax[1].scatter(df["coeff_velocity"], df["delta"]/df["Ric_adim"], c=df["Ric_adim"], marker='+', cmap=cmap)
        #ax[2].scatter(df["coeff_velocity"], df["sum_phi"], c=df["Ric_adim"], marker='+', cmap=cmap)
        ax[0].set_yscale("log")
        ax[0].set_xscale("log")
        ax[1].set_yscale("log")
        ax[1].set_xscale("log")
        #ax[2].set_yscale("log")
        #ax[2].set_xscale("log")
        ax[0].set_ylim([1e-2, 1e4])
        ax[1].set_ylim([1e-4, 1e2])
        #ax[2].set_ylim([1e-3, 0.5])
        ax[0].set_xlim([0.2e-3, 1.1e-1])
        ax[0].set_ylabel("$\delta$")
        ax[1].set_ylabel("$\delta$/$R_{ic}$")
        #ax[2].set_ylabel("<$\phi$>")
        #ax[2].set_xlabel("Growth rate")
        #fig.subplots_adjust(right=0.8)
        #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cbar = plt.colorbar(sc, ax=ax.ravel().tolist())
        cbar.set_label("$R_{ic}$")
        #plt.tight_layout()
        plt.savefig("scaling_laws_2lines.pdf")


def scaling_laws(file="data.csv"):
        df = pd.read_csv(file)

        small_v = df[df["coeff_velocity"]<3e-2]
        large_v = df[df["coeff_velocity"]>1e0]
        fig, ax = plt.subplots(1,2)
        ax[0].scatter(small_v["coeff_velocity"], small_v["delta"])
        ax[1].scatter(large_v["coeff_velocity"], large_v["delta"]/large_v["Ric_adim"])
        ax[0].set_yscale("log")
        ax[0].set_xscale("log")
        ax[1].set_yscale("log")
        ax[1].set_xscale("log")
        ax[0].set_xlim([1e-3, 1e-1])
        ax[0].set_ylim([1e-2, 1e0])
        ax[1].set_xlim([1e0, 1e2])
        ax[1].set_ylim([1e0, 1e2])

        powerlaw = lambda x, amp, index: amp * (x**index)

        def power_law_fit(xdata, ydata):
            logx = np.log10(xdata)
            logy = np.log10(ydata)
            logyerr = 0.2
            # define our (line) fitting function
            fitfunc = lambda p, x: p[0] + p[1] * x
            errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
            pinit = [1.0, -1.0]
            out = optimize.leastsq(errfunc, pinit,
                                args=(logx, logy, logyerr), full_output=1)
            pfinal = out[0]
            covar = out[1]

            index = pfinal[1]
            amp = 10.0**pfinal[0]

            return pfinal, covar, index, amp

        #large = power_law_fit(large_v["coeff_velocity"].values, large_v["delta"].values)
        small = power_law_fit(small_v["coeff_velocity"].values, small_v["delta"].values)

        fig, ax = plt.subplots(2, 1)
        #xdata = large_v["coeff_velocity"].values
        #ydata = large_v["delta"]/large_v["Ric_adim"]
        #ydata = ydata.values
        #pfinal, covar, index, amp = power_law_fit(xdata, ydata)
        #ax[0].loglog(xdata, powerlaw(xdata, amp, index))     # Fit
        #ax[0].loglog(xdata, ydata, 'x')     # Fit
        #print("exposant for >1e0: {}".format(index))
        xdata = small_v["coeff_velocity"].values
        ydata = small_v["delta"].values
        pfinal, covar, index, amp = power_law_fit(xdata, ydata)
        ax[1].loglog(xdata, powerlaw(xdata, amp, index))     # Fit
        ax[1].loglog(xdata, ydata, 'x')     # Fit
        print("exposant for <1e-1: {}".format(index))

        #print(df)


if __name__ == "__main__":

    folder = "/home/marine/ownCloud/Research/Projets/output_mush/compaction/presentation"
    #all_figures(folder+"/")
    
    list_folder = os.listdir(folder)

    for name in list_folder:
       if name[:3] == "exp":
           print(folder+name)
           all_figures(folder+'/'+name)
           plt.close("all")

    
    fig_porosity_thickness(folder)
    scaling_laws()
    #fig_thickness(folder)
    #fig_porosity(folder)
    plt.show()
    # snippet for ordering dictionnary and print values.
    # from operator import itemgetter
    # for key, value in sorted(d.items(), key = itemgetter(1)):
    #     print(key, value)
