import re
import sys
import os
from subprocess import Popen, PIPE
import multiprocessing
from string import Template
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot as plt


n_bins=90
regex = r"traj_comp_([0-9]?.?[0-9])_([0-1]).xtc"
#regex = r"traj_merged.xtc"
vector_template = """["name CA and (bynum $atomNumbers)", "bynum $atomNumbers and resnum 12"]"""
reference_template = """["name CA and (bynum $atomNumbers)",[50.51,51.66,0]]"""
args = sys.argv[1:]
commands_to_run = []
outputs = {
    "window": [],
    "index": [],
    "output": []
}

for arg in args:
    res = re.search(regex, arg)
    window = res.group(1)
    index = res.group(2)
    trajectory = arg
    topology = f"w_{window}/prod_0/traj_comp_ca.gro"
    output = f"temp/w_{window}_{index}.xvg"
    if index == "0":
        template_values = dict(atomNumbers="1:23")
    else:
        template_values = dict(atomNumbers="24:46")
    vector = Template(vector_template).substitute(template_values)
    reference = Template(reference_template).substitute(template_values)
    outputs["window"].append(window)
    outputs["index"].append(f"{index}")
    outputs["output"].append(f"{output}")
    command = ["python", "./theta/analysis/theta.py"]
    command += [f"--traj={trajectory}"]
    command += [f"--top={topology}"]
    command += [f"--output={output}"]
    command += [f"--vector={vector}"]
    command += [f"--reference={reference}"]
    commands_to_run.append(command)


def work(cmd):
    print(cmd)
    #proc = Popen(cmd)
    print(proc)
    return proc

def print_window_angle_ditribution(data):
    plt.figure()
    plt.errorbar(data["window"], data["theta_max_index"],data["gauss_sigma"])
    plt.plot(data["window"], data["theta_max_index"])
    plt.ylim(0,360)
    plt.savefig("window_angle_distribution.png")

def gauss_for(x0):
    return lambda x, a, sigma : a*exp(-(x-x0)**2/(2*sigma**2))

def read_outputs():
    output_df = pd.DataFrame.from_dict(outputs)

    for index, row in output_df.iterrows():
        fileName = row["output"]
        print(fileName)
        df = pd.read_csv(fileName,sep="\t", names=["t", "theta"])
        hist, bins = np.histogram(df["theta"], bins=n_bins)
        hist_df = pd.DataFrame({"hist": hist})
        max_index = hist_df.idxmax().values[0]
        fun = gauss_for(bins[max_index])
        popt, popc = curve_fit(fun, bins[:-1], hist)

        output_df.at[index,"max"] = np.max(hist)
        output_df.at[index,"theta_max_index"] = max_index
        output_df.at[index,"gauss_a"] = popt[0]
        output_df.at[index,"gauss_sigma"] = popt[1]
        
        plt.figure()
        plt.plot(bins, fun(bins, popt[0], popt[1]))
        plt.hist(df["theta"], n_bins)
        plt.savefig(f"w_{row['window']}_{row['index']}.png")
        plt.close()
       # print(df)
    output_df = output_df.sort_values("window")
    zeros_df = output_df[output_df["index"] == "0"]
    
    print_window_angle_ditribution(zeros_df)

    print(zeros_df)


count = multiprocessing.cpu_count()
print(count)
pool = multiprocessing.Pool(processes=count)
print(pool)
r = pool.map_async(work, commands_to_run)
print(r)
r.wait()
print("Analyzing outputs")
read_outputs()
print("Finished")
