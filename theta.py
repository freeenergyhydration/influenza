#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys, getopt
import MDAnalysis
import numpy as np
from theta_analysis import ThetaAnalysis as analysis
import json

class POPC:
    def __init__(self, name, atoms):
        self.Name = name
        self.Atoms = atoms


def main(argv):
    topology_file = ''
    trajectory_file = ''
    output_file = ''
    reference_vector = ''
    vector = ''

    try:
        opts, args = getopt.getopt(argv, "h", ["traj=", "top=", "output=", "reference=", "vector="])
    except getopt.GetoptError:
        handle_error()

    if len(opts) != 5:
        print("Invalid number of parameters")
        handle_error()

    for opt, arg in opts:
        if opt == '-h':
            print_help()
            sys.exit()
        elif opt == "--traj":
            trajectory_file = arg
        elif opt == "--top":
            topology_file = arg
        elif opt == "--output":
            output_file = arg
        elif opt == "--vector":
            print("Decoding: " + arg)
            vector = json.loads(arg)
        elif opt == "--reference":
            print("Decoding: " + arg)
            reference_vector = json.loads(arg)

    analyze_theta(topology_file, trajectory_file, output_file, reference_vector, vector)


def handle_error():
    print('Error!')
    print_usage()
    sys.exit(2)


def analyze_theta(topology, trajectory, output, reference_vector, vector):
    """
    Analyze misalignment in cholesterol surface
    :param topology: topology to read
    :param trajectory: trajectory to read
    :param output: output for index with misaligned frames
    :param vector: target vector
    :return:
    """
    print(topology)
    print(trajectory)
    u = MDAnalysis.Universe(topology, trajectory)
    print("Reading finished")
    verification = analysis(u, reference_vector, vector,  output_file_name=output)
    print("Analysis starting")
    verification.run()
    print("Analysis finished")
    print_result(output)


def print_result(output_path):
    print("Generated output file: " + output_path)


def print_help():
    print_usage()
    print_description()


def print_usage():
    print("Usage: python " + sys.argv[
        0] + " --traj=trajectory_file_path --top=topology_file_path --output=output_file_path --vector=[x,y,z] --reference=[x,y,z]")
    print("Or: python " + sys.argv[
        0] + " --traj=trajectory_file_path --top=topology_file_path --output=output_file_path --vector=[start_selector, end_selector] --reference=[start_selector, end_selector]")


def print_description():
    print(
        "Program will generate file with angles between vectors based on selectors or constant vectors for whole trajectory."
        "Selectors can be atoms selectors compatible with MDAnalysis or constant values in format [x,y,z], eg. --vector=[[1,2,3],[4,5,6]], --vector=['protein', [1,2,3]]")

main(sys.argv[1:])
