from console_progressbar import ProgressBar
import matplotlib.pyplot as plt
import MDAnalysis.analysis.base as anal
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.ndimage import gaussian_filter
import matplotlib
import MDAnalysis.lib.mdamath as mdmath
import types

matplotlib.use('Agg')


class ThetaAnalysis(anal.AnalysisBase):
    def __init__(self, universe, reference_vector, vector, output_file_name, **kwargs):
        super(ThetaAnalysis, self).__init__(universe.trajectory, **kwargs)
        print("Analyzing angle between vector %s and %s" % (reference_vector,
                                                            vector))
        self.reference_is_const = len(reference_vector) == 3
        self.vector_is_const = len(vector) == 3
        self.reference = reference_vector
        self.vector = vector
        self.u = universe
        self.output_file_name = output_file_name

    def _prepare(self):
        # Initialize
        print("prepare")
        self.angles = {
            "Frame": [],
            "Theta": []
        }
        self.pb = ProgressBar(len(self.u.trajectory), length=30)

    def __prepare_vector__(self, vector, is_const):
        if is_const:
            return vector

        if type(vector[0]) is unicode or type(vector[0]) is str:
            start = self.u.select_atoms(vector[0]).center_of_geometry()
        else:
            start = vector[0]

        if type(vector[1]) is unicode or type(vector[1]) is str:
            end = self.u.select_atoms(vector[1]).center_of_geometry()
        else:
            end = vector[1]

        result = end - start
        return result

    def __angle__(self, reference, target):
        v1 = [reference[0], reference[1], 0]
        v2 = [target[0], target[1], 0]

        angle = mdmath.angle(v1, v2)
        return np.degrees(angle)

    def _single_frame(self):
        # process single frame
        reference_vector = self.__prepare_vector__(
            self.reference, self.reference_is_const)
        target_vector = self.__prepare_vector__(
            self.vector, self.vector_is_const)

        self.angles["Frame"].append(self._ts.frame)
        self.angles["Theta"].append(
            self.__angle__(reference_vector, target_vector))

        self.pb.print_progress_bar(self._ts.frame)

    def _conclude(self):
        # Volume in each radial shell
        self.pb.print_progress_bar(len(self.u.trajectory))

        result = pd.DataFrame.from_dict(self.angles)

        result.to_csv(self.output_file_name, sep="\t",
                      header=False, index=False)
        result.to_pickle(self.output_file_name + ".pkl")
        print("conclude")
