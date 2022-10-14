import pandas as pd
import numpy as np


class GaussianKernel:
    """compute gaussian kernel weighting between y1 and y2"""

    def __init__(self, bandwidth=1.0):
        """
        bandwith (float): bandwith argument to kernel, default = 1.0
        """
        self.bandwidth = bandwidth

    def values(self, y1, y2):
        """
        compute gaussian kernel between y1 and y2
        """
        return np.exp(-((y1 - y2) ** 2) / 2 / self.bandwidth)


class BoxKernel:
    """compute box kernel between y1 and y2"""

    def __init__(self, bandwidth=1.0):
        """
        bandwith (float): bandwith argument to kernel, default = 1.0
        """
        self.bandwidth = bandwidth

    def values(self, y1, y2):
        """
        compute box kernel between y1 and y2
        """
        return 1.0 * (np.abs(y1 - y2) <= self.bandwidth / 2)
