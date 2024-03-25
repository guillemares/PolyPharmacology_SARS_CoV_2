#!/usr/bin/env python
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

class superimposition(object):
    print('Hello, I am superimposition class')

"""
---------------------------
---------------------------
---------------------------
"""

if __name__ == '__main__':
    import doctest
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--method',
                        help='Choose the preferred superimposition method.\n'
                             '1vs1: Perform superimposition between pairs of structures.\n'
                             'allvsall: Perform superimposition between all possible pairs of structures.',
                        choices=['1vs1','allvsall'])
    parser.add_argument('--test',
                        help='Test the code',
                        action='store_true')
    parser.add_argument('--dir',
                        help='Provide the directory with multiple pdb files to be superimposed')


    args = parser.parse_args()