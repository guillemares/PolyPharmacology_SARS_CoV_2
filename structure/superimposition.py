#!/usr/bin/env python
import sys
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
    def __init__(self, method, dir, test=False):
        self.method = method
        self.dir = dir
        self.test = test

        self._get_files()

    def _get_files(self):
        """
        Get all the pdb files in the directory and create a list
        with the file names.

        -----------

        Returns:
        files: list
            List of pdb files in the directory

        >>> superimpobj = superimposition(method='1vs1', dir='../test/superimposition/', test=True)
        >>> superimpobj._get_files()
        0
        """

        os.chdir(self.dir)
        os.system("""find . -type d -exec sh -c 'cd "{}" && find . -type f -name "*.pdb" > pdb_filename.txt' \;""")

        with open('pdb_filename.txt', 'r') as f:
            models_pdb = f.readlines()
        models_pdb = [x.strip() for x in models_pdb]
        self.models_pdb = models_pdb

        if self.test:
            return models_pdb
        else:
            return 0

    def _1_vs_1(self):
        """
        Perform superimposition between pairs of structures to
        create a matrix of RMSD and TM scores.

        -----------
        Returns:
        tm_score: numpy array
            Matrix with TM scores between pairs of structures
        rmsd: numpy array
            Matrix with RMSD values between pairs of structures

        # >>> superimpobj = superimposition('1vs1', 'test')
        # >>> superimpobj._1_vs_1()
        #
        """

        self.tm_score = np.zeros((len(self.models_pdb), len(self.models_pdb)))
        self.rmsd = np.zeros((len(self.models_pdb), len(self.models_pdb)))

        for i, pdb1 in enumerate(self.models_pdb):
            self.tm_score[i, i] = 1
            self.rmsd[i, i] = 0
            pdb1_filename = os.path.basename(pdb1)
            for j, pdb2 in enumerate(self.models_pdb):
                pdb2_filename = os.path.basename(pdb2)
                if i <= j:
                    os.system(f"../../USalign {pdb1} {pdb2} -mol RNA > 'align_{pdb1_filename}_{pdb2_filename}.txt'")
                    # os.system(f"USalign {pdb1} {pdb2} -mol RNA -rasmol sup{pdb1_filename}_{pdb2_filename} > 'align_{pdb1_filename}_{pdb2_filename}.txt'")
                    rmsd_result = os.popen(f"grep 'RMSD=' align_{pdb1_filename}_{pdb2_filename}.txt").read()
                    rmsd_line = rmsd_result.splitlines()[0]
                    rmsd_value = rmsd_line.split("RMSD=")[1].split(",")[0].strip()
                    self.rmsd[i, j] = rmsd_value
                    self.rmsd[j, i] = rmsd_value

                    tm_score_result = os.popen(f"grep 'TM-score=' align_{pdb1_filename}_{pdb2_filename}.txt").read()
                    tm_score_line = tm_score_result.splitlines()[0]
                    tm_score_value = tm_score_line.split("TM-score=")[1].split(",")[0].strip().split(" ")[0]
                    self.tm_score[i, j] = tm_score_value
                    self.tm_score[j, i] = tm_score_value

        if self.test:
            return self.tm_score, self.rmsd
        else:
            return 0

    def compress_files_1vs1(self):
        """
        Compress the 1vs1 pdb files and the alignment files in a .tar.gz file.
        """
        os.system("tar -czvf alignment_1vs1.tar.gz align_*")
        os.system("rm align_*")


    def df_1vs1(self):
        """
        Create a dataframe with the RMSD and TM scores between pairs of structures.
        """
        self.data_rmsd = pd.DataFrame(self.rmsd, columns=self.models_pdb, index=self.models_pdb)
        self.data_tm = pd.DataFrame(self.tm_score, columns=self.models_pdb, index=self.models_pdb)

        return 0

    def _all_vs_all(self):
        """
        Perform superimposition between all structures simultaneously.

        -----------

        Returns:
        rmsd: float
            Average RMSD value
        tm_score: float
            Average TM-score value

        # >>> superimpobj = superimposition('allvsall', 'test')
        # >>> superimpobj._all_vs_all()
        #
        """
        os.system(""" ../../USalign -dir ./ pdb_filename.txt -mm 4 -rasmol sup > 'alignment.txt' """)

        rmsd_result = os.popen(""" grep 'RMSD=' alignment.txt""").read()
        rmsd_line = rmsd_result.splitlines()[0]
        rmsd_value = rmsd_line.split("RMSD=")[1].split(",")[0].strip()
        self.rmsd = rmsd_value

        tm_score_result = os.popen(""" grep 'TM-score=' alignment.txt""").read()
        tm_score_line = tm_score_result.splitlines()[0]
        tm_score_value = tm_score_line.split("TM-score=")[1].split(",")[0].strip().split(" ")[0]
        self.tm_score = tm_score_value

        if self.test:
            return self.rmsd, self.tm_score
        else:
            return 0

    def compress_files_allvsall(self):
        """
        Compress the allvsall pdb files and the alignment files in a .tar.gz file.
        """
        os.system("tar -czvf alignment_allvsall.tar.gz sup*")
        os.system("rm sup*")

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

    if args.test:
        doctest.testmod(
            optionflags=doctest.ELLIPSIS | doctest.REPORT_ONLY_FIRST_FAILURE)
        sys.exit()

    superimpobj = superimposition(args.method, args.dir)

    if args.method == '1vs1':
        superimpobj._1_vs_1()
        superimpobj.compress_files_1vs1()
        superimpobj.df_1vs1()
    elif args.method == 'allvsall':
        superimpobj._all_vs_all()
        superimpobj.compress_files_allvsall()

