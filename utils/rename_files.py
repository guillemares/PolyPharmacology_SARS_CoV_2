import os
import glob

def change_name(path, old_name, new_name):
    """
    Function to change the name of files in a directory
    to a new name with a number at the end.

    ----------------
    Example:
    SL1_3dRNA_pred1.pdb -> 3dRNA_1.pdb

    """
    os.chdir(path)
    i = 1
    for file in glob.glob(old_name):
        file_name, file_extension = os.path.splitext(file)
        os.rename(file, new_name + str(i) + file_extension)
        i += 1
    return 0

