import os
import glob

carpeta_structures = "../Experimental/RNA_predicted/3dRNA/SL1234"

def reemplazar_nt_terminal(pdb_file): 

    r""" En els fitxers de 3dRNA, els nt terminals 
         estan etiquetats com a G5, G3, C5, C3, A5, A3, U5, U3
         i això no és compatible amb biotite, que espera G, C, A, U
    """
    
    with open(pdb_file, "r") as file: 
        lines = file.readlines()      

    for i, line in enumerate(lines):
        if (line.startswith("ATOM") or line.startswith("TER")) and (line[18:20] == "G5" or line[18:20] == "G3"):
            lines[i] = line[:18] + " G" + line[20:]
        elif (line.startswith("ATOM") or line.startswith("TER")) and (line[18:20] == "C5" or line[18:20] == "C3"):
            lines[i] = line[:18] + " C" + line[20:]
        elif (line.startswith("ATOM") or line.startswith("TER")) and (line[18:20] == "A5" or line[18:20] == "A3"):
            lines[i] = line[:18] + " A" + line[20:]
        elif (line.startswith("ATOM") or line.startswith("TER")) and (line[18:20] == "U5" or line[18:20] == "U3"):
            lines[i] = line[:18] + " U" + line[20:]

    with open(pdb_file, "w") as file:
        file.writelines(lines)

# os.system("cd " + carpeta_structures)
archivos = os.path.join(carpeta_structures, "*.pdb")

for pdb_file in glob.glob(archivos):
    reemplazar_nt_terminal(pdb_file)

print("Done!")