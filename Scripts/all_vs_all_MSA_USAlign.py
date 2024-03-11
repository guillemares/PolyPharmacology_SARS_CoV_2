import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

# Canviar a la carpeta on es tenen les estrucures a alinear
ruta_input = "/home/guillem_ares/OneDrive_UB/ASSIGNATURES/TFM/Experimental/RNA_predicted/3dRNA/SL1234" 

os.chdir(ruta_input)

# Si hi ha arxius en la carpeta que comencin amb l'extensió
# sup*, esborrar-los (són els resultats de l'últim alineament)
# i tenen extensió .pdb. Si no s'esborren, el programa els
# agafarà com a fitxers d'entrada i els alinearà amb tots

os.system("rm sup*.pdb")

# Aquesta comanda serveix per, a la carpeta 3dRNA, trRosetta i DRFold, generar una llista amb TOTS els arxius *.pdb que hi hagi
# aquesta llista serveix més tard per executar un alineament amb US align
# cal assegurar-se que les carpetes no contenen .zip (fer unzip) o .tar (tar -xz)
# i que NOMÉS contenen els arxius .pdb que es vol alinear (per tant,
# si s'ha executat 1_vs_1_MSA_USAlign.py, s'ha d'esborrar els altres fitxers)

os.system("""find . -type d -exec sh -c 'cd "{}" && find . -type f -name "*.pdb" > lista_pdb.txt' \;""")


with open("lista_pdb.txt", "r") as f:
    models_pdb = f.readlines()
models_pdb = [x.strip() for x in models_pdb]

# Per alinear tots els arxius que hi ha al directori (comanda -mm 4) 
# i a més generar els arxius 3d que es poden veure amb el programa rasmol 
# (comanda -rasmol sup); l'output que s'imprimeix per pantalla es guarda 
# a alignment.txt (els valors de RMSD, TM-score... entre d'altres)

os.system(""" USalign -dir ./ lista_pdb.txt -mm 4 -rasmol sup > 'alignment.txt' """)

# Per si es vol comprimir els arxius del rasmol
# os.system("tar -cvzf sup*")

# Buscar els valors promig de RMSD i TM-score a l'arxiu alignment.txt
rmsd_result = os.popen(""" grep 'RMSD=' alignment.txt""").read()
rmsd_line = rmsd_result.splitlines()[0]
rmsd_value = rmsd_line.split("RMSD=")[1].split(",")[0].strip()
print("Average RMSD value: ", rmsd_value)

# El càlcul de TM-score que es mostra serà unicament per estrutures
# amb exactament el mateix nombre de nucleòtids.

tm_score_result = os.popen(""" grep 'TM-score=' alignment.txt""").read()
tm_score_line = tm_score_result.splitlines()[0]
tm_score_value = tm_score_line.split("TM-score=")[1].split(",")[0].strip().split(" ")[0]
print("Average TM-score value: ", tm_score_value)