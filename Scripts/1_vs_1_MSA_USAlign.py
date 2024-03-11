import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

# Canviar a la ruta desitjada
ruta_archivo = "/home/guillem_ares/OneDrive_UB/ASSIGNATURES/TFM/Experimental/RNA_predicted/HBonds/Proves/" 
os.chdir(ruta_archivo)

# Dins de les carpetes amb les estructures (3dRNA, trRosetta i DRFold), generar una 
# llista amb TOTS els arxius *.pdb que hi hagi.
# Aquesta llista serveix més tard per executar un alineament amb US align
# Cal assegurar-se que les carpetes no contenen pdbs que no es vulguin alinear,
# incloses estructures dins de .zip (fer unzip), .tar (tar -xz), etc.
os.system("""find . -type d -exec sh -c 'cd "{}" && find . -type f -name "*.pdb" > lista_pdb.txt' \;""")

with open("lista_pdb.txt", "r") as f:
    models_pdb = f.readlines()
models_pdb = [x.strip() for x in models_pdb]

tm_score = np.zeros((len(models_pdb), len(models_pdb)))
rmsd = np.zeros((len(models_pdb), len(models_pdb)))

# Aquest bucle executa l'alineament amb US align entre tots els arxius 
# pdb de la llista 1 vs 1 i guarda els valors de RMSD i TM-score en 
# dues matrius

for i, pdb1 in enumerate(models_pdb):
    tm_score[i, i] = 1
    rmsd[i, i] = 0
    pdb1_filename = os.path.basename(pdb1)
    for j, pdb2 in enumerate(models_pdb):
        pdb2_filename = os.path.basename(pdb2)
        if i <= j:
            os.system(f"USalign {pdb1} {pdb2} -mol RNA > 'align_{pdb1_filename}_{pdb2_filename}.txt'")
            # SI ES VOLEN TOTS ELS OUTPUTS EXECUTAR AQUESTA COMANDA
            # os.system(f"USalign {pdb1} {pdb2} -mol RNA -rasmol sup{pdb1_filename}_{pdb2_filename} > 'align_{pdb1_filename}_{pdb2_filename}.txt'")
            rmsd_result = os.popen(f"grep 'RMSD=' align_{pdb1_filename}_{pdb2_filename}.txt").read()
            rmsd_line = rmsd_result.splitlines()[0]
            rmsd_value = rmsd_line.split("RMSD=")[1].split(",")[0].strip()
            rmsd[i, j] = rmsd_value
            rmsd[j, i] = rmsd_value

            tm_score_result = os.popen(f"grep 'TM-score=' align_{pdb1_filename}_{pdb2_filename}.txt").read()
            tm_score_line = tm_score_result.splitlines()[0]
            tm_score_value = tm_score_line.split("TM-score=")[1].split(",")[0].strip().split(" ")[0]
            tm_score[i, j] = tm_score_value
            tm_score[j, i] = tm_score_value

print("RMSD values")
print(rmsd)
print("TM-score values")
print(tm_score)

# Aquesta part del codi serveix per generar un heatmap amb els valors de RMSD i TM-score
data_rmsd = pd.DataFrame(rmsd, columns=models_pdb, index=models_pdb)
data_tm_score = pd.DataFrame(tm_score, columns=models_pdb, index=models_pdb)

sns.set_theme()
plt.figure(figsize=(10, 10))
sns.heatmap(data_rmsd, annot=True, cmap="YlGnBu", fmt=".2f", linewidths=0.5, linecolor="white")
plt.title("RMSD heatmap")
plt.savefig("rmsd_heatmap.png")

sns.set_theme()
plt.figure(figsize=(10, 10))
sns.heatmap(data_tm_score, annot=True, cmap="YlGnBu", fmt=".2f", linewidths=0.5, linecolor="white")
plt.title("TM-score heatmap")
plt.savefig("tm_score_heatmap.png")

# Aquesta part del codi serveix per generar clusters amb els valors de RMSD i TM-score
sns.set_theme()
sns.clustermap(data_rmsd, cmap="YlGnBu", linewidths=0.5, linecolor="white", figsize=(10, 10))
plt.title("RMSD clustermap")
plt.savefig("rmsd_clustermap.png")

sns.set_theme()
sns.clustermap(data_tm_score, cmap="YlGnBu", linewidths=0.5, linecolor="white", figsize=(10, 10))
plt.title("TM-score clustermap")
plt.savefig("tm_score_clustermap.png")