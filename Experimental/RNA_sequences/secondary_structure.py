# Secondary Structure Consensus
import glob
from collections import Counter
import os

os.chdir('../data/sequences')

rna_seq = {}

for filename in glob.glob('*_rna.fasta'):
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            else:
                rna_seq[filename] = line.strip()

os.chdir('../SS')

posicion_final = 480

rna_struct = {}
rna_segmentos = {}

segmentos = {
    'SL1': (7, 33), 
    'SL2': (45, 60), 
    'SL23': (45, 75),
    'SL3': (61, 75), 
    'SL4': (84, 125),
    'SL1234': (7, 127), 
} 

for filename in glob.glob('*.ss'):
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('nt'):
                inicio, fin = map(int, line.split()[1].split('-'))
                secuencia = ['*'] * (inicio - 1)
            else:
                secuencia.extend(list(line.strip()))
                secuencia.extend(['*'] * (posicion_final - len(secuencia)))
                rna_struct[filename] = secuencia
                rna_segmentos[filename] = {nombre: secuencia[inicio-1:fin] for nombre, (inicio, fin) in segmentos.items()}

with open('rna_struct.txt', 'w') as f:
    for filename, secuencia in rna_struct.items():
        f.write(f'{filename}:\n\tfull seq: {"".join(secuencia)}:\n')
        for nombre, secuencias in rna_segmentos[filename].items():
            f.write(f'\t{nombre}:\n\t{"".join(secuencias)}\n')

contadores = {nombre: [] for nombre in segmentos.keys()}

for filename, segs in rna_segmentos.items():
    for nombre, secuencia in segs.items():
        if not contadores[nombre]:
            contadores[nombre] = [Counter() for _ in secuencia]
        for i, base in enumerate(secuencia):
            contadores[nombre][i][base] += 1

secuencias_consenso = {nombre: ''.join(counter.most_common(1)[0][0] for counter in contadores_nombre) for nombre, contadores_nombre in contadores.items()}

with open('consenso.txt', 'w') as f:
    for nombre, secuencia in secuencias_consenso.items():
        inicio, fin = segmentos[nombre]
        f.write(f'{nombre} ({inicio},{fin}):\n\t{secuencia}\n')