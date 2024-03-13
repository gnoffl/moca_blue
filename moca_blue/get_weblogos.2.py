from weblogo import *
import re
import numpy as np

# Read the file and extract the PWM motifs
with open('rdf5_ArthD6_cwm-motifs.jaspar', 'r') as file:
    data = file.read()

motifs = re.findall(r'>.*?[\s\S]*?(?=>|$)', data)

for motif in motifs:
    lines = motif.strip().split('\n')
    title = lines[0][1:]
    matrix = {line.split()[0]: [int(x) for x in re.findall(r'\d+', line)] for line in lines[1:]}

    # Transpose the matrix
    counts = np.array([matrix[nuc] for nuc in ['A', 'C', 'G', 'T']]).T

    # Create the logo options and format the data
    options = LogoOptions()
    data = LogoData.from_counts('ACGT', counts)
    options.title = title
    options.show_errorbars = False
    options.fineprint = ''
    format = LogoFormat(data, options)

    # Create the logo
    eps = eps_formatter(data, format)
    with open(f'{title}.eps', 'wb') as f:
        f.write(eps)
