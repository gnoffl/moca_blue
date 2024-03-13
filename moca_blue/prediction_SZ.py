import os
import pandas as pd
import numpy as np
from tensorflow.keras.models import load_model
from utils2 import onehot
from pyfaidx import Fasta
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
pd.options.display.width = 0


def prepare_seqs(genome, df_gtf, upstream=1000, downstream=500):
    chroms = [f"{i}" if i < 10 else f"{i}" for i in range(1, 13)]
    gene_model_df = df_gtf
    gene_model_df = gene_model_df[gene_model_df[0].isin(chroms)]

    x_test, gene_ids = [], []
    for chrom, start, end, strand, gene_id in gene_model_df.values:
        if strand == '+':
            prom_start, prom_end = start - upstream, start + downstream
            term_start, term_end = end - downstream, end + upstream
            if prom_start > 0 and term_start > 0:
                encoded_seq = np.concatenate([onehot(genome[chrom][prom_start:prom_end]),
                                              np.zeros(shape=(20, 4)),
                                              onehot(genome[chrom][term_start:term_end])])
                if encoded_seq.shape[0] == 2 * (upstream + downstream) + 20:
                    x_test.append(encoded_seq)
                    gene_ids.append(gene_id)

        else:
            prom_start, prom_end = end - downstream, end + upstream
            term_start, term_end = start - upstream, start + downstream
            if prom_start > 0 and term_start > 0:
                encoded_seq = np.concatenate([onehot(genome[chrom][prom_start:prom_end])[::-1, ::-1],
                                              np.zeros(shape=(20, 4)),
                                              onehot(genome[chrom][term_start:term_end])[::-1, ::-1]])

                if encoded_seq.shape[0] == 2 * (upstream + downstream) + 20:
                    x_test.append(encoded_seq)
                    gene_ids.append(gene_id)

    x_test = np.array(x_test)
    x_test[:, upstream:upstream+3, :] = 0
    x_test[:, upstream+1017:upstream+1020, :] = 0
    return x_test, gene_ids


path_genome = './genomes/Alyrata.v.1.0.dna.toplevel.fa'
path_annot = './gene_models/Alyrata_v.1.0.57.gtf'

fasta = Fasta(path_genome, as_raw=True, sequence_always_upper=True, read_ahead=10000)
gene_models = pd.read_csv(path_annot, sep='\t', comment='#', header=None)
gene_models = gene_models[gene_models[2] == 'gene']
gene_models[9] = [x.split('Name=')[-1].split(';')[0] for x in gene_models[8]]
gene_models = gene_models[[0, 3, 4, 6, 9]]

seqs, genes = prepare_seqs(genome=fasta, df_gtf=gene_models)
print(seqs.shape)
print(len(genes))
print(genes[:6])

model_prefix = 'Alyr-6h_model_'
output_folder = './saved_models/predictions/'

matching_files = [file_name for file_name in os.listdir('saved_models') if model_prefix in file_name]

if not matching_files:
    print(f"No matching files found with the prefix '{model_prefix}'.")
else:
    all_predictions_df = pd.DataFrame({'gene_id': genes})
    
    for file_name in matching_files:
        print(file_name)
        saved_model = load_model(f'saved_models/{file_name}')
        predictions = saved_model.predict(seqs)
        model_number = file_name.replace(model_prefix, '').replace('_promoter_terminator.h5', '')
        model_predictions_df = pd.DataFrame({'gene_id': genes, f'Alyr_{model_number}': np.mean(predictions, axis=1)})
        
        all_predictions_df = pd.merge(all_predictions_df, model_predictions_df, on='gene_id', how='outer')
    
    all_predictions_df.to_csv(f'{output_folder}{model_prefix}all_predictions.csv', index=False)
    
    print(all_predictions_df.head())