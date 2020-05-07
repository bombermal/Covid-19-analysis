#Imports
import pandas as pd

#My imports
import Codes.function as fc
import Codes.graph as gf
import Codes.read as rdF
#%% 1 - Load files

#Load trimmed aln
path = "Read/"
file = "3_Human_cov30032020_sequences.aln.trimmed.fasta"
loaded_files = rdF.read_aligned_files(path, file)
#Put Sequences in a dataframe
raw_df_aln = rdF.seq_to_df(loaded_files)

#%% Processo de contagem

#Df vazio
counted_df = pd.DataFrame()
#função que faz a transposta de Seq e conta as ocorrencias
fc.transpose_seq_and_count(raw_df_aln.Seq, counted_df)
#Limpa os Nan e corrige a Pos
counted_df = counted_df.fillna(0).reset_index().rename(columns={"index" : "Pos"})
counted_df.Pos = counted_df.Pos.apply(lambda x: x+1)
#Salva o trabalho
counted_df.to_csv("Saved/3_Counting_df_Covid19_1302S.csv", index=False)
#%% REading saved counted DF

counted_df = pd.read_csv("Saved/3_Counting_df_Covid19_1302S.csv")

#%% Ploting

conditions = [.5,1,2]
list_of_filtered_dfs = []

for ii in conditions:
  list_of_filtered_dfs.append(fc.filter_criteria(counted_df, 1964, ii))
  
gf.three_plots(conditions, list_of_filtered_dfs, "3_Covid19_1964S", 250, True)
gf.three_plots(conditions, list_of_filtered_dfs, "3_Covid19_1964S", 250)

#%% Vamos pensar
