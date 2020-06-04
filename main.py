#%% 0 - Imports
#Imports
#import pandas as pd
#rom datetime import datetime
#import re
#My imports
#import Codes.function as fc
#import Codes.graph as gf
#import Codes.read as rdF
import Codes.workFlow as wf

#%% 1 - Lê arquivos

#Load trimmed aln
path = "Read/4_Human_cov03042020_sequences.aln.trimmed.fasta"
raw_df_aln = wf.read_raw_files(path)

#%% 2 - Processo de contagem

path = "Saved/4_Counting_df_Covid19_2340.csv"

# True: lê csv anteriormente criado, False: Lê fasta e cria arquivo csv
counted_df = wf.read_Or_create(path=path,  ler=True)#False, df=raw_df_aln)

#%% 3 - Ploting
conditions = [.5,1,2]

# Descomente para rodar
#wf.plot(conditions, len(raw_df_aln), counted_df, "4_Covid19_2340" )

#%% 4 - Ler PDB

# Read PDB's
pdbsNames = [ "6vyo", "6wji"]
pdbsNodesFiles, pdbsEdgesFiles, pdbsModifiedFiles = wf.read_pdbs(pdbsNames)

#%% 5 - Tratar Raw data, separa da tabela de alinhamento o pdb das amostras


csvPath = "Saved/pdb_6vyo.csv"
dfPath = "Saved/ProcessedData_6vyo.csv"
path_list = ["Read/PDBs RIN/align-refs/6vyo_A-gs.fasta.aln", "Read/PDBs RIN/align-refs/6wji_A-gs.fasta.aln"]
# True: lê csv anteriormente criado, False: Lê fasta.aln e cria arquivo csv
#Salva o resultado em um dicionário com key = "nome do pdb" e value = [series(informações do pdb), df(informações das amostras)]
processed = {}
for ii in path_list:
    auxPdb, auxProcessedData = wf.pdbRead_or_Create(alnPath=ii, ler=False)
    processed[str.upper(auxPdb.PdbName)] = [auxPdb, auxProcessedData] 
"""
Bug conhecido - O arquivo pdb salvo não mantém a primeira coluna como index ( salvar series é diferente?)
"""
#%% 6 - Transformar sequencias em nós
"""
Filtrar baseado na entrada, verifico quais cadeias estão presentes na entrada
e a partir disso, filtro a parte dos pdbs com as mesmas cadeias, para evitar
futuros problemas com alinhamento
"""
filtered_dfs = wf.filter_df_on_pdbsDict(processed, pdbsNodesFiles)

#2 - converter as sequencias em lista de nós
#Dicionário com pdbs Filtrados na forma de nós, prontos para linhar contra pdb do multifasta
nodes_list_from_pdbs = {}
for ii in filtered_dfs.keys():
    nodes_list_from_pdbs[ii] = wf.pdb_file_to_list_of_nodes(filtered_dfs[ii])

nodes_list_from_sample = {}
for ii in processed.keys():
    nodes_list_from_sample[ii] = wf.simple_seq_to_list_of_nodes(processed[ii][0])

#%%
#3 - alinhar

"""
    Verificar o alinhamento - alinhando apenas um por vez. Ampliar
"""
import Codes.smithWaterman as sw

obj = sw.smithWaterman()
seqIds, seqPos, cover = obj.constructor(2, -1, -1, nodes_list_from_sample["6VYO"], nodes_list_from_pdbs["6VYO"], True, False)
#%%
#4 salvar resultado em df
#Lista de Ids dos aminoácidos alinhados
pdb["AlnResult"] = "|".join(seqIds)
#Lista de Degrees criadas a partir da lista de Ids
pdb["AlnDgree"] = "|".join([str(filtered_dfs[filtered_dfs.NodeId == x].Degree.values[0]) for x in seqIds])

#%% 7 - Transpor resultados
import matplotlib.pyplot as plot
import seaborn as sns

#temp = [float(x) for x in pdb.AlnDgree.split("|")]
#sns.distplot(temp, kde=False)
processedData["Degree"] = 0
temp = processedData.loc[11278]

def range_converter(rng):
    """
    Converte PDBRange '1:124' em um range(0:124)

    Parameters
    ----------
    rng : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    temp = [int(x) for x in rng.split(":")]

    #return range(temp[0]-1, temp[1])
    return temp[0]-1, temp[1]
aux = range_converter(temp.PDBAlign)
temp.Degree = "|".join(pdb.AlnDgree.split("|")[aux[0]:aux[1]])
# %%
