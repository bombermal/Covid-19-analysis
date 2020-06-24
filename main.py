#%% 0 - Imports
#Imports

#My imports
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

#%% 4 - Ler PDB and Betwenness

# Read PDB's
pdbsNames = [ "6vyo", "6wji"]
pdbsNodesFiles, pdbsEdgesFiles, pdbsModifiedFiles = wf.read_pdbs(pdbsNames)
# Load Betwenness
path = "Read/Betweness e Coef. Clusterização/Nucleocapsid_Protein_NodesResult.csv"
betwenness_dict = wf.read_betwenness(path)
# Merge Dataframes 
wf.join_pdb_betw(pdbsNodesFiles, betwenness_dict)

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

#1 - Separa os pdbs em dicionários
filtered_dfs = wf.filter_df_on_pdbsDict(processed, pdbsNodesFiles)

#2 - converter as sequencias em lista de nós
nodes_list_from_pdbs, nodes_list_from_sample = wf.sequences_to_nodes_list(filtered_dfs, processed)

#3 - alinhar

"""
    Alinhas os pdbs com as sequências iniciais dos arquivos multifasta
"""
wf.align_head_pdb(filtered_dfs, processed, nodes_list_from_pdbs, nodes_list_from_sample)

#%% 7 - Transpor resultados
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

from tqdm import tqdm
def resultString(key, indexName):
    return processed[key][0][indexName].split("|")

for key in processed.keys():
    indexes = ["AlnDgree", "AlnClust", "AlnBetw"]
    columns = [ 'Degrees', 'ClusteringCoef', 'BetweennessWeighted']

    processed[key][1]["Degrees"] = 0
    processed[key][1]["ClusteringCoef"] = 0
    processed[key][1]["BetweennessWeighted"] = 0
    for idx, row in tqdm(processed[key][1].iterrows(), total=len(processed[key][1])):
        aux = range_converter(row.PDBAlign)
        for indx, col in zip(indexes, columns):
            processed[key][1][col].loc[idx] = "|".join(resultString(key, indx)[aux[0]:aux[1]])
#%%
for keys, name in zip(processed.keys(), ["Saved/6VYO_Processed.csv", "Saved/6WJI_Processed.csv"]):
    processed[key][1].to_csv(name)

            
# %%
import matplotlib.pyplot as plt
import seaborn as sns

### To Do
# 1 - Gerar Plots para todas os pdbs
for key in processed.keys():
    x = []#[ float(x) for x in processed[key][0].AlnDgree.split("|")]
    for idx, row in tqdm(processed[key][1].iterrows(), total=len(processed[key][1])):
        x += [ float(x) for x in row.Degrees.split("|")]

    gf = sns.distplot(x)
    fig = gf.fig

    fig.suptitle("teste")
# 2 - Gerar plot mostrando apenas as posições polimórficas, pegar a tabela filtrada modificando o
# algoritmo de filtragem



# %%
