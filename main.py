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
    return temp[0], temp[1]

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
            processed[key][1][col].loc[idx] = "|".join(resultString(key, indx)[aux[0]-1:aux[1]])

for keys, name in zip(processed.keys(), ["Saved/6VYO_Processed.csv", "Saved/6WJI_Processed.csv"]):
    processed[key][1].to_csv(name)

            
#%% 8
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import Codes.function as fc

### To Do
# 1 - Criar lista com posições polimórficas
def snps_list(df, condition, size):
    return fc.filter_criteria(df, size, condition).dropna().Pos.to_list()

#counted_df[counted_df.Pos.isin(templist)]
templist = snps_list(counted_df, .5, len(raw_df_aln))
#2 - selecionar degrees da amostra 6VYO

#%% 9 
### Problema com a identificação dos valores no DF baseados no indice das sequencias do PDB


def tempplot(column, title):
    temp = [ x.split("|") for x in column.to_list()]
    temp2 = [float(el) for sublist in temp for el in sublist]
    sns.distplot(temp2).set_title(title)
    plt.show()

for col in ["Degrees", "ClusteringCoef", "BetweennessWeighted"]:
    for key in processed.keys():
        tempplot(processed[key][1][col], "Distribuição "+col+" "+key)


# %%
from tqdm import tqdm
tidy = {}

#Preparo a lista de filtagem, aqui possuo a lista com todos os SNPS
preFilterList =  [float(x) for x in ":".join(list(processed["6WJI"][1].SeqAlign)).split(":")]
snpsRange = (min(preFilterList), max(preFilterList))
filteredRange = list(filter(lambda x : x >= snpsRange[0] and x<= snpsRange[1], templist))

for key in processed.keys():
    temptidy = pd.DataFrame(columns=["Id", "Pos", "Amino", "AlnPos", "SeqPos", "Degr", "Clust", "Betw"])
    df = processed[key][1] 
    for idx, row in tqdm(df.iterrows(), total=len(df)):
        id = row.Id
        cleanSeq = row.Seq.replace("-","") 
        pos = range(1, len(cleanSeq)+1)
        amino = list(cleanSeq)
        degre = [float(x) for x in row.Degrees.split("|")]
        clust = [float(x) for x in row.ClusteringCoef.split("|")]
        betwe = [float(x) for x in row.BetweennessWeighted.split("|")]
        pdbAlnRange = range_converter(row.PDBAlign)
        pdbSeqRange = range_converter(row.SeqAlign)
        pdbAlnRange = list(range(pdbAlnRange[0],pdbAlnRange[1]+1))
        pdbSeqRange = list(range(pdbSeqRange[0], pdbSeqRange[1]+1))+[0]
        rangeControl = 0
        for p, a, pa, d, c, b in zip(pos, amino, pdbAlnRange, degre, clust, betwe):
            for _ in range(3):
                tpRange = pdbSeqRange[rangeControl]
                # aqui, adiciono ao DF apenas as posições que são SNPS, ganhando performance 
                if tpRange in filteredRange: 
                    temptidy.loc[len(temptidy)] = [id,p,a,pa,tpRange,d,c,b]
                rangeControl +=1
    tidy[key] = temptidy
    temptidy.to_csv("Saved/TempTidy_"+key+".csv")


# %%

for colP, colT in zip(["Degrees", "ClusteringCoef", "BetweennessWeighted"],['Degr', 'Clust', 'Betw']):
    for keyP, keyT in zip(processed.keys(), tidy.keys()):
        fig, ax = plt.subplots()
        temp = [float(x) for x in "|".join(processed[keyP][1][colP].to_list()).split("|") ]
        sns.distplot(temp, ax=ax, label="Norm").set_title("Distribuição "+colP+" "+keyP)
        sns.distplot(tidy[keyT][colT], ax=ax, label="SNPs")
        ax.legend()

# for col in ['Degr', 'Clust', 'Betw']:
#     fig, ax = plt.subplots()
#     for key in tidy.keys():
#         sns.distplot(tidy[key][col], ax=ax).set_title("Distribuição "+col+" "+key)

# %%
