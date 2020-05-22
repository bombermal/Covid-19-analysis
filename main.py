#%% 0 - Imports
#Imports
import pandas as pd
from datetime import datetime
import re
#My imports
import Codes.function as fc
import Codes.graph as gf
import Codes.read as rdF
import Codes.workFlow as wf

#%% 1 - Lê arquivos

#Load trimmed aln
path = "Read/4_Human_cov03042020_sequences.aln.trimmed.fasta"
raw_df_aln = wf.read_files(path)

#%% Processo de contagem

path = "Saved/4_Counting_df_Covid19_2340.csv"

def read_Or_create(path, ler = True, df = ""):
    if ler:
        return pd.read_csv(path)
    else:
        #Df vazio
        counted_df = pd.DataFrame()
        #função que faz a transposta de Seq e conta as ocorrencias
        fc.transpose_seq_and_count(df.Seq, counted_df)
        #Limpa os Nan e corrige a Pos
        counted_df = counted_df.fillna(0).reset_index().rename(columns={"index" : "Pos"})
        counted_df.Pos = counted_df.Pos.apply(lambda x: x+1)
        #Salva o trabalho
        counted_df.to_csv(path, index=False)
        
        return counted_df

# True: lê csv anteriormente criado, False: Lê fasta e cria arquivo csv
counted_df = read_Or_create(path=path,  ler=True)#False, df=raw_df_aln)

#%% Ploting
def plot(conditions, counted_df, name_id, target="Img/", dpi=50):
    
    list_of_filtered_dfs = []
    
    for ii in conditions:
      list_of_filtered_dfs.append(fc.filter_criteria(counted_df, len(raw_df_aln), ii))
      
    gf.three_plots(conditions, list_of_filtered_dfs, name_id, target=target, dpi=dpi)

conditions = [.5,1,2]

# Descomente para rodar
#plot(conditions, counted_df, "4_Covid19_2340" )

#%% Ler PDB

def read_pdbs(pdbsNames):
    """
        rd.readPDBs(fileNames, option, discotop value)
        fileName = list with pdbs names
        option = 'node', 'edges', 'discotop'
        discotop value = -3.7 or -7.7
    """
    pdbsNodesFiles = rdF.readPDBs(pdbsNames, ".pdb.nodes")
    pdbsEdgesFiles = rdF.readPDBs(pdbsNames, ".pdb.edges")

    pdbsModifiedFiles = rdF.readPDBs(pdbsNames, ".pdb_modified", head=None)
    for id, pdbItr in pdbsModifiedFiles.items():
        pdbItr.columns = ["record_name", "atom_number", "atom_name", "residue_name", "chain", 
                    "residue_number", "x", "y", "z", "occupancy", "temperature_factor", 
                    "element_symbol"]

    return pdbsNodesFiles, pdbsEdgesFiles, pdbsModifiedFiles

# Read PDB's
pdbsNames = [ "6vyo", "6wji"]
pdbsNodesFiles, pdbsEdgesFiles, pdbsModifiedFiles = read_pdbs(pdbsNames)

#%% Transformar sequencias em nós

#Dividir PDBs em listas de listas
pdbs_tuple_list = fc.prepare_PDBs(pdbsNodesFiles)
#Dicinário de PDBs com listas de cadeias e nós
pdbs_list_of_tuples = fc.pdbs_ids_to_nodes(pdbs_tuple_list)
#Lista de amostras em forma de nós
samples_node_list = fc.df_to_node_list(raw_df_aln)
#%% Tratar Raw data, separa da tabela de alinhamento o pdb das amostras
def fix_pdb(pdb):
    aux = re.split(" |:|[._]", pdb.Id)
    pdb["PdbRange"] = ":".join(aux[-2:])
    pdb.Id = aux[0].replace(">", "")
    pdb["PdabName"] = aux[1]
    pdb["Chain"] = aux[3]
    
    return pdb

def fix_raw_data_one(row):
    temp = row.Id.split(" ")
    raw_id = re.split(":|/|[|]", temp[0])
    pdb_pos = temp[-2]
    seq_pos = temp[-1]
    id = raw_id[0].replace(">", "")
    virus = raw_id[1]
    country = raw_id[2]
    lab_id = raw_id[3] 
    acc_id = raw_id[5]
    collect_data = raw_id[-1]
    
    return [id, row.Seq, pdb_pos, seq_pos, virus, country, lab_id, acc_id, collect_data]


def read_files(path):
    aux = pd.read_csv(path, header=None, sep="\t", names=["Caos"])
    treated = pd.DataFrame(columns=["Id", "Seq"])
    
    seq = ""
    id = ""
    for ii in aux.iterrows():
        if ">" in ii[1].Caos:
            #salva
            treated.loc[len(treated)] = [id, seq]
            #limpa
            seq = ""
            #começa
            id = ii[1].Caos
        else:
            seq += ii[1].Caos
            
    treated.loc[len(treated)] = [id, seq]
        
    return treated.drop(treated.index[0]).reset_index().drop(columns=["index"])

path = "Read/PDBs RIN/align-refs/6vyo_A-gs.fasta.aln"
raw_data = read_files(path)
pdb = raw_data.loc[0]
raw_data.drop(raw_data.index[0], inplace=True)
pdb = fix_pdb(pdb)

processedData = pd.concat([raw_data, pd.DataFrame(columns=["PDBAlign", "SeqAlign", "Virus", "Country", "LabId", "AccessionID", "CollectionData"])], axis=1, sort=False)
for id, row in processedData.iterrows():
    processedData.loc[id] = fix_raw_data_one(row)
processedData.to_csv("Saved/ProcessedData_6vyo.csv", index=False)
pdb.to_csv("Saved/pdb_6vyo.csv")
#pdb = pd.read_csv("Saved/pdb_6vyo.csv")
#processedData = pd.read_csv("Saved/ProcessedData_6vyo.csv")
#%% Atribuir valores de degree ao pdb
