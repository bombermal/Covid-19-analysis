# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:54:54 2020

@author: ivan
"""
#System Imports
import sys
import argparse
import pandas as pd
#from datetime import datetime
from tqdm import tqdm
import re
from Bio.SeqUtils import seq1

#My Imports
import Codes.read as rdF
import Codes.function as fc
import Codes.graph as gf
import Codes.node as nd
import Codes.smithWaterman as sw

def read_raw_files(path):#
    """
    Lê os arquivos de entrada

    Parameters
    ----------
    path : str
        Caminho do arquivo.

    Returns
    -------
    df_aln : DataFrame
        DataFrame com o alinhamento.

    """
    aux = path.split("/")
    path = "/".join(aux[:-1])+"/"
    name = aux[-1]
    raw_aln = rdF.read_aligned_files(path, name)
 	#Converto alinhado para um DF
    df_aln = rdF.seq_to_df(raw_aln)
    
    return df_aln

def read_Or_create(path, ler = True, df = ""):#
    """
    Comanda a leitura ou a criação do arquivo CSV contendo a contagem dos aminoácidos por posição

    Parameters
    ----------
    path : str
        string contendo a posição do arquivo que será lido, no caso de um CSV anteriormente criado
    ler : bool, optional
        True: Lê o CSV anteriormente criado, False: Lê o fasta e cria o CSV, by default True
    df : str, optional
        DataFrame utilizado para a contabilização dos SNPs, by default ""

    Returns
    -------
    DataFrame
        Contém a contagem de aminoácidos por posição     
    """
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

def plot(conditions, samplesSize, counted_df, name_id, target="Img/", dpi=50):
    """
    Recebe o dataframe de contagem de posições, filtra e aciona a o plot dos gráficos

    Parameters
    ----------
    conditions : list
        Condições para o filtro
    samplesSize : int
        número de amostras
    counted_df : DataFrame
        DataFrame fonte dos dados
    name_id : str
        nome utilizado como ID para a definição dos nomes dos arquivos que serão salvos
    target : str, optional
        caminho onde será salvo os arquivos criados, by default "Img/"
    dpi : int, optional
        Dpi da imagem, by default 50
    """
    
    list_of_filtered_dfs = []
    
    for ii in conditions:
      list_of_filtered_dfs.append(fc.filter_criteria(counted_df, samplesSize, ii))
      
    gf.three_plots(conditions, list_of_filtered_dfs, name_id, target=target, dpi=dpi)

def read_pdbs(pdbsNames):
    """
    Recebe uma lista de com os arquivos dos pdbs que estão salvos na pasta read/PDBs RIN/

    Parameters
    ----------
    pdbsNames : list
        Lista com os nomes dos pdbs

    Returns
    -------
    Dict
        Retorna 3 dicionários com os pdbs carregados em DataFrames
    """
    pdbsNodesFiles = rdF.readPDBs(pdbsNames, ".pdb.nodes")
    pdbsEdgesFiles = rdF.readPDBs(pdbsNames, ".pdb.edges")

    pdbsModifiedFiles = rdF.readPDBs(pdbsNames, ".pdb_modified", head=None)
    for id, pdbItr in pdbsModifiedFiles.items():
        pdbItr.columns = ["record_name", "atom_number", "atom_name", "residue_name", "chain", 
                    "residue_number", "x", "y", "z", "occupancy", "temperature_factor", 
                    "element_symbol"]

    return pdbsNodesFiles, pdbsEdgesFiles, pdbsModifiedFiles

def read_betwenness(path):
    betwenness_raw = pd.read_csv(path)
    betwenness_dict = {}
    for ii in betwenness_raw.filename.unique():
        key = str.upper(ii.split(".")[0])
        betwenness_dict[key] = betwenness_raw[betwenness_raw.filename == ii][["node", "aminoAcid", "clusteringCoef", "betweennessWeighted"]]
        betwenness_dict[key].rename(columns = {"node": "NodeId", "aminoAcid": "Residue", "clusteringCoef": "ClusteringCoef",
        "betweennessWeighted": "Betweennessweighted"}, inplace=True)

    return betwenness_dict

def join_pdb_betw(pdbDict, betwDict):
    for key in pdbDict.keys():
        pdbDict[key] = pd.merge(pdbDict[key], betwDict[key], on=['NodeId', "Residue"])

def fix_pdb(pdb):
    """
    Separa as partes do arquivo pdb utilizando REGEX

    Parameters
    ----------
    pdb : DataFrama
        DataFrame com os dados que serão organizados

    Returns
    -------
    DataFrame
        DataFrame processado e limpo
    """
    aux = re.split(" |:|[._]", pdb.Id)
    pdb["PdbRange"] = ":".join(aux[-2:])
    pdb.Id = aux[0].replace(">", "")
    pdb["PdbName"] = aux[1]
    pdb["Chain"] = aux[3]
    
    return pdb

def fix_raw_data_one(row):
    """
    Utilização de REGEX para a separação das informações das amostras, uma linha por vez

    Parameters
    ----------
    row : DataFrame
        linha do DataFrame que será analizada e limpa

    Returns
    -------
    DataFrame
        Linha pós processada e organizada
    """
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
    """
    Lê arquivo multifasta e separa as amostras

    Parameters
    ----------
    path : str
        caminho do arquivo multifasta

    Returns
    -------
    DataFrame
        DataFrame com amostras separadas
    """
    aux = pd.read_csv(path, header=None, sep="\t", names=["Caos"])
    treated = pd.DataFrame(columns=["Id", "Seq"])
    
    seq = ""
    id = ""
    for ii in tqdm(aux.iterrows(), total=len(aux)):
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

def pdbRead_or_Create(alnPath, csvPath="", dfPath="", ler=True):
    """
    Comanda a leitura de arquivos CSVs anteriormente criados ou o novo processamento dos arquivos fasta 
    para a criação de um CSV.

    Parameters
    ----------
    alnPath : str
        Caminho do arquivo .fasta.aln que será utilizado para a criação do CSV
    csvPath : str
        Caminho do arquivo CSV anteriormente criado com informações apenas do PDB
    dfPath : str
        Caminho do arquivo CSV anteriormente criado com informações das amostras
    ler : bool, optional
        True: lê csv anteriormente criado, False: Lê fasta.aln e cria arquivo csv, by default True

    Returns
    -------
    DataFrame
        Dois dataframes com informações dos pdbs e das amostras
    """
    if ler:
        pdb = pd.read_csv(csvPath)
        processedData = pd.read_csv(dfPath, skiprows=0)
    else:
        raw_data = read_files(alnPath)
        pdb = raw_data.loc[0]
        raw_data.drop(raw_data.index[0], inplace=True)
        pdb = fix_pdb(pdb)

        processedData = pd.concat([raw_data, pd.DataFrame(columns=["PDBAlign", "SeqAlign", "Virus", "Country", "LabId", "AccessionID", "CollectionData"])], axis=1, sort=False)
        for id, row in tqdm(processedData.iterrows(), total=len(processedData)):
            processedData.loc[id] = fix_raw_data_one(row)

        processedData.to_csv("Saved/ProcessedData_6vyo.csv", index=False)
        pdb.to_csv("Saved/pdb_6vyo.csv", index=True)

    return pdb, processedData

def filter_df_on_pdbsDict(dataDict, pdbDict):
    """
    Recebe uma um dicionário com os pdbs+amostras alinhados já lidos e um outro dicionário com pdbs apenas e faz uma filtragem
    das informações dos PDBs baseado no tipo de cadeia que foi ultilizada no alinhamento dos pdbs+amostras

    Parameters
    ----------
    dataDict : dict
        Dicionário com os pdbs+amostras dos arquivos multifasta
    pdbDict : dict
        Dicionário com os pdbs puros, de onde serão selecionados as partes com as cadeias desejadas

    Returns
    -------
    dict
        Dicionário com DataFrames já filtrados
    """
    temp = {}
    for key, value in dataDict.items():
        tempChain = value[0].Chain
        upperName = str.upper(value[0].PdbName)
        tempDict = pdbDict[upperName]
        temp[upperName] = tempDict[tempDict.Chain == tempChain]
    
    return temp

def pdb_file_to_list_of_nodes(df):
    """
    Converte PDB em uma lista de nós, contendo todos os dados que serão usados posteriormente no alinhamento

    Parameters
    ----------
    df : DataFrame
        DataFrame contendo os dados do pdb

    Returns
    -------
    list
        lista de nós com os valores da tabela agragados a cada resíduo(Degree, Id, nome, posição, posição do alvo)
    """
    temp = []
    for num, row in df.iterrows():
        temp.append(nd.node(seq1(row.Residue), row.Degree, row.NodeId, row.Position, -1, row.ClusteringCoef, row.Betweennessweighted))

    return temp

def simple_seq_to_list_of_nodes(pdb):
    """
    Transforma uma sequencia de residuos em um nó, por vir de um arquivo simples, não contem valores como:
    degree e id.

    Parameters
    ----------
    pdb : DataFrame
        DataFrame contendo as dados da amostra

    Returns
    -------
    list
        lista de nós contendo os valores agregados a cada resíduo(residuo e posição)
    """
    temp = []
    for num, val in enumerate(pdb.Seq):
        temp.append(nd.node(val, -1, "sampleNoId", num, -1, 0, 0))
    
    return temp

def sequences_to_nodes_list(filtered_dfs, processed):
    #Dicionário com pdbs Filtrados na forma de nós, prontos para alinhar contra pdb do multifasta
    nodes_list_from_pdbs = {}
    for ii in filtered_dfs.keys():
        nodes_list_from_pdbs[ii] = pdb_file_to_list_of_nodes(filtered_dfs[ii])

    nodes_list_from_sample = {}
    for ii in processed.keys():
        nodes_list_from_sample[ii] = simple_seq_to_list_of_nodes(processed[ii][0])

    return nodes_list_from_pdbs, nodes_list_from_sample

def align_head_pdb(filtered_dfs, processed, nodes_list_from_pdbs, nodes_list_from_sample):
    """
    Alinha a sequencia do multifasta simples(sem degree) com a sequencia do pdb

    Parameters
    ----------
    filtered_dfs : dict
         Dicionário com DataFrames já filtados
    processed : dict

    nodes_list_from_pdbs : dict
        [description]
    nodes_list_from_sample : dict
        [description]

    """
    for key in nodes_list_from_pdbs.keys():
        obj = sw.smithWaterman()
        seqIds, seqPos, cover = obj.constructor(2, -1, -1, nodes_list_from_sample[key], nodes_list_from_pdbs[key], False, False)
        #4 salvar resultado em df
        #Lista de Ids dos aminoácidos alinhados
        processed[key][0]["AlnResult"] = "|".join(seqIds)
        #Lista de Degrees/Clusterin/Betweness criadas a partir da lista de Ids
        for name, ii in zip(["AlnDgree", "AlnClust", "AlnBetw"], ["Degree", "ClusteringCoef", "Betweennessweighted"]):
            processed[key][0][name] = "|".join([str(filtered_dfs[key][filtered_dfs[key].NodeId == x][ii].values[0]) for x in seqIds])
  