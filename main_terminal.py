# -*- coding: utf-8 -*-
"""
Created on Tue May  5 13:01:54 2020

@author: ivan
"""
#System Imports
import sys
import argparse
import pandas as pd
from datetime import datetime

#My Imports
import Codes.read as rdF
import Codes.function as fc
import Codes.graph as gf

def read_files(path):
    """
    Lê os arquivos de entrada

    Parameters
    ----------
    path : TYPE
        DESCRIPTION.

    Returns
    -------
    df_aln : TYPE
        DESCRIPTION.

    """
    aux = path.split("/")
    path = "/".join(aux[:-1])+"/"
    name = aux[-1]
    raw_aln = rdF.read_aligned_files(path, name)
	#Converto alinhado para um DF
    df_aln = rdF.seq_to_df(raw_aln)
    
    return df_aln

def default_name(source=""):
    """
    Gera um nome padrtão para as amostras, caso não seja dado um nome

    Parameters
    ----------
    plot : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # Pega apenas o nome do arquivo passando no caminho: exemplo/pasta/arquivo.csv -> arquivo
    source = source.split("/")[-1].split(".")[0]+"_"
    # Cria um nome único usando a data atual > Resultado_01-01-01_12-12-12
    date = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
    result = "Resultado_"
    """
        Preciso revisar
    """
    if source == "":
        return result+date+".csv"
    else:
        return result+source+date


def create(df_aln, source, target=default_name(), read=True):
    """
    Lê arquivo de contagem já pronto ou cria um

    Parameters
    ----------
    df_aln : TYPE
        DESCRIPTION.
    target : TYPE, optional
        DESCRIPTION. The default is default_name().
    read : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    counted_df : TYPE
        DESCRIPTION.

    """
    if read:
        #lê salvo
        counted_df = pd.read_csv(source)
    else:
        #Df vazio
        counted_df = pd.DataFrame()
        #função que faz a transposta de Seq e conta as ocorrencias
        fc.transpose_seq_and_count(df_aln.Seq, counted_df)
        #Limpa os Nan e corrige a Pos
        counted_df = counted_df.fillna(0).reset_index().rename(columns={"index" : "Pos"})
        counted_df.Pos = counted_df.Pos.apply(lambda x: x+1)
        #Salva o trabalho
        counted_df.to_csv(target+default_name(), index=False)
    
    return counted_df


def docker_flow():
    """
    Work flow do trabalho preparado para terminal e Docker

    Returns
    -------
    int
        DESCRIPTION.

    """
    #Parser
    parser = argparse.ArgumentParser(description="Plot de SNPs")
    parser.add_argument("--option", "-o", required=True, help="ler/criar - se igual a 'ler', recebe um csv já processado para gerar o grafico, se igual a 'criar' recebe fasta, e retorna um csv e um gráfico de SNPs")
    parser.add_argument("--source", "-s", required=True, help="Endereço do arquivo fonte")
    parser.add_argument("--target", "-t", required=True, help="Endereço onde os arquivos gerados serão salvos")                 
    parser.add_argument("--filter", "-f", required=True, nargs='+', type=float, help="Valor float de 0 a 100. Representa em quantos % da população um SNP precisa aparecer para ser contabilizado")
    parser.add_argument("--dpi", "-d", type=int, default=50, help="DPI da imagem gerada")
    
    args = parser.parse_args()
    
    if args.option == "ler":
        counted_df = pd.read_csv(args.source)
    else:
        #Ler aqrquivos
        df_samples = read_files(args.source)
        counted_df = create(df_samples, target=args.target, read=False)
    
    #Plotar
    ## Melhorar o filtro para receber quantidade variável de valores
    ## Atual: 3 valores float
    filter = args.filter
    list_of_filtered_dfs = []
    
    for ii in conditions:
        list_of_filtered_dfs.append(fc.filter_criteria(counted_df, counted_df.loc[0].sum()-1, ii))
    	    
    gf.three_plots(filter, list_of_filtered_dfs, default_name(args.source), target=args.target, dpi=args.dpi)
        
    return 0

sys.exit(docker_flow())