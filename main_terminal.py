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

def default_name(source, csv = False):
    """
    Define um nome para os arquivos gerados, tanto para o csv quanto para o png

    Parameters
    ----------
    source : str
        Path do arquivo fonte, usado para o processamento
    csv : bool, optional
        Determina se será criado o nome para um png ou para um scv. The default is False.

    Returns
    -------
    str
        o nome para o arquivo que será salvo. 

    """
    # Pega apenas o nome do arquivo passando no caminho: exemplo/pasta/arquivo.csv -> arquivo
    source = source.split("/")[-1].split(".")[0]+"_"
    # Cria um nome único usando a data atual > Resultado_01-01-01_12-12-12
    date = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
    if csv:
        return "Counted_df_"+source+date+".csv"
    else:
        return "Resultado_"+source+date


def create(df, source, target):
    """
    Fata a contagem dos SNPs nas sequencias e salva em csv

    Parameters
    ----------
    df : Dataframe
        DataFrame com as sequências para a análise 
    source : str
        nome que será utilizado no arquivo csv
    target : str
        caminho onde o csv será salvo

    Returns
    -------
    counted_df : DataFrame
        Resultado da contagem de SNPs

    """
    #Df vazio
    counted_df = pd.DataFrame()
    #função que faz a transposta de Seq e conta as ocorrencias
    fc.transpose_seq_and_count(df.Seq, counted_df)
    #Limpa os Nan e corrige a Pos
    counted_df = counted_df.fillna(0).reset_index().rename(columns={"index" : "Pos"})
    counted_df.Pos = counted_df.Pos.apply(lambda x: x+1)
    #Salva o trabalho
    counted_df.to_csv(target+default_name(source, True), index=False)
    
    return counted_df


def docker_flow():
    """
    Cordena todo o processo de leitura do multifasta, criação da tabela de SNPs e png do gráfico

    Returns
    -------
    int
        0, para finalizar o programa

    """
    
    # Argparse
    parser = argparse.ArgumentParser(description="Plot de SNPs")
    parser.add_argument("--option", "-o", required=True, help="ler/criar - se igual a 'ler', recebe um csv já processado para gerar o grafico, se igual a 'criar' recebe fasta, e retorna um csv e um gráfico de SNPs")
    parser.add_argument("--source", "-s", required=True, help="Endereço do arquivo fonte")
    parser.add_argument("--target", "-t", required=True, help="Endereço onde os arquivos gerados serão salvos")                 
    parser.add_argument("--filter", "-f", required=True, nargs='+', type=float, help="Valor float de 0 a 100. Representa em quantos % da população um SNP precisa aparecer para ser contabilizado")
    parser.add_argument("--dpi", "-d", type=int, default=50, help="DPI da imagem gerada")
    
    args = parser.parse_args()
    
    if args.option == "ler":
        #Ler CSV já processado
        counted_df = pd.read_csv(args.source)
    else:
        #Ler multifasta
        df_samples = read_files(args.source)
        # Cria DF com as posições contadas
        counted_df = create(df_samples, source=args.source, target=args.target)
        
        
    # Processar
    # valores usados como filtro, passados no terminal
    filter = args.filter
    #Lista onde será armazenado os DataFrames com os valores filtrados
    list_of_filtered_dfs = []

    for ii in filter:
        list_of_filtered_dfs.append(fc.filter_criteria(counted_df, counted_df.loc[0].sum()-1, ii))
 
    # Gerar gráficos
    gf.three_plots(filter, list_of_filtered_dfs, default_name(args.source), target=args.target, dpi=args.dpi)
    
    return 0

sys.exit(docker_flow())