from tqdm import tqdm
import pandas as pd
import numpy as np
from collections import Counter 

def percentage_to_quantity(total, prc):
    """
    Converts percentage to a integer value

    Parameters
    ----------
    total : integer
        Total amount of samples
    prc : integer
        choosed percentage

    Returns
    -------
    integer
        Value equivalent of given percentage in the samples amount

    """
    return (prc*total)/100

def filter_criteria(df, total, condition=1):
    """
    Filtering data using pandas.

    Parameters
    ----------
    df : DataFrame
        DataFrame with the info
    total : integer
        Total of samples
    condition : integer, optional
        Criteria for the filter. The default is 1.

    Returns
    -------
    temp : DataFrame
        Filtered data

    """
    temp = df[df < total - percentage_to_quantity(total, condition)]
    temp.Pos = df.Pos
      
    return temp

def transpose_seq_and_count(column, counted_df):
    """
    Algoritmo para a contagem de nucletotídeos por porsição
    
    ex ilustrativo:
    Posição 1  2  3    4    5  6 ... n    
            A  B  C    D    E  E ... A
            A  B  D    E    E  C ... A
    Total   2A 2B 1C1D 1D1E 2E 1E1C  2A

    Parameters
    ----------
    column : list 
        coluna de sequencias vinda do DataFrame de amostras
    counted_df : DataFrame
        DataFrame com o resultado da contagem

    Returns
    -------
    None.

    """
    #Deixa todas as listas em Upper Case
    aux = column.apply(str.upper)
    #Transformei cada seq.row em lista
    aux = aux.apply(list)
    #Tranformei a coluna seq em uma lista de listas e fiz a Transposta
    aux = np.transpose(list(aux))
    
    #Paraca cada posição(lista) contar Counter
    for idx, row in tqdm(enumerate(aux), total=len(aux)):
        row = Counter(row)
        for nuc, count in row.items():
            if nuc == "U":
                nuc = "T"
            counted_df.loc[idx, nuc] = count