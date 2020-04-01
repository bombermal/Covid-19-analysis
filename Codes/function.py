from tqdm import tqdm
import pandas as pd
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing as mp

def workDivide(dfWork):
    """
    Divide the DataFrame in a list of rows, to use in multiprocess

    Parameters
    ----------
    dfWork : pd.DataFrame
        Receive a DataFrame as input.

    Returns
    -------
    listWork : list()
        A list that receaves all rows of the DataFrame as a single element

    """

    listWork = []
    for i in tqdm(range(len(dfWork))):
        listWork.append(dfWork.iloc[i])
    return listWork

def single_nuc_count(sample):
    """
    Receave a string "XCDAE" and return a dataframe with the position of each
    charcater
    
    Pos, Nuc
    1 - X
    2 - C
    3 - D
    4 - A
    5 - E
    ....

    Parameters
    ----------
    sample : String
        A string that will be divided

    Returns
    -------
    aux : DataFrame
        Dataframe cointaning the position of each character in the string

    """
    aux = pd.DataFrame(columns=["Pos", "Nuc"])

    for idx, value in enumerate(sample.Seq):
        aux.loc[len(aux)] = [idx+1, value]  
    return aux

def df_to_list(df):
    """
    Converts the DataFrame in a list of rows

    Parameters
    ----------
    df : DataFrame
        A pandas Dataframe

    Returns
    -------
    temp : list
        each element is a DataFrame row

    """
    temp = []
    for idx, row in df.iterrows():
        temp.append(row)
        
    return temp


def process_4(partial_list_of_samples, num):
    """
    Will convert the list of rows in a dataframe that contains two columns
    Pos, Nuc -> Position and Nucleotide
    
    Save the work in a csv

    Parameters
    ----------
    partial_list_of_samples : list 
        List of rows from a DataFrame
    num : Integer
        Id to identify the saved csv file

    Returns
    -------
    None.

    """
    tidy_count_df2 = pd.DataFrame(columns=["Pos", "Nuc"])
    poolSize = mp.cpu_count()
    pool = ThreadPool(poolSize)
    
    for ii in tqdm(pool.imap_unordered(single_nuc_count, partial_list_of_samples), total=len(partial_list_of_samples)):
        tidy_count_df2 = pd.concat([tidy_count_df2, ii])
    pool.close()
    pool.join()
    
    tidy_count_df2.to_csv("df_pos_nuc"+str(num)+".csv", index=False)
    
def work_on_specific_files(divided_df_file_names):
    """
    Give a list of filenames and they'll be loaded and sent to process_4
    
    Parameters
    ----------
    divided_df_file_names : List
        List of filenames, that will be used to read a csv and and call
        process_4
    
    Returns
    -------
    None.
    
    """
    for key, value in enumerate(divided_df_file_names):
        df = pd.read_csv(value)
        
        process_4(df_to_list(df), key+1)
        
def list_of_poli_positions(df, samples_total, prc):
    """
    OLD --- Crate a list with the positions that correspond to a criteria

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    samples_total : TYPE
        DESCRIPTION.
    prc : TYPE
        DESCRIPTION.

    Returns
    -------
    aux : TYPE
        DESCRIPTION.

    """
    prc_val = (prc*samples_total)/100
    aux = pd.DataFrame(columns=["Poli_Pos"])
    
    for ii in tqdm(range(df.Pos.max())):
        temp = df[df.Pos == ii].Nuc.value_counts()
        condition = samples_total - temp.max()
        if condition > prc_val:
            aux.loc[len(aux)] = ii
            
    return aux

def plus_one(df, row):
    """
    Add +1 e a DataFrame cell

    Parameters
    ----------
    df : DataFrame
        Receave a DataFrame to work
    row : Series
        use the row info to find the target cell

    Returns
    -------
    None.

    """
    df.loc[row.Pos, row.Nuc] += 1
    
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