# -*- coding: utf-8 -*-
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt

def three_plots(conditions, list_of_filtered_dfs, name_id, target="Img/", dpi=100):
    """
    Plota o gráfico das análises

    Parameters
    ----------
    conditions : list
        valores usados como filtro para a análise, valor float entre 0 e 100
        Representa a porcentagem da população que apresenta uma modificação para 
        que ela seja contada como SNP (0.5%, 1%, 2%, ...)
    list_of_filtered_dfs : list
        Lista de DataFrames já pre-filtrados usados para o Plot
    name_id : str
        Nome único para o arquivo que será salvo após o processamento
    target : str, optional
        Caminho onde será salvo o arquivo. O padrão é "Img/".
    dpi : int, optional
        DPI da imagem criada, o padrão é 100.
    full_or_simple : boolean, optional
        Booleano que define se será mostrado no plot apenas nucleotídeos ou aminoácidos também.
        O padrão é False, mostra apenas nucleotídeos.

    Returns
    -------
    None.

    """
    fig, (ax1, ax2, ax3 ) = plt.subplots(len(conditions),1, figsize=(90, 30), dpi=dpi)#, sharex=True)
    aux = [ax1, ax2, ax3]
    alpha = .7
    plt.rcParams.update({'font.size': 20})
    
    for ii, jj, tb in zip(aux, conditions, list_of_filtered_dfs):
      ii.set_title("SNPs "+str(jj)+"%")#, fontweight="bold", size=20)
      ii.xaxis.set_major_locator(ticker.FixedLocator(range(1,29400, 500)))
           
      for col in tb.columns[1:]:
        ii.plot('Pos', col, data=tb, alpha=alpha)
        
      ii.legend(loc=2)
    
    plt.savefig(target+str(name_id)+"_"+'-'.join(str(x) for x in conditions)+".png", dpi=dpi)

    #plt.plot()