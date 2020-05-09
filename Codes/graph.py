# -*- coding: utf-8 -*-
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt

def three_plots(conditions, list_of_filtered_dfs, name_id, target="Img/", dpi=320, full_or_simple=False):
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
        DPI da imagem criada, o padrão é 320.
    full_or_simple : boolean, optional
        Booleano que define se será mostrado no plot apenas nucleotídeos ou aminoácidos também.
        O padrão é False, mostra apenas nucleotídeos.

    Returns
    -------
    None.

    """
    fig, (ax1, ax2, ax3 ) = plt.subplots(3,1, figsize=(90, 30), dpi=dpi)#, sharex=True)
    aux = [ax1, ax2, ax3]
    alpha = .7
    plt.rcParams.update({'font.size': 20})
    
    for ii, jj, tb in zip(aux, conditions, list_of_filtered_dfs):
      ii.set_title("SNPs "+str(jj)+"%")#, fontweight="bold", size=20)
      ii.xaxis.set_major_locator(ticker.FixedLocator(range(1,29400, 500)))
      
      # ii.plot('Pos', 'A', data=tb, alpha=alpha)
      # ii.plot('Pos', 'C', data=tb, alpha=alpha)
      # ii.plot('Pos', 'T', data=tb, alpha=alpha)
      # ii.plot('Pos', 'G', data=tb, alpha=alpha)
      # ii.plot('Pos', '-', data=tb, alpha=alpha)
      # ii.legend(loc=2)
      # if full_or_simple:
      #   ii.plot('Pos', 'N', data=tb, alpha=alpha)
      #   ii.plot('Pos', 'W', data=tb, alpha=alpha)
      #   ii.plot('Pos', 'R', data=tb, alpha=alpha)
      #   ii.plot('Pos', 'Y', data=tb, alpha=alpha)
      #   ii.plot('Pos', 'S', data=tb, alpha=alpha)
      #   ii.plot('Pos', 'M', data=tb, alpha=alpha)
      #   ii.plot('Pos', 'K', data=tb, alpha=alpha)
      #   ii.plot('Pos', 'V', data=tb, alpha=alpha)
      #   ii.plot('Pos', 'H', data=tb, alpha=alpha)
      
      for col in tb.columns[1:]:
        ii.plot('Pos', col, data=tb, alpha=alpha)
        
      ii.legend(loc=2)
    
    if full_or_simple:
        plt.savefig(target+"Full_"+str(name_id)+"_"+'-'.join(str(x) for x in conditions)+".png", dpi=dpi)
    else:
        plt.savefig(target+"Semi_"+str(name_id)+"_"+'-'.join(str(x) for x in conditions)+".png", dpi=dpi)
    
    #plt.plot()