# -*- coding: utf-8 -*-
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt

def three_plots(conditions, list_of_filtered_dfs, name_id, target="Img/", dpi=100):
    """
    Plota gráfico da análise de SNPs    

    Parameters
    ----------
    conditions : list
        valores usados como filtro nos DataFrames
    list_of_filtered_dfs : list
        Lista de DataFrames
    name_id : str
        Nome que será dado ao arquivo criado
    target : str, optional
        Caminho onde será salvoo gráfico criado. The default is "Img/".
    dpi : int, optional
        Valor de dpi do gráfico, quanto maior, mais qualidade e tamanho. The default is 100.

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