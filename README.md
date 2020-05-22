# Covid-19-analysis
### Com a chegada do novo Covid19 se fez necessário que todos os cientístas do pais, que podem ajudar, cedam um pouco do seu tempo e conhecimento para pesquisar e encontrar novos dados que possam ajudar no combate dessa pândemia. 
#### O intuito dessa pesquisa é ultilizar RINs( Residue Iteration Network ) para identificar quais regiões do genoma do vírus são mais propensos a mudanças.

## Pt 1
 ![fluxo 1](Uml/Uml/Covid - Terminal.png)
 
 - [X] Read aligned Fasta files
    - [X] Load Fasta on DataFrame
 - [X] Counting SNPs
    - [X] Save counted DF as csv
 - [X] Plot SNPs positions
 
## Pt 2
 ![fluxo 2(partial)](Uml/Uml/Covid - Total.png)
 - [X] Read PDBs
    - [X] Load PDBs as DataFrames
    - [X] Convert PDBs list to nodes
    - [ ] Convert samples list to nodes
 - [ ] Align PDBs with samples
 - [ ] Determine values of degree and betwenness
 - [ ] Plot degree and betwenness dristribution
