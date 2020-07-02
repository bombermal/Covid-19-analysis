# Covid-19-analysis
### Com a chegada do novo Covid19 se fez necessário que todos os cientístas do pais, que podem ajudar, cedam um pouco do seu tempo e conhecimento para pesquisar e encontrar novos dados que possam ajudar no combate dessa pândemia. 
#### O intuito dessa pesquisa é ultilizar RINs( Residue Iteration Network ) para identificar quais regiões do genoma do vírus são mais propensos a mudanças.

![fluxo 2(Under revision)](https://github.com/bombermal/Covid-19-analysis/blob/master/Uml/Covid%20-%20Total.png)
## Pt 1
 - [X] Read aligned Fasta files
    - [X] Load Fasta on DataFrame
 - [X] Counting SNPs
    - [X] Save counted DF as csv
 - [X] Plot SNPs positions
 
## Pt 2
 
 - [X] Read PDBs
    - [X] Load PDBs as DataFrames
    - [X] Convert PDBs list to nodes
    - [X] Convert samples list to nodes
 - [X] Align PDBs with samples
 - [X] Determine values of degree
 - [X] Determine values of betwenness
 - [X] Determine values of clustering coeficient
 - [X] Plot degree dristribution
 - [X] Plot betwenness dristribution
 - [X] Plot clustering coeficient
 - [ ] Find a way to correlact SNPs with found data
 - [ ] Plot degree dristribution with SNPs
 - [ ] Plot betwenness dristribution with SNPs
 - [ ] Plot clustering coeficient with SNPs