#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np

############## Se produce una matriz para todos los clusters de la busqueda ##############################

#Se crea el archivo order_out_mci que contiene todos los genes en una sola linea
#este conserva el orden original de los clusters, de mayor a menor
#get_ipython().system('tr "\\n" "\\t" < data_out_mci > order_out_mci')
#get_ipython().system("sed -i 's/\\t$//' order_out_mci")

#Se leen los archivos usando pandas
genes1 = pd.read_csv('order_out_mci', sep='\t')
data1 = pd.read_csv('evalue_filter2', sep='\t', lineterminator='\n', header=None)

#Se crea una matriz vacia que respeta el orden original de los genes
cols1=rows1=[x for x in genes1]
matrix1=pd.DataFrame(99, columns=cols1, index=rows1, dtype='float')

#Se rellena la matriz con los valores de la tercer columna de data1
for i_index, i_row in data1.iterrows():
        
    valor_actual = matrix1[i_row[0]][i_row[1]] 
    valor_transpuesto = matrix1[i_row[1]][i_row[0]]
    valor_final=min(valor_actual, valor_transpuesto)
    
    if valor_final > i_row[2]:
        matrix1[i_row[0]][i_row[1]] = i_row[2]
        matrix1[i_row[1]][i_row[0]] = i_row[2]
    
    #if valor_actual > i_row[2]: 
    #    matrix1[i_row[0]][i_row[1]] = i_row[2]
    #if valor_transpuesto > i_row[2]:
    #    matrix1[i_row[1]][i_row[0]] = i_row[2]

#Se transforman los valores de la matris usando -log10(x + 1e-200)
def minus_log(x = float):
    if x != 99:
        return(-(np.log(x+10**(-200))/np.log(10)))
    else:
        return(0)
log_matrix1=matrix1.applymap(minus_log)

#########################################################################################################

import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sc 

######## Se crea el heatmap a partir de la matriz que contiene todolos clusters de la busqueda ##########
sns.set(rc={'figure.figsize':(100,100)})

heat1 = sns.heatmap(log_matrix1, annot=False)

fig = heat1.get_figure()
fig.savefig("All_clusters_heat_map.png") 
plt.close()

###################################### Se obtiene una lista de indices #########################################

index_list_temp=open('data_index').read().split('\n')
ind_list=[]
for i in index_list_temp:
    if i != '':
        ind_list.append(int(i))

###########################   Se realizan dendrogramas+heatmaps para cada cluster   ############################

#Obtengo el archivo 'gene_number' indicando la cantidad de genes por cluster:
#get_ipython().system("awk '{print NF}' data_out_mci > gene_number")
#get_ipython().system('tr "\\n" "\\t" < gene_number > A')
#get_ipython().system('mv A gene_number')

#### Funcion de StackOverflow para convertir dendrogramas de Scipy a formato Newick ####

def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

#### - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - ####


#Guardo la cantidad en la lista 'gene_num' la cantidad dde genes por cluster
gene_num=open('gene_number')
gene_num=list(gene_num)[0].split('\t')
if gene_num[len(gene_num)-1] == '':
    gene_num.pop(len(gene_num)-1)
#print("gene_num:")
#print(gene_num)
#A partir del df log_matrix1, selecciono las submatrices correspondientes a
#cada cluster y elaboro un dendograma
prev_num=0
num=-1
index=-1
for number in gene_num:
    num+=int(number)
    index+=1
    #print(num, prev_num)
    
    ## Se descartan clusters dos o menos elementos
    if num - prev_num > 0: 
    
        print("Creando drendrograma para el cluster: "+str(ind_list[index])+", "+str(num-prev_num + 1)+" genes")
        sub_log_matrix=log_matrix1.iloc[prev_num:num+1,prev_num:num+1]
        prev_num = num+1
    
    ### Heatmap y Dendrograma ###
        
        clustermap = sns.clustermap(sub_log_matrix, xticklabels=False, yticklabels=True,
                         cbar_pos=(0.8, 0.8, 0.05, 0.18) ,square=True,
                         dendrogram_ratio=(0.00001,0.2), vmax=100, annot=True)
        #plt.show()
        plt.savefig('Cluster_'+str(ind_list[index])+'_Heatmap_n_Dendrogram')
        plt.close()
        
    ### Dendrograma ###
        
        link = clustermap.dendrogram_row.linkage
        genes = clustermap.dendrogram_row.data.index

        plt.figure(figsize=(25*1.5, 10*1.5))
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('Genes')
        plt.ylabel('Distance')

        sc.set_link_color_palette(None)  
        sc.dendrogram(
               link,
               leaf_rotation=90.,  # rotates the x axis labels
               leaf_font_size=8.,  # font size for the x axis labels
               color_threshold=None,
               labels=genes
        )
        
        plt.savefig("Cluster_" + str(ind_list[index])+"_Dendrogram")
        plt.close()

    ### Guardar el dendrograma en formato Newick ###
    
        
        tree = sc.to_tree(link, False)
        newick = get_newick(tree, tree.dist, genes)

        text_file = open("Cluster_" + str(ind_list[index]) + "_Dendrogram.newick", "w")
        n = text_file.write(newick)
        text_file.close()
        
    else:
        print("Salteandose Dendrograma y Heatmap para el cluster pequeño Nº: "+ str(ind_list[index]))    
        prev_num = num+1

    #plt.imshow(sub_log_matrix, cmap='hot', interpolation='nearest')
    #plt.show()


### Cosas que quedan por hacer:
### -Detectar y descartar outliers (Opcionalmente)
### -Marcar de ser posible la seleccion
### -Marcar con un mapa de colores cada cluster 
### -Marcar con colores los genes de cada organismo
### -Hacer mas legible los genes dependiendo del caso




