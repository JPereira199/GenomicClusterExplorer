#!/usr/bin/env python
# coding: utf-8

##Histograma del los encabezados
##Input:
## file1 y file2: archivos en formato clust que contengan una serie de 
## valores numericos
## index_file: Archivo con el indice de clusters
## xlabel: Nombre del ejex
## fig_name: Nombre de salida de las imagenes
## gpb: Genes per bin (1, 2, 3, ...) 

## Agregar opcion que permita sustituir gpb por total bins

##Output: 
##- Histogramas con los datos de file1 y file2

import matplotlib.pyplot as plt
import sys

file1=sys.argv[1] 	#'clust_largos_whole'
file2=sys.argv[2] 	#'clust_largos_selection'
index_file=sys.argv[3] 	#'data_index'
xlabel=sys.argv[4] 	#'Largo en AA'
fig_name=sys.argv[5] 	#'largos'
gpb=int(sys.argv[6])	# integrer


#Se carga el archivo con los indices
ind_file=open(index_file)
read_indx=ind_file.read()
clust_index=read_indx.split("\n") 
clust_index.pop(len(clust_index)-1)

#Se cargan los archivos de clusters
clust_file_whole=open(file1)
clust_file_selection=open(file2)

#Se leen los archivos
read1=clust_file_whole.read()
read2=clust_file_selection.read()

#Se divide cada cluster en un elemento de una lista
clust_list1=read1.split('\n')
clust_list2=read2.split('\n')

#Se quita el ultimo elemento de la lista 
#(Se asume que es '')
clust_list1.pop(len(clust_list1)-1)
clust_list2.pop(len(clust_list2)-1)

##Se grafica un par de histogramas para cada cluster
for i in range(len(clust_list1)):
    
    #Se separan los elementos de cada cluster
    values_list1=clust_list1[i].split('\t')
    values_list2=clust_list2[i].split('\t')
    
    #Se convierten los elementos de str a int
    values1=[]
    for j in range(len(values_list1)):
        values1.append(float(values_list1[j]))
    
    values2=[]
    for j in range(len(values_list2)):
        values2.append(float(values_list2[j]))
    
    values1.sort()
    plt.hist(values1, bins=len(values1)//gpb+1, range=(values1[0], values1[len(values1)-1]))
    plt.hist(values2, bins=len(values1)//gpb+1, range=(values1[0], values1[len(values1)-1]))
    plt.title('Cluster '+ str(clust_index[i]))
    plt.ylabel('Genes')
    plt.xlabel(xlabel);

    plt.savefig('Cluster_'+str(clust_index[i]+'_'+fig_name))
    plt.close()
    

##Se grafica un histograma considerando el total de los genes
all_clust1_tmp='\t'.join(clust_list1)
all_clust2_tmp='\t'.join(clust_list2)

all_clust1=all_clust1_tmp.split('\t')
all_clust2=all_clust2_tmp.split('\t')

all_values1=[]
for i in all_clust1:
    all_values1.append(float(i))

all_values2=[]
for i in all_clust2:
    all_values2.append(float(i))

plt.hist(all_values1, bins=len(all_values1)//gpb+1, range=(min(all_values1), max(all_values1)), label='Cluster total')
plt.hist(all_values2, bins=len(all_values1)//gpb+1, range=(min(all_values1), max(all_values1)), label='Seleccion')
plt.title('All Clusters')
plt.ylabel('Genes')
plt.xlabel(xlabel)
plt.legend()

plt.savefig('All_clusters_'+fig_name)
plt.close()


