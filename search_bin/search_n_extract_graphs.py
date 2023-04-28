#!/usr/bin/env python
# coding: utf-8

# In[14]:

#######
##Objetivos
## Crear graficos del tipo
## Barplot de atributos en la tabla:
##  -Que tome en cuenta todos los tipos de atributos que aparecen en la busqueda
##  -Que el usuario puede decidir sobre que atributos del cuerpo hacer histogramas
##  -Opcion para decidir si usar nombres cajas especificos para un atributo, se debe incluir un box donde entre
##   todo lo que no coincidio con la caja
##  -Induvidual para cada cluster, que incluya la seleccion y el total
##  -Histograma que incluya, seleccionados + complementos
## 
## Histograma del encabezado (data_cl, jceegee, total de genes)
##  -Seleccionados, totales en un archivo


def atribute_barplot ( att_file = str, index_file = str ,xlabel = '', img_name = str  ):
    
    #att_file='cluster_chroms'
    #xlabel='Cromosomas'
    #clust_index=['44','45']
    
    ind_file=open(index_file)
    read_indx=ind_file.read()
    clust_index=read_indx.split("\n")
    clust_index.pop(len(clust_index)-1)

    
    clust_file_whole=open(att_file+'_whole')
    clust_file_selection=open(att_file+'_selection')

    read1=clust_file_whole.read()
    read2=clust_file_selection.read()

    ##Lista con un cluster por elemento
    clust_list1=read1.split("\n") 
    clust_list2=read2.split("\n")

    all_clust1_tmp=' '.join(clust_list1)  
    all_clust2_tmp=' '.join(clust_list2)

    ## Cada elemento de la busqueda es un elemento de la lista
    all_clust1=all_clust1_tmp.split(" ") 
    all_clust2=all_clust2_tmp.split(" ") 

    ############################# Todos los nombres y elementos de la busqueda #####################################

    ## Se crea un diccionario que contiene los nombres de todos los elementos
    ## que aparecieron en la busqueda. A cada elemento se le asigna un valor que 
    ## indica la cantidad de veces que se encuentra repetido
    
    clust_dict1 = {}
    clust_dict2 = {}

    for i in all_clust1:  #Nc_11 Nc_11I Nc_11 Nc_11 
        if i in clust_dict1:
            clust_dict1[i] += 1
        else:
            clust_dict1[i] = 1

    for i in all_clust2:
        if i in clust_dict2:
            clust_dict2[i] += 1
        else:
            clust_dict2[i] = 1

    clust_dict1.pop('')
    clust_dict2.pop('')

    ### Se agregan los indices faltantes de clust_dict1 a clust_dict2
    for i in clust_dict1:
        if i not in clust_dict2:
            clust_dict2[i] = 0
        


    ## Posteriormente obtiene un vector ordenado (keys) con los nombres del diccionario
    ## Y otro vector (values) que contiene los valores correspondiente a cada elemento de keys
    values1=[]
    values2=[]
    keys=sorted(clust_dict1)
    for i in keys:
        values1.append(clust_dict1[i])
        values2.append(clust_dict2[i])
    

    import matplotlib.pyplot as plt

    plt.bar(keys,values1,alpha=0.7)
    plt.title('Cluster ')
    plt.ylabel('N. de Genes')
    plt.xlabel( xlabel );


    plt.bar(keys,values2,color='red',alpha=0.7)
    plt.title('Cluster ')
    plt.ylabel('N. de Genes')
    plt.xlabel( xlabel);
    
    plt.savefig('All_clusters_' + img_name)
    plt.close()

    ####################################### Cantidad de elementos en cada cluster #####################################
  

    ## Se dejan en 0 todos los valores del diccionario
    for i in clust_dict1:
        clust_dict1[i]=0
        clust_dict2[i]=0
    
    if '' in clust_list1:
        clust_list1.remove('')
    if '' in clust_list2:
        clust_list2.remove('')
    
    indx=0
    for i in range(len(clust_list1)):

        ## Numero de elementos en cada cluster
        elements = clust_list1[i].split(" ")
        elements.remove('')
        for j in elements:
            if j in clust_dict1:
                clust_dict1[j] += 1
            else:
                print("ERROR unespected element: ", j)
    
        elements = clust_list2[i].split(" ")
        #elements.remove('')
        for j in elements:
            if j in clust_dict2:
                clust_dict2[j] += 1
            else:
                print("ERROR unespected element: ", j)  
 
    
        values1=[]
        values2=[]  
        partial = True
        if partial == True:
            keys=[]
            for i in clust_dict1:
                if  clust_dict1[i] != 0:
                    values1.append(clust_dict1[i])
                    values2.append(clust_dict2[i])
                    keys.append(i)
                
        else:
            for i in clust_dict1:
                values1.append(clust_dict1[i])
                values2.append(clust_dict2[i])
        
    
        
        ##Se hace un grafico de barras para cada cluster
        plt.bar(keys,values1, label="Cluster Total")
        plt.title('Cluster '+clust_index[indx])
        plt.ylabel('N. de Genes')
        plt.xlabel(xlabel) 
        
        plt.bar(keys,values2, label="Seleccion")
        plt.title('Cluster '+clust_index[indx])
        plt.ylabel('N. de Genes')
        plt.xlabel(xlabel)
        plt.legend()
        
        plt.savefig('Cluster_' +clust_index[indx] +'_'+ img_name)
        plt.close()
            
        indx+=1
            
        ## Se dejan en 0 todos los valores del diccionario
        for i in clust_dict1:
            clust_dict1[i] = 0
            clust_dict2[i] = 0


# In[23]:

import sys

#sys.argv[1]='cluster_chroms'
#sys.argv[2]='data_index'
#sys.argv[3]='Cromosomas'
#sys.argv[4]='chroms'

atribute_barplot(att_file=sys.argv[1], index_file=sys.argv[2] ,xlabel=sys.argv[3], img_name=sys.argv[4])


# In[148]:


##Objetivos
## Crear graficos del tipo
## Barplot de atributos en la tabla:
##  -Que tome en cuenta todos los tipos de atributos que aparecen en la busqueda
##  -Que el usuario puede decidir sobre que atributos del cuerpo hacer histogramas (A)
##  -Opcion para decidir si usar nombres cajas especificos para un atributo, se debe incluir un box donde entre
##   todo lo que no coincidio con la caja (X1)
##  -Induvidual para cada cluster, que incluya la seleccion y el total (B)
##  -Barplot que incluya, seleccionados + complementos (C)
## 
## Histograma del encabezado (data_cl, jceegee, total de genes)
##  -Seleccionados, totales en un archivo


# In[ ]:


#(A): Por defecto se trabajan sobre todos los tipos de atributos, si el
####  usuario lo desea puede pasar por la terminal una lista de atributos de interes
####  en caso que no de equivocacion debe aparecer un cuadro dialogo similar al que aparece
####  en feature_extractor
#(B): Hecho
#(C): Hecho
# Por hacer:
# -Agregar data_index
# -Opcion para quedarse solo con el grafico comparativo, el del cluster total o la seleccion
# -Crear archivos con los graficos producidos
# -Asociar la funcion a la linea de comando

#X1: Por el momento el usuario solo puede decidir si quiere que los box abarquen todos los tipos
#### de elementos encontrados en toda la busqueda, o si solo desean que abarquen los presentes 
#### en el cluster. 
#### En un futuro se puede agregar una opcion que incluya una lista de nombres de elementos y una lista
#### de regex asociados a esos elementos. En este caso se debe agregar un box 'Otros', que incluya todos 
#### los elementos no asociados al la lista de regexp


# In[ ]:


#### Histograma por hacer:
#### -Pasar los archivos de largos y cg a formato cluster
#### -Moverlos al directorio de trabajo
#### -Graficar a ambos superpuestos
#### -Producir archivos 
#### -Asociarlo a la linea de comando de la terminal

