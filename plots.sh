#!/bin/bash

function search_graphs
{

echo "Plotting Histogrmas and Barplots..."
#pwd

if [[ -d ./search_graphs ]]
then
        rm -r search_graphs
        mkdir search_graphs
else
        mkdir search_graphs

fi

cd ./search_graphs
## Asumiendo que estoy en scripts/
cp ./../search_bin/search_n_extract_graphs.py .
cp ./../../Datos/searcher_dir/data_index .
cp ./../../Datos/searcher_dir/crude_results/total_cluster_directorio_whole/cluster_chroms .
mv cluster_chroms cluster_chroms_whole
cp ./../../Datos/searcher_dir/crude_results/total_cluster_directorio_selection/cluster_chroms .
mv cluster_chroms cluster_chroms_selection

python3 search_n_extract_graphs.py "cluster_chroms" "data_index" "Cromosomas" "chroms"
rm cluster_chroms_selection cluster_chroms_whole

############################## HISTOGRAMA DE LARGOS ###############################
echo "Histograma de Largos"
#pwd

cp ./../../Datos/searcher_dir/crude_results/data_directory_whole/largos_directorio/largo_* .
paste -s largo_* > clust_largos_whole

cp ./../../Datos/searcher_dir/crude_results/data_directory_selection/largos_directorio/largo_* .
paste -s largo_* > clust_largos_selection

cp ./../search_bin/search_n_extract_hist.py .

python3 search_n_extract_hist.py 'clust_largos_whole' 'clust_largos_selection' 'data_index' 'Largo en AA' 'largos' 1
rm clust_largos_whole clust_largos_selection largo_*

############################### HISTOGRAMA DE CG ####################################
echo "Histograma de GC"

#Estoy en searcher_graphs
cp ./../../Datos/searcher_dir/crude_results/jceegee_directorio_whole/clusters_cg/cluster_cg_* .
sed -i 's/ .*//' cluster_cg_*
sed -i 's/,/\./' cluster_cg_*
paste -s cluster_cg_* > clust_cg_whole

cp ./../../Datos/searcher_dir/crude_results/jceegee_directorio_selection/clusters_cg/cluster_cg_* .
sed -i 's/ .*//' cluster_cg_*
sed -i 's/,/\./' cluster_cg_*
paste -s cluster_cg_* > clust_cg_selection

python3 search_n_extract_hist.py 'clust_cg_whole' 'clust_cg_selection' 'data_index' 'Fraccion de GC' 'gc' 5
rm clust_cg_whole clust_cg_selection cluster_cg_*

################################## HISTOGRAMA DE CG TOTAL ############################

echo "Histograma de GC total"
#pwd

cp ./../../Datos/searcher_dir/crude_results/jceegee_directorio_whole/clusters_cg/cluster_cg_* .
sed -i 's/ *$//' cluster_cg_*
sed -i 's/.* //' cluster_cg_*
paste -s cluster_cg_* > clust_total_cg_whole

cp ./../../Datos/searcher_dir/crude_results/jceegee_directorio_selection/clusters_cg/cluster_cg_* .
sed -i 's/ *$//' cluster_cg_*
sed -i 's/.* //' cluster_cg_*
paste -s cluster_cg_* > clust_total_cg_selection

python3 search_n_extract_hist.py 'clust_total_cg_selection' 'clust_total_cg_selection' 'data_index' 'Total de GC' 'total_gc' 1
rm clust_total_cg_whole clust_total_cg_selection cluster_cg_*

################################## HeatMaps and Dendrograms ###########################

echo "Plotting Heatmaps and Dendrograms..."
#pwd

cp ./../search_bin/abc_filter ./../search_bin/heatmap_dendrogram.py .
cp ./../../Datos/searcher_dir/data_index ./../../Datos/searcher_dir/out_mci_whole .
## Usar out_mci_whole,  y señalar seleccionados #Anteriormente se uasba data_out_mci
cp ./../../Datos/blastp_evalue.abc  .

source ./../search_bin/abc_filter
abc_filter_func out_mci_whole

tr "\\n" "\\t" < out_mci_whole > order_out_mci
sed -i 's/\\t$//' order_out_mci

awk '{print NF}' out_mci_whole > gene_number
tr "\\n" "\\t" < gene_number > A
mv A gene_number
python3 heatmap_dendrogram.py

################################### ORDENANDO FIGURAS GENERADAS ########################

#pwd
for clust in $(cat data_index)
do
        mkdir cluster_$clust
        mv Cluster_"$clust"_* cluster_$clust
done

## Se remueven otros archivos usados para graficar
#rm data_index search_n_extract_graphs.py search_n_extract_hist.py

echo "Ploteo de histogramas y mapas de calor finalizado"

## Agregar opcion para modificar el numero de bins
## Despues hacer un barplot generico para los elementos del cuerpo (opcional)
## Opcion para poder usar o no los graficos
## Averiguar sobre los heat maps

cd ..
}

##### Integrar el resto de los scripts producidos en R



