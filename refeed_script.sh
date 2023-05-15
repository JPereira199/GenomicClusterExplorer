#!/bin/bash

#Script Refeed

#1
#Obtener un grupo genes de out_mci
#Crear un faa restrigido a ese grupo de genes

#La idea es obtenerlos de la carpeta deesde la salida de plotting_function
#Se usa el samtools

function refeed_start
{

  #$1: out_mci 
  #$2: faa
  #$3: Mfasta
  #$4: gff
  #$5: sufix_name 

  mkdir refeed_dir_"$5" 

  cp ./../$1 out_mci
  cp ./../$2 faa
  cp ./../$3 Mfasta
  cp ./../$4 gff

  mv  out_mci faa Mfasta gff ./refeed_dir_"$5"
  cd ./refeed_dir_"$5"

}


function refeed_mini
{

  tr "\t" "\n" < out_mci > names.temp
  sed -i '/^$/d ' names.temp

  ## Se crean mini faa y Mfasta con las secuencias nombradas en names.temp
  xargs samtools faidx faa < names.temp > mini_faa_"$1"
  xargs samtools faidx Mfasta < names.temp > mini_Mfasta_"$1"

  rm names.temp

}

## Agrega el atributo Org_name a los genes y cromosomas de gff
function refeed_orgname
{

  ## $1 suffix_name
  ## $2 org_name

  ## Se particiona el archivo out_mci. Con el objetivo de evitar el error 
  ## "La lista de argumentos es demasiado larga"

  mkdir split_mci
  cp out_mci ./split_mci
  cd ./split_mci
 
  tr "\t" "\n" < out_mci > temp		#Se agrupan todos losgenes en una sola columna
  sed -i '/^[[:space:]]*$/d' temp	#Se borran posible lineas vacias
  split -l 3000 temp			#Se subdivide temp en archivos de 3000 lineas
  rm out_mci temp
  
  for file in $( ls )
  do
 
    ## Se crea una expresión regular para agregar el atributo Org_name
    ## a todos los genes de una sola vez
    sed_regexp=$(tr "\n" "\t" < $file |\
          	  sed 's/^\t\|\t$//g' |\
        	  sed 's/^/\\(ID=/; s/$/;\\)/' |\
        	  sed 's/\t/;\\|ID=/g')

    sed -i "s/$sed_regexp/&Org_name=$1;/g" ./../gff 

  done
 
  cd ..  
  rm -rf split_mci
  
  ## Se agrega el atributo Org_name a todos los contigs/cromosomas
  ## encontrados al inicio del archivo

  sed -i "s/##sequence.*/& $1/" gff

}


function refeed_recursive
{

## Creado de carpeta
if [ -d refeed_dir ]
then
        rm -r refeed_dir
fi
mkdir refeed_dir
cd ./refeed_dir

## Convirtiendo los inputs en vectores
local -n mci_vect=($1) &> temp		# Vector de genes a agrupar
local -n faa_vect=($2) &> temp		# Vector de faa a agrupar 
local -n Mfasta_vect=($3) &> temp	# Vector de fasta a agraupar
local -n gff_vect=($4) &> temp		# Vector de gff a agrupar
local -n name_vect=($5) &> temp		# Vector de nombres
local -n org_vect=($6) &> temp		# Vector de nombre de organismos
rm temp


#Numero de elementos de cada vector
mci_var=${#mci_vect[@]}
faa_var=${#faa_vect[@]}
Mfasta_var=${#Mfasta_vect[@]}
gff_var=${#gff_vect[@]}
name_var=${#name_vect[@]}
org_var=${#org_vect[@]}


## Checkeando que todos los vectores tengan el mismo tamaño
if (( $mci_var != $faa_var || $mci_var != $gff_var || $mci_var != $name_var || $mci_var != $Mfasta_var || $mci_var != $org_var ))
then
	echo "ERROR: Las entradas deben tener el mismo número de elementos"
	echo "Clusters:		$mci_var"
	echo "Fastas de AA:	$faa_var"
	echo "Multifastas:	$Mfasta_var"
	echo "Archivos GFF:	$gff_var"
        echo "Sufijos:		$name_var"  	
	echo "Organismos: 	$org_var"
	
	exit
fi


## Aplicando los programas refeed_start, refeed_mini_faa y refeed_orgname a cada elemento de los vectores
touch all_faa
touch all_Mfasta
touch all_mci
touch all_gff
for i in ${!mci_vect[@]} #La variable 'i' toma valores 0 1 2 . . 
do

	#echo "Index de mci_vect: $i"
	echo "Procesando genoma Nº$i "

	mci_var=${mci_vect[$i]}
	faa_var=${faa_vect[$i]}
	Mfasta_var=${Mfasta_vect[$i]}
	gff_var=${gff_vect[$i]}
	name_var=${name_vect[$i]}
	org_var=${org_vect[$i]}
	
	refeed_start $mci_var $faa_var $Mfasta_var $gff_var $name_var
	refeed_mini $name_var
	refeed_orgname $org_var

	cd ..
	
	cat all_faa ./refeed_dir_"$name_var"/mini_faa_"$name_var" > A
	cat all_Mfasta ./refeed_dir_"$name_var"/mini_Mfasta_"$name_var" > B  
	cat all_mci ./refeed_dir_"$name_var"/out_mci > C  
	cat all_gff ./refeed_dir_"$name_var"/gff > D
	
	mv A all_faa
       	mv B all_Mfasta
       	mv C all_mci
	mv D all_gff
	
done

}

##Falta revisar que esten todos los archivos de entrada

#refeed_recursive "Ncaninum.mci Tgondii.mci" "Ncaninum_LIV.faa TgondiiME49_Genome.faa" "Ncaninum_LIV.Mfasta TgondiiME49.Mfasta" "Ncaninum_LIV.gff3 TgondiiME49.gff" "NcSRS TgSRS" "Ncaninum Tgondii"

#refeed_recursive "./single_mci/out_mci_whole_NcSRS ./single_mci/out_mci_whole_TgSRS" "Ncaninum_LIV.faa TgondiiME49_Genome.faa" "Ncaninum_LIV.Mfasta TgondiiME49.Mfasta" "Ncaninum_LIV.gff3 TgondiiME49.gff" "NcSRS TgSRS" "Ncaninum Tgondii"
#refeed_recursive "./single_mci/out_mci_whole_NcMIC ./single_mci/out_mci_whole_TgMIC" "Ncaninum_LIV.faa TgondiiME49_Genome.faa" "Ncaninum_LIV.Mfasta TgondiiME49.Mfasta" "Ncaninum_LIV.gff3 TgondiiME49.gff" "NcMIC TgMIC" "Ncaninum Tgondii"

#refeed_recursive "./single_mci/out_mci_whole_NcGRA ./single_mci/out_mci_whole_TgGRA" "Ncaninum_LIV.faa TgondiiME49_Genome.faa" "Ncaninum_LIV.Mfasta TgondiiME49.Mfasta" "Ncaninum_LIV.gff3 TgondiiME49.gff" "NcGRA TgGRA" "Ncaninum Tgondii"

#refeed_recursive "out_mci_whole_Neo out_mci_whole_Toxo" "Ncaninum_LIV.faa TgondiiME49.faa" "Ncaninum_LIV.Mfasta TgondiiME49.Mfasta" "Ncaninum_LIV.gff3 TgondiiME49.gff" "Neo Toxo" "Ncaninum Tgondii"



#2 
#Repetir uno para otros grupos de genes de interes
#Concatenar los faa entre si y sus respectivos gff entre si

#3
#Realizar correr el blastp_mcl_protocolo desde el blastp
