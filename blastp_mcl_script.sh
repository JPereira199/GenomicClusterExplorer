#!/bin/bash


function mcl_blastp_function
{

# $1: Archivo .fasta 
# $2: Archivo con anotaciones correspondientes .gff3 

#Se crean los directorios de trabajo
mkdir mcl_blastp_directorio

if [ -z  $re_feed  ]
then
	echo "Iniciando rutina por defecto"

	cd ./mcl_blastp_directorio 
	cp ./../$1 fasta 
	cp ./../$2 gff

	#Crea un string que cambia la terminacion de $1 de fasta a faa
	#### Evaluar la posibilidad de que el usuario decida el nombre de salida de los archivos, de forma de evitar errorres
	FAA=$(echo $1 | sed 's!.*\/!!g' | sed 's/fasta$/faa/' ) 
	MFASTA=$(echo $1 | sed 's!.*\/!!g'  | sed 's/fasta$/Mfasta/' )

	#:Creacion de un multifasta traducido (-y) y sin traducir (-x) a partir de del genome.fasta y el .gff

	echo "Creacion de $FAA y $MFASTA"
	gffread -y $FAA -g fasta gff 
	gffread -x $MFASTA -g fasta gff

	echo "Realizando blastp de  $FAA sobre si mismo"
	blastp -query $FAA -subject $FAA -outfmt "6 qseqid sseqid evalue bitscore" > blastp_tabla
else
	
	echo "Iniciando rutina re-feed"
	echo "Ejecutando BLASTP"

	cd ./mcl_blastp_directorio
	cp ./../$re_feed faa
       	#cp ./../$1 gff  #Entrada ineecesaria con esta opción
	
	blastp -query faa -subject faa -outfmt "6 qseqid sseqid evalue bitscore" > blastp_tabla

fi

cut -f 1,2,3 blastp_tabla > blastp_evalue.abc
cut -f 1,2,4 blastp_tabla > blastp_bitscore.abc ## Sirve para en el futuro , una vez creado los cluster por mcl, observar el bitscore de cada cluster mediante un heatmap

#Agregar opcion para controlar el 'ceil' de mcxload
echo "Ejecutando MCL"
mcxload -abc blastp_evalue.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o blastp_evalue.mci -write-tab blastp_evalue.tab

mcl blastp_evalue.mci -I 4.5 -use-tab blastp_evalue.tab

}


function mcl_blastp_function2
{
if [ -d ./mcl_blastp_directorio ]
then
 echo -n "mcl_blastp_directorio ya existe, desea removerlo? (y|n): "
 read input
 case $input in
         y|yes|YES) rm -r mcl_blastp_directorio
                    mcl_blastp_function $1 $2;;
         *)echo "total_cl se ha detenido";;
esac
else
 mcl_blastp_function $1 $2
fi
}

#mcl_blastp_function2 TgondiiME49_Genome.fasta TgondiiME49.gff 


######################################################################


#function master_input2 { ### Nombre antiguo de la función
function clustering_blastp_mcl {

function usage()
{

echo "Usage: blastp_mcl_rutine 	[-R | --Re-Feed FILE]
                           	[-g | --gff FILE]
                           	[-f | --fasta FILE]"
exit 2

}
## -a: Modo alternativo. Permite opciones largas con un solo guion (-)
## -n: Indica a getopt el nombre del programa en ejecucion. Usado al momento de  devolver errores
PARSED_ARGUMENTS=$( getopt -a -n blastp_mcl_rutine -o R:f:g: --long Re-Fee: -- "$@" )
VALID_ARGUMENTS=$?

if [ "$VALID_ARGUMENTS" != "0"  ]; then
        usage
fi

echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"

while :
do
        case "$1" in
                -R | --Re-Feed) re_feed="$2" ; shift 2;;
                -f | --fasta) fasta="$2" ; shift 2;;
                -g | --gff) gff="$2"; shift 2;;

               # -P | --plot) plot="TRUE"; shift 1;;
                --) shift 1; break ;;
                *) echo "Unexpected option: $1"; usage
        esac
done

echo "Final arguments: $@"

################ Ejecutando el resto del Script ###############

mcl_blastp_function2 "$@"

}

##NEW NAME: clustering_blastp_mcl 

#master_input2 -R Full_TgNc.faa 
#master_input2 -R ./refeed_inputs/refeed_dir_MIC/all_faa refeed_inputs/refeed_dir_MIC/all_gff ; cd .. ; mv mcl_blastp_directorio mcl_blastp_MIC
#master_input2 -R ./refeed_inputs/refeed_dir_GRA/all_faa refeed_inputs/refeed_dir_GRA/all_gff ; cd ..  ; mv mcl_blastp_directorio mcl_blastp_GRA
#master_input2 -R ./refeed_inputs/refeed_dir_SRS/all_faa refeed_inputs/refeed_dir_SRS/all_gff ; cd ..  ; mv mcl_blastp_directorio mcl_blastp_SRS
#master_input2 -R all_faa all_gff 
#master_input2 -plot Ncaninum_LIV.gff3 out.blastp_evalue.mci.I45 "product ortholog_cluster" Ncaninum_LIV.faa Ncaninum_LIV.Mfasta "SRS|srs|SAG"




