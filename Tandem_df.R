## A partir de las posiciones de cada gen, se indican si estos
## forman un tandem. También crea una tabla que indica cuantos
## tandems hay por cromosoma
## Usa como entrada los DF generados previamente en lafunciíon 
## gene_datafreme

## Requiere el paquete "dplyr"
if (!require(dplyr)) {
  cat("El paquete dplyr no está instalado \n")
} 

## ./../Datos/searcher_dir/df_dir/
Arg <- commandArgs(trailingOnly = TRUE)
Path <- Arg[1]
TanDist <- as.numeric(Arg[2])
Input <- Arg[3]
Output <- Arg[4]

cat("Directorio:", Path, "\n")
cat("Distancia de corte de Tandem:", Arg[2], " bases\n")
cat("Entrada:", Arg[3], "\n")
cat("Salida:", Arg[4], "\n")

setwd(as.character(Path))
df_coords <- read.csv( file = Input, header = TRUE, sep = "\t", dec = ".", stringsAsFactors = TRUE )

#Se arreglan los nombres 
#Quitar en la versiòn final
df_coords[grep("Ncaninum", df_coords$genID), "Org_name"] <- "Ncaninum"
df_coords[grep("TG", df_coords$genID), "Org_name"] <- "Tgondii"


## Se mide la distancia entre genes
{df_order <- df_coords[order(df_coords$Chrom, df_coords$Start),]
  dist <- c()
  chroms <- as.character(df_order$Chrom)
  for (i in 1:(length(df_coords$Chrom) - 1) ) {
    if (chroms[i] == chroms[i+1]) {
      dist <-append(dist, df_order$Start[i+1] - df_order$End[i]   )
    } else {
      dist <- append(dist, 99999999)
    }
  }; dist <- append(dist, 99999999)
  
  df_order$Dist <- dist}

{i=1
  t=0
  tand_size = 8000
  tandem <- c()
  while ( i < ( length(dist) ) ) {
    
    if (dist[i] < tand_size) {
      t = t + 1
      tandem <- append(tandem, t)
      i = i + 1
      
      while (dist[i-1] < tand_size) {
        tandem <- append(tandem, t)
        i = i + 1
      }
    }
    
    
    if ( dist[i] > tand_size) {
      tandem <- append(tandem, 0)
      i = i + 1
    }
    
  }
}
df_order$Tandem <- tandem

# Se indica a que tandem pertenece cada gen. Los genes con número de tandem
# igual a cero refieren a genes que no estan incluidos en ningùn tandem
tandem_list <- unique(df_order[order(df_order$Tandem),c("Tandem")])
clust_per_tandem <- c()
for ( i  in tandem_list) {
  
  clust_per_tandem_i <- df_order[ df_order$Tandem == i,  c("clust_num")] %>%
    unique() %>% length() 
  clust_per_tandem <- append(clust_per_tandem, clust_per_tandem_i)
}

tand_df <- data.frame(
  Tandem_num = tandem_list,
  clustXtandem = clust_per_tandem
)
write.table(x = df_order, file = Output, sep = "\t", col.names = T, quote = F, row.names = F)
df_tandems <- read.csv2(file = Output, header = T, sep = "\t", dec = ".", stringsAsFactors = T)


## Obtiene distintos datos por cromosoma de la familia multigénica buscada
## Números de genes por cromosoma, número de tandems, promedio de largo
## aminoácidico y GC de esos genes 
ChromData <- function(df) {
  
  numXChrom <- c()
  TandXChrom <- c()
  cgXChrom <- c()
  longXChrom <- c()
  
  for ( chrom in unique(df$Chrom) ) {
    
    numXChrom <- append(numXChrom, length(df[ df$Chrom == chrom, c("genID") ]))
    TandXChrom <- append(TandXChrom, length(unique(df[ df$Chrom == chrom & df$Tandem != 0, c("Tandem") ])))
    cgXChrom <- append(cgXChrom, median(as.numeric(df[ df$Chrom == chrom, c("genes_cg") ])))
    longXChrom <- append(longXChrom, median(as.numeric(df[ df$Chrom == chrom, c("genes_longs") ])))
    
  }
  
  new_df <- data.frame( Chrom = unique(df$Chrom),
                        numXChrom = numXChrom,
                        TandXChrom = TandXChrom,
                        cgXChrom = cgXChrom,
                        longXChrom  = longXChrom  )
  
  return(new_df)
  
}

## Función sin uso. Al igual que ChromData, Tandem_count obtiene distintos
## datos por cromosoma de la familia multigénica buscada. Además separa
## esos datos según los genes pertenescan a un tandem o no
## P. ej: num0XChrom y num1XChrom indican respectivamente el número de genes
## aislados y en tandem para la familia multigénica buscada 
Tandem_count <- function(df)
{
  num0XChrom <- c()
  num1XChrom <- c()
  TandXChrom <- c()
  cg0XChrom <- c()
  cg1XChrom <- c()
  long0XChrom <- c()
  long1XChrom <- c()
  
  for ( chrom in unique(df$Chrom) ) {
    
    num0XChrom <- append(num0XChrom, length(df[ df$Chrom == chrom & df$Tandem == 0, c("genID") ]))
    num1XChrom <- append(num1XChrom, length(df[ df$Chrom == chrom & df$Tandem != 0, c("genID") ]))
    TandXChrom <- append(TandXChrom, length(unique(df[ df$Chrom == chrom & df$Tandem != 0, c("Tandem") ])))
    cg0XChrom <- append(cg0XChrom, median(as.numeric(df[ df$Chrom == chrom & df$Tandem == 0, c("genes_cg") ])))
    cg1XChrom <- append(cg1XChrom, median(as.numeric(df[ df$Chrom == chrom & df$Tandem != 0, c("genes_cg") ])))
    long0XChrom <- append(long0XChrom, median(as.numeric(df[ df$Chrom == chrom & df$Tandem == 0, c("genes_longs") ])))
    long1XChrom <- append(long1XChrom, median(as.numeric(df[ df$Chrom == chrom & df$Tandem != 0, c("genes_longs") ])))
    
  }
  
  new_df <- data.frame( Chrom = unique(df$Chrom), 
                        num0XChrom = num0XChrom,
                        num1XChrom = num1XChrom,
                        TandXChrom = TandXChrom,
                        cg0XChrom = cg0XChrom,
                        cg1XChrom = cg1XChrom,
                        long0XChrom = long0XChrom,
                        long1XChrom  = long1XChrom  )
  
  return(new_df)
}


Chrom_tandem_df <- ChromData(df_tandems)
write.table(x = Chrom_tandem_df, file = paste(Output,"chrom", sep="_"), sep = "\t", col.names = T, quote = F, row.names = F)

#Tandem_count(df_tandems)

#Tandem_count(df_tand_Tg)
#Tandem_count(df_tand_Nc)
#ChromData(df_tand_Tg)
#ChromData(df_tand_Nc)

