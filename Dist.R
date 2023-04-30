## Medir la distancia entre genes de un mismo cromosoma
## usando un dataframe como entrada
## requiere el paquete biogenerics
getwd()
#setwd("/home/jpereira/Escritorio/vacation_work/final_work/Full_TgNc/R_work/R_chrom/")
setwd("/home/juaco/Desktop/vacation_work/final_work/Full_TgNc/R_work/R_chrom/")
df_coords <- read.csv( file = "df_coords", header = TRUE, sep = "\t", dec = ".", stringsAsFactors = TRUE )

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
#TGME49_224790 y TGME49_224780 son genes con una superposición de -266 nc
df_order[299,]

#dist <- c( rep( 1000, 4), rep(9999, 4), rep(1000,4), 9999 )

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

df_order[df_order$clust_num == 7, c("clust_num","Tandem", "Dist")]


install.packages("dplyr") 
library("dplyr")

#Lista con todos los tandems ordenados
tandem_list <- unique(df_order[order(df_order$Tandem),c("Tandem")])

tand_per_clust <- c()
for ( i  in tandem_list) {

  tand_per_clust_i <- df_order[ df_order$Tandem == i,  c("clust_num")] %>%
    unique() %>% length() 
  tand_per_clust <- append(tand_per_clust, tand_per_clust_i)
}



order(tand_per_clust)
df_order[ df_order$Tandem == 7,  c("clust_num")] %>% unique() %>% length()
