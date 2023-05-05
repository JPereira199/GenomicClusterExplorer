#setwd("/home/jpereira/Escritorio/respaldo.Inflation_45/refeed_Toxo_Neo/R_chroms/")
setwd("/home/jpereira/Escritorio/vacation_work/final_work/Full_TgNc/R_work/R_chrom/")

#BiocManager::install("DECIPHER")
#BiocManager::install("GRanges")
#BiocManager::install("GenomicRanges")
#install.packages("XML")
#install.packages("openssl")
R_chrom.R

library("dplyr")
library("chromPlot")
library("GenomicRanges")

df_SRS <- read.csv(file = "df_coords", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
df_coords <- df_SRS
chr_lens <- read.csv(file = "chrom_longs", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

## Se 'Arreglan' los nombres de los cromosomas
df_coords$Chrom <- gsub("^I", "Tg I", df_coords$Chrom)
df_coords$Chrom <- gsub("^X", "Tg X", df_coords$Chrom)
df_coords$Chrom <- gsub("^V", "Tg V", df_coords$Chrom)
df_coords$Chrom <- gsub("_", " ", df_coords$Chrom)

## Se agregan los orgnames
### Se arreglan los Org Names
df_coords$Org_name <- as.character(df_coords$Org_name)
df_coords[grep( "Ncaninum", df_coords$genID ), "Org_name"] <- "Ncaninum"
df_coords[grep( "TG", df_coords$genID ), "Org_name"]  <- "Tgondii"
df_coords$Org_name <- as.factor(df_coords$Org_name)

## Se separan los datos de Ncaninum y Tgondii
df_coords_Nc <- df_coords[ grep("Ncaninum", df_coords$Org_name), ] 
df_coords_Tg <- df_coords[ grep("Tgondii", df_coords$Org_name), ] 
chr_Nc <- chr_lens[grep("Nc", chr_lens$V1),]
chr_Tg <- chr_lens[grep("Tg", chr_lens$V1),]

## ChromPlot para Ncaninum
{df_bands <-data.frame(Chrom = df_coords_Nc$Chrom, 
                      Start = df_coords_Nc$Start,
                      End = df_coords_Nc$End) 
df_bands$Colors <- "red"

GR <- GRanges( seqnames = chr_Nc$V1,  
               ranges = IRanges(rep(0, length(chr_Nc$V2)), end = chr_Nc$V2, names =  NULL ),
               score = rep(1, length(chr_Nc$V2)) )
#png(paste("Chrom_map_clust_", as.character(clust), ".png"),  width = 590, height = 625)}
chromPlot(gaps = GR,
          annot1 =  df_bands,
          figCols = 8,
          bin=1.0e5)
}

## Ploteo para Tgondii
{df_bands <-data.frame(Chrom = df_coords_Tg$Chrom, 
                      Start = df_coords_Tg$Start,
                      End = df_coords_Tg$End) 
df_bands$Colors <- "red"

GR <- GRanges( seqnames = chr_Tg$V1,  
               ranges = IRanges(rep(0, length(chr_Tg$V2)), end = chr_Tg$V2, names =  NULL ),
               score = rep(1, length(chr_Tg$V2)) )
#png(paste("Chrom_map_clust_", as.character(clust), ".png"),  width = 590, height = 625)
chromPlot(gaps = GR,
          annot1 =  df_bands,
          figCols = 8,
          bin=1e5)
}



###########################

df_merge <- df_coords
df_merge["product"][df_merge["product"] == "-----"] <- NA
df_merge$product <- as.character(df_merge$product)
df_merge$description <- as.character(df_merge$description)

df_merge <- within(df_merge, product <- ifelse(is.na(product), description, product))
df_merge$product <- as.factor(df_merge$product)
df_merge$description <- as.factor(df_merge$description)


df_SRS <- df_merge
summary(df_SRS$product)
df_SRS$product <- gsub("hypothetical protein.*", "hypothetical protein", df_SRS$product)
df_SRS$product <- gsub("putative SRS.*", "putative SRS", df_SRS$product)
df_SRS$product <- gsub(".*SAG-related sequence SRS.*", "SRS", df_SRS$product)
df_SRS$product <- gsub("srs domain-containing protein", "SRS", df_SRS$product)
df_SRS$product <- gsub("SRS domain-containing protein", "SRS", df_SRS$product)
df_SRS$product <- gsub("SAG1-related sequence 2, related", "putative SRS", df_SRS$product)
df_SRS$product <- gsub("SRS domain containing protein, putative", "putative SRS", df_SRS$product)
df_SRS$product <- gsub("SAG3 protein, related", "putative SRS", df_SRS$product)
df_SRS$product <- gsub("SAG2 related antigen SAG2D, related", "putative SRS", df_SRS$product)

df_SRS$product <- as.factor(df_SRS$product)
summary(df_SRS$product)

df_SRS[grep("hypothetical protein", df_SRS$product),]$clust_num

#### Hipoteticas 
## Tgondii
#Clust 2  321460 Imagen ilegible
#Clust 38 238530 E-value 41 y 50, tamaño de 104 pb ||  Yes
#Clust 43 259280 E-value muy bajo, todos los miembros menores a 9
#Clust 94 207015 E-value Excelente con vario miembros, aa 575 || yes
#Clust 374 271145 E-value 200 a un miembro de Ncaninum: Ncaninum_LIV_000034300.1, putative SRS, tamaño de 5645 aa, cg 61.7 || Maybe
#Clust 374 229300 E-value muy bajo
#Clust 374 217875 E-value muy bajo
#Clust 605 224160 E-value excelente con un miembro de Ncaninum (187600), este gen de Ncaninum tiene excelente E-value con un gen relacionado a SAG en Tgondii || Maybe
#TGME49_224170 y TGME49_224160 parecen formar a 187600
#Clust 1782 271040 E-value excelente con un miembro anotado en Ncaninum, casi mismo GC 52.9 y mismo len 204 || Yess Small
#Clust 2024 328700 E-value 50 con un miembro SRS, diferencia de tamaños (89 Hip vs 389 en Annt ), pertenece a un contig pequeño KE139801 || Small

df_SRS[df_SRS$clust_num == 605,]
df_SRS[df_SRS$clust_num == 1782,]
df_SRS[df_SRS$clust_num == 2024,]

## Ncaninum 
# Clust 621700 3942 E-value 79 con un miembro anotado de Tgondii, gene long 191 - 192, CG Hip 55.5 vs Annt51.4    || Small
# Clust 571600 4333 E-value 110 con un miembro anotado de Tgondii, gene long 244 - 294, CG Hip 61.1 vs Annt 60.0  || Small
# Clust 571600 4333 E-value 80 con un miembro anotado de Tgondii, gene long 214 - 211, CG Hip 53.4 vs Annt 51.8   || Small
df_SRS[df_SRS$clust_num == 3942,]
df_SRS[df_SRS$clust_num == 4333,]
df_SRS[df_SRS$clust_num == 5337,]

############################################################################################
## Chromosome Map Para cada cluster
clusters <- unique(df_coords$clust_num)
for (clust in clusters){
  df_clust <- df_coords[df_coords$clust_num == clust,]

df_bands <-data.frame(Chrom = df_clust$Chrom, 
                      Start = df_clust$Start,
                      End = df_clust$End) 
df_bands$Colors <- "red"
#chr_lens$V1 <- gsub("_", " ̲ ", chr_lens$V1 )
#df_bands$Chrom <-  gsub("_", " ̲ ", df_bands$Chrom ) 
head(df_bands)

# Si se usa chr[num] chromPlot solo utiliza [num] 
# chromPlot no acepta "_" en los nombres de los cromosomas 
GR <- GRanges( seqnames = chr_lens$V1,  
               ranges = IRanges(rep(0, length(chr_lens$V2)), end = chr_lens$V2, names =  NULL ),
               score = rep(1, length(chr_lens$V2)) )
GR
#png(paste("Chrom_map_clust_", as.character(clust), ".png"),  width = 590, height = 625)
chromPlot(gaps = GR,
          annot1 =  df_bands,
          figCols = 8,
          bin=1e5)
dev.off()
}

data(hg_cytoBandIdeo)
data(hg_gap)
chromPlot(bands=hg_cytoBandIdeo, gaps=hg_gap)

########################################################

#TGME49_239090, Clust5 110 (Cual es su Chrom) Tg_VI
#Ncaninum_LIV_000771900, Clust5 109 (Esta en el cluster con el resto?) No, todos estan en Nc_V y estte esta en Nc_XIII 
#TGME49_238530, Clust7 173 (P. Hip, esta en el clust?) Si, está en Tg_VI con los demas
#TGME49_207015, Clust10 202 (P. Hip, esta en el clust?) Si, está en Tg_X con los demas

df_coords[ df_coords$genID == "TGME49_239090-t26_1",]       # Clust5  Tg_VI   457703    459579
df_coords[ df_coords$genID == "Ncaninum_LIV_000771900.1",]  # Clust5  Nc_XIII 1010962   1011864 
df_coords[ df_coords$genID == "TGME49_238530-t26_1",]       # Clust7  Tg_VI   267452    268679
df_coords[ df_coords$genID == "TGME49_207015-t26_1",]       # Clust10 Tg_X    7338883   7341116
############################################################################
############################################################################

IRanges(rep(0, length(chr_lens$V1)), end = chr_lens$V2, names = chr_lens$V1 )
IRanges(rep(1, ), end = chr_lens$V2[1:4], names = NULL )

        
data("hg_gap")
hg_gap
chromPlot(gaps=hg_gap)
?chromPlot
chromPlot

gr <- GRanges( seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
               ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
               strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
               score = 1:10,
               GC = seq(1, 0, length=10))
gr

gr2 <- GRanges( seqnames = c("chr1", "chr2", "chr3", "chr4"),
               ranges = IRanges(c(0,0,0,0), end = c(1100000000,12000000,13000000,1400000), names = head(letters, 4)),
               score = 1:4,
               GC = seq(1, 0, length=4))

gr2
chromPlot(gaps=gr2)

gr3 <- GRanges( seqnames = c("chr1", "chr2", "chr3", "chr4"),
               ranges = IRanges(rep(0, length(chr_lens$V2[c(1,2,3,4)])), end = chr_lens$V2[c(1,2,3,4)], names = head(letters, 4)),
               )

gr3
chromPlot(gaps=gr3)
