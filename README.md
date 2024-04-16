# cryptosporidium_parvum_evolution
Repository containing scripts related to the research paper 'Comparative Genomics of Cryptosporidium parvum: Uncovering the Emergence of an Outbreak-Associated Population in Europe and Its Transmission to the USA.' https://www.biorxiv.org/content/10.1101/2023.09.19.558430v1


## Assess multiplicity of infection
R script for assessing the multiplicity of infection across the samples
```{R}
library(SeqArray)
library(moimix)

# Load VCF and covert to GDS
seqVCF2GDS ("input_moimix.vcf" , "input_moimix.gds")
# Import newly created GDS
isolates <- seqOpen("input_moimix.gds")
# Check name of the isolates
seqSummary(isolates)
sample.id <- seqGetData(isolates, "sample.id")
coords <- getCoordinates(isolates)
head(coords)

# Compute B allele frequency matrix
isolate_baf <- bafMatrix(isolates)

# Obtain a plot for a single sample
plot(isolate_baf, "CN11707") 
class(isolate_baf)
str(isolate_baf)

par(mfrow=c(2,2))
pdf(name_plot)
for (i in sample.id){
  
  plot(isolate_baf, i, main = i)
}
dev.off()


# Compute Fws values for each sample in VCF file
fws_all <- data.frame(getFws(isolates))
write.csv(fws_all, "moimixout.tsv", row.names = TRUE, quote = FALSE)
```


## SNP density across the chromosomes
Compute SNP density per kilobase
```{bash}
vcftools --vcf input.vcf --SNPdensity 1000 --out SNPs_densityPerKb
```
R script to obtain Supplemental Figure S3
```{R}
library(ggplot2)
library(dplyr)
library(plyr)

# Remove scientific notation for the final image
options(scipen=999)

# Load the input file
d <- read.table("SNPs_densityPerKb", stringsAsFactors=T, header=T) 
### change SNPs_densityPerKb with the name of the input file 

# Replace with chromosome name (in wanted, can also be skipped)
d$CHROM <- revalue(d$CHROM, c("CP044422.1"="1"))
d$CHROM <- revalue(d$CHROM, c("CP044421.1"="2"))
d$CHROM <- revalue(d$CHROM, c("CP044420.1"="3"))
d$CHROM <- revalue(d$CHROM, c("CP044419.1"="4"))
d$CHROM <- revalue(d$CHROM, c("CP044418.1"="5"))
d$CHROM <- revalue(d$CHROM, c("CP044417.1"="6"))
d$CHROM <- revalue(d$CHROM, c("CP044416.1"="7"))
d$CHROM <- revalue(d$CHROM, c("CP044415.1"="8"))


# order the result
d["SNP_COUNT"][d["SNP_COUNT"] == 0] <- NA
d <- d %>% arrange(desc(CHROM))


# Plot using ggplot
p <- ggplot(data=d, aes(x=BIN_START, y=-1)) +
  facet_grid(CHROM ~ ., switch='y') +
  geom_tile(aes(fill=SNP_COUNT)) +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle=180),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#EDECEA"),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradientn(colours=c("white", "#FFDB43","#fcac08ff","#fc7b05ff", "#ff2603ff"), na.value="white")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), labels=scales::unit_format(unit="kb", scale=1e-3)) +
  labs(title="Number of SNPs within 10kb window size" , fill="")
```
## Heatmap SNP distances
R script used to obtain Figure 2 
```{R}
library(reshape2)
library(ggplot2)
library(dplyr)
library(ape)  
library(dendextend)
library(gplots)
library(autoimage)

heatmap_data <- read.csv("matrix.tsv", sep='\t', check.names=FALSE, stringsAsFactors = FALSE)
h_d <- as.matrix(heatmap_data[,-1])
row.names(heatmap_data) <- heatmap_data[,1] 
heatmap_data[,1] <- NULL
heatmap_data[heatmap_data==0]=NA
h_d <- as.matrix(heatmap_data)

### define the colors within 4 zones
breaks = seq(0, max(h_d), length.out=1000)
gradient1 = colorpanel( sum( breaks[-1]>0 & breaks[-1]<=50 ), "red", "orange" )
gradient2 = colorpanel( sum( breaks[-1]>50 & breaks[-1]<=1000 ), "orange", "green" )
gradient3 = colorpanel( sum( breaks[-1]>1000 & breaks[-1]<=8000 ), "green", "violet" )
hm.colors = c(gradient1, gradient2, gradient3)
graphics.off()


#with included legend
a <- heatmap.2(h_d, trace="none", col=hm.colors, breaks = breaks, na.color = "red",
               key=TRUE, key.title = "", key.xlab = "", symkey=FALSE, symbreaks=FALSE, keysize = 10, key.par=list(mar=c(3,1,2,2)),#( "bottom.margin", "left.margin", "top.margin", "right.margin" )
               sepcolor = "white", sepwidth = c(0.005, 0.005), colsep = 1:ncol(h_d), rowsep=1:nrow(h_d), cexCol = 0.001, cexRow = 0.001, lhei=c(0.75,5), lwid=c(0.75,3.5))
```


## Recombinations
To convert the VCF file into a human readble excel file
```{python}
python3 vcf2xlsx.py used_sites.tsv fasta_file.fa chrom_name
````
used_sites.tsv and fasta_file.fa are the output of https://github.com/edgardomortiz/vcf2phylip 

Each SNP has been annotated with the corresponding gene 
```{python}
python3 add_annotation_to_SNP.py GFF_file_reference Chr_file.xlsx
````
Chr_file.xlsx is the result of the previous python script


