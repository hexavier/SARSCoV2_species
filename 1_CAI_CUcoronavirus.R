library(tAI)
library(ggplot2)
library(ggpubr)

extract_cod <- function (trnas, anticod){
  output = data.frame(row.names = anticod)
  trnas_acod = sapply(rownames(trnas), function(x) substr(x,nchar(x)-2,nchar(x)))
  for (s in colnames(trnas)){
    output[,s] = sapply(anticod, function(x) if(any(trnas_acod==x)){mean(trnas[trnas_acod==x,s])}else{0})
  }
  return(output)
}

AAnormalize <- function(data,codons){
  # Compute relative data
  aa = sapply(codons,function(x) substr(x,1,nchar(x)-3))
  uniqueaa = unique(aa)
  outdata = numeric(length=length(data))
  for (n in uniqueaa){
    idx = (aa %in% n)
    total = max(data[idx],na.rm=T)
    outdata[idx] = data[idx]/total
    if (total %in% 0){
      outdata[idx] = 1.0/sum(idx)
    }
  }
  return(outdata)
}

## Load trna and weighted CU
# Coronavirus strains
sars2 = "NC_045512.2"

# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# Codon frequencies
codonfreq = read.csv("data/species_codontable.csv",row.names = 1)
rownames(codonfreq) = paste0(codons[rownames(codonfreq),"AA"],rownames(codonfreq))
species = c("Homo_sapiens","Gallus_gallus","Anas_platyrhynchos","Mus_musculus","Rattus_norvegicus","Sus_scrofa","Canis_lupus_familiaris","Mustela_putorius_furo","Felis_catus")
# Normalize by the maximim codon of each amino acid family
relcodus = apply(codonfreq,2,AAnormalize,rownames(codonfreq)); rownames(relcodus) = rownames(codonfreq)
CUweights = extract_cod(relcodus,rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

# Genomic codon usage
codus = read.csv("data/CU_coronaviruses.tsv",sep="\t", row.names = 1)
codus = codus[codus$Accession %in% sars2,]
# Keep only columns with codon info
codus_clean = data.frame(sapply(rownames(codus),function(x) as.numeric(codus[x,11:ncol(codus)])), row.names = colnames(codus)[11:ncol(codus)])
rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))

## Calculate tAI for genomic CU
codon = extract_cod(codus_clean,rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

CAI = data.frame(matrix(ncol = ncol(CUweights), nrow = ncol(codon)),row.names = colnames(codon)); colnames(CAI) = colnames(CUweights)
# Calculate tAI
for (sample in colnames(CUweights)){
  # Calculate relative adaptiveness values (ws)
  sample.ws = as.numeric(CUweights[,sample])
  # Calculate tAI for all CUs
  sample.cai <- get.tai(t(codon), sample.ws)
  CAI[,sample] = sample.cai
}
CAI[,c("Species","description")] = codus[,c("Species","description")]

## Plot
# Create dataset structure
dataset = c()
for (l in species){
  dataset_temp = data.frame(row.names = rownames(CAI))
  dataset_temp$CAI = as.numeric(CAI[,l])
  dataset_temp$virus = as.character(CAI$Species)
  dataset_temp$species = l
  dataset_temp$protein = as.character(CAI$description)
  dataset = rbind(dataset,dataset_temp)
}
dataset$species <- factor(dataset$species,levels = species)


ggplot(dataset, aes(x=species, y=CAI, fill=virus)) +  
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2) +
  geom_jitter(na.rm=T , size=2.5, alpha = 0.75 ,aes(color = protein), width = 0.15) +
  scale_color_manual(values=c("#29338B","#496CB4","#54bfc4","#d049a3","#707272","#CCCCCC",
                              "#803082","#e6e600","#53c653","#2d862d","#D49F9F","#E91F24")) +
  scale_fill_manual(values=c("white")) +
  theme_classic() +
  theme(axis.text.x=element_text(size = rel(1), angle=45, hjust=1, vjust =  1)) +
  labs(title="Coronaviruses",y = "CAI")

speciessig = compare_means(CAI ~ species, method = "wilcox.test",paired =T, data=dataset)
