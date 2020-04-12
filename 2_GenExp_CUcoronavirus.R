library(ggplot2)

# Load data
rnaseq = read.csv("data/species_ACE2expression.csv",row.names = 1)
species = c("Homo_sapiens","Gallus_gallus","Anas_platyrhynchos","Mus_musculus","Rattus_norvegicus","Sus_scrofa","Canis_lupus_familiaris","Mustela_putorius_furo","Felis_catus")

## Plot
# Create dataset structure
dataset = data.frame(row.names = species)
dataset$ACE2 = sapply(species,function(x) median(rnaseq[rnaseq$species %in% x,"norm_expression"]))
dataset$species = factor(species,levels = species)

ggplot(dataset, aes(x=species, y=ACE2)) +  
  geom_bar(stat="identity",color="black",fill="grey") +
  theme_classic() +
  theme(axis.text.x=element_text(size = rel(1), angle=45, hjust=1, vjust =  1)) +
  labs(title="Coronaviruses",y = "ACE2")
