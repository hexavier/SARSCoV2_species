
# Find trna data
species = dir("data/stAI_stadium_codon_weights")
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)
# Create output
output = data.frame(row.names = rownames(codons))
# Analyze copy numbers per species
for (s in species){
  weights = read.csv(sprintf("data/stAI_stadium_codon_weights/%s",s),row.names = 2)
  abbr = strsplit(s,"\\.")[[1]][1]
  output[,abbr] = 0
  output[rownames(weights),abbr] = as.numeric(weights$weight)
}

write.csv(output,"data/stAI_weights.csv")
