
# 3 Sequencing Technologies
# 4 Data Management, Reproducibility & Integrity
# 5 Effective Population Size

# 6 Management Of Small Populations
# 7 Natural Selection
# 8 Landscape Genetics

# 9 Lineage Divergence

# 11 From Genes to Kin: Dissecting Relatedness & Kinship
# 12 Genetic Structure
# 13 Combining Genomic Resources



dir.create(
  './sessions/'
)

for (i in 11:13) {
  file.copy(from = './session04.qmd',
            to = './sessions/session04.qmd')
            
  file.rename(from = './sessions/session04.qmd',
              to = paste0('./sessions/session', i,'.qmd'))
          
}

