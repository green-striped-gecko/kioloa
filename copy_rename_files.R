dir.create(
  './sessions/'
)

fromfile <- list.files(pattern = '_code')
tofile <- sub('_code', '', fromfile)

for (i in 1:length(fromfile)) {
  file.copy(from = paste0('./',fromfile[i]),
            to = paste0('./',tofile[i]), overwrite = T)
  
  # file.rename(from = './sessions/session04.qmd',
  #             to = paste0('./sessions/session', i,'.qmd'))
  
}


