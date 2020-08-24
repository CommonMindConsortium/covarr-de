# This just an implementation detail: x and rotation are complimentary.  Since we ran prcomp with genes as rows, the rotation matrix stores the relevant PCs:
Y = matrix(rnorm(10000), ncol=50)
rownames(Y) = paste0('gene_', 1:nrow(Y))
colnames(Y) = paste0('ID_', 1:ncol(Y))

dcmp = prcomp( Y )

# Components for each Gene
dcmp$x[1:4, 1:4]

# Components for each Individual
dcmp$rotation[1:4, 1:4]


# PCA on the transpose
dcmp2 = prcomp( t(Y) )

# Components for each individual
dcmp2$x[1:4, 1:4]

# Components for each gene
dcmp2$rotation[1:4, 1:4]