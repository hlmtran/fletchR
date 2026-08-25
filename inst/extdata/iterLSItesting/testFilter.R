roxygen2::roxygenise()
c = filterAndGetMat(SCEsmall,useMatrix = "TileMatrix")
d = filterAndGetMat(SCEsmall,useMatrix = "TileMatrix",subsetLSI = T)
identical(a,c)
identical(b,d)

e = filterAndGetMat(SCEsmall,useMatrix = "TileMatrix",replaceZeros = F)
f = filterAndGetMat(SCEsmall,useMatrix = "TileMatrix",subsetLSI = T,replaceZeros = T)

identical(c,e)
identical(f[rownames(d),],d)


z = f[rownames(d),]
