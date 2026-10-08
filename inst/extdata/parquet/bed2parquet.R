
# testing
bed <- "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_31726.bed.gz"
tbi <- paste(bed, "tbi", sep=".") 

# if a tabix index exists, the file is sorted:
stopifnot(file.exists(tbi))

library(rtracklayer)
system.time(gr <- import(bed))
#    user  system elapsed
#  53.852   1.048  55.275

# we can also consult the tabix headers to at least guess the genome used
library(Rsamtools)
hdr <- headerTabix(tbi) 
contigs <- hdr$seqnames
# there's a trick to obtain the chrom.sizes which I forget at the moment
# anyways

system.time(tbxbed <- TabixFile(bed))
#    user  system elapsed 
#   0.001   0.000   0.001 
#
# The above is obviously not the same as computing overlaps 

# current fragment BEDs are actually a weird BED5-ish format 
# e.g. 
# GL000008.2      71      369     ATGGGAAC_TTGCCTAA_H05-17        1
# GL000008.2      109     200     TGGATCTG_CAGCAACG_G05-17        2

# BED4 columns: chrom, chromStart, chromEnd, name
# BED5 columns: chrom, chromStart, chromEnd, name, score
system.time(
  bed_data <- 
    read.table(bed, 
      sep="\t", header=FALSE,
      col.names=c("chrom", "chromStart", "chromEnd", "name", "score"),
      colClasses=c("character","integer","integer","character","numeric")
    )
)

# screw that, do it in parallel 
library(vroom)
bed5_col_types <- list(chrom = "c",
		      chromStart = "i",
		      chromEnd = "i",
    		      name = "c", 
	    	      score = "i")
system.time(
  bed_tbl <- 
    vroom::vroom(
		 bed,
		 delim = "\t",
		 col_names = names(bed5_col_types), 
		 col_types = bed5_col_types, 
		 col_select = 1:4 # for now
		 )
)
#    user  system elapsed                                                       
#   8.953   3.702  12.038 

# vroom can also automatically read from remote compressed files 

# can we coerce this directly to a GRanges?
# Of course:
system.time(bed_gr2  <- as(bed_tbl, "GRanges"))

#    user  system elapsed 
#  27.290   1.511  28.957 
#
# observe that reading with vroom and then coercing to a GR is faster
# than read.table or similar (!) 

library(GenomeInfoDb)
bed_genome <- "hg38" 
seqinfo <- SeqinfoForUCSCGenome(bed_genome)
shared <- intersect(seqlevels(seqinfo), hdr$seqnames)
tiles <- tileGenome(seqinfo[shared], tilewidth=50000, cut.last=TRUE)
system.time(tiles$frags <- countOverlaps(tiles, bed_gr2))
#    user  system elapsed 
#   2.890   0.237   3.142 


# obviously need to do this by cell barcode, which implies selection, so...

system.time(cells <- unique(bed_tbl$name))
#    user  system elapsed 
#   2.655   0.107   2.775 

length(cells)
# [1] 81729

# this would be a COUNT operation in duckdb obvs
frags_per_cell <- function(cell, frags) length(which(frags$name == cell))

# this is slow AF 
system.time(fpc <- sapply(cells, frags_per_cell, frags=bed_tbl))

# write to parquet
library(nanoparquet) 
stub <- sub("\\.bed\\.gz$", "", bed) 
bed_pqt <- paste(stub, "parquet", sep=".")
system.time(write_parquet(bed_tbl, bed_pqt))
#    user  system elapsed 
#  11.648   0.166  12.081 

# open with duckdb
# which takes FOREVER to install btw 
library(duckdb)

# see https://bwlewis.github.io/duckdb_and_r/ranges/ranges_redux.html
