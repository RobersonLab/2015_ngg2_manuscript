load_bioc_libraries = function( libName=NA )
  {
  if ( !require( libName, character.only = TRUE ) ) {
    source( "http://www.bioconductor.org/biocLite.R" )
    biocLite( libName )
    if ( !require( libName, character.only = TRUE ) ) {
      stop( "Couldn't install", libName, "from Bioconductor. Please install manually." )
    }
  }
}

############
# packages #
############
load_bioc_libraries( "plyr" )
load_bioc_libraries( "dplyr" )
load_bioc_libraries( "magrittr" )
load_bioc_libraries( "tidyr" )
load_bioc_libraries( "stringr" )
load_bioc_libraries( "GenomicRanges" )

##################
# parallel ddply #
##################
# this section should be modified if you want to run the analysis with parallel processing
# original analysis on ubuntu
# may or may not transfer to other platforms
# doMC can probably be replaced by your parallel backend of choice
# but be sure to modify the deregister at the end if necessary
PARALLEL = FALSE
CPUS = 6

if ( PARALLEL ) {
  load_bioc_libraries( "doMC" )
  registerDoMC( CPUS )
}

#############
# functions #
#############
readFaiFile = function( faiFilename )
  {
  fai = read.table( faiFilename, header=FALSE, sep="\t", stringsAsFactors=FALSE )
  colnames( fai ) = c( "name", "size", "offset", "basesPerLine", "bytesPerLine" )
  return( fai )
  }

readGrnaFile = function( grnaFilename )
  {
  # grna file
  grna = read.csv( grnaFilename, header=TRUE, stringsAsFactors=FALSE ) %>%
    mutate( Strand = ifelse( Strand=="sense", "+", "-" ) )
  
  return( grna )
  }

grepAndReturnRecord = function( values, grepVal )
  {
  index = grep( grepVal, values )
  
  if ( length( index ) == 0 ) {
    return( NA )
  } else {
    
    values[index] %>%
      str_replace_all( pattern=paste( grepVal, " ", sep="" ), "" ) %>%
      return
  }
}

gtfInfoRecordToNameList = function( s )
  {
  vals = s %>%
    str_replace_all( pattern='"', "" ) %>%
    str_replace_all( "; ", ";" ) %>%
    str_split( pattern=";" ) %>%
    unlist
  
  gene_id = grepAndReturnRecord( vals, "gene_id" )
  gene_name = grepAndReturnRecord( vals, "gene_name" )
  gene_biotype = grepAndReturnRecord( vals, "gene_biotype" )
  
  return( c( "gene_id"=gene_id, "gene_name"=gene_name, "gene_biotype"=gene_biotype ) )
}

readGtfFile = function( gtfFilename, faiObj )
  {
  # gtf file
  dat = read.table( gtfFilename, header=FALSE, sep="\t", stringsAsFactors=FALSE )
  dat = subset( dat, V3 == "gene" )
 
  if ( PARALLEL ) {
    otherInfo = ldply( dat$V9, .fun=gtfInfoRecordToNameList, .parallel=TRUE )
  } else {
    otherInfo = ldply( dat$V9, .fun=gtfInfoRecordToNameList )
  }
  
  dat = dat[,-c(3,6,8,9)]
  colnames( dat ) = c( "Contig", "Source", "Start", "End", "Strand" )
  
  dat = cbind( dat, otherInfo )

  dat %>% subset( Contig %in% faiObj$name ) %>% return
  }

gtf2Ranges = function( gtfObj, faiObj, genomeName )
  {
  # seq info from fai file
  mySeqInfo = Seqinfo( seqnames = faiObj$name, seqlengths = faiObj$size, isCircular = rep(FALSE, dim(faiObj)[1]), genome = genomeName )
  
  gtfObj = subset( gtfObj, Contig %in% faiObj$name )
  
  # GRanges object
  gr = GRanges( seqnames = as.factor( gtfObj$Contig ), ranges = IRanges( start=gtfObj$Start, end=gtfObj$End ), strand = as.factor( gtfObj$Strand ), mcols = gtfObj[,c("Source", "gene_id", "gene_name", "gene_biotype")], seqinfo = mySeqInfo )
  return( gr )
}

grna2Ranges = function( grnaObj, faiObj, genomeName )
  {
  # seq info from fai file
  mySeqInfo = Seqinfo( seqnames = faiObj$name, seqlengths = faiObj$size, isCircular = rep(FALSE, dim(faiObj)[1]), genome = genomeName )
  
  # GRanges object
  gr = GRanges( seqnames = grnaObj$Contig, ranges = IRanges( start=grnaObj$Start, end=grnaObj$End), strand = grnaObj$Strand, mcols = grnaObj[,-c(1,2,3,6)], seqinfo = mySeqInfo )
  
  return( gr )
}

collapseUsingIndexes = function( inputDataframe, grnaGranges, gtfGranges )
  {
    
    gRNA = inputDataframe[1,"queryHits"]
    
    gene_ids = inputDataframe[,"subjectHits"] %>%
      mcols( gtfGranges )[.,"mcols.gene_id"] %>%
      paste( collapse="|" )
    
    gene_names = inputDataframe[,"subjectHits"] %>%
      mcols( gtfGranges )[.,"mcols.gene_name"]
    
    if ( all( is.na( gene_names ) )  ) {
      return( c( gRNA_Index = gRNA, gene_id = gene_ids, gene_name = "" ) )
    } else {
      gene_names[is.na( gene_names ) ] = "N/A"
      gene_names = paste( gene_names, collapse="|" )
      
      return( c( gRNA_Index = gRNA, gene_id=gene_ids, gene_name=gene_names ) )
    }
  }

########################
# data capture objects #
########################
species = c( 'saccharomyces_cerevisiae', 'caenorhabditis_elegans', 'drosophila_melanogaster', 'danio_rerio', 'mus_musculus', 'homo_sapiens' )

countColumns = c( "total_gRNAs", "unique_gRNAs", "gStarts", "uniqGstarts", "MbGenome" )
countInfo = matrix( nrow=length( species ), ncol=length( countColumns ), dimnames=list( species, countColumns ) )

pamSites = c( "AGG", "CGG", "GGG", "TGG" )
pamInfo = matrix( nrow=length( species ), ncol=length( pamSites ), dimnames=list( species, pamSites ) )

gStartPamInfo = pamInfo

strandTypes = c( "+", "-", "*" )
strandInfo = matrix( nrow=length(species), ncol=length(strandTypes), dimnames=list( species, strandTypes) )

gStartStrandInfo = strandInfo

geneBiotypeAllgRNA = c()
geneBiotypeUniquegRNA = c()
geneBiotypeGstartgRNA = c()
geneBiotypeUniqueGstartgRNA = c()

#################################################
# analyze each species, accumulating group data #
#################################################
for ( sp.index in 1:length( species ) ) {
  curr.species = species[sp.index]
  ################
  # import files #
  ################
  fai = readFaiFile( paste( curr.species, ".fa.fai", sep="" ) )
  gtf = readGtfFile( paste( curr.species, ".gtf", sep="" ), fai )
  grna = readGrnaFile( paste( curr.species, "_grnas.csv", sep="" ) )
  
  ###########################
  # count gRNA occurences   #
  # identify unique cutters #
  ###########################
  grnaSeqCounts = grna$gRNA_Seq %>% table %>% data.frame
  colnames( grnaSeqCounts ) = c( "gRNA_Seq", "NumPerfectGenomeMatches" )
  grnaSeqCounts = grnaSeqCounts %>%
    mutate( gRNA_Seq = as.character( gRNA_Seq ) )
  rownames( grnaSeqCounts ) = as.character( grnaSeqCounts$gRNA_Seq )
  
  singleCutters = grnaSeqCounts %>% subset( NumPerfectGenomeMatches == 1 )
  
  listOfGstarts = grna[ which( grna$G_start == "Yes" ), "gRNA_Seq" ] %>% unique
  
  #####################
  # convert to ranges #
  #####################
  gr.gtf = gtf2Ranges( gtf, fai, curr.species )
  
  gr.grna = grna2Ranges( grna, fai, curr.species )
  
  rm( grna )

  ############
  # overlaps #
  ############
  olap = findOverlaps( gr.grna, gr.gtf, ignore.strand=TRUE )
  
  olap.names = data.frame(
    gRNA_Seq = mcols( gr.grna )[queryHits(olap),"mcols.gRNA_Seq"],
    gene_id = mcols( gr.gtf )[subjectHits(olap),"mcols.gene_id"],
    gene_name = mcols( gr.gtf )[subjectHits(olap),"mcols.gene_name"] )
  
  olap.names.na.index = which( is.na( olap.names$gene_name ) )
  
  olap.names = olap.names %>%
    mutate( gRNA_Seq = as.character( gRNA_Seq ) ) %>%
    mutate( gene_id = as.character( gene_id ) ) %>%
    mutate( gene_name = as.character( gene_name ) )
  
  olap.names$gene_name[olap.names.na.index] = ""
  
  olap.names.unique = olap.names %>%
    subset( gRNA_Seq %in% singleCutters$gRNA_Seq )
 
  olap.names.gstart = olap.names %>%
    subset( gRNA_Seq %in% listOfGstarts )

  olap.names.gstart.unique = olap.names %>%
    subset( gRNA_Seq %in% singleCutters$gRNA_Seq & gRNA_Seq %in% listOfGstarts )

  #######################
  # find gRNAs per gene #
  #######################
  # I use the gene_ids and gRNA_Seq for this
  # We're interested in total and unique cutters.
  # That is derived from the actual seq, not just the site
  grnasPerGene = olap.names$gene_id %>% table %>% data.frame
  colnames( grnasPerGene ) = c( "gene_id", "gRNAs" )
  
  # unique separate from all
  uniqueGrnasPerGene = olap.names.unique$gene_id %>% table %>% data.frame
  colnames( uniqueGrnasPerGene ) = c("gene_id", "unique_gRNAs" )

  # gstarts per gene
  gstartGrnaPerGene = olap.names.gstart$gene_id %>% table %>% data.frame
  colnames( gstartGrnaPerGene ) = c( "gene_id", "gStart_gRNAs" )

  # uniq gstarts per gene
  gstartUniqGrnaPerGene = olap.names.gstart.unique$gene_id %>% table %>% data.frame
  colnames( gstartUniqGrnaPerGene ) = c( "gene_id", "gStart_uniq_gRNAs" )
  
  # merge with gene data
  genesWithCounts = merge( gtf, grnasPerGene, all.x=TRUE ) %>%
    merge( uniqueGrnasPerGene, all.x=TRUE ) %>%
    merge( gstartGrnaPerGene, all.x=TRUE ) %>%
    merge( gstartUniqGrnaPerGene, all.x=TRUE ) %>% 
    arrange( Contig, Start, End ) %>%
    mutate( gRNAs = ifelse( is.na( gRNAs ), 0, gRNAs ) ) %>%
    mutate( unique_gRNAs = ifelse( is.na( unique_gRNAs ), 0, unique_gRNAs ) ) %>%
    mutate( gStart_gRNAs = ifelse( is.na( gStart_gRNAs ), 0, gStart_gRNAs ) ) %>%
    mutate( gStart_uniq_gRNAs = ifelse( is.na( gStart_uniq_gRNAs ), 0, gStart_uniq_gRNAs ) ) %>%
    mutate( gene_name = ifelse( is.na( gene_name ), "", gene_name ) )
  
  rm( gtf ) 

  genesWithCounts = genesWithCounts[, c("Contig", "Source", "Start", "End", "Strand", "gene_id", "gene_name", "gene_biotype", "gRNAs", "unique_gRNAs", "gStart_gRNAs", "gStart_uniq_gRNAs" )]
  
  # write gene info per organism
  write.csv( genesWithCounts, paste( curr.species, "_cutsPerGene.csv", sep=""), row.names=FALSE, quote=FALSE )
  
  #####################
  # instances of gRNA #
  #####################
  # I switch back to using gRNA and gtf index.
  # The gene overlapped by a specific gRNA depends on *site*, not sequence.
  
  grnaNumGeneOverlaps = olap %>% queryHits %>% table %>% data.frame
  colnames( grnaNumGeneOverlaps ) = c( "gRNA_Index", "Freq" )
  
  singleGeneOverlaps = subset( grnaNumGeneOverlaps, Freq == 1 )
  multiGeneOverlaps = subset( grnaNumGeneOverlaps, Freq > 1 )
 
  if ( PARALLEL ) {
    grnaGenes = olap %>%
      as.data.frame %>%
      subset( queryHits %in% multiGeneOverlaps$gRNA_Index ) %>%
      ddply( .variables="queryHits", .fun=collapseUsingIndexes, gr.grna, gr.gtf, .parallel=TRUE ) %>%
      select( -queryHits )

  } else { 
    grnaGenes = olap %>%
      as.data.frame %>%
      subset( queryHits %in% multiGeneOverlaps$gRNA_Index ) %>%
      ddply( .variables="queryHits", .fun=collapseUsingIndexes, gr.grna, gr.gtf ) %>%
      select( -queryHits )
  }
  
  singleIndexes = olap %>%
    as.data.frame %>%
    subset( queryHits %in% singleGeneOverlaps$gRNA_Index )
  
  singleGenes = data.frame(
    gRNA_Index = singleIndexes$queryHits,
    gene_id = mcols( gr.gtf )[ singleIndexes$subjectHits, "mcols.gene_id" ],
    gene_name = mcols( gr.gtf )[ singleIndexes$subjectHits,"mcols.gene_name" ]
    )
  
  singleGenes$gene_id = as.character( singleGenes$gene_id )
  singleGenes.na.index = which( is.na( singleGenes$gene_name ) )
  singleGenes$gene_name = as.character( singleGenes$gene_name )
  singleGenes$gene_name[singleGenes.na.index] = ""
  
  grnaGenes = rbind( grnaGenes, singleGenes )
  
  grIndex = seq( 1, length( gr.grna ), by=1 )
  
  grnasFromRanges = data.frame(
    gRNA_Index = grIndex,
    Contig = as.character( seqnames( gr.grna[grIndex] ) ),
    Start = start( ranges( gr.grna[grIndex] ) ),
    End = end( ranges( gr.grna[grIndex] ) ),
    Strand = as.character( strand( gr.grna[grIndex] ) ),
    mcols( gr.grna )[grIndex,]
    )
  
  colnames( grnasFromRanges ) = gsub( "mcols.", "", colnames( grnasFromRanges ) )
  
  grnasFromRanges = merge( grnasFromRanges, grnaGenes, all.x=TRUE ) %>%
    mutate( gene_id = ifelse( is.na( gene_id ), "", gene_id ) ) %>%
    mutate( gene_name = ifelse( is.na( gene_name ), "", gene_name ) )
  
  grnasFromRanges = cbind( grnasFromRanges, PerfMatches = grnaSeqCounts[grnasFromRanges$gRNA_Seq,"NumPerfectGenomeMatches"] ) %>%
    mutate( UniqueSite = ifelse( PerfMatches == 1, "Yes", "No" ) ) %>%
    arrange( Contig, Start, Strand ) %>%
    select( -gRNA_Index )
  
  # write per species
  write.csv( grnasFromRanges, paste( curr.species, "_gRNA_gene_annotation.csv", sep=""), row.names=FALSE, quote=FALSE )
  
  ################################
  # Accumulate stats for species #
  ################################
  
  # cut gene biotypes
  geneBiotypeAllgRNA  = genesWithCounts[,c("gene_biotype", "gRNAs")] %>%
    mutate( gRNAs=ifelse( gRNAs>0, "Cut", "Uncut" ) ) %>%
    xtabs( ~gRNAs+gene_biotype, data=. ) %>%
    as.data.frame %>%
    data.frame( species=curr.species, .) %>%
    rbind( geneBiotypeAllgRNA )
  
  geneBiotypeUniquegRNA  = genesWithCounts[,c("gene_biotype", "unique_gRNAs")] %>%
    mutate( gRNAs=ifelse( unique_gRNAs>0, "Cut", "Uncut" ) ) %>%
    xtabs( ~gRNAs+gene_biotype, data=. ) %>%
    as.data.frame %>%
    data.frame( species=curr.species, .) %>%
    rbind( geneBiotypeUniquegRNA  )
 
  geneBiotypeGstartgRNA  = genesWithCounts[,c("gene_biotype", "gStart_gRNAs")] %>%
    mutate( gRNAs=ifelse( gStart_gRNAs>0, "Cut", "Uncut" ) ) %>%
    xtabs( ~gRNAs+gene_biotype, data=. ) %>%
    as.data.frame %>%
    data.frame( species=curr.species, .) %>%
    rbind( geneBiotypeGstartgRNA  )
 
  geneBiotypeUniqueGstartgRNA  = genesWithCounts[,c("gene_biotype", "gStart_uniq_gRNAs")] %>%
    mutate( gRNAs=ifelse( gStart_uniq_gRNAs>0, "Cut", "Uncut" ) ) %>%
    xtabs( ~gRNAs+gene_biotype, data=. ) %>%
    as.data.frame %>%
    data.frame( species=curr.species, .) %>%
    rbind( geneBiotypeUniqueGstartgRNA  )

  # total gRNA
  countInfo[curr.species,"total_gRNAs"] = length( gr.grna )
  
  # unique gRNA
  countInfo[curr.species,"unique_gRNAs"] = dim( singleCutters )[1]
  
  # gStarts
  countInfo[curr.species, "gStarts"] = length( listOfGstarts )

  # uniqGstarts
  uniqGstartsCount = grnasFromRanges %>%
    subset( G_start == "Yes" & PerfMatches == 1 ) %>%
    dim(.)
  countInfo[curr.species, "uniqGstarts"] = uniqGstartsCount[1]
  
  # Mb genome
  countInfo[curr.species, "MbGenome"] = sum( as.numeric(fai$size) ) / 1E6
  
  # PAM sites
  pamTab = mcols( gr.grna )[,"mcols.PAM"] %>% table
  pamInfo[curr.species,] = pamTab[colnames(pamInfo)]
  
  # sense / antisense
  strandTab = strand( gr.grna ) %>% table
  strandInfo[curr.species,] = strandTab[colnames( strandInfo )]

  # gstart mcol index
  mcolGstartIndex = which( mcols( gr.grna )[,"mcols.G_start"] == "Yes" )

  # PAM - canonical
  pamTab = mcols( gr.grna )[mcolGstartIndex,"mcols.PAM"] %>% table
  gStartPamInfo[curr.species,] = pamTab[colnames(gStartPamInfo)]

  # strand - canonical
  strandTab = strand( gr.grna )[mcolGstartIndex] %>% table
  gStartStrandInfo[curr.species,] = strandTab[ colnames( gStartStrandInfo ) ]
}

###########################
# write accumulated stats #
###########################
write.csv( countInfo, file="all_species_counts.csv", quote=FALSE )
write.csv( pamInfo, file="all_species_pam.csv", quote=FALSE )
write.csv( strandInfo, file="all_species_strand.csv", quote=FALSE )
write.csv( gStartPamInfo, file="all_species_gstart_pam.csv", quote=FALSE )
write.csv( gStartStrandInfo, file="all_species_gstart_strand.csv", quote=FALSE )
write.csv( geneBiotypeAllgRNA, file="all_species_all_gRNA_biotypes.csv", quote=FALSE, row.names=FALSE )
write.csv( geneBiotypeUniquegRNA, file="all_species_unique_gRNA_biotypes.csv", quote=FALSE, row.names=FALSE )
write.csv( geneBiotypeGstartgRNA, file="all_species_gstart_gRNA_biotypes.csv", quote=FALSE, row.names=FALSE )
write.csv( geneBiotypeUniqueGstartgRNA, file="all_species_gstart_unique_gRNA_biotypes.csv", quote=FALSE, row.names=FALSE )

if ( PARALLEL ) {
  registerDoSEQ()
}

