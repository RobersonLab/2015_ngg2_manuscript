load_bioc_libraries = function( libName=NA )
{
	if ( !require( libName, character.only = TRUE ) ) {
		source( "http://www.bioconductor.org/biocLite.R" )
		biocLite( libName, character.only = TRUE, suppressUpdates=TRUE  )
		
		if ( !require( libName, character.only = TRUE ) )
		{
			stop( "Couldn't install", libName, "from Bioconductor. Please install manually." )
		}
	}
}

load_r_libraries = function( libName=NA )
{
	if ( !require( libName, character.only = TRUE ) ) {
		install.packages( libName )
		if ( !require( libName, character.only = TRUE ) ) {
			stop( "Couldn't install", libName, "from R packages. Please install manually." )
		}
	}
}

############
# packages #
############
load_bioc_libraries( "GenomicRanges" )
load_bioc_libraries( "GenomicFeatures" )
load_bioc_libraries( "plyr" )
load_bioc_libraries( "dplyr" )
load_bioc_libraries( "magrittr" )
load_bioc_libraries( "tidyr" )
load_bioc_libraries( "stringr" )
load_r_libraries( "knitr" )

##################
# parallel ddply #
##################
# this section should be modified if you want to run the analysis with parallel processing
# original analysis on ubuntu
# may or may not transfer to other platforms
# doMC can probably be replaced by your parallel backend of choice
# but be sure to modify the deregister at the end if necessary
PARALLEL = TRUE
NUM_CPUS = scan( "num_rscript_cpus", quiet=TRUE )[1] %>% as.integer

if ( NUM_CPUS > 1 )
{
	PARALLEL = TRUE
} else
{
	PARALLEL = FALSE
}

if ( PARALLEL )
{
	load_bioc_libraries( "doMC" )
	registerDoMC( NUM_CPUS )
}

#############
# functions #
#############
read_fasta_index = function( faiFilename )
{
	fai = read.table( faiFilename, header=FALSE, sep="\t", stringsAsFactors=FALSE )
	colnames( fai ) = c( "name", "size", "offset", "basesPerLine", "bytesPerLine" )
	return( fai )
}

read_grna_csv = function( grnaFilename )
{
	# grna file
	grna = read.csv( grnaFilename, header=TRUE, stringsAsFactors=FALSE )
  
	return( grna )
}

grep_and_return_record = function( values, grep_value )
{
	index = grep( grep_value, values )
	
	if ( length( index ) == 0 )
	{
		return( NA )
	} else
	{
		values[index] %>%
			str_replace_all( pattern=paste( grep_value, " ", sep="" ), "" ) %>%
			return
	}
}

gtf_info_record_to_name_list = function( s )
{
	vals = s %>%
		str_replace_all( pattern='"', "" ) %>%
		str_replace_all( "; ", ";" ) %>%
		str_split( pattern=";" ) %>%
		unlist
  
	gene_id = grep_and_return_record( vals, "gene_id" )
	gene_name = grep_and_return_record( vals, "gene_name" )
	gene_biotype = grep_and_return_record( vals, "gene_biotype" )
	
	return( c( "gene_id"=gene_id, "gene_name"=gene_name, "gene_biotype"=gene_biotype ) )
}

read_gtf_file = function( gtfFilename, faiObj, in_parallel=FALSE )
{
	# gtf file
	dat = read.table( gtfFilename, header=FALSE, sep="\t", stringsAsFactors=FALSE ) %>%
		filter( V3 == "gene" )
	
	other_info = ldply( dat$V9, .fun=gtf_info_record_to_name_list, .parallel=in_parallel )
	
	dat = dat[,c(1,2,4,5,7)]
	colnames( dat ) = c( "Contig", "Source", "Start", "End", "Strand" )
	
	dat = cbind( dat, other_info )
	rownames( dat ) = dat$gene_id %>% as.character
	
	dat %>%
		filter( Contig %in% faiObj$name & Source == 'ensembl' ) %>%
		return
}

grna_to_GRanges = function( grna_obj, fai_obj, genome_name )
{
	circ = rep( FALSE, dim( fai_obj )[1] )
	
	mitoIndex = which( fai_obj$name %in% c( "M", "chrM", "MT", "chrMT" ) )
	
	if ( length( mitoIndex ) > 0 )
	{
		circ[mitoIndex] = TRUE
	}
	
	# seq info from fai file
	mySeqInfo = Seqinfo( seqnames = fai_obj$name, seqlengths = fai_obj$size, isCircular = circ, genome = genome_name )
	
	# GRanges object
	gr = GRanges( seqnames = grna_obj$Contig, ranges = IRanges( start=grna_obj$Start, end=grna_obj$End), strand = grna_obj$Strand, mcols = grna_obj[,-c(1,2,3,6)], seqinfo = mySeqInfo )
	
	return( gr )
}

collapse_using_indexes = function( input_data_frame, gtf_df, exon_GRanges )
{
	gRNA = input_data_frame[1,"queryHits"]
	
	gene_id_names = input_data_frame[,"subjectHits"] %>%
		names( exon_GRanges )[.]
	
	gene_ids = gene_id_names %>%
		paste( collapse="|" )
	
	gene_names = gtf_df[gene_id_names, "gene_names"]
	
	if ( all( is.na( gene_names ) )  )
	{
		return( c( gRNA_Index = gRNA, gene_id = gene_ids, gene_name = "" ) )
	} else
	{
		gene_names[is.na( gene_names ) ] = "N/A"
		gene_names = paste( gene_names, collapse="|" )
		
		return( c( gRNA_Index=gRNA, gene_id=gene_ids, gene_name=gene_names ) )
	}
}

########################
# data capture objects #
########################
species = scan( 'species', 'character', quiet=TRUE )
pamSites = c( "AGG", "CGG", "GGG", "TGG" )
strandTypes = c( "+", "-", "*" )

countColumns = c( "total_gRNAs", "unique_gRNAs", "gStarts", "uniqGstarts", "MbGenome" )
countInfo = matrix( nrow=length( species ), ncol=length( countColumns ), dimnames=list( species, countColumns ) )

pamInfo = matrix( nrow=length( species ), ncol=length( pamSites ), dimnames=list( species, pamSites ) )
gStartPamInfo = pamInfo

strandInfo = matrix( nrow=length(species), ncol=length(strandTypes), dimnames=list( species, strandTypes) )
gStartStrandInfo = strandInfo

geneBiotypeAllgRNA = c()
geneBiotypeUniquegRNA = c()
geneBiotypeGstartgRNA = c()
geneBiotypeUniqueGstartgRNA = c()

#################################################
# analyze each species, accumulating group data #
#################################################
for ( species_index in 1:length( species ) )
{
	curr_species = species[species_index]
	################
	# import files #
	################
	fai = read_fasta_index( paste( "fasta/", curr_species, ".fa.fai", sep="" ) )
	
	gtf = read_gtf_file( paste( "gtfs/", curr_species, ".gtf", sep="" ), fai, in_parallel = PARALLEL )
	
	grna = paste( "grna_out/", curr_species, "_full_grnas.csv", sep="" ) %>%
		read.csv( header=TRUE, stringsAsFactors=FALSE ) %>%
		arrange( Contig, Start, End )
	
	txdb = makeTxDbFromGFF( paste( "gtfs/", curr_species, ".gtf", sep="" ), format="gtf" )
	exons = exonsBy( txdb, by="gene" )
	
	###########################
	# count gRNA occurences   #
	# identify unique cutters #
	###########################
	singleCutters = grna %>%
		subset( Unique == "Yes" ) %>%
		select( gRNA_Seq )
	singleCutters = singleCutters$gRNA_Seq
	
	listOfGstarts = grna %>%
		subset( G_start == "True" ) %>%
		select( gRNA_Seq )
	listOfGstarts = listOfGstarts$gRNA_Seq %>%
		unique
	
	#####################
	# convert to ranges #
	#####################
	gr.grna = grna_to_GRanges( grna, fai, curr_species )
	rm( grna ) # this is mostly just to save memory for the species with lots of grna sites

	############
	# overlaps #
	############
	olap = findOverlaps( gr.grna, exons, ignore.strand=TRUE )
	
	olap.names = data.frame(
		gRNA_Seq = mcols( gr.grna )[queryHits(olap),"mcols.gRNA_Seq"],
		gene_id = names(exons)[subjectHits(olap)] ) %>%
		mutate( gRNA_Seq = as.character( gRNA_Seq ) ) %>%
		mutate( gene_id = as.character( gene_id ) )
	
	olap.names.unique = olap.names %>%
		filter( gRNA_Seq %in% singleCutters )
 
	olap.names.gstart = olap.names %>%
		filter( gRNA_Seq %in% listOfGstarts )

	olap.names.gstart.unique = olap.names %>%
		filter( gRNA_Seq %in% singleCutters & gRNA_Seq %in% listOfGstarts )

	#######################
	# find gRNAs per gene #
	#######################
	# I use the gene_ids and gRNA_Seq for this
	# We're interested in total and unique cutters.
	# That is derived from the grna_seqs and not the PAM site
	grnasPerGene = olap.names$gene_id %>%
		table %>%
		data.frame
	colnames( grnasPerGene ) = c( "gene_id", "gRNAs" )
	
	# unique separate from all
	uniqueGrnasPerGene = olap.names.unique$gene_id %>%
		table %>%
		data.frame
	colnames( uniqueGrnasPerGene ) = c("gene_id", "unique_gRNAs" )

	# gstarts per gene
	gstartGrnaPerGene = olap.names.gstart$gene_id %>%
		table %>%
		data.frame
	colnames( gstartGrnaPerGene ) = c( "gene_id", "gStart_gRNAs" )

	# uniq gstarts per gene
	gstartUniqGrnaPerGene = olap.names.gstart.unique$gene_id %>%
		table %>%
		data.frame
	colnames( gstartUniqGrnaPerGene ) = c( "gene_id", "gStart_uniq_gRNAs" )
	
	# merge with gene data
	genesWithCounts = merge( gtf, grnasPerGene, by="gene_id", all.x=TRUE ) %>%
		merge( uniqueGrnasPerGene, by="gene_id", all.x=TRUE ) %>%
		merge( gstartGrnaPerGene, by="gene_id", all.x=TRUE ) %>%
		merge( gstartUniqGrnaPerGene, by="gene_id", all.x=TRUE ) %>% 
		mutate( gRNAs = ifelse( is.na( gRNAs ), 0, gRNAs ) ) %>%
		mutate( unique_gRNAs = ifelse( is.na( unique_gRNAs ), 0, unique_gRNAs ) ) %>%
		mutate( gStart_gRNAs = ifelse( is.na( gStart_gRNAs ), 0, gStart_gRNAs ) ) %>%
		mutate( gStart_uniq_gRNAs = ifelse( is.na( gStart_uniq_gRNAs ), 0, gStart_uniq_gRNAs ) ) %>%
		mutate( gene_name = ifelse( is.na( gene_name ), "", gene_name ) ) %>%
		select( Contig, Source, Start, End, Strand, gene_id, gene_name, gene_biotype, gRNAs, unique_gRNAs, gStart_gRNAs, gStart_uniq_gRNAs ) %>%
		arrange( Contig, Start, End, Strand, gene_id )
	
	########################
	# quick validity check #
	########################
	if (! all( genesWithCounts$unique_gRNAs <= genesWithCounts$gRNAs, na.rm=TRUE ) )
	{
		stop( "More unique gRNAs than total gRNAs for some sites!!!" )
	}
	
	if (! all( genesWithCounts$gStart_uniq_gRNAs <= genesWithCounts$gStart_gRNAs, na.rm=TRUE ) )
	{
		stop( "More unique G-start gRNAs that total G-start gRNAs for some sites!!!" )
	}
	
	# write gene info per organism
	write.csv( genesWithCounts, paste( "R_output/", curr_species, "_cutsPerGene.csv", sep=""), row.names=FALSE, quote=FALSE )
	
	#####################
	# instances of gRNA #
	#####################
	# I switch back to using gRNA and gtf index.
	# The gene overlapped by a specific gRNA depends on *site*, not sequence.
	
	grnaNumGeneOverlaps = olap %>%
		queryHits %>%
		table %>%
		data.frame
	colnames( grnaNumGeneOverlaps ) = c( "gRNA_Index", "Freq" )
	
	singleGeneOverlaps = subset( grnaNumGeneOverlaps, Freq == 1 )
	multiGeneOverlaps = subset( grnaNumGeneOverlaps, Freq > 1 )
 
	grnaGenes = olap %>%
		as.data.frame %>%
		filter( queryHits %in% multiGeneOverlaps$gRNA_Index ) %>%
		ddply( .variables="queryHits", .parallel=PARALLEL, .fun=collapse_using_indexes, gtf, exons ) %>%
		dplyr::select( -queryHits )
		
	if ( dim( grnaGenes )[1] != dim( multiGeneOverlaps )[1] )
	{
		stop( "The gRNAs that overlap multiple genes did not collapse into a dataframe correctly!"  )
	}
	
	singleIndexes = olap %>%
		as.data.frame %>%
		filter( queryHits %in% singleGeneOverlaps$gRNA_Index )

	grnaGenes = gtf %>%
		dplyr::select( gene_id, gene_name ) %>%
		merge( 
			data.frame(gRNA_Index = singleIndexes$queryHits, gene_id = names( exons )[ singleIndexes$subjectHits ] ),
			by="gene_id", all.y=TRUE ) %>%
		select( gRNA_Index, gene_id, gene_name ) %>%
		mutate( gene_id = as.character( gene_id ) ) %>%
		mutate( gene_id = ifelse( is.na( gene_id ), "", gene_id ) ) %>%
		mutate( gene_name = as.character( gene_name ) ) %>%
		mutate( gene_name = ifelse( is.na( gene_name ), "", gene_name ) ) %>%
		rbind( grnaGenes )
	
	if ( dim( singleIndexes )[1] != dim( grnaGenes )[1] - dim( multiGeneOverlaps )[1] )
	{
		stop( "The gRNAs that overlap single genes did not correctly bind correctly to the multiple gene overlaps!" )
	}
	
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
	
	grnasFromRanges = grnasFromRanges %>%
		merge( grnaGenes, by="gRNA_Index", all.x=TRUE ) %>%
		mutate( gene_id = ifelse( is.na( gene_id ), "", gene_id ) ) %>%
		mutate( gene_name = ifelse( is.na( gene_name ), "", gene_name ) ) %>%
		select( Contig, Start, End, Strand, gRNA_Seq, PAM, G_start, Unique ) %>%
		arrange( Contig, Start, Strand )
	
	# write per species
	write.csv( grnasFromRanges, paste( "R_output/", curr_species, "_gRNA_gene_annotation.csv", sep=""), row.names=FALSE, quote=FALSE )
	
	################################
	# Accumulate stats for species #
	################################
	
	# cut gene biotypes
	geneBiotypeAllgRNA  = genesWithCounts %>%
		select( gene_biotype, gRNAs ) %>%
		mutate( gRNAs=ifelse( gRNAs>0, "Cut", "Uncut" ) ) %>%
		xtabs( ~gRNAs+gene_biotype, data=. ) %>%
		as.data.frame %>%
		data.frame( species=curr_species, .) %>%
		rbind( geneBiotypeAllgRNA )
	
	geneBiotypeUniquegRNA  = genesWithCounts %>%
		select( gene_biotype, unique_gRNAs ) %>%
		mutate( gRNAs=ifelse( unique_gRNAs>0, "Cut", "Uncut" ) ) %>%
		xtabs( ~gRNAs+gene_biotype, data=. ) %>%
		as.data.frame %>%
		data.frame( species=curr_species, .) %>%
		rbind( geneBiotypeUniquegRNA  )
 
	geneBiotypeGstartgRNA  = genesWithCounts %>%
		select( gene_biotype, gStart_gRNAs ) %>%
		mutate( gRNAs=ifelse( gStart_gRNAs>0, "Cut", "Uncut" ) ) %>%
		xtabs( ~gRNAs+gene_biotype, data=. ) %>%
		as.data.frame %>%
		data.frame( species=curr_species, .) %>%
		rbind( geneBiotypeGstartgRNA  )
 
	geneBiotypeUniqueGstartgRNA  = genesWithCounts %>%
		select( gene_biotype, gStart_uniq_gRNAs ) %>%
		mutate( gRNAs=ifelse( gStart_uniq_gRNAs>0, "Cut", "Uncut" ) ) %>%
		xtabs( ~gRNAs+gene_biotype, data=. ) %>%
		as.data.frame %>%
		data.frame( species=curr_species, .) %>%
		rbind( geneBiotypeUniqueGstartgRNA  )

	#######################################################################
	# Reality check:                                                      #
	# The analysis does *not* use only G-start sites                      #
	# Total sites should therefore not be identical to only G-start sites #
	# Also it is not limited to unique sites. Therefore all sites should  #
	# not be identical to only unique sites                               #
	#######################################################################
	if ( all( geneBiotypeAllgRNA == geneBiotypeGstartgRNA ) )
	{
		stop( "There is not difference in all gRNA sites and only G-start sites!" )
	}
	
	if ( all( geneBiotypeAllgRNA == geneBiotypeUniquegRNA ) )
	{
		stop( "There is no difference in all gRNA sites and unique sites only!" )
	}
		
	# total gRNA
	countInfo[curr_species,"total_gRNAs"] = length( gr.grna )
	
	# unique gRNA
	countInfo[curr_species,"unique_gRNAs"] = length( singleCutters )
	
	# gStarts
	countInfo[curr_species, "gStarts"] = length( listOfGstarts )

	# uniqGstarts
	uniqGstartsCount = grnasFromRanges %>%
		filter( G_start == "True" & Unique == "Yes" ) %>%
		dim
	countInfo[curr_species, "uniqGstarts"] = uniqGstartsCount[1]
	
	# Mb genome
	countInfo[curr_species, "MbGenome"] = sum( as.numeric(fai$size) ) / 1E6
	
	# PAM sites
	pamTab = mcols( gr.grna )[,"mcols.PAM"] %>% table
	pamInfo[curr_species,] = pamTab[colnames(pamInfo)]
	
	# sense / antisense
	strandTab = strand( gr.grna ) %>% table
	strandInfo[curr_species,] = strandTab[colnames( strandInfo )]

	# gstart mcol index
	mcolGstartIndex = mcols( gr.grna ) %>%
		as.data.frame %>%
		filter( mcols.G_start == "True" ) %>%
		rownames %>%
		as.integer

	# PAM - canonical
	pamTab = mcols( gr.grna )[mcolGstartIndex,"mcols.PAM"] %>% table
	gStartPamInfo[curr_species,] = pamTab[colnames(gStartPamInfo)]

	# strand - canonical
	strandTab = strand( gr.grna )[mcolGstartIndex] %>% table
	gStartStrandInfo[curr_species,] = strandTab[ colnames( gStartStrandInfo ) ]
}

###########################
# write accumulated stats #
###########################
write.csv( countInfo, file="R_output/all_species_counts.csv", quote=FALSE )
write.csv( pamInfo, file="R_output/all_species_pam.csv", quote=FALSE )
write.csv( strandInfo, file="R_output/all_species_strand.csv", quote=FALSE )

write.csv( gStartPamInfo, file="R_output/all_species_gstart_pam.csv", quote=FALSE )
write.csv( gStartStrandInfo, file="R_output/all_species_gstart_strand.csv", quote=FALSE )

write.csv( geneBiotypeAllgRNA, file="R_output/all_species_all_gRNA_biotypes.csv", quote=FALSE, row.names=FALSE )
write.csv( geneBiotypeUniquegRNA, file="R_output/all_species_unique_gRNA_biotypes.csv", quote=FALSE, row.names=FALSE )
write.csv( geneBiotypeGstartgRNA, file="R_output/all_species_gstart_gRNA_biotypes.csv", quote=FALSE, row.names=FALSE )
write.csv( geneBiotypeUniqueGstartgRNA, file="R_output/all_species_gstart_unique_gRNA_biotypes.csv", quote=FALSE, row.names=FALSE )

if ( PARALLEL ) {
  registerDoSEQ()
}
