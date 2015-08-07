#!/usr/bin/python

##########
# Import #
##########
import argparse
import re
import sys
import time

####################
# Version and name #
####################
SCRIPT_PATH = sys.argv[0]
SCRIPT_NAME = SCRIPT_PATH.split( '/' )[-1].split( '\\' )[-1]
VERSION = '1.0.0'

########
# fxns #
########
def log_msg( msg, printTime=True ):
	if printTime:
		print time.strftime( '%H:%M %z on %b %d, %Y' )
	print msg
	sys.stdout.flush()
	
def error_msg( msg ):
	log_msg( "Error: %s\n" % ( msg.rstrip() ) )
	sys.exit( 1 )
	
###########
# classes #
###########
class ContigIndex:
	def __init__( self, string ):
		vals = string.split( '\t' )
		
		self.name = vals[0]
		self.size = int( vals[1] )
		self.byteStart = int( vals[2] )
		self.basesPerLine = int( vals[3] )
		self.bytesPerLine = int( vals[4] )
		self.extraBytesPerLine = self.bytesPerLine - self.basesPerLine
	
	def getBytePosition( self, position ):
		if position < 1 or position > self.size:
			return None
		
		generalOffset = self.byteStart + position - 1
		sequenceLineWraps = int( ( position - 1 ) / self.basesPerLine )
		whiteSpaceChars = self.extraBytesPerLine * sequenceLineWraps
		
		return generalOffset + whiteSpaceChars

if __name__ == '__main__':
	############
	# argparse #
	############
	parser = argparse.ArgumentParser( prog=SCRIPT_NAME, epilog="%s v%s" % ( SCRIPT_NAME, VERSION ) )
	
	parser.add_argument( 'fastaFile', help="Path to indexed FASTA file to be scanned for gRNA sites" )
	parser.add_argument( '--outputFile', help="Defaults to gRNAs.csv", default="gRNAs.csv" )
	parser.add_argument( '--fastaIndex', help="Expects the index to be the full FASTA name plus '.fai'. A different name can be specified here.", default=None )
	parser.add_argument( "--quiet", default=True, action='store_false', dest='verbose' )

	args = parser.parse_args()
	
	if args.fastaIndex == None:
		args.fastaIndex = "%s.fai" % ( args.fastaFile )
	
	options = "%s v%s\n\nOptions\n=======\n" % ( SCRIPT_NAME, VERSION )
	options += "FASTA: %s\n" % ( args.fastaFile )
	options += "FASTA Index: %s\n" % ( args.fastaIndex )
	options += "Output file: %s\n" % ( args.outputFile )
	options += "Verbose: %s\n" % ( str( args.verbose ) )

	if args.verbose:
		log_msg( options )

	#################
	# Compile RegEx #
	#################
	whiteSpaceOrN = re.compile( '\n|\r|N' )
	
	########################
	# Read the FASTA index #
	########################
	contigIndexDict = {}
	
	try:
		FAIFH = open( args.fastaIndex, 'r' )
	except:
		error_msg( "Couldn't open FASTA index [%s]" % ( args.fastaIndex ) )
		
	for line in FAIFH:
		line = line.rstrip()
		
		if len( line ) == 0 or line[0] == '#':
			continue
		else:
			contigData = ContigIndex( line )
			
			if contigData.name in contigIndexDict:
				error_msg( "Contig [%s] listed twice in FASTA index!" % ( contigData.name ) )
			else:
				contigIndexDict[ contigData.name ] = contigData
	
	FAIFH.close()
	
	###################
	# Open FASTA file #
	###################
	try:
		FAFH = open( args.fastaFile, 'r' )
	except:
		error_msg( "Could not open FASTA file [%s]" % ( args.fastaFile ) )
	
	###############
	# Open output #
	###############
	try:
		OUTFH = open( args.outputFile, 'w' )
	except:
		error_msg( "Could not open output file [%s] for writing" % ( args.outputFile ) )
	
	##################################
	# Operate over specified regions #
	##################################
	goodBases = 0
	gcBases = 0
	
	for contig in contigIndexDict:
		start = 1
		end = contigIndexDict[contig].size
	
		#################
		# Pull sequence #
		#################
		startBytes = contigIndexDict[contig].getBytePosition( start )
		endBytes = contigIndexDict[contig].getBytePosition( end )
		
		totalBytes = endBytes - startBytes + 1
		
		FAFH.seek( startBytes, 0 )
		
		sequence = FAFH.read( totalBytes )
		sequence = whiteSpaceOrN.sub( '', sequence )
		
		goodBases += len( sequence )
		gcBases += sequence.count( "G" ) + sequence.count( "g" ) + sequence.count( "C" ) + sequence.count( "c" )
		
	OUTFH.write( "%.4f" % ( float( gcBases ) / float( goodBases ) ) )
		
	OUTFH.close()
	FAFH.close()
