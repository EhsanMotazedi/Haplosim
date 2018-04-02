#!/usr/bin/env python
"""Takes for input two files containing blocks of variant haplotypes, as produced by haplogenertor, HapCompass (Aguiar and Istrail 2013, 
HapTree (Berger et al. 2014) and SDhaP (Das and Vikalo 2015); and computes for each corresponding block pair the perfect solution rate 
(Berger et al. 2014), as well as the vector error rate (Berger et al. 2014) and allelic correlation score defined as the probability of 
having the same allele at the same locus after having matched the haploytpes to the most similar using MEC criterion. In -t, --true mode, 
hapcompare assumes the first input file contains the true haplotypes (given as a single continuous block) and takes the average of the
scores if gaps have been introduced in estimation, i.e. in the second input file. In that case, the number of gaps per SNP is also given.
In a similar manner, VER is normalized by the number of SNPs and homologues in the true mode.
Written by Ehsan Motazedi, 08-10-2015, Wageningen UR
Last updated, 19-07-2016"""
import sys
import argparse
import copy
import functools
import os 
import os.path
import re 
import tempfile
import traceback
from collections import Counter
from collections import OrderedDict
from itertools import permutations
from math import isnan, exp
from scipy.special import gammaln

class CustomException(Exception):
	def __init__(self, value):
		self.parameter = value
	def __str__(self):
		return repr(self.parameter)

def nCr(n,r):
	return(int(round(exp(gammaln(n+1)-gammaln(r+1)-gammaln(n-r+1)))))

def get_match_score(hap1,hap2):
	"""Calculates the similarity ratio between the two input strings"""
	try:
		if not ((isinstance(hap1,tuple) or isinstance(hap1,str) or isinstance(hap1,list) ) 
				and (isinstance(hap2,tuple) or isinstance(hap2,str) or isinstance(hap2,list))):
			raise TypeError('Please provide string, list or tuple inputs!')
		elif len(hap1)!=len(hap2):
			raise CustomException("The input strings must have equal lengths!")
		else:
			_score=0
			for _n, _letter in enumerate(hap1):
				if _letter==hap2[_n]:
					_score+=1
			return(round(float(_score)/len(hap1),3))
	except CustomException as instance:
		print(instance)
		return(float('NaN'))
	except TypeError as err:
		print(err)
		print "Error in get_match_score()"
		return(float('NaN'))
	except:
		print "Unexpected error in get_match_score():"
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback,
									  limit=2, file=sys.stdout)
		return(float('NaN'))

def get_best_match(ref,**haps):
	"""Chooses the closest matching haplotype, compared to the first input, 
	from a set of input haplotypes; according to the similarity ratio from
	get_match_score()"""
	try:
		if ref and haps:
			scores={}
			for _name, _string in haps.iteritems():
				scores[_name]=get_match_score(ref,_string)
				if isnan(scores[_name]):
					raise RuntimeError('NaN values produced!')
				else:
					pass
			best_match=max(scores, key=scores.get)
			return(best_match, scores.get(best_match))
		else:
			raise CustomException("A reference haplotype, as well as matching candidate haplo's is \
 needed!")
	except CustomException as instance:
		print(instance)
		print "Error in get_best_match()"
		return(None, float('NaN'))
	except:
		print "Unexpected error in get_best_match():"
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback,
									  limit=2, file=sys.stdout)
		return(None,float('NaN'))

def get_correlation(block1, block2):
	"""Calculates the correlation score (0<r<1) between the two input blocks,
	suitable for small ploidy as of order O(k!)."""
	try:
		if isinstance(block1, list) and isinstance(block2, list):
			if len(block1) == len(block2):	
				n_total=len(block1)
				_hapnames=list('hap'+str(_i) for _i in range(1, 2*n_total+1))
				block1 = {_hapnames[x]: block1[x] for x in range(n_total)} 					 # First haplotype block
				block2 = {_hapnames[x]: block2[x-n_total] for x in range(n_total,2*n_total)} # Second haplotype block
				corr=[]
				match=[]
				for _i, _haps in enumerate(permutations(block1.keys())):
					match.append({})
					_corr=0
					for _n, _refs in enumerate(_haps):
						_corr+=get_match_score(block1[_refs], block2.values()[_n])# Each haplotype in the second block is matched to \ 
						if isnan(_corr):
							raise RuntimeError('NaN values produced!')
						else:
							match[_i].update({_refs:block2.keys()[_n]})			  # the corresponding haplotype in the permutation of block 1.
					corr.append(round(_corr/float(n_total),3))	
				return(match[corr.index(max(corr))],max(corr))
			else:
				raise CustomException("Two haplotype blocks of the same ploidy needed!")
		else:
			raise CustomException("Two haplotype blocks must be given as lists of haplotypes!")
	except CustomException as instance:
		print(instance)
		print "Error in get_correlation()"
		return(None, float('NaN'))
	except:
		print "Unexpected error in get_correlation():"
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback,
									  limit=2, file=sys.stdout)
		return(None, float('NaN'))

def getVER(A,B):
	g = homologue_group([A,B])
	mg = g.get_minimalgroups()
	return len(mg)-g.ploidy

def quick_get_correlation(block1, block2):
	"""Calculates a lower bound for the correlation score (0<r<1) between the two input blocks,
	suitable for large ploidy, e.g. k>5."""
	try:
		if isinstance(block1, list) and isinstance(block2, list):
			if len(block1) == len(block2):
				n_total=len(block1)
				_hapnames=list('hap'+str(_i) for _i in range(1, 2*n_total+1))
				block1 = {_hapnames[x]: block1[x] for x in range(n_total)} 					 # First haplotype block
				block2 = {_hapnames[x]: block2[x-n_total] for x in range(n_total,2*n_total)} # Second haplotype block
				corr=0
				match={}
				_tmp=copy.deepcopy(block2)
				for _refs in (block1.keys()):
					match[_refs], _dummy=get_best_match(ref=block1[_refs],**_tmp)
					if isnan(_dummy):
						raise RuntimeError('NaN values produced!')
					else:
						_tmp.pop(match[_refs])
						corr+=_dummy
				corr=round(corr/float(n_total),3)	
				return(match,corr)
			else:
				raise CustomException("Two haplotype blocks of the same ploidy needed!")					
		else:
			raise CustomException("Two haplotype blocks must be given as lists of haplotypes!")
	except CustomException as instance:
		print(instance)
		print "Error in quick_get_correlation()"
		return(None, float('NaN'))	
	except:
		print "Unexpected error in quick_get_correlation():"
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback,
									  limit=2, file=sys.stdout)
		return(None, float('NaN'))		

def HapLiner(header, line):			   	# Function to extract the columns containing variant alleles from each \
	_snp_cols=re.split('\s*|\t', line.rstrip())	# line of the variant haplotype file. 
	new_block=False
	if len(_snp_cols)==1:	# => That is, if the line is a blank line				
		out_snp=''
	elif _snp_cols[0].upper()=='BLOCK':  # => That is, if the line is a header line
		out_snp=''
	else: 
		if 'HapTree' in header.items()[0][1]:						# For each heterozygous site, i.e. line, the alleles are written to a list \									
			if header.values()[0][2]<header.values()[0][3]:
				out_snp=[str(header.values()[0][2])]+_snp_cols[1:]	# with the variant number is also added to the beginning of the list.
				header[header.items()[0][0]]=('HapTree',header.values()[0][1],
										header.values()[0][2]+1,header.values()[0][3])
			elif header.values()[0][2]==header.values()[0][3]:
				out_snp=[str(header.values()[0][2])]+_snp_cols[1:]
				header={x:header[x] for x in header if x!=header.keys()[0]}   # When the block variants are all read,\ 
				header=OrderedDict(sorted(header.items(), key=lambda t: t[0]))# the block is omitted from the header
				new_block=True
			else:
				raise CustomException("ERROR: Wrong header detected in the input file!")
		elif 'SDhaP' in header.items()[0][1]:			# SDhaP haplotypes must be in ConvertAlleleSDhaP format (ConvertAlelleSDhaP.py)
			if header.values()[0][1]>1:
				out_snp=_snp_cols
				header[header.items()[0][0]]=('SDhaP',header.values()[0][1]-1)
			elif header.values()[0][1]==1:
				out_snp=_snp_cols
				header={x:header[x] for x in header if x!=header.keys()[0]}   # When the block variants are all read,\ 
				header=OrderedDict(sorted(header.items(), key=lambda t: t[0]))# the block is omitted from the header
				new_block=True
			else:
				raise CustomException("ERROR: Wrong header detected in the input file!")
		elif 'HapComp' in header.items()[0][1]:	
			if header.values()[0][2]<header.values()[0][3]:
				out_snp=[_snp_cols[2],_snp_cols[1]]+_snp_cols[3:]
				header[header.items()[0][0]]=('HapComp',header.values()[0][1],
											header.values()[0][2]+1,header.values()[0][3])
			elif header.values()[0][2]==header.values()[0][3]:
				out_snp=[_snp_cols[2],_snp_cols[1]]+_snp_cols[3:]
				header={x:header[x] for x in header if x!=header.keys()[0]}   # When the block variants are all read,\ 
				header=OrderedDict(sorted(header.items(), key=lambda t: t[0]))# the block is omitted from the header
				new_block=True
			else:
				raise CustomException("ERROR: Wrong header detected in the input file!")
		elif 'HapGen' in header.items()[0][1]:	
			if header.values()[0][1]>1:
				out_snp=[_snp_cols[0],_snp_cols[2]]+_snp_cols[5:]
				header[header.items()[0][0]]=('HapGen',header.values()[0][1]-1)
			elif header.values()[0][1]==1:
				out_snp=[_snp_cols[0],_snp_cols[2]]+_snp_cols[5:]
				header={x:header[x] for x in header if x!=header.keys()[0]}   # When the block variants are all read,\ 
				header=OrderedDict(sorted(header.items(), key=lambda t: t[0]))# the block is omitted from the header
				new_block=True
			else:
				raise CustomException("ERROR: Wrong header detected in the input file!")
	return(out_snp, header, new_block)

class homologue_group(object):
	""" Detect the segments of no-switch between the two phasings, divide the phasings into the segments and obtain the minimal set of common segments between the two phasings."""
	def __init__(self, phasings):
		self.size = len(phasings[0][0]) 	# Size of the genome, i.e. number of SNPs
		self.ploidy = len(phasings[0])  	# The ploidy level
		self.homologues = [[phasings[i][k] for i in range(len(phasings))] for k in range(self.ploidy)] # A list whose elements are lists. Each list element contains the two phasings (from the two genomes) for one of the homologues.
		self.phasings = copy.deepcopy(phasings) # A list containing the two genomes phasings.
		self.segments = []
		self.groups = []
	
	def put_segs(self):
		""" Determine the no-switch segments for the phasings and put it in 'segments' field."""
		segments = []
		_start = 0
		new_segment=True
		while _start < self.size:
			if new_segment:
				_genotypes=[[[] for _x in range(0,self.ploidy)] for _y in range(len(self.phasings))]
				segments.append([])
			for _genome in range(len(self.phasings)):
				for _homolo in range(0, self.ploidy):
					_genotypes[_genome][_homolo].append(self.phasings[_genome][_homolo][_start]) # Extract the two genotype vectors for SNP pisition _start and append them to the current segment.
			new_segment=not sorted(_genotypes[0])==sorted(_genotypes[1]) # If no switch has occured, the two segments must be the same except maybe for the order of their homologues.
			if not new_segment:
				segments[-1].append(_start)
				_start+=1
		self.segments = segments

	def put_seggroups(self):
		""" Divide the phasings into the detected segments."""
		if not self.segments:
			self.put_segs()
		groups = [[[[] for _n in range(len(self.segments))] for _k in range(self.ploidy)] for _n in range(len(self.phasings))]
	        for _n, _seg in enumerate(self.segments):
			for _genome in range(len(self.phasings)):
				for _homolo in range(self.ploidy):
					groups[_genome][_homolo][_n].extend(self.homologues[_homolo][_genome][_x] for _x in _seg) 	
		self.groups = groups

	def get_minimalgroups(self):
		""" Obtain the minimal set of common segments between the two phasings."""
		def _hash(optim, *args):
			_nlsts=[dict() for _lst in args]
			for _n, _lst in enumerate(args): 
				for _m, _y in enumerate(_lst):
					_h=''
					for _x in _y:
						if _x!=[]:
							_h+=''.join(str(_z) for _z in _x)
						else:
							_h+='.'
					if optim: 
						_nlsts[_n][_m]=(''.join(str(_c) for _c in [ord(_x)^ord(_y) for _x,_y in zip(_h, '1'*len(_h))]), ''.join(bin(int(_x, 16))[2:].zfill(2)[0:2] for _x in toHex(str(int(toHex(_h), 16)% 2**(len(_h)-1)))).zfill(8)[0:8])
					else:
						_nlsts[_n][_m] = (_h, 0)
			return _nlsts 		
		if not self.groups:
			self.put_seggroups()
		merged_segs, _verified = [], 0
		_optim=-1
		while _verified < self.ploidy*len(self.segments):
			_optim+=1
			_homolo, _n, new_merge = 0, 0, True
			hashes=_hash(_optim<self.ploidy*len(self.segments), *self.groups)
			for _genome in range(len(self.phasings)):
				self.groups[_genome]=sorted(self.groups[_genome], key=lambda x: hashes[_genome][self.groups[_genome].index(x)])
			while True:
				if _homolo<self.ploidy:
					if _n<len(self.segments):
						if new_merge:
							merged_segs.append([])
						new_merge = not self.groups[0][_homolo][_n] == self.groups[1][_homolo][_n]
						if new_merge:
							_n=0
							_homolo+=1
						else:
							if self.groups[0][_homolo][_n]:
								merged_segs[-1].extend(self.groups[0][_homolo][_n])
								for _genome in range(len(self.phasings)):
									self.groups[_genome][_homolo][_n]=[]
								_verified+=1
							_n+=1
					else:
						new_merge=True
						_n = 0
						_homolo+=1
				else:
					break
		return [_mseg for _mseg in merged_segs if _mseg!=[]] 

def toHex(string):
    """Convert a string to its hexadecimal representation. Digits keep their integer values and '.' has zero value. Other characters are represented by their unicode equivalent."""
    lst = []
    for _chr in string:
        if _chr=='.':
		continue
	if _chr not in [str(_x) for _x in range(10)]: 
		hex_chr = hex(ord(_chr)).replace('0x', '')
        else:
		hex_chr = hex(int(_chr)).replace('0x', '')
	if len(hex_chr) == 1:
            hex_chr = '0'+hex_chr
        lst.append(hex_chr)
    if len(lst)==0:
	return('00')
    return functools.reduce(lambda x,y: x+y, lst)

if __name__ == "__main__":					    # <<The __main__ module starts here!>>
	try:
		exit_msg=0		
		hap1=None       # File stream for haplotypes 1
		hap2=None       # File stream for haplotypes 2
		tmp=()	        # Temporary file for extracting haplotype variants 
		parser = argparse.ArgumentParser(description='Compares two variant haplotype files and\
		reports the perfect solution rate, as well as the vector error rate between them. The first file is\
		assumed to contain the true haplotypes by definition. In case the comparison measures are reported for\
		multiple blocks, their order will be specified by the position of the first variant within the blocks.', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog='Notes: If -t, --true is not specified, the two input files must contain equal number of blocks, separated by header lines.\n\
If -t, --true is specified, the first file must contain ONE SINGLE haplotype block!')
		parser.add_argument('file1', metavar='hap_file1', type=str, nargs=1,
						   help='the first (true) variant haplotype file.',default=None)
		parser.add_argument('file2', metavar='hap_file2', type=str, nargs=1,
						   help='the second variant haplotype file.', default=None)
		parser.add_argument('-v', '--verbose', action='store_true',
						   help='verbose flag to print warning messages.', default=False)
		parser.add_argument('-q', '--quick', action='store_true',
						   help='flag to calculate a lower bound of the correlation score between two\
 blocks instead of the MEC based correlation score. Suitable for high ploidies, e.g. >8.', default=False)
		parser.add_argument('-t', '--true', action='store_true',
						   help='flag indicating that the first file contains the original haplotypes, i.e. not the estimated ones.',
						   default=False)
		try:
			if not len(sys.argv)>1:
				parser.parse_args('-h'.split())
			else:
				args = parser.parse_args() 
		except SystemExit:
			exit_msg=3
			raise
		verbose=args.verbose
		quick=args.quick
		true_hap=args.true
		hapfiles=(args.file1[0],args.file2[0])
		header=[dict(),dict()]
		block_count=[0,0]	
		for _i, _hapfile in enumerate(hapfiles):
			with open(_hapfile,'rU') as _hap:	# First read the header of each haplotype block
				for _header_line in _hap:								 # The type of input the varinat hapltype files: SDhaP/HapTree/HapComp/HapGen \
					_header_cols=re.split('\s*|\t',_header_line.rstrip())# representing SDhaP_Converted, HapTree, HapCompass and haplogenerator formats					
					_header_count=len(_header_cols)
					if (_header_count>8 and _header_cols[0]=='Block' and _header_cols[1]=='length' and _header_cols[3]=='var_id' and _header_cols[4]=='contig' and 
					_header_cols[5]=='varpos' and _header_cols[6]=='ref_allele' and _header_cols[7]=='alt_allele'):
						block_count[_i]+=1
						header[_i].update({block_count[_i]:('HapGen',int(_header_cols[2]))}) # (file type, length of the block)
					elif _header_count==14 and _header_cols[0]=='Block':
						block_count[_i]+=1
						header[_i].update({block_count[_i]:('SDhaP',int(_header_cols[6]))})	# (filetype, length of the block)
					elif _header_count==4 and _header_cols[0]=='BLOCK': 	 # Property of a header line, as distinguished from the allele lines
						block_count[_i]+=1
						header[_i].update({block_count[_i]:('HapTree',int(_header_cols[1]),
														int(_header_cols[2]),int(_header_cols[3]))}) # (filetype, first variant pos in the block, first and last variant numbers in the block)
					elif _header_count==7 and _header_cols[0]=='BLOCK':
						block_count[_i]+=1
						header[_i].update({block_count[_i]:('HapComp',int(_header_cols[1]), # (filetype, first variant pos in the block, first and last variant numbers in the block)
							int(_header_cols[3]),int(_header_cols[4]))})					
					else:
						pass
			_hap=None
			_header_line=None
			_header_cols=None
			_header_count=None
		del(block_count, _header_line, _header_count, _header_cols)
		_i=None
		_hapfile=None
		if len(header[0])==0 or len(header[1])==0:
			raise CustomException("ERROR: Input format must be of HapTree/HapCompass/haplogenerator/\
SDhaP. The header is compatible with none!")
		elif len(header[0])!=len(header[1]) and true_hap==False:
			raise CustomException("ERROR: The two inputs contain different number of haplotype blocks!\
Comparison is not possible unless -t, --true flag is on!")
		elif len(header[0])!=1 and true_hap==True:
			raise CustomException("ERROR: The first haplotype file is assumed to contain the true haplotypes!\
 It should therefore contain just one haplotype block!")
		else:
			for _i in (0,1):
				header[_i]= OrderedDict(sorted(header[_i].items(), key=lambda t: t[0]))
			n_blocks=len(header[1])			 # Number of haplotype blocks present in the estimation file
		tmp=['','']								 # Extract the alleles for each heterozygous site\  
		for _i,_hapfile in enumerate(hapfiles):	 # and writing them to tmp files.
			with open(_hapfile,'rU') as _hap, tempfile.NamedTemporaryFile(delete=False) as snps_hap:	
				_block=1						 # Block id to be added to the variant info
				for line in _hap:
						if not (header[_i] or line=='\n'):
							raise CustomException("ERROR: The input file or its header info is corrupt!")
						_snp_out,header[_i], new_block=HapLiner(header[_i],line)
						if _snp_out:	#_snp_out contains [variant number, variant VCF position, allele1, allele2, ...]
							snps_hap.write('['+'Block'+str(_block)+','+','.join(_snp_out)+'];')
						else:			# Block id (Block1, Block2,...) is added to the beginning of _snp_out
							snps_hap.write('['+','.join(_snp_out)+'];')
						if new_block:
							_block+=1
				tmp[_i]=snps_hap.name
			snps_hap=None
			line=None
			_snp_out=None
		_i=None
		_block=None
		snps=[[],[]]
		__snps__=None
		for _n,_tmp in enumerate(tmp):				# Put all the variants (from each variant haplotype file) into a list
			with open(_tmp,'rU') as snps_hap:		# The two lists are in turn put into a list, which is called 'snps'. Comparisons will be performed 'snps'.
				__snps__=list(_snp[1:-1].split(',') for _snp in snps_hap.readline().rstrip(' ;\n').split(';') if _snp!='[]')
				snps[_n]=list([_snp[0]]+map(int,_snp[1:]) for _snp in __snps__)
				__snps__=None
			snps_hap=None
		del(__snps__)
		Correct_haplos=list(0 for x in range (0,n_blocks)) # The proportion of correctly estimated homologues within each block
		PSR=list(0 for x in range (0,n_blocks))	# The Perfect Solution Rate PSR (also called Phasing Accuracy Rate) for each block
		Total_PSR=list(0 for x in range (0,n_blocks))	# PSR for each block including variants with discordant dosages
		VER=list(0 for x in range (0,n_blocks))	# The Vector Error Rate VER for each block
		corr=list(0 for x in range (0,n_blocks))# Allelic Correlation AC for each block pair
		Total_corr=list(0 for x in range (0,n_blocks))# AC for each block pair for each block including variants with discordant dosages 
		Missing_Variants_Dosages=list(0 for x in range (0,n_blocks)) # Number of variant loci within each block, whose dosages differ between the two block
		Missing_Variants_Uncommon=list(0 for x in range (0,n_blocks))# Number of variant loci within each block, which have been present in just one input file.
		snps_in_block=list([[],[]] for x in range(0,n_blocks))	# Variants sitting in each haplotype block that are common between the two inputs
		_snps_in_block=list([[],[]] for x in range(0,n_blocks))	# Temporary list of block variants
		if not true_hap:
			for _n, __snps__ in enumerate(snps):	# The two files contain equal number of blocks, which are matched in the coming order
				for _snp in __snps__:				# Group the corresponding haplotype blocks into a list
					for _block in range(1, n_blocks+1):
						if _snp[0]==('Block'+str(_block)):
							_snps_in_block[_block-1][_n].append(_snp[1:])
							break
				sorted_blocks_input=sorted((_blk[_n] for _blk in _snps_in_block), key= lambda x: min(t[1] for t in x) if len(x)>=2 else x[0][1] if len(x)==1 else -1) # Blocks are sorted according to their start variant position
				for _block in range(1, n_blocks+1):													# This is important in order to later match the blocks for comparison
					_snps_in_block[_block-1][_n]=copy.deepcopy(sorted_blocks_input[_block-1])
				sorted_blocks_input, _blk = (None, None)
		else:							# The first file contains just one haplotype block with true_hap
			Total_var_true=len(snps[0]) # Total number of true variants
			Total_miss_uncommon_true=0  # Initial number of missing variants, increment by 1 as missing true variants are observed
			Total_spurious_estimate=0   # Initial number of spuriously detected variants, increment by 1 as false variants are observed
			spurious_indx=[]			# The index of spurious variants
			if Total_var_true==0:
				raise CustomException("ERROR: the original haplotype contains no variants! Nothing to compare!")
			else:
				for _i, _snp in enumerate(snps[1]): # The SNPs from the estimated haplotypes that have no counterpart in the original haplotypes are counted as "spurious" 	
					if _snp[2] in (_x[2] for _x in snps[0]):
						pass
					else:
						Total_spurious_estimate+=1
						spurious_indx.append(_i) 
				for _x in reversed(spurious_indx):
					del snps[1][_x]
			for _snp in snps[1]:	# Separate the blocks of the estimated haplotypes file	
				for _block in range(1, n_blocks+1):
					if _snp[0]==('Block'+str(_block)):
						_snps_in_block[_block-1][1].append(_snp[1:])
						break
			sorted_blocks_input=sorted((_blk[1] for _blk in _snps_in_block), key= lambda x: min(t[1] for t in x) if len(x)>=2 else x[0][1] if len(x)==1 else -1) # Blocks are sorted according to their start variant position
			for _block in range(1, n_blocks+1):													# This is important in order to later match the blocks for comparison
				_snps_in_block[_block-1][1]=copy.deepcopy(sorted_blocks_input[_block-1])
			sorted_blocks_input, _blk = (None, None)
			for _snp in snps[0]: # Break the large block of the true haplotype file into blocks, corresponding to the blocks of the estimated haplotypes file
				for _block in range(1, n_blocks+1): # Correspondance is checked using the VCF position for each variant
					if _snp[2] in (_x[1] for _x in _snps_in_block[_block-1][1]):
						_snps_in_block[_block-1][0].append(_snp[1:])
						break
			sorted_blocks_input=sorted((_blk[0] for _blk in _snps_in_block), key= lambda x: min(t[1] for t in x) if len(x)>=2 else x[0][1] if len(x)==1 else -1) # Blocks are sorted according to their start variant position
			for _block in range(1, n_blocks+1):													# This is important in order to later match the blocks for comparison
				_snps_in_block[_block-1][0]=copy.deepcopy(sorted_blocks_input[_block-1])
			sorted_blocks_input, _blk = (None, None)
			for _snp in snps[0]: # The SNPs from the original haplotype that have no counterpart in the estimated haplotypes are counted as "missing" 	
				already_included=False
				for _block in range(1, n_blocks+1):
					if _snp[2] in (_x[1] for _x in _snps_in_block[_block-1][0]):
						already_included=True
						break
				if not already_included:
					Total_miss_uncommon_true+=1				
			already_included=None
		_snp, __snps__=(None,None)				
		_block, _n, snps, ploidy1, ploidy2 = (None, None, None, 0, 0)
		Total_variants=[]				# Total number of variants within each block
		for y in _snps_in_block:	# Calculate the ploidy level from the number of columns in the input haplotype blocks 
			if y[0]!=[]: # Escape empty blocks
				ploidy1=len(zip(*list(x[2:] for x in y[0])))
				break
		for y in _snps_in_block:
			if y[1]!=[]:
				ploidy2=len(zip(*list(x[2:] for x in y[1])))
				break
		for _ploidy in [ploidy1, ploidy2]:
			if not _ploidy:
				raise CustomException("ERROR: At least One of the input files did not contain any haplotypes or did not have the required format!")
		if ploidy2==ploidy1:
			ploidy=ploidy1
			del ploidy1, ploidy2
		else:
			raise CustomException("ERROR: Different ploidy levels detected for the two input files! Ploidies must be the same\
 for meaningful comparisons!")
		if verbose:
			print("The detected ploidy level is {} ...".format(ploidy))
		for _block, snps in enumerate(_snps_in_block):    # Throw away the variants which are uncommun between the two input files
			_warn=False	
			Total_variants.append(len(snps[0]+snps[1]))  
			_all_snps=Counter(_snp[1] for _snp in snps[0]+snps[1]) # Accumulate all the variant positions
			for key, value in _all_snps.iteritems():
				if value==2:
					snps_in_block[_block][0].extend(list(_snp for _snp in snps[0] if _snp[1]==key))
					snps_in_block[_block][1].extend(list(_snp for _snp in snps[1] if _snp[1]==key))
				elif value>2:
					raise CustomException("ERROR: Same variant detected two times or more within Block "+
						str(_block+1)+"! The haplotypes are not valid!")
				else:
					_warn=True
					Missing_Variants_Uncommon[_block]+=1
			if verbose and _warn:
				print("WARNING: one or more loci of Block %s were missing in one of the inputs!") %(_block+1) 		
			snps_in_block[_block][0]=sorted(snps_in_block[_block][0], key=lambda t: t[0]) # Sort the variants according to their attached id number
			snps_in_block[_block][1]=sorted(snps_in_block[_block][1], key=lambda t: t[0])
		_block, snps, _snps_in_block = (None, None, list([[],[]] for x in range(0,n_blocks)))
		total_correcthaplos_count=0
		for _block, snps in enumerate(snps_in_block): 
			true_haplos_count=0				# The count of the true whole block phasings
			true_phase_count=0				# The count of the true pairwise phasings 
			n_snps=len(snps[0])				# Total number of SNPs within the block
			haplos=(zip(*(snps[0][i][2:] for i in range(0,n_snps))), 
					zip(*(snps[1][j][2:] for j in range(0,n_snps))))
			for _haplo in haplos[0]:
				if _haplo in haplos[1]:
					true_haplos_count+=1
					haplos[1].remove(_haplo)
			for i in range(0,n_snps-1):										# Total PSR is defined the percentage of variant site pairs\				
				for j in range(i+1,n_snps): 							  	# whose phase are concordant between the corresponding blocks\		 
					pair_phase_true=zip(*(snps[0][i][2:],snps[0][j][2:])) 	# of two haplotype files, including variants with discordant dosages
					pair_phase_estim=zip(*(snps[1][i][2:],snps[1][j][2:]))
					if Counter(pair_phase_true)==Counter(pair_phase_estim): # This condition means the two phasings agree
						true_phase_count+=1
					else:
						pass
			if n_snps==0 or n_snps==1:
				Total_PSR[_block]=float(0)
				Total_corr[_block]=float(0)
				Correct_haplos[_block]=float('NaN')
				total_correcthaplos_count=float('NaN')
			else:	# Total PSR and the probability of correct estimation of the haplotypes with three decimals precision for each block
				Correct_haplos[_block]=round(float(true_haplos_count)/ploidy*1000)/1000. 
				Total_PSR[_block]=round(float(true_phase_count)/(n_snps*(n_snps-1)/2)*1000)/1000.
				total_correcthaplos_count+=true_haplos_count
				if quick and verbose:
					print("WARNING: the total correlation score is a lower bound of the MEC correlation score!")
				if quick:
					Total_corr[_block] = quick_get_correlation(zip(*list(x[2:] for x in snps[0])),
												zip(*list(x[2:] for x in snps[1])))[1]
				else:
					Total_corr[_block] = get_correlation(zip(*list(x[2:] for x in snps[0])),
												zip(*list(x[2:] for x in snps[1])))[1]	
			true_haplos_count, true_phase_count, pair_phase_estim, pair_phase_true = (None, None, None, None)		
		_haplo, _block, snps = (None, None, None)
		Total_common_variants=[] 			# Total number of common variants within each block
		for _block, snps in enumerate(snps_in_block):      # Throw away variants whose dosages are not the same between the two blocks
			_warn=False
			Total_common_variants.append(len(snps[0]+snps[1])) 
			for _n, _snp in enumerate(snps[0]):
				if sum(_snp[2:]) == sum (snps[1][_n][2:]): # Check if the dosages are the same at the same locus
					_snps_in_block[_block][0].append(_snp)
					_snps_in_block[_block][1].append(snps[1][_n])
				else:
					_warn=True
					Missing_Variants_Dosages[_block]+=2
			if verbose and _warn:
				print("WARNING: one or more variants of Block %s did not have the same dosage within the two inputs!") %(_block+1)				
		snps_in_block=copy.deepcopy(_snps_in_block)		
		_block, snps, _n, _snps_in_block, _warn = (None, None, None, None, None)
		for _block, snps in enumerate(snps_in_block): 
			true_phase_count=0
			n_snps=len(snps[0])			    # PSR is defined similar to Total PSR, but it is calculated with 
			for i in range(0,n_snps-1):	    # only those SNPs whose dosages have been correctly estimated.
				for j in range(i+1,n_snps):
					pair_phase_true=zip(*(snps[0][i][2:],snps[0][j][2:]))
					pair_phase_estim=zip(*(snps[1][i][2:],snps[1][j][2:]))
					if Counter(pair_phase_true)==Counter(pair_phase_estim):   # This condition means the two phasings agree
						true_phase_count+=1
					else:
						pass
			if n_snps==0 or n_snps==1:
				PSR[_block]=float(0)
				corr[_block]=float(0)
			else:							# PSR with three decimals precision 
				PSR[_block]=round(float(true_phase_count)/(n_snps*(n_snps-1)/2)*1000)/1000. 			
				if quick and verbose:									
					print("WARNING: the correlation score is a lower bound of the MEC correlation score!")
				if quick: # Calculate allelic correlation after dropping the SNPs whose dosages have been wrongly estimated
					corr[_block] = quick_get_correlation(zip(*list(x[2:] for x in snps[0])),
												zip(*list(x[2:] for x in snps[1])))[1]
				else:
					corr[_block] = get_correlation(zip(*list(x[2:] for x in snps[0])),
												zip(*list(x[2:] for x in snps[1])))[1]	 
			true_phase_count, pair_phase_estim, pair_phase_true=(None,None,None)
			if n_snps==0 or n_snps==1:
				VER[_block]=float("inf")
			else:
				phase_true=tuple(zip(*list(x[2:] for x in snps[0])))    # Calculate the Vector Error Rate (VER)
				phase_estim=tuple(zip(*list(x[2:] for x in snps[1])))   # VER is defined as the minimum number of preserved \
				VER[_block]=getVER(phase_estim, phase_true)	        # segments between the true and estimated homolgues, minus the ploidy; corresponding to the number of switches \
		if not true_hap:         # between the true homologues resulting in the estimated homologues.
			if verbose:
				print >>sys.stderr, "WARNING: the measures will be given for each block, in the order specified by the position of the first variant within the block!"
			print "Allelic correlation between the corresponding blocks (including variants with discordant dosages):{0}\
\nPerfect Solution Rates for each block pair (including variants with discordant dosages): {1}".format(' ,'.join(map(str,Total_corr)),' ,'.join(map(str,Total_PSR))) 
			print "Vector Error Rates for common variants of each block pair: {0} \nPerfect Solution Rates for each block pair (excluding variants with discordant dosages): {1}".format(' ,'.join(map(str,VER)),' ,'.join(map(str,PSR))) 
			print "Allelic correlation between the corresponding blocks (excluding variants with discordant dosages): {0}".format(' ,'.join(map(str,corr)))
			print "Proportion of loci omitted due to uncommonness, for each block pair: {0}".format(' ,'.join(map(str,list(round(x/(y+1e-10),3) for x,y in zip(Missing_Variants_Uncommon, Total_variants)))))
			print "Proportion of common loci omitted due to discordant dosages for each block pair: {0}".format(' ,'.join(map(str,list(round(x/(y+1e-10),3) for x,y in zip(Missing_Variants_Dosages,Total_common_variants)))))
			print "Total proportion of omitted loci for each block pair: {0}".format(' ,'.join(map(str,list(round(x/(y+1e-10),3) for x,y in zip(list(a+b for a,b in	zip(Missing_Variants_Dosages, Missing_Variants_Uncommon)), Total_variants)))))	
			print "Proportion of correctly estimated homologues within each haplotype block: {0}".format(' ,'.join(map(str,Correct_haplos)))
		else:
			print "Allelic correlation for common variants (including variants with wrong dosages):{0}\
\nPerfect Solution Rate for common variants (including variants with wrong dosages): {1}".format(round(sum(x*y for x,y in zip(Total_corr,Total_common_variants))/(1e-10+sum(z for z in Total_common_variants)),3),
				round(sum(x*nCr(y/2,2) for x,y in zip(Total_PSR,Total_common_variants))/(1e-10+sum(nCr(y/2,2) for y in Total_common_variants)),3))# Division bt 2 (y/2) as Total_common_variants is double the number of variant sites 
			avg_VER=float("Inf") if float("Inf") in VER else round(sum(x*y for x,y in zip(VER,Total_common_variants))/((1e-20+sum(z for z in Total_common_variants))*(1e-20+Total_var_true)*ploidy),3)
			print "Average Vector Error Rate per 'variant-homologue': {0} \nPerfect Solution Rate for common variants (excluding variants with wrong dosages): {1}".format(avg_VER, 
				round(sum(x*nCr((y-z)/2,2) for x,y,z in zip(PSR,Total_common_variants,Missing_Variants_Dosages))/(1e-10+sum(nCr((y-z)/2,2) for y,z in zip(Total_common_variants,Missing_Variants_Dosages))),3))# Division by 2 as Total_common_variants and Missing_Variants_Dosages are double the number of sites involved
			print "Allelic correlation for common variant (excluding variants with wrong dosages): {0}".format(round(sum(x*(y-z) for x,y,z in zip(corr,Total_common_variants,Missing_Variants_Dosages))/(1e-10+sum(y-z for y,z in zip(Total_common_variants,Missing_Variants_Dosages))),3))
			print "Proportion of loci missing in the estimation: {0}".format(round(Total_miss_uncommon_true/(Total_var_true+1e-10),3))
			print "Ratio of the spurious variants in the estimate to the original variants: {0}".format(round(Total_spurious_estimate/(Total_var_true+1e-10),3))
			print "Proportion of common loci omitted due to wrong dosage: {0}".format(round(sum(x for x in Missing_Variants_Dosages)/(2*(Total_var_true-Total_miss_uncommon_true)+1e-10),3))
			print "Positive predictive value of the estimated variants: {0}".format(round((Total_var_true-Total_miss_uncommon_true)/(1e-10+Total_var_true-Total_miss_uncommon_true+Total_spurious_estimate),3))
			print "Number of gaps introduced in the original block during estimation, per variant: {0}".format(round((n_blocks-1)/(Total_var_true+1e-10),5))			
			print "Total proportion of correctly estimated homologues averaged over the haplotype blocks: {0}".format(round(float(total_correcthaplos_count)/(n_blocks*ploidy+1e-10),3))
	except CustomException as err:
		print >> sys.stderr, err
		exit_msg=2
	except IOError as e:											
		print >> sys.stderr, "I/O error({0}): {1} {2}".format(e.errno, e.strerror, e.filename)
		exit_msg=2
	except:
		if exit_msg!=3:
			print >> sys.stderr, 'Unexpected error:'
			exc_type, exc_value, exc_traceback = sys.exc_info()
			traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
			traceback.print_exception(exc_type, exc_value, exc_traceback,
									  limit=2, file=sys.stderr)
			exit_msg=2
			raise
		else:
			pass
	finally:
		for _tmp in tmp:
			if os.path.isfile(_tmp):
				os.remove(_tmp)
		sys.exit(exit_msg)
