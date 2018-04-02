#!/usr/bin/env python
"""Converts the allele coding in SDhaP output, i.e. 1,2,3..., into the conventional VCF coding \
0,1,2,... In addition, variant id's are read from the corresponding vcf file and added as a column\
to the output. Written by Ehsan Motazedi, 19-08-2015, Wageningen UR"""
import sys
import getopt
import os 
import os.path
import tempfile
import re
from shutil import move

class CustomException(Exception):
	def __init__(self, value):
		self.parameter = value
	def __str__(self):
		return repr(self.parameter)

def getinput(argv):
	SDhaPphasefile = ''
	outputfile = ''
	vcffile=''
	try:
		opts, args = getopt.getopt(argv,"hp:o:v:",["help","phase=","output=","vcf="])
		if not (opts or args):
			opts=[('-h','')]
		elif args:
			raise getopt.GetoptError("Input arguments not in the proper format!", None)
		else:
			pass
	except getopt.GetoptError as err:
		print str(err)
		print 'Command not correct!Try -h, --help!'
		sys.exit(2)
	else: 
		for opt, arg in opts:
			if opt in ('-h','--help'):
				print 'Convert the 1,2,3,... allele coding in the SDhaP phased solution to 0,1,2,...'
				print "In addition, variant id's are read from the corresponding vcf file and added as a column to the output." 
				print 'Usage: ConvertAllelesSDhaP.py -p, --phase <SDhaP phased solution> -o, --output <output file> -v, --vcf <vcf file>'
				sys.exit()
			elif opt in ("-p", "--phase"):
				SDhaPphasefile = arg
			elif opt in ("-o", "--output"):
				outputfile = arg
			elif opt in ("-v", "--vcf"):
				vcffile = arg
		return (SDhaPphasefile, outputfile, vcffile)        
if __name__ == "__main__":   
	args=getinput(sys.argv[1:])
	try:
		outfile=None                                  # Output file stream
		infile=None                                   # Input file stream
		vcffile=None
		ploidy=0
		for i in args:
			if i=='':
				print 'Usage: ConvertAllelesSDhaP.py -p, --phase <SDhaP phased solution> -o, --output <output file> -v, --vcf <vcf file>'
				raise IOError("output and/or input and/or vcf file not specified!")
			else:
				pass  
		SDhaPphasefile=args[0]
		outputfile=args[1]
		vcffile=args[2] 
		if not os.path.isfile(SDhaPphasefile):
			try:
				raise IOError(2, 'No such file or directory',SDhaPphasefile)
			except IOError as err:
				if not err.args: 
					err.args=('',)
				err.args = err.args + ("Could not find SDhaP's original phasing solution!",)
				raise 
		elif not os.path.isfile(vcffile):
			try:
				raise IOError(2, 'No such file or directory',vcffile)
			except IOError as err:
				if not err.args: 
					err.args=('',)
				err.args = err.args + ("Could not find the vcf file!",)
				raise
		else:   											
			varpos_dic=dict()
			ploidy=-1
			with tempfile.NamedTemporaryFile(delete=False) as tmpvcf, open(vcffile, 'rU') as myvcf: # VCF file is first changed to be compatible with\
				for _line in myvcf:	                                    # extractHAIR, which has produced SDhaP fragment file. myvcf should have been from the VCF file given to extractHAIR
					if '#' not in _line:
						genotype = re.split('/|\|', _line.rstrip().split()[9].split(':')[0])
						if ploidy==-1 and '.' not in set(genotype): # The ploidy is first determined from the VCF file, as it is needed later to remove probable homozygous variants. 
							ploidy=len(genotype)
						if '.' not in set(genotype) and len(set(genotype))>1: # Homozygous variants are omitted
							tmpvcf.write(_line)
			if ploidy<2:
				raise CustomException("ERROR: Unexpected genotype detected in the VCF file or no called genotypes detected! Provide the correct format!")
			with open(tmpvcf.name,'rU') as tmpvcf2:                     # The variant position is linked to the variant number
				_i=0
				for _line in tmpvcf2:	
					_i+=1
					varpos_dic[str(_i)]=_line.rstrip().split()[1]              
			os.remove(tmpvcf.name)
			tmpvcf, _i = None, None
		with tempfile.NamedTemporaryFile(delete=False) as outfile, open(SDhaPphasefile, 'rU') as infile: # A temporary output is first made and the results are written to the output only \
			for line in infile:     # if everything goes smooth! This will protect the already existing files in case of error!                       
				if (line=="\n"):                     
					outfile.write(line)     # Skip blank lines between haplotype blocks              
				else:
					columns=line.rstrip().split()
					if (columns[0]=="Block"):
						outfile.write(line)    
					else:
						l_col=len(columns)
						if (l_col>=3):
							try:
								var_num=str(int(columns[0])+1)      # Get the variant number
								if ploidy>2:			# Do NOT convert alleles for diploids as the SDhaP diploid output is already 0,1,...	
									for col_num in range(1,l_col):      # Change allele coding 1,2,3,... to 0,1,2,... for polyploid output
										srr=str(int(columns[col_num])-1)  # Temporary string...
										if int(srr)<0:
											raise ValueError 
										else:
											columns[col_num]=srr
								columns=(var_num,)+(varpos_dic.get(var_num),)+tuple(columns[1:]) # Add the corresponding VCF position to the variant number and dosage
								outfile.write('\t'.join(columns)+'\n')
							except ValueError as err:
								if not err.args: 
									err.args=('',)
								err.args =err.args + ("Error in the index file! Allele coding should be 1,2,3...\
\nzero or invalid allele values observed!",)
								raise 
						else:
							raise CustomException("The SDhap phase file is corrupt!")
	except IOError as e:
		print >> sys.stderr, "I/O error({0}): {1} {2}".format(e.errno, e.strerror, e.filename)
		print >> sys.stderr, str(e.args[-1])
		if outfile:
			os.remove(outfile.name)
		sys.exit(2)
	except ValueError as e:
		print >> sys.stderr, "ERROR: "+str(type(e))+ '\n' + '\n'.join(e.args)
		if outfile:
			os.remove(outfile.name)
		sys.exit(2)
	except CustomException as instance:
		print >> sys.stderr, instance
		if outfile:
			os.remove(outfile.name) 
		sys.exit(2)        
	except:
		print >> sys.stderr, "Unexpected Error: %s" %sys.exc_info()[0]
		if outfile:
			os.remove(outfile.name)
		raise
	else:
		move(outfile.name, outputfile)
		print "Succesfully converted the alleles from 1,2,3... to 0,1,2... into {0}" \
			.format(outputfile) 
		sys.exit()
