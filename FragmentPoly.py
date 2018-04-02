#!/usr/bin/env python
"""Converts the HapCut input, i.e. the fargment matrix produced by extractHAIRS from .bam and .vcf files, 
into HapTree or SDhaP input formats. Written by Ehsan Motazedi, 19-08-2015, Wageningen UR """
import sys
import getopt
import re 
import os 
import os.path
import tempfile
from shutil import move
class CustomException(Exception):
	def __init__(self, value):
		self.parameter = value
	def __str__(self):
		return repr(self.parameter)
def getinput(argv):
	fragmentfile = ''
	outputfile = ''
	outformat=''
	try:
		opts, args = getopt.getopt(argv,"hf:o:x:",["help","fragment=","output=","outformat="])
		if args:
			raise getopt.GetoptError("Input arguments not in the proper format!", None)
		else:
			pass
	except getopt.GetoptError as err:
		print str(err)
		print 'Command not correct! Try -h, --help!' 
		sys.exit(2)
	else:    
		for opt, arg in opts:
			if opt in ('-h', '--help'):
				print 'FragmentPoly.py -f, --fragment <fragment> -o, --output <outputfile>, \
-x, --outformat <outputformat: SDhaP or HapTree>'
				sys.exit()
			elif opt in ("-f", "--fragment"):
				fragmentfile = arg
			elif opt in ("-o", "--output"):
				outputfile = arg
			elif opt in ("-x", "--outformat"):
				outformat = arg
		return (fragmentfile, outputfile,outformat)        
if __name__ == "__main__":									# <<The __main__ module starts here!>>
	args=getinput(sys.argv[1:])
	try:
		infile=None   										# Input stream for the fragment matrix
		outfile=None										# Output stream for HapTree/SDhaP fragment file
		writefile=None
		for i in args:
			if i=='':
				print 'All parameters must be specified! \n \
Usage: FragmentPoly.py -f, --fragment <fragment> -o, \
--output <outputfile>, -x, --outformat <outputformat> (SDhaP or HapTree)' 
				raise CustomException("output and/or input and/or format not specified!")
			else:
				pass  
		fragmentfile=args[0]
		outputfile=args[1]
		outformat=args[2]
		if not (outformat=='SDhaP' or outformat=='HapTree'):
			print 'Output format must be either SDhaP or HapTree!'
			raise CustomException("Output format not recognized!")
		else:
			pass
		num_reads=0
		num_snps=0                               			  								   
		with tempfile.NamedTemporaryFile(delete=False) as outfile, open(fragmentfile, 'rU') as infile: # A temporary output is first made and the results\
			for line in infile:	   												 					   # are written to the specified output only if everything goes smooth!
				columns=re.split('\s*|\t',line.rstrip())											   # This will protect the already existing files in case of error!
				l_col=len(columns)-1
				if ((l_col>=4) and (l_col%2==0)):
					reads_start=columns[2:l_col:2]
					hetro_sites=columns[3:l_col:2]
					hetro_dict={}
					for i in range(0,len(reads_start)):
						start_pos=int(reads_start[i])
						for j in range(0,len(hetro_sites[i])):
							hetro_dict[start_pos-1]=int(hetro_sites[i][j])
							start_pos+=1
					if outformat=='HapTree':
						if (len(hetro_dict)>1):
							num_reads+=1
							num_snps=max(num_snps,max(hetro_dict.keys()))
							outfile.write(str(hetro_dict))
							outfile.write("\n")
						else:
							pass
					else:                                               # if the format is SDhaP
						if (len(hetro_dict)>1):
							num_reads+=1
							num_snps=max(num_snps,max(hetro_dict.keys()))
							for col in range(3,l_col,2):                # change allele coding 0,1,2,... to 1,2,3,...
								srr=""                                  # temporary string to convert VCF coding to SDhap coding
								for allele in range(0,len(columns[col])):
									srr+=str(int(columns[col][allele])+1) 
								columns[col]=srr
							outfile.write(" ".join(columns))
							outfile.write("\n")
						else:
							pass
				else:
					raise CustomException("The fragment matrix is corrupt!")
		if os.path.getsize(fragmentfile)>0:
			num_snps+=1 
		else:
			pass        
		num_snps=str(num_snps)
		num_reads=str(num_reads)
		if outformat=='SDhaP':
			seqs=""
			with open(outfile.name, 'rU') as readfile, tempfile.NamedTemporaryFile(delete=False) as writefile:
				writefile.write(num_reads+"\n"+num_snps+"\n")
				for seqs in readfile:
					writefile.write(seqs)
		else:
			pass
	except IOError as e:
		print >> sys.stderr, "I/O error({0}): {1} {2}".format(e.errno, e.strerror, e.filename)
		if outfile:
			os.remove(outfile.name)
		if writefile:
			os.remove(outfile.name)
		sys.exit(2)
	except ValueError:
		print >> sys.stderr, "Error in the fragment matrix! Invalid values observed!"
		if outfile:
			os.remove(outfile.name)
		if writefile:
			os.remove(writefile.name)
		sys.exit(2)
	except CustomException as instance:
		print >> sys.stderr, instance
		if outfile:
			os.remove(outfile.name)
		if writefile:
			os.remove(writefile.name) 
		sys.exit(2)        
	except:
		print >> sys.stderr, "Unexpected Error: %s" %sys.exc_info()[0]
		if outfile:
			os.remove(outfile.name)
		if writefile:
			os.remove(writefile.name) 
		raise
	else:
		if outformat=='HapTree':
			move(outfile.name, outputfile)
		else:
			os.remove(outfile.name)
			move(writefile.name, outputfile)
		print "%s SNPs written from %s reads." % (num_snps, num_reads)
		print "Output successfully written to {0} fragment file.".format(outputfile)
		sys.exit()
