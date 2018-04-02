#!/usr/bin/env python
"""The simulation pipeline: from generating the (polyploid) genome by haplogenerator.py to producing raw NGS 
reads by ART and PBSIM, mapping the reads back to the reference, indexing the bam files, sorting the bam files 
(all by samtools), calling the variants by FreeBayes and finally estimating the haplotypes by HapTree, SDhaP and 
HapCompass. The output results summarize the quality of estimation using hapcompare.py. Iterations over the input
fasta, i.e. choosing different regions of the same length, and multiple random genomes, i.e. mutiple calls to 
haplogenerator for the same target region are also possible.
Written by Ehsan Motazedi, 10-12-2015, Wageningen UR
Last updated 02-04-2018"""
import argparse
import copy
import gc
import logging
import os 
import os.path
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import random as rnd
import resource
import signal
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import zipfile
from argparse import ArgumentParser, PARSER, REMAINDER, ArgumentError
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from collections import deque, namedtuple
from distutils.errors import DistutilsExecError
from functools import partial, wraps
from gettext import gettext as _, ngettext
from math import exp, log
from matplotlib.backends.backend_pdf import PdfPages
from re import sub as substitute
from resource import struct_rusage
from shutil import copyfile, rmtree, move
from tempfile import NamedTemporaryFile

""" Method that allows a free intermix of positional and optional argument strings
Its API is the same as for parse_known_args. 
obtained from http://bugs.python.org/file30204/test_intermixed.py """

# from warnings import warn
logging.basicConfig(format='%(levelname)s: Intermix: %(message)s')

def parse_intermixed_args(self, args=None, namespace=None):
	args, argv = self.parse_known_intermixed_args(args, namespace)
	if argv:
		msg = _('unrecognized arguments: %s')
		self.error(msg % ' '.join(argv))
	return args

def parse_known_intermixed_args(self, args=None, namespace=None, _fallback=None):
	# self - argparse parser
	# simplified from http://bugs.python.org/file30204/test_intermixed.py 
	# args, namespace - as used by parse_known_args
	# _fallback - action to take if it can't handle this parser's arguments (default raises an error)
	# returns a namespace and list of extras
	# positional can be freely intermixed with optionals
	# optionals are first parsed with all positional arguments deactivated
	# the 'extras' are then parsed
	# positionals 'deactivated' by setting nargs=0
	try:
		positionals = self._get_positional_actions()
		a = [action for action in positionals if action.nargs in [PARSER, REMAINDER]]
		if a:
			if _fallback is None:
				a = a[0]
				err = ArgumentError(a, 'parse_intermixed_args: positional arg with nargs=%s'%a.nargs)
				self.error(str(err))
			else:
				return _fallback(args, namespace)
		if [action.dest for group in self._mutually_exclusive_groups
			for action in group._group_actions if action in positionals]:
			if _fallback is None:
				self.error('parse_intermixed_args: positional in mutuallyExclusiveGroup')
			else:
				return _fallback(args, namespace)
		save_usage = self.usage
		if self.usage is None:
			self.usage = self.format_usage()[7:] # capture the full usage for use in error messages
		for action in positionals: 		 # deactivate positionals
			action.save_nargs = action.nargs # save nargs to be restored later
			action.nargs = 0		 # put nargs to zero to temporarily deactivate the positionals
		try:
			namespace, remaining_args = self.parse_known_args(args, namespace)
			for action in positionals: 		 # remove the empty positional values from namespace
				if hasattr(namespace, action.dest):
					delattr(namespace, action.dest)
		except SystemExit:
			warn('error from the optional input arguments!')
			raise
		finally:
			for action in positionals:		 # restore nargs and usage message before exiting
				action.nargs = action.save_nargs	
				# self.usage = save_usage
		logging.info('1st: %s,%s'%(namespace, remaining_args))
		optionals = self._get_optional_actions()		
		for action in optionals:			# perform the actions for the optionals
			action.save_required = action.required
			action.required = False
		for group in self._mutually_exclusive_groups:
			group.save_required = group.required
			group.required = False
		try:
			namespace, extras = self.parse_known_args(remaining_args, namespace) # parse positional arguments
		except SystemExit:
			warn('error from the positional input arguments!')
			raise
		finally:
			for action in optionals: 		# restore parser values before exiting
				action.required = action.save_required
			for group in self._mutually_exclusive_groups:
				group.required = group.save_required
	finally:
		self.usage = save_usage
	return namespace, extras

def warn(message):
	logging.warning(message)

ArgumentParser.parse_intermixed_args = parse_intermixed_args
ArgumentParser.parse_known_intermixed_args = parse_known_intermixed_args

def check_and_store():
	""" Checks if the sequencing method and library type are correctly specified. In case the values are switched for these two, they will be put back\
 at their proper places."""
	class Check_And_Store(argparse.Action):
		def __call__(self, parser, args, value, option_string=None):
			S = value
			T = getattr(args, 'pe')
			if T and T not in ['SE','PE','MP','GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS']:
					raise argparse.ArgumentError(self, "Unallowed sequencing method or library type '{}'!".format(T))
			if S and S not in ['SE','PE','MP','GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS']:
					raise argparse.ArgumentError(self, "Unallowed sequencing method or library type '{}'!".format(S))
			if (S in ['SE','PE','MP']) and (T in ['GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS']):
				setattr(args, self.dest, T)
				setattr(args, 'pe', S)
			elif (S and S not in ['SE','PE','MP']) and (T in ['GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS']):
				raise argparse.ArgumentError(self, "Only 'SE', 'PE' and 'MP' are allowed for the library type! '{}' is not allowed!".format(S))
			elif (S in ['SE','PE','MP']) and (T and T not in ['GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS']):
				raise argparse.ArgumentError(self, "Only 'GA1','GA2','HS10','HS20','HS25','MS','CLR' and 'CCS' are allowed for the sequencing method! '{}' is not allowed!".format(T))
			elif (not S) and (T in ['GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS']):
				setattr(args, self.dest, T)
				setattr(args, 'pe', '')
			elif (not T) and (S in ['SE','PE','MP']):
				setattr(args, self.dest, '')
				setattr(args, 'pe', S)
			else:
				setattr(args, self.dest, S)
	return Check_And_Store

def check_file():
	""" Checks if the given file path exists."""
	class Check_File(argparse.Action):
		def __call__(self, parser, args, values, option_string=None):
			if values is None:
				setattr(args, self.dest, '')
			elif False in [os.path.exists(_val) for _val in values]:
				raise argparse.ArgumentError(self, "no file called '{}' was found!".format(_val))
			else:
				setattr(args, self.dest, values)
	return Check_File

def check_freq(value):
	""" Checks if the passed argument is a valid frequency, i.e. real number between zero and 1."""
	rvalue = None
	rvalue = value if value in ['PE','SE','MP', 'GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS'] else rvalue
	if rvalue:
		return rvalue
	try:
		rvalue = round(float(value),5)
	except ValueError:
		raise argparse.ArgumentTypeError("%s is not a real number!" % value)
	except TypeError:
		raise argparse.ArgumentTypeError("%s is not a real number!" % value)
	else:
		if rvalue<0 or rvalue>1:
			raise argparse.ArgumentTypeError("%s is not between zero and 1!" % value)
		else:
			pass
	return rvalue

def check_pid(pid):
	""" Checks if a unix process exists or not."""
	try:
		os.kill(pid,0)
	except OSError:
		return False
	else:
		return True

def check_positive(value):
	""" Checks if the passed argument is a positive integer."""
	ivalue = None
	ivalue = value if value in ['PE','SE','MP','GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS'] else ivalue
	if ivalue:
		return ivalue
	try:
		ivalue = int(float(value))
	except ValueError:
		raise argparse.ArgumentTypeError("%s is not a positive integer value!" % value)
	except TypeError:
		raise argparse.ArgumentTypeError("%s is not a positive integer value!" % value)
	else:
		if ivalue <= 0:
			raise argparse.ArgumentTypeError("%s is not a positive integer value!" % value)
		else:
			pass
	return ivalue

def check_positive_float(value):
	""" Checks if the passed argument is a positive real number."""
	fvalue=None
	fvalue = value if value in ['PE','SE','MP', 'GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS'] else fvalue
	if fvalue:
		return fvalue
	try:
		fvalue = float(value)
	except ValueError:
		raise argparse.ArgumentTypeError("%s is not a positive real value!" % value)
	except TypeError:
		raise argparse.ArgumentTypeError("%s is not a positive real value!" % value)
	else:
		if fvalue <= 0:
			raise argparse.ArgumentTypeError("%s is not a positive real value!" % value)
		else:
			pass
	return fvalue

def clean_up_files(file_lst):
	""" Cleans up the files in its input list"""
	for _file in file_lst:
		if os.path.exists(_file):
			os.remove(_file)

def clean_up_folders(dir_lst):
	""" Cleans up the directories in its input list"""
	for _dir in dir_lst:
		if os.path.isdir(_dir):
			rmtree(_dir)

class Conditional_decorator(object):
	""" Generate a decorator function class object with a conditional __call__ method."""
	timem = False
	@classmethod
	def set_timem(cls, timem):
		""" Set the value of the static variable timem to apply the decorator to __call__'s argument if timem is true"""
		cls.timem = timem

	@classmethod
	def get_timem(cls):
		""" Get the value of the static variable timem"""
		return cls.timem

	def __init__(self, dec):
		self.decorator = dec

	def __call__(self, func):
		if not self.get_timem(): # Return the function unchanged, i.e. not decorated.
			return func
		return self.decorator(func)

class Command(Conditional_decorator):
	""" Handles the spawning of UNIX command line processes."""
	
	verbose = False # Verbose is just here to make debugging easy. It is not used in normal calls to run.

	@classmethod
	def set_verbose(cls, verbose=False):
		""" Set the value of the static variable verbose """
		cls.verbose = verbose

	@classmethod
	def get_verbose(cls):
		""" Get the value of the static variable verbose """
		return cls.verbose

	def __init__(self, cmd, shell=True, outfile=None, keep_error_file=False):	
		#self.cmd ="ulimit -v 2097152;"+cmd # Set a memory limit of 2GBYTES on the spawned process
		self.cmd = cmd
		self.keep_error_file = keep_error_file
		self.outfile = outfile
		self.process = None
		self.std_out = None
		self.std_err = None
		self.shell = shell
		self.mem = None
		self.time_elapsed = None

	def _fn_timer(function): 
		@wraps(function)
		def function_timer(self, *args, **kwargs): 
			""" Return the wall-clock time spent during a call to the fucntion """
			t0 = os.times()
			try:
				result = function(self, *args, **kwargs)
			except:
				raise
			finally:
				t1 = os.times()
				self.time_elapsed = round(sum(y-x for y, x in zip(t1[:-1], t0[:-1])), 4)
			return result
		return function_timer

	@Conditional_decorator(_fn_timer)
	def run(self, timeout): 
		""" Spawn the subprocess by Popen method and wait for it to finish until at most the timeout. Report the total 
 execution time, as well as the maximum resident set size used by the subprocess. """
		verbose = self.get_verbose()
		
		def reader(iostream, outlist, chunk=64):
			"""function to read the I/O byte stream and append it to outlist in blocks of size chunk.
			   The I/O pipe is closed when read."""
			read_block = partial(iostream.read, chunk) #at most 64 bytes are read and stored in each block
			for _block in iter(read_block, b''):
				outlist.append(_block)
			iostream.close()
		
		def target():
			try:
				sys.stdout.flush()
				sys.stderr.flush()
				self.std_err, self.std_out = [], []
				self.mem = struct_rusage(['0']*16) # a struct rusage object with all zero values
				options = os.WUNTRACED | os.WCONTINUED
				self.process = subprocess.Popen(self.cmd, shell=self.shell, close_fds=True,
								preexec_fn = preexec, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				reader1 = threading.Thread(target=reader, args=(self.process.stdout, self.std_out))
				reader2 = threading.Thread(target=reader, args=(self.process.stderr, self.std_err))
				for _reader in (reader1, reader2):
					_reader.start()
				for _reader in (reader1, reader2):
					_reader.join()
				while 1:
					_pid, _status, self.mem = os.wait4(self.process.pid, options)
					if os.WIFSTOPPED(_status): # The child has been stopped
						continue
					elif os.WIFSIGNALED(_status):
						if os.WCOREDUMP(_status):
							raise DistutilsExecError("ERROR: segmentation fault (core dumped) by '{0}'. Check the core file!".format(self.cmd if self.shell==True else ' '.join(map(str, self.cmd))))
						raise DistutilsExecError("ERROR: '{0}' terminated by signal '{1}'".format(self.cmd if self.shell==True else ' '.join(map(str,self.cmd)), s.WTERMSIG(_status)))
					elif os.WIFEXITED(_status):
						self.process.returncode = os.WEXITSTATUS(_status)
						if self.process.returncode == 0:
							break   # The child terminated successfully
						else:
							raise DistutilsExecError("ERROR: '{0}' failed with exit status '{1}'".format(self.cmd if self.shell==True else ' '.join(map(str,self.cmd)), self.process.returncode))
					else:
						raise DistutilsExecError("ERROR: unknown error while executing '{0}': termination status {1}".format(self.cmd if self.shell==True else ' '.join(map(str,self.cmd)), _status))
			except DistutilsExecError as e:
				self.process.returncode = -24 # Arbitrary number for DistutilsExecError
				self.std_err.append(str.encode('\n'+str(e), 'utf-8'))
			except OSError as e:
				self.process.returncode = -64 # Arbitrary number for OSError
				self.std_err.append(str.encode("\nOSError ({0}): {1}".format(e.errno, e.strerror), 'utf-8'))
			except:
				self.process.returncode = -85 # Arbitrary number for unknown exception	
				self.std_err.append(str.encode("\nUnknown error caused forking to fail for '{0}'!".format(self.cmd if self.shell==True else ' '.join(map(str,self.cmd))), 'utf-8'))
			finally:
				if self.process and self.process.returncode!=0:
					if check_pid(self.process.pid):
						os.killpg(self.process.pid, signal.SIGTERM)
				self.std_err = b''.join(self.std_err).decode('utf-8')
				self.std_out = b''.join(self.std_out).decode('utf-8')
				if self.outfile is not None:
					if not self.keep_error_file:
						with open(self.outfile,'w') as fileobject:
							pass # delete the contents of the error file if it already exists
					with open(self.outfile,'a') as fileobject:
						fileobject.write(self.std_out+'\n*************\n')
						fileobject.write(self.std_err+'\n*************\n')
				elif verbose:
					print >> sys.stderr, self.std_err
					print >> sys.stdout, self.std_out
				else:
					pass
				if self.process and ("HapTree_multiblock.py" in self.cmd) and self.mem != struct_rusage(['0']*16):
					_spawned_mem=struct_rusage([0]*16)
					_get_mem = False
					for _line in self.std_out.split(chr(10)):
						if _get_mem:
							dummy=_line.strip().split('\t')
							if dummy==['ru_'+_x for _x in 'utime,stime,maxrss,ixrss,idrss,isrss,minflt,majflt,\
nswap,inblock,oublock,msgsnd,msgrcv,nsignals,nvcsw,nivcsw'.split(',')]:
								pass
							elif dummy!=['']:
								_spawned_mem = struct_rusage([sum(map(float, _x)) if _fld!=2 else max(map(float, _x)) for _fld, _x in enumerate(zip(_spawned_mem[:], dummy[:]))])
						if _line.strip()=='RESOURCE USAGE REPORT:':
							_get_mem=True
					self.mem=struct_rusage(self.mem[0:2]+(max(self.mem[2], _spawned_mem[2]),)+
									tuple(_x+_y for _x,_y in zip(_spawned_mem[3:], self.mem[3:])))
		try:
			thread = threading.Thread(target=target)
			thread.start()
			thread.join(timeout)
		except:
			thread.join()
			raise subprocess.CalledProcessError(1, "Unknown error occured while threading for '{0}'!".format(self.cmd if self.shell==True else ' '.join(map(str,self.cmd))), None)
		
		if thread.is_alive():
			try:
				os.killpg(self.process.pid, signal.SIGTERM)
			except OSError as e:
				if e.errno==12:
					raise subprocess.CalledProcessError(1, "OSError ({0}): {1} occured for '{2}'".format(e.errorno, e.strerror, self.cmd if self.shell==True else ' '.join(map(str,self.cmd))), None)
				else:	
					pass # The process finished between the 'is_alive()' and 'terminate()'
			except AttributeError as e:
				pass # The process finished between the 'is_alive()' and 'terminate()'
			else:
				print >> sys.stderr, 'Process #{0:<5d} was not responding! Timeout termination occured after {1:<5.3f} seconds!'.format(self.process.pid, round(float(timeout), 3))
				self.process.returncode=-9 # Arbitrary number indicating timeout error
			finally:
				thread.join() # Be sure the thread is terminated
		if (not self.process) or self.process.returncode != 0:
			raise subprocess.CalledProcessError(1, self.cmd if self.shell==True else ' '.join(map(str,self.cmd)), None)
		elif self.process:
			return(self.process.returncode)
			
class CostumeFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
	"""Help message formatter based on ArgumentDefaultsHelpFormatter from argparse standard module"""
	def _get_help_string(self, action):
		help = action.help
		return help

class CustomException(Exception):
	def __init__(self, value):
		self.parameter = value
	def __str__(self):
		return self.parameter
	
def fixed_snps(seq):
	""" Determine the fixed SNP positions for its input seq"""
	if len(seq)<=2:
		Mu=[1,len(seq)]
	elif len(seq)>2 and len(seq)<=1000: # Fixed mutation sites at the start, end and middle of the reference
		Mu=[1, (len(seq)+1)/2,len(seq)]
	elif len(seq)>1000 and len(seq)<=10000:
		Mu=range(1,len(seq)-500)[::500]+[len(seq)]
	elif len(seq)>int(1e4) and len(seq)<=int(1e5):
		Mu=range(1,9500)[::500]+range(int(1e4),len(seq)-int(1e4))[::int(1e4)]+[len(seq)]
	elif len(seq)>int(1e5) and len(seq)<=int(1e6):
		Mu=(range(1,9500)[::500]+range(int(1e4),int(1e5)-int(1e4))[::int(1e4)]+
			range(int(1e5),len(seq)-int(1e5))[::int(1e5)]+[len(seq)])
	else:
		Mu=(range(1,9500)[::500]+range(int(1e4),int(1e5)-int(1e4))[::int(1e4)]+
			range(int(1e5),int(1e6)-int(1e5))[::int(1e5)]+range(int(1e6),len(seq)-int(1e6))[::int(1e6)]+
			[len(seq)])
	return Mu

def make_pophap(ref, size, output, emp_filename, lognormal_dist, lognormal_sd, freq, timeout, verbose=False):
	"""makes a number of haplotypes from which the simulated haplotypes will be selected by sampling with replacement.""" 
	with NamedTemporaryFile(delete=False) as tmphap, NamedTemporaryFile(delete=False) as subprocess_err:
		if verbose:
			print >>sys.stdout, "The whole population haplotype set is being generated with the given size {0}...".format(size)
	try: 
		zip_pophap=zipfile.ZipFile(output,'w')
		rnd_seed=np.random.uniform(low=9000, high=95000 , size=(size,))
		for _pophap in range(0, size):
			if emp_filename:
				cmd=Command("haplogenerator.py -f "+ref+" -o "+tmphap.name+" -p 1"+
					' -s "'+str([1,0,0])+'"'+" -m "+'"'+"{'A':'C','C':'A','G':'T','T':'G'}"+'"'+ 
					' --empirical '+str(emp_filename)+' --genomic --rndseed '+str(rnd_seed[_pophap]), outfile=subprocess_err.name)  # Generate the artificial genome with the lognormal model
			elif lognormal_dist[0]>0:
				cmd=Command("haplogenerator.py -f "+ref+" -o "+tmphap.name+" -p 1"+
					' -s "'+str(lognormal_dist)+'"'+" -m "+'"'+"{'A':'C','C':'A','G':'T','T':'G'}"+'"'+ 
					' --model lognormal --sdlog "'+str(lognormal_sd)+'" --genomic --rndseed '+str(rnd_seed[_pophap]), outfile=subprocess_err.name)  # Generate the artificial genome with the lognormal model
			else:
				cmd=Command("haplogenerator.py -f "+ref+" -o "+tmphap.name+" -p 1"+
					' -s "'+str(freq)+'"'+" -m "+'"'+"{'A':'C','C':'A','G':'T','T':'G'}"+'"'+ 
					" --model poisson --genomic --rndseed "+str(rnd_seed[_pophap]), outfile=subprocess_err.name)  # Generate the artificial genome with the (default) Poisson stochastic model
			genome=tmphap.name+"_hap1.fa"
			genome_indx=tmphap.name+"_hap1.fa.fai"
			genomehap=tmphap.name+"_varianthaplos.txt"
			haplogen_err_msg=''
			try:
				exit_msg=cmd.run(timeout=timeout)
			except subprocess.CalledProcessError:
				with open(subprocess_err.name, 'rU') as subprocess_err:
					subprocess_err_msg=[_x if ('*' not in _x and _x!='\n') else '' for _x in subprocess_err.readlines()]
					subprocess_err_msg='\n'+'Original error message:'+'\t'+''.join(''.join(str(_y) for _y in _x) for _x in subprocess_err_msg)+'\n'+"-------------------------------"+'\n'
				raise CustomException("ERROR: Failed to make the artificial genome! Check the help for haplogenerator!"+subprocess_err_msg)
			else:
				for _file in [genome, genome_indx, genomehap]:
					zip_pophap.write(_file, arcname=substitute(tmphap.name, 'hap'+str(_pophap), _file))
	except:
		print >>sys.stdout, "ERROR: failed to construct the whole population haplotype set! Check the help for haplogenerator!"
		raise
	else:
		if verbose:
			print >>sys.stdout, "The whole population haplotype set was successfully constructed!"
	finally:
		for _file in [genome, genome_indx, genomehap, subprocess_err.name]:
			if os.path.isfile(_file):
				os.remove(_file)
		if not zip_pophap.testzip():
			archive=zip_pophap.filename
			zip_pophap.close()
			return(archive)
		return(None)		

def preexec():
	signal.signal(signal.SIGPIPE, signal.SIG_DFL)
	os.setsid()

def required_length(nmin, ex1, ex2, nmax=-1):
	""" Checks if the number of passed arguments falls within the valid range."""
	class RequiredLength(argparse.Action):
		def __call__(self, parser, args, values, option_string=None):
			if set(['GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS']) & set(values):
				S=list(set(['GA1','GA2','HS10','HS20','HS25','MS','CLR','CCS']) & set(values))
				values=[_x for _x in values if _x not in S]
				ex1.extend(S)
			if set(['PE','SE','MP']) & set(values):
				S=list(set(['PE','SE','MP']) & set(values))
				values=[_x for _x in values if _x not in S]
				ex2.extend(S)
			if not (nmin<=len(values)<=nmax or (nmin<=len(values) and nmax==-1)): # nmax=-1 means that there is no upper limit on the required length
				msg='argument "{f}" requires between {nmin} and {nmax} arguments! {num} specified.'.format(
					f=self.dest,nmin=nmin,nmax=nmax, num=len(values))
				raise argparse.ArgumentError(self, msg)
			setattr(args, self.dest, values)
	return RequiredLength

def write_list(lst, openfile):
	""" Maps its input list to the columns of its input file object.""" 
	for _x in lst[:-1]:
		openfile.write(_x+'\t')
	openfile.write(lst[-1]+'\n')
	return None

if __name__ == "__main__":					    # <<The __main__ module starts here!>>
	start_time_cpu=os.times()
	gc.enable()
	gc.set_threshold(20,10,10)
	try:	   	   	 
		contig_found=False
		dontrun=None # do not run SDhaP or HapTree if their input format could not be provided
		estimate=('HapCompass', 'HapTree', 'SDhaP')     # The name of the three haplotyping algorithms, in the order they are applied
		exit_msg=0
		_extra1=[]
		_extra2=[]
		z=None
		zip_arc=''
		report_aln=None
		pophapname=None
		pophap=None
		ARTSPath=os.getenv("ARTSPath")#1
		blasrPath=os.getenv("blasrPath")#2
		bowtiePath=os.getenv("bowtiePath")#3
		bwaPath=os.getenv("bwaPath")#4
		FreebayesPath=os.getenv("FreebayesPath")#5
		HapCompassPath=os.getenv("HapCompassPath")#6
		HapCutPath=os.getenv("HapCutPath")#7
		HapTreePath=os.getenv("HapTreePath")#8
		HOME=os.getenv("HOME").rstrip('/')+'/'
		PBSIMPath=os.getenv("PBSIMPath")#9
		PicardtoolsPath=os.getenv("PicardtoolsPath")#10
		samtoolsPath=os.getenv('samtoolsPath')#11
		SDhaPPath=os.getenv("SDhaPPath")#12
		if not ARTSPath:#1
			ARTSPath=HOME+"art_bin_ChocolateCherryCake"
			if os.path.isdir(ARTSPath):
				ARTSPath+='/'
			else:
				ARTSPath=''
			os.environ['ARTSPath']=ARTSPath.rstrip('/')
		else:
			ARTSPath=ARTSPath.rstrip('/')+'/'
		if not blasrPath:#2
			blasrPath=HOME+"blasr"
			if os.path.isdir(blasrPath):
				blasrPath+='/'
			else:
				blasrPath=''
			os.environ['blasrPath']=blasrPath.rstrip('/')
		else:
			bowtiePath=bowtiePath.rstrip('/')+'/'
		if not bowtiePath:#3
			bowtiePath=HOME+"bowtie2-2.2.4"
			if os.path.isdir(bowtiePath):
				bowtiePath+='/'
			else:
				bowtiePath=''
			os.environ['bowtiePath']=bowtiePath.rstrip('/')
		else:
			bowtiePath=bowtiePath.rstrip('/')+'/'	
		if not bwaPath:#4
			bwaPath=HOME+"bwa"
			if os.path.isdir(bwaPath):
				bwaPath+='/'
			else:
				bwaPath=''
			os.environ['bwaPath']=bwaPath.rstrip('/')
		else:
			bwaPath=bwaPath.rstrip('/')+'/'
		if not FreebayesPath:#5
			FreebayesPath=HOME+"freebayes/bin"
			if os.path.isdir(FreebayesPath):
				FreebayesPath+='/'
			else:
				FreebayesPath=''
			os.environ['FreebayesPath']=FreebayesPath.rstrip('/')
		else:
			FreebayesPath=FreebayesPath.rstrip('/')+'/'
		if not HapCompassPath:#6
			HapCompassPath=HOME+"hapcompass_v0.8.2"
			if os.path.isdir(HapCompassPath):
				HapCompassPath+='/'
			else:
				HapCompassPath=''
			os.environ['HapCompassPath']=HapCompassPath.rstrip('/')
		else:
			HapCompassPath=HapCompassPath.rstrip('/')+'/'
		if not HapCutPath:#7
			HapCutPath=HOME+"hapcut"
			if os.path.isdir(HapCutPath):
				HapCutPath+='/'
			else:
				HapCutPath=''
			os.environ['HapCutPath']=HapCutPath.rstrip('/')
		else:
			HapCutPath=HapCutPath.rstrip('/')+'/'
		if not HapTreePath:#8
			HapTreePath=HOME+"HapTree_v0"
			if os.path.isdir(HapTreePath):
				HapTreePath+='/'
			else:
				HapTreePath=''
			os.environ['HapTreePath']=HapTreePath.rstrip('/')
		else:
			HapTreePath=HapTreePath.rstrip('/')+'/'
		if not PBSIMPath:#9
			PBSIMPath=HOME+"pbsim-1.0.3"
			if os.path.isdir(PBSIMPath):
				PBSIMPath+='/'
			else:
				PBSIMPath=''
			os.environ['PBSIMPath']=PBSIMPath.rstrip('/')
		else:
			PBSIMPath=PBSIMPath.rstrip('/')+'/'
		if not PicardtoolsPath:#10
			PicardtoolsPath=HOME+"picard-tools-2.9.0"
			if os.path.isdir(PicardtoolsPath):
				PicardtoolsPath+='/'
			else:
				PicardtoolsPath=''
			os.environ['PicardtoolsPath']=PicardtoolsPath.rstrip('/')
		else:
			PicardtoolsPath=PicardtoolsPath.rstrip('/')+'/'
		if not samtoolsPath:#11
			samtoolsPath=HOME+"samtools-1.4"
			if os.path.isdir(samtoolsPath):
				samtoolsPath+='/'
			else:
				samtoolsPath=''
			os.environ['samtoolsPath']=samtoolsPath.rstrip('/')
		else:
			samtoolsPath=samtoolsPath.rstrip('/')+'/'
		if not SDhaPPath:#12
			SDhaPPath=HOME+"SDhaP/SDhaP_poly"
			if os.path.isdir(SDhaPPath):
				SDhaPPath+='/'
			else:
				SDhaPPath=''
			os.environ['SDhaPPath']=SDhaPPath.rstrip('/')
		else:
			SDhaPPath=SDhaPPath.rstrip('/')+'/'
		HapCompass_HapTree_SDhaP_errors=[0,0,0]				# Number of times that each haplotyping algorithm fails
		HapCompass_HapTree_SDhaP_hapcompare_errors=[0,0,0]	# Number of times hapcompare failed to evaluate each haplotyping algorithm
		method_report=[]									# The final output files, having 'output' as base name and _HapCompass, _HapTree and _SDhaP as suffix
		observed_dosage=[]									# The number of observed simplex, duplex, ..., ploidy-plex alleles
		pdffile=None
		pdfpages=None
		plotting=False 									# Flag to plot variant density histograms
		target_record = SeqIO.SeqRecord(seq='', id= '')
		target_str=''
		tmpdirs=[]										# List of the temporary folders needed
		tmpfiles=[]										# List of the temporary files needed
		tmp_ref=None
		timem=False 									# Flag to report time and memory use of haplotyping algorithms in a separate dat file  
		verbose=False 									# Flag to print warnings and messages
		errors_in_haplotyping=None
		parser = argparse.ArgumentParser(description="The simulation pipeline: from generating the (polyploid) genome by haplogenerator.py to producing NGS\
 reads by ART and PBSIM, mapping the reads back to the reference, indexing the bam files, sorting the bam\
 files (all by samtools), calling the variants by FreeBayes and finally estimating the haplotypes by HapTree,\
 SDhaP and HapCompass. The output results summarize the quality of estimation using hapcompare.py.", 
		formatter_class=CostumeFormatter, epilog='WARNING: If uncalled bases, specified with "N", exist in the input fasta file,\
 haplosim \ndiscards them should they occur in its chosen reference. However, if more than \n20% of the chosen reference happens\
 to be "N", haplosim throws an exception and stops \n(if a fixed start point is given) OR tries to choose another reference from a\
 new \nrandom start point on the input fasta sequence (if the random mode is used). In \nthe latter case, haplosim performs up to 15 unsuccessful\
 attempts before throwing\nan exception.')
		parser.add_argument('fastafile', metavar='fasta file', action=check_file(), type=str, nargs=1, help='the fasta file to select a targeted reference region from.')
		parser.add_argument('output', metavar='output file', action='store', type=str, nargs=1, help='the output files base name to report comparison statistics and log(s).')
		parser.add_argument('pe', metavar='Single, paired-end or mate-pair', action='store', nargs='?', type=str, 
			help='Choose SE, PE or MP sequencing protocol for the Illumina sequencers (default = PE). Ignored when the sequencer is PacBio ("CLR" or "CCS").',
			default='')
		parser.add_argument('method', metavar='Sequencing method', action=check_and_store(), type=str, nargs='?', help='the sequencing method (technology) used to generate the reads (default = HS25). Available:\n\
	1) Illumina: GA1, GA2, HS10, HS20, HS25, MS. Single read lengths are fixed with these technologies to 35bp, 55bp, 75bp, 100bp, 125bp and 250bp, respectively. Paired read length, i.e. insert-size, can be specified using --insert, --insert-size option.\n\
	2) PacBio: CLR, CCS. Read lengths for CLR and CCS can be specified using --insert, --insert-size option.', default='')
		parser.add_argument('-p','--ploidy', dest='ploidy', action='store', type=check_positive, nargs=1, metavar="int+",
			help='the ploidy level of the artificial genome (default = 2)',default='2')
		parser.add_argument('--insert', '--insert-size', dest='insert', action='store', metavar="int+",
			type=check_positive, nargs=1, help='the insert-size in base pairs used with Illumina PE technologies\
 (default chosen from (70,100,150,250,300,550) for (GA1,GA2,HS10,HS20,HS25,MS), respectively). Determines the read length instead with "CLR" and "CCS" sequencing methods, as these methods do not use paired libraries (default 1000 and 5000 for CCS and CLR, respectively).',
			default=0)
		parser.add_argument('-c','--contig', dest='contig_name', action='store', type=str, nargs=1, metavar="str",
			help='the name of the contig in the fasta file from which the reference is selected\
 (default: the first contig in the fasta file)', default='')
		parser.add_argument('-v', '--verbose', action='store_true',
		   help='flag to print warning and progress messages. (default = False)',
		   default=False)
		parser.add_argument('--remove-duplicate', action='store_true',
		   help='if set, exact duplicate reads will be removed from the alignment using Picardtools, before the indexing and variant calling\
		   steps.', default=False)
		parser.add_argument('--bwa-mem', action='store_true',
		   help='if set, bwa-mem is used to align the reads instead of the default blasr for PacBio and Bowtie2 for Illumina.',
		   default=False)
		parser.add_argument('--report-aln', action='store_true',
		   help='if set, the reference, sam, bam, vcf and original read files are reported as my{ref,sam,bam,vcf,reads}xy.{fa,sam,bam,vcf,fastq} files\
 in a single zip archive, where x is the region id and y\
 is the iteration number. All of the generated and estimated haplotypes will also be written to the log file.',
		   default=False)
		parser.add_argument('--plothist','--hist','--plot-hist', dest='plot', action='store_true',
		   help='if set, the histograms of the distances between the consecutive variants will be saved in output.pdf file for all of the\
 simulated genomes.', default=False)
		parser.add_argument('--timem-use','--timem','-timem-use', '-timem', action='store_true',
			help='option to report the total CPU time and the maximum physical memory consumed by each haplotyping\
 algorithm: HapTree, SDhaP and HapCompass, in a separate dat file called output_timem.dat.', default=False)
		parser.add_argument('--no-missing', action='store_true',
		   help='if set, the output of a haplotyping algorithm for a generated genome will be evaluated, even if the other haplotyping algorithms fail on the same genome.\
 Otherwise, the haplotype estimates are evaluated only for the genomes that are successfully phased by all the three haplotyping algorithms.',
		   default=False)
		parser.add_argument('-l','--length', dest='ref_length', action='store', type=check_positive, nargs=1, metavar="int+",
		   help='the length of the reference. (default: the length of the contig from the start position)',
		   default=0)
		process = parser.add_mutually_exclusive_group()
		process.add_argument('-m','--mufreq', dest='freq', action='store', type=check_freq, nargs=1, metavar="freq",
			help='mutation frequency passed to haplogenerator, used in the mixed poisson process to\
 generate genomic variants (default = 0.01). NOT to be used with --lognormal, --log-normal and --empirical!',
			default='0.01')
		process.add_argument('--lognormal', '--log-normal',metavar="REAL+",
		   help='if specified, a lognormal stochastic process will be used to generate mutations by haplogenerator.\
 Otherwise a mixed poisson process is applied. Values of the mean and of the standard deviation of the log\
 distance between the variants (in bp) could be passed to this argument (default = 1, 1). NOT to be used with\
 -m, --mufreq and --empirical!', nargs='*', type=check_positive_float, action=required_length(0, _extra1, _extra2, 2), default=[-1,-1])
		process.add_argument('--empirical', action=check_file(), nargs=1, metavar="str", type=str,
			help='the name of the file containing the distances between neighbor SNPs, from which the simulated distances will be sampled with replacement.\
 NOT to be used with --lognormal, --log-normal, -m  and --mufreq!')
		group_dosage = parser.add_mutually_exclusive_group()
		group_dosage.add_argument('--dosages', '--dosage', dest='dosages', metavar="freq",
		   help='the proportion of simplex, duplex, ..., ploidy-plex alleles in the simulated genomes, specified respectively. Must\
 be a valid number between 0 and 1 (default = 1/p for every type).\
 INCOMPATIBLE with -maxhap, --maxhap!', nargs='*', type=check_freq, action=required_length(0, _extra1, _extra2), 
		default=[-1,-1])
		parser.add_argument('-t','--timeout', action='store', type=check_positive_float, nargs=1, metavar="Real+",
		   help='the timeout limit in seconds for individual processes, i.e. samtools, haplogenerator, HapCompass, etc.\
(default: 500)', default='500')
		parser.add_argument('-x','--coverage', dest='coverage', action='store', type=check_positive, nargs=1, metavar="int+",
		   help='sequencing coverage per homologue. (default = 10)', default='10')
		group = parser.add_mutually_exclusive_group()
		group.add_argument('-s','--start', metavar='int', dest='start_pos', type=check_positive, nargs='*', action=required_length(0, _extra1, _extra2, 2),
		   help="the start position for the reference upon the contig (1'st value) and the number of random mutations on\
 it (default = 1 1). NOT to be combined with -r, --random!", default=[1,1])
		group.add_argument('-r', '--random', metavar='int', dest='random_str_mu', type=check_positive, nargs='*', action=required_length(0, _extra1, _extra2, 2),
			help="set the number of random start points (1'st value) and the number of random mutations for each of these\
 regions (2'nd value) (default = 1 1). NOT to be combined with -s, --start, which is the default mode!.", default=[0,1])
		group_dosage.add_argument('-maxhap', '--maxhap', metavar='int', dest='maxhap', type=check_positive, nargs=1, action='store',
			help="maximum number of haplotypes in the simulated population. This maximum haplotype set is first simulated by the\
 chosen model, and the haplotypes of each genome are selected from this set by sampling with replacement using random\
 weights for each haplotype (default to no maximum). INCOMPATIBLE with --dosages, --dosage!", default=0)
		try:
			if not len(sys.argv)>1:
				parser.parse_args('-h'.split())
			elif ('-h' in sys.argv) or ('--help' in sys.argv):
				parser.parse_args('-h'.split())
			else:
				args = vars(parser.parse_intermixed_args()) 
		except SystemExit:
			exit_msg=3
			raise
		fastafile=args['fastafile'][0] 		# The fasta file to extract the references from it
		ploidy=args['ploidy'] if isinstance(args['ploidy'], int) else args['ploidy'][0] 					# The ploidy level, passed to haplogenerator
		insert=args['insert'] if isinstance(args['insert'], int) else args['insert'][0] 					# The insert-size for PE reads, i.e. the distance in bp between the two adaptors, or the read-length for PacBio CCS and CLR reads.
		output=args['output'] if isinstance(args['output'], str) else args['output'][0] 		    	    # The name of the output file
		contig_name=args['contig_name'] if isinstance(args['contig_name'],str) else args['contig_name'][0] 	# The name of the contig to extract the reference region, first contig in the fasta used if not specified
		if len(args['start_pos'])==2:						# If -s, --start is specified with two arguments or is not specified at all (default)
			start_pos=args['start_pos'][0]					# The start position of the target region on the contig (to be used as the reference) NOT compatible with random
			num_mu=args['start_pos'][1]						# Number of random mutations for the target sequence
		elif len(args['start_pos'])==1:						# If -s, --start is specified with just one argument
			start_pos=args['start_pos'][0]
			num_mu=1
		else:												# If -s, --start is specified with no arguments
			start_pos=1										# The default values should be passed in this case
			num_mu=1
		ref_length=args['ref_length'] if isinstance(args['ref_length'],int) else args['ref_length'][0]	# The length of the reference
		timeout=args['timeout'] if isinstance(args['timeout'],float) else args['timeout'][0]			# The timeout limit for individual processes
		freq=[args['freq'],0,0]	if isinstance(args['freq'],float) else [args['freq'][0],0,0]            # Mutation frequency passed to haplogenerator
		emp_filename='' if args['empirical'] is None else args['empirical'][0] 
		if len(args['lognormal'])==2:
			lognormal_dist=[args['lognormal'][0],0,0]       # Mean of the log mutation distance passed to haplogenerator with lognormal model
			lognormal_sd=[args['lognormal'][1],1,1] 		# sd of the log mutation distance passed to haplogenerator with lognormal model
		elif len(args['lognormal'])==1:
			lognormal_dist=[args['lognormal'][0],0,0]
			lognormal_sd=[1,1,1]
		else:
			lognormal_dist=[1,0,0]
			lognormal_sd=[1,1,1]
		if len(args['random_str_mu'])==2:				# If -r, --random is specified with two arguments or is not specified at all (default)
			num_random_str=args['random_str_mu'][0]		# Number of random regions (NOT compatible with start_pos)
			num_random_mu=args['random_str_mu'][1]		# Number of random mutations on each random region (NOT compatible with start_pos)
		elif len(args['random_str_mu'])==1:				# If -r, --random is specified with just one argument
			num_random_str=args['random_str_mu'][0]	
			num_random_mu=1
		else:											# If -r, --random is specified with no arguments
			num_random_str=1							# Default random values are passed in this case
			num_random_mu=1	
		bwa_mem=args['bwa_mem']					        # flag to use bwa mem for NGS read alignment
		no_missing=args['no_missing']					# argument to save all of the simulation results including those simulations that fail for one or two of the haplotyping algorithms 
		remove_duplicate=args['remove_duplicate']       # flag to remove duplicates by Picardtools
		report_aln=args['report_aln']					# flag to save alignment and vcf files for each simulation in a zip archive
		timem=args['timem_use']							# flag to report the time and memory consumption of each haplotyping method in a separate dat file
		verbose=args['verbose']							# verbosity flag to print warning and progress messages to stderr
		coverage=args['coverage'] if isinstance(args['coverage'],int) else int(args['coverage'][0])# Sequencing coverage per homologue to generate the NGS reads
		method=args['method'] if isinstance(args['method'],str) else args['method'][0]  # Sequencing technology to generate the NGS reads 
		pe=args['pe'] if isinstance(args['pe'],str) else args['pe'][0]   				# Single, paired end or mate pair reads
		maxhap=args['maxhap'] if isinstance(args['maxhap'],int) else args['maxhap'][0]  
		if len(_extra1)>1:
			raise CustomException('ERROR: too many library type arguments: '+'\t'.join(_x for _x in _extra1))
		elif _extra1:
			method=_extra1[0]
		if len(_extra2)>1:
			raise CustomException('ERROR: too many library type arguments: '+'\t'.joint(_x for _x in _extra2))
		elif _extra2:
			pe=_extra2[0]
		if not pe:	# Assign the true default values for the library type and sequencing method here
			pe='PE'
		if not method:
			method='HS25'
		plotting=args['plot']								
		if args['dosages']==[-1,-1]: 
			if verbose:
				print >>sys.stderr, "WARNING: no dosage distribution has been specified! All markers will have the same probability {0:.5f}...".format(1./ploidy)
			dosage_dist=[round(1./ploidy, 5) for _x in range(0, ploidy)]
		else:
			dosage_dist=args['dosages']
		if len(dosage_dist)>ploidy:
			if verbose:
				print >> sys.stderr, ("WARNING: You have specified probabilities for dosages which are larger than\
 the ploidy level!\nThese will be ignored!")
			dosage_dist=dosage_dist[0:ploidy]
		dosage_sum=sum(dosage_dist)
		if round(dosage_sum,4)>1 or (round(dosage_sum,4)<1 and len(dosage_dist)==ploidy):
			raise CustomException("ERROR: dosage probabilities from simplex to ploidy-plex must sum up to 1!")
		elif round(dosage_sum,4)<1:
			if verbose:
				print >> sys.stderr, "WARNING: You have not specified probabilities for dosages more than %i and the sum of probabilities is less than 1!\n\
Missing probabilities will be set equally, so that the total sum of dosage probabilities becomes one!" %len(dosage_dist)
			_perc=(1-dosage_sum)/float(ploidy-len(dosage_dist))
			for _i in range(len(dosage_dist), ploidy):
				dosage_dist.append(_perc)
		elif len(dosage_dist)<ploidy:
			if verbose:
				print >> sys.stderr, "WARNING: You have not specified probabilities for dosages more than %i!\
 Missing probabilities will be set to zero!" %len(dosage_dist)
			for _i in range(len(dosage_dist), ploidy):
				dosage_dist.append(0)
		else:
			pass
		for _i in range(0, ploidy):
			observed_dosage.append(0)
		if plotting:
			with NamedTemporaryFile(delete=True) as pdffile:
				pass
			axis_font = {'size':'5', 'color':'black', 'weight':'normal'}
			suptitle_font = {'size':'8', 'color':'black', 'weight':'bold','verticalalignment':'bottom'}
			title_font = {'size':'7', 'color':'black', 'weight':'normal','verticalalignment':'bottom'} # Bottom vertical alignment for more space
			gs=gridspec.GridSpec(3, 3, wspace=0.25, hspace=0.49) 	# grid to plot the histograms of SNP distances in a pdf file, 
			next_plot=1						# with pages containing 3x3 subplots.
			pdfpages=PdfPages(pdffile.name)				# the pdf pages to save the figures
			fig=plt.figure()					# initial figure to be added to pdf
			fig.suptitle('Histograms of the distances between neighbor SNPs', **suptitle_font)
		if verbose:
			sequencer_log=output+'_'+'_'.join(time.strftime("%c").split()[0:3])+'_sequencer.log'
			with open(sequencer_log,'w') as _seq_log:
				pass
			if lognormal_dist[0]>0:
				print >>sys.stderr, ("WARNING: lognormal process will be applied to generate random genomes!")
		if timem:
			timem_dat=output+'_timem.dat'
			with open(timem_dat,'w') as _timem_dat:
				_timem_dat.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format("MAX_RSS", "US_TIME", "METHOD", "ID_REGION", "ID_MUTATION", "RUN_TIME_ERROR"))
		if report_aln:
			zip_arc=output+'_'+'_'.join(time.strftime("%c").split()[0:3])+'_ALLFILES.zip'
			try:
				import zlib
				compression = zipfile.ZIP_DEFLATED
			except:
				compression = zipfile.ZIP_STORED			
			z=zipfile.ZipFile(zip_arc,'w',compression)
			if verbose and compression==zipfile.ZIP_STORED:
				print >>sys.stderr, ("WARNING: the files in the zip archive {} will not be compressed as zlib library was not available!").format(zip_arc)
		if ref_length<2:
			raise CustomException("ERROR: minimum required reference length is 2! A smaller value was given!")
		with NamedTemporaryFile(delete=False) as errors_in_haplotyping: 	# The log file to write haplotyping and evaluation errors,
			errors_in_haplotyping.write("Haplosim log report, "+		# erroneous haplotypes, etc.
			time.strftime("%a %d %b %Y")+" "+time.strftime("%X")+".\n")
		if not num_random_str: # The default is to use the first base of the contig as the start of the target region and do just one random mutation on it
			N=1			 # Number of regions
			M=num_mu 	 # Number of random mutations (default 1)
			random=False # Random mode is Off
		else:
			if verbose:
				print >>sys.stderr, ("WARNING: random start position(s) will be used for the selecting the reference!")
			N=num_random_str # Number of random regions to be generated
			M=num_random_mu  # Number of random mutations on each random region
			random=True      # Random mode is On
		if maxhap>0:	     # Start building the haplotype set for the whole population.
			with NamedTemporaryFile(delete=True) as zip_pophap: # Name of the zip archive to put the population haplotypes in it.
				if verbose:
					print >>sys.stderr, ("WARNING: the total number of haplotypes will be limited to maxhap {0} for each reference genome!\n\
For each of {1} references, {2} random genomes will be simulated by sampling with replacement from the {0} haplotypes of each reference!").format(maxhap, N, M)
		with NamedTemporaryFile(delete=True) as subprocess_err: # tmp file to store and display the eventual error messages from the subprocesses
			tmpfiles.append(subprocess_err.name)
		with open(fastafile,'rU') as handle: # Extract the target record from the input fasta for selecting the references 
			if not contig_name: # if contig name is not specified, use the first contig in the fasta
				if verbose:
					print >>sys.stderr, ("WARNING: you have provided no contig name for the reference! The first contig of the fasta file will be used!")
				try:
					target_record_original = copy.deepcopy(SeqIO.parse(handle, "fasta").next())
					contig_name = target_record_original.__getattribute__('id')
				except StopIteration:
					raise CustomException("ERROR: The provided fasta file is empty!")
			else:		        # if the contig name is given, search for it within the fasta 
				for record in SeqIO.parse(handle,"fasta"):
					if not contig_found:
						_id=record.__getattribute__('id')
						if _id==contig_name:
							target_record_original = copy.deepcopy(record)
							contig_found = True		
					else:
						break
				if not contig_found:
					raise CustomException("ERROR: fasta file did not contain the specified contig!")
		if not str(target_record_original.__getattribute__('seq')):
			raise CustomException("ERROR: Empty contig was chosen! Choose another contig!")
		elif random and ref_length and len(str(target_record_original.__getattribute__('seq')))<ref_length:
			raise CustomException("ERROR: The given reference length is larger than the whole target contig!\n\
Either choose another contig or choose a smaller reference length!")
		elif random and (not ref_length) and len(str(target_record_original.__getattribute__('seq')))<2:
			raise CustomException("ERROR: Contig length is smaller than 2bp! Contig length should be at least 2bp to contain a valid reference!")
		elif (not random) and ref_length and (len(str(target_record_original.__getattribute__('seq')))<ref_length+start_pos):
				raise CustomException("ERROR: Given reference length and start position exceed the contig range! Choose another length, start position or contig!")
		elif (not random) and (not ref_length) and len(str(target_record_original.__getattribute__('seq')))<start_pos+2:
			raise CustomException("ERROR: Contig length from the given start point is smaller than 2bp! At least 2bp is neaded from the start point to contain a valid reference!")
		elif verbose:
			print("Contig selected from the fasta file to extract the reference region...")
		else:
			pass
		for _n in range(0,N):
			with NamedTemporaryFile(delete=False) as tmp_ref:
				_selection_attempts=0
				_selection_done=False
				while (not _selection_done) and _selection_attempts<15:	
					if random:         # write the target region to the (temporary) fasta file
						if ref_length:
							start_pos=rnd.sample(range(0,len(str(target_record_original.__getattribute__('seq')))-ref_length+1),1)[0]
							target_str=str(target_record_original.__getattribute__('seq'))[start_pos:start_pos+ref_length]
						else:
							start_pos=rnd.sample(range(0,len(str(target_record_original.__getattribute__('seq')))-1),1)[0]
							target_str=str(target_record_original.__getattribute__('seq'))[start_pos:len(str(target_record_original.__getattribute__('seq')))]
							if verbose and _n==0:
								print >>sys.stderr, ("WARNING: no reference length specified! The reference(s) will stretch from the start position to the end of the contig!")
					else:
						start_pos-=1   # Change the genomic coordinates to python coordinates 
						if ref_length:
							target_str=str(target_record_original.__getattribute__('seq'))[start_pos:start_pos+ref_length]
						else:
							target_str=str(target_record_original.__getattribute__('seq'))[start_pos:len(str(target_record_original.__getattribute__('seq')))]
							if verbose and _n==0:
								print >>sys.stderr, ("WARNING: no reference length specified! The reference(s) will stretch from the start position to the end of the contig!")
					target_record=copy.deepcopy(target_record_original)
					target_record.seq=Seq(''.join(_base for _base in target_str if _base!='N'), generic_dna) # Remove the 'N's from target seq
					if len(str(target_record.__getattribute__('seq')))<0.8*len(target_str):					 # This missing of the reference length should not be more than 20%	
						_selection_attempts+=1
						if _selection_attempts==15 or (not random):
							raise CustomException("ERROR: the target contig contains too much 'N' values! Either choose a different start point or a different contig and try again!\n\
	'N' values should not compose more than 20% of the reference! They will be omitted from the reference in case <20%!")
					else:
						_selection_done=True
				_write=SeqIO.write(target_record, tmp_ref, "fasta") # Write the target seq on a fasta file
				_write, _base, target_record = (None, None, None)
			if verbose:
				print("Reference region %s extracted from the specified contig..." % str(_n+1))
			cmd=Command(samtoolsPath+"samtools faidx "+tmp_ref.name, outfile=subprocess_err.name)   # index the (temporary) reference file
			try:
				exit_msg=cmd.run(timeout=timeout)
			except subprocess.CalledProcessError:
				with open(subprocess_err.name, 'rU') as subprocess_err:
					subprocess_err_msg=[_x if ('*' not in _x and _x!='\n') else '' for _x in subprocess_err.readlines()]
					subprocess_err_msg='\n'+'Original error message:'+'\t'+''.join(''.join(str(_y) for _y in _x) for _x in subprocess_err_msg)+'\n'+"-------------------------------"+'\n'
				raise CustomException("Failed to make fasta index! Check the help for samtools faidx!"+subprocess_err_msg)
			with NamedTemporaryFile() as tmp_hap:
				pass
			#Mu=fixed_snps(target_str) # Fixed SNP positions at the beginning, end and along the target_str.
			target_str=None
			if maxhap>0:
				if _n==0:
					import numpy.random as nprnd
				pophapname = make_pophap(tmp_ref.name, maxhap, zip_pophap.name, emp_filename, lognormal_dist, lognormal_sd, freq, timeout, verbose)
				poplist=['hap'+str(_x) for _x in range(0, maxhap)]
				tmpfiles.extend(tmp_ref.name+"_for_"+_haploid+"_hap1.fa.fai" for _haploid in poplist)
				tmpfiles.extend(tmp_ref.name+"_for_"+_haploid+"_hap1.fa" for _haploid in poplist)
				tmpfiles.extend(tmp_ref.name+"_for_"+_haploid+"_varianthaplos.txt" for _haploid in poplist)
				pophap=zipfile.ZipFile(pophapname, 'r')
			for _m in range(0,M): # The reference has been constructed and indexed by now. The next step is to make mutations.
				if maxhap>0: 	  # If maxhap is set to some finite positive integer
					_homozygous=True
					while _homozygous:  # Continue the selection of the haplotypes until a heterozygous genome is generated
						_x=0
						dosage_hap = [] 	# Random dosage for each to be selected haplotype
						shuffled_haps = []  # The list to contain the haplotypes with their respective dosages
						rnd_seed = np.random.uniform(low=9000, high=95000 , size=(ploidy+1,))  # Maximum number of seeds needed to select the haplotypes
						for _y in range(0,ploidy):
							dosage_hap.append(1+int((64*(2**10)-1)*rnd.random())%ploidy)       # Generate maximum number of random dosages needed for the selected haplotypes, an integer >=1 and <=ploidy
						while _x<ploidy:
							r = nprnd.RandomState(int(rnd_seed[_x]))
							tmp_shuffled_haps = r.choice(poplist, size=1).tolist()*dosage_hap[_x] # Select a random haplotype (Sampling With Replacement). 
							_x+=1
							shuffled_haps.extend(tmp_shuffled_haps)
						r = nprnd.RandomState(int(rnd_seed[_x]))
						shuffled_haps = r.choice(shuffled_haps, size=ploidy, replace=True).tolist()  # Throw away extra haplotypes, i.e. more than the ploidy, randomly
						if len(set(shuffled_haps))>1:
							_homozygous=False
					for _haploid in shuffled_haps:
						for _file in [_haploid+"_hap1.fa.fai", _haploid+"_hap1.fa", _haploid+"_varianthaplos.txt"]:
							with open(tmp_ref.name+"_for_"+_file,'w') as _extractee, pophap.open(_file, 'rU') as _zipop:
								for _line in _zipop:
									_extractee.write(_line)
					cmd=Command(('allowrapper -p '+str(ploidy)+' -f '+tmp_ref.name+' '+' '.join(tmp_ref.name+"_for_"+_haploid for _haploid in shuffled_haps)+
								' '+tmp_hap.name).split(' '), shell=False, outfile=subprocess_err.name) #Use allowrapper to combine the seleceted haplotypes into a single polyploid genome
				elif emp_filename:
					cmd=Command("haplogenerator.py -f "+tmp_ref.name+" -o "+tmp_hap.name+" -p "+str(ploidy)+
						' -s "'+str([1,0,0])+'"'+" -m "+'"'+"{'A':'C','C':'A','G':'T','T':'G'}"+'"'+ 
						' --dosage "'+str(dosage_dist)+'"'+#' -S "'+str(int(ploidy/2)*[Mu,[]])+'"'+ 
						' --empirical '+str(emp_filename)+' --genomic', outfile=subprocess_err.name)  # Generate the artificial genome with the lognormal model
				elif lognormal_dist[0]>0:
					cmd=Command("haplogenerator.py -f "+tmp_ref.name+" -o "+tmp_hap.name+" -p "+str(ploidy)+
						' -s "'+str(lognormal_dist)+'"'+" -m "+'"'+"{'A':'C','C':'A','G':'T','T':'G'}"+'"'+ 
						' --dosage "'+str(dosage_dist)+'"'+#' -S "'+str(int(ploidy/2)*[Mu,[]])+'"'+ 
						' --model lognormal --sdlog "'+str(lognormal_sd)+'" --genomic', outfile=subprocess_err.name)  # Generate the artificial genome with the lognormal model
				else:
					cmd=Command("haplogenerator.py -f "+tmp_ref.name+" -o "+tmp_hap.name+" -p "+str(ploidy)+
						' -s "'+str(freq)+'"'+" -m "+'"'+"{'A':'C','C':'A','G':'T','T':'G'}"+'"'+ 
						' --dosage "'+str(dosage_dist)+'"'+#' -S "'+str(int(ploidy/2)*[Mu,[]])+'"'+ 
						" --model poisson --genomic", outfile=subprocess_err.name)  # Generate the artificial genome with the (default) Poisson stochastic model
				homolos=list(tmp_hap.name+"_hap"+str(_i)+".fa" for _i in range(1, ploidy+1))
				tmpfiles.extend(_x+'.fai' for _x in homolos)
				tmpfiles.extend(homolos)
				tmpfiles.append(tmp_hap.name+"_varianthaplos.txt")
				haplogen_err_msg=''
				try:
					exit_msg=cmd.run(timeout=timeout)
				except subprocess.CalledProcessError:
					with open(subprocess_err.name, 'rU') as subprocess_err:
						subprocess_err_msg=[_x if ('*' not in _x and _x!='\n') else '' for _x in subprocess_err.readlines()]
						subprocess_err_msg='\n'+'Original error message:'+'\t'+''.join(''.join(str(_y) for _y in _x) for _x in subprocess_err_msg)+'\n'+"-------------------------------"+'\n'
					raise CustomException("ERROR: Failed to make the artificial genome! Check the help for haplogenerator!"+subprocess_err_msg)
				else:
					if verbose:
						print("Successfully produced the artificial genome for region %i, mutation %i...") %(_n+1,_m+1)
				tmp_reads=list("tmp_reads_hap"+str(_i) for _i in range(1, ploidy+1))
				for _i, _homolo in enumerate(homolos):
					with NamedTemporaryFile() as tmp_reads[_i]:
						tmpfiles.extend(tmp_reads[_i].name+_fq_suffix for _fq_suffix in ("1.fq","2.fq","fq","_0001.fastq","_0001.maf","_0001.ref"))
				with NamedTemporaryFile() as tmp_seq_log, NamedTemporaryFile() as removed_homozygous: # 1: Temporary log for the sequencer simulator
					tmpfiles.extend([tmp_seq_log.name, removed_homozygous.name])					  # 2: Temporary file for the filtered variant haplotypes
				if plotting:	
					variants=deque()
				with open(tmp_hap.name+"_varianthaplos.txt", 'rU') as variants_get, open(removed_homozygous.name, 'w') as new_haplovariants:
					new_haplovariants.write(variants_get.readline()) # Write the header to new_haplovariants and skip it
					for _line in variants_get:
						_alleles=[int(_x) for _x in _line.split()[5:]]
						observed_dosage[np.sum(np.asanyarray(_alleles)>0)-1]+=1 # Update the counts of the variant types 
						if np.sum(np.asanyarray(_alleles)>0)!=ploidy: # Throw away homozygous variants from the new variant haplotypes
							new_haplovariants.write(_line)
						if plotting:
							variants.append(int(_line.split()[2])) 	 # The variant (including homozygous) positions needed to draw density plots
				move(removed_homozygous.name, tmp_hap.name+"_varianthaplos.txt")
				_alleles=None
				if plotting:
					variants_distance=deque()
					try:
						variants_distance.append(variants.popleft())
						last=variants_distance.pop()
						variants_distance.append(last)
					except:
						raise
					while variants:
						annex=variants.popleft()
						variants_distance.append(annex-last)
						last=annex
					del variants
					variants_distance=list(variants_distance)
					variants_distance.sort()
					_1stQuartile=variants_distance[max(0, int(round(0.25*len(variants_distance)))-1)]
					_3rdQuartile=variants_distance[min(len(variants_distance)-1, int(round(0.75*len(variants_distance)))+1)]
					num_bins=max(1, int((variants_distance[-1]-variants_distance[0])*exp(log(len(variants_distance))/3)/(2*(_3rdQuartile-_1stQuartile))))#Choosing the number of bins according to Freedman and Diaconis (1981)
					if next_plot>9:
						pdfpages.savefig()  # save the current page and go to the next page
						plt.close(fig)			# close the current figure 
						fig=plt.figure()	# start the plots in the next page
						fig.suptitle('Histograms of the distances between neighbor SNPs', **suptitle_font)
						next_plot=1			# start plot numbering for the new current page
					hist, bin_edges=np.histogram(variants_distance, bins=num_bins, density=True)
					histo=fig.add_subplot(gs[int(next_plot-0.001)/3, next_plot%3-1])
					histo.bar(bin_edges[:-1], hist, width = (max(variants_distance)-min(variants_distance))/float(num_bins), 
												color="b", edgecolor="b")
					histo.set_xlabel('Distance', **axis_font)   # add xlabel to the subplots
					if next_plot%3==1:  # add xlabel to the plots on the first column
						histo.set_ylabel('Probability', **axis_font)
					else:
						pass
					histo.set_title('Region {0}, Mutation {1}'.format(_n+1, _m+1),
									**title_font)
					histo.set_xlim(0, int(1.05*np.amax(bin_edges)))
					histo.set_ylim(0.55*np.amin(hist), min(1, 1.19*np.amax(hist)))
					xticks = [int(round(_x+(_y-_x)/2)) for _x, _y in zip(bin_edges[:-1], bin_edges[1:])][::max(1, int(len(bin_edges)-0.001)/10+1)]
					histo.xaxis.set_ticks(xticks)
					histo.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
					for _tick in histo.xaxis.get_major_ticks():
						_tick.label.set_fontsize(4) 
					for _tick in histo.yaxis.get_major_ticks():
						_tick.label.set_fontsize(4)
					next_plot+=1
					histo, _tick, xticks, hist, bin_edges = None, None, None, None, None
				if verbose:
					print("NGS reads being generated from the artificial genome...")
				rnd_seed=np.random.uniform(low=1847, high=7919 , size=(ploidy,)) # Generate a random seed (for each homologue) to generate NGS reads for it
				if method in ['GA1','GA2','HS10','HS20','HS25','MS']:
					read_lengths=dict(zip(['GA1','GA2','HS10','HS20','HS25','MS'],[35,55,75,100,125,250]))
					if not insert:
						frag_lengths=dict(zip(['GA1','GA2','HS10','HS20','HS25','MS'],[70,100,150,250,300,550]))
					else:
						frag_lengths=dict(zip(['GA1','GA2','HS10','HS20','HS25','MS'],list(insert for _x in [35,55,75,100,125,250])))
					for _i, _homolo in enumerate(homolos):
						if pe=='PE':	
							_method_id=0
							cmd=Command((ARTSPath+"art_illumina -na -i "+_homolo+" -p -l "+
								str(read_lengths[method])+" -ss "+method+" --fcov "+str(coverage)+" --mflen "+
								str(frag_lengths[method])+" --sdev 10 -o "+tmp_reads[_i].name+" --id ART_PE_"+method+"-"+str(_i)+"-"+
								" --rndSeed "+str(rnd_seed[_i])).split(' '), shell=False, outfile=tmp_seq_log.name, keep_error_file=True)
						elif pe=='SE':
							if insert and verbose and _i==0:
								print >>sys.stderr, ('WARNING: the specified insert size will be ignored as the reads are set to be single!')
							_method_id=1
							cmd=Command((ARTSPath+"art_illumina -na -i "+_homolo+" -l "+str(read_lengths[method])+" -ss "+method+
							" --fcov "+ str(coverage)+ " -o "+tmp_reads[_i].name+" --id ART_SINGLE_"+method+"-"+str(_i)+"-"+
							" --rndSeed "+str(rnd_seed[_i])).split(' '), shell=False, outfile=tmp_seq_log.name, keep_error_file=True)
						else:
							_method_id=2
							cmd=Command((ARTSPath+"art_illumina -na -i "+_homolo+" -mp -l "+str(read_lengths[method])+" -ss "+method+
							" --fcov "+str(coverage)+ " --mflen "+str(frag_lengths[method])+" --sdev 50 -o "+tmp_reads[_i].name+" --id ART_MP_"+method+"-"+str(_i)+"-"+ 
							" --rndSeed "+str(rnd_seed[_i])).split(' '), shell=False, outfile=tmp_seq_log.name, keep_error_file=True)
						try:
							exit_msg = cmd.run(timeout=timeout)
						except subprocess.CalledProcessError:
							raise CustomException("Failed to generate {} Illumina reads! Check the sequencing log file (with verbose set to true) and the help for ART!".format(("paired-end", "single-end", "mate-pair")[_method_id])) 
						finally:
							if verbose and os.path.exists(tmp_seq_log.name):
								with open(sequencer_log,'a') as _seq_log, open(tmp_seq_log.name,'rU') as _tmp_log:
									for _line in _tmp_log:
										_seq_log.write(_line)
								os.remove(tmp_seq_log.name)
				else:	# method == 'CCS' or method=='CLR', PacBio "Circular Concensus Sequencing" and "Continuous Long Reads", respectively. 
					min_accuracy={'CCS':str(0.989),'CLR':str(0.75)}
					max_accuracy={'CCS':str(1),'CLR':str(0.89)}
					accuracy_mean={'CCS':str(0.995),'CLR':str(.82)}
					length_max={'CCS':str(15000),'CLR':str(75000)}
					if not insert:
						length_mean={'CCS':str(1000),'CLR':str(5000)}
					else:
						length_mean={'CCS':str(insert),'CLR':str(insert)}
					for _i, _homolo in enumerate(homolos):
						cmd=Command((PBSIMPath+"Linux-amd64/bin/pbsim "+_homolo+" --prefix "+tmp_reads[_i].name+
						" --depth "+ str(coverage)+ " --data-type "+ method +" --model_qc "+ 
						PBSIMPath+"data/model_qc_"+method.lower()+ " --seed "+str(rnd_seed[_i])+" --accuracy-min " + min_accuracy[method] +
						" --accuracy-max " + max_accuracy[method] + " --length-mean " + length_mean[method] + 
						" --accuracy-mean " + accuracy_mean[method]+ " --length-max "+ length_max[method]).split(' '), shell=False, outfile=tmp_seq_log.name, keep_error_file=True)
						try:
							exit_msg=cmd.run(timeout=timeout)
						except subprocess.CalledProcessError:
							raise CustomException("Failed to generate PacBio reads! Check the sequencing log file (with verbose set to true) and the help for PBSIM!")
						finally:
							if verbose and os.path.exists(tmp_seq_log.name):
								with open(sequencer_log,'a') as _seq_log, open(tmp_seq_log.name,'rU') as _tmp_log:
									for _line in _tmp_log:
										_seq_log.write(_line)
								os.remove(tmp_seq_log.name)
				for _homolo in homolos:	# Clean-up, delete the separate fasta files and indexes of each homologue
					os.remove(_homolo)	# They are not needed anymore as reads have been generated from them
					os.remove(_homolo+'.fai')
				if os.path.exists(tmp_seq_log.name):
					os.remove(tmp_seq_log.name)
				homolos, _homolo = (None, None)
				with NamedTemporaryFile(delete=False) as raw_reads1, NamedTemporaryFile(delete=False) as raw_reads2: # Put the reads from all of the homologues in a single file
					tmpfiles.extend([raw_reads1.name,raw_reads2.name])
					for _i, _reads in enumerate(tmp_reads):
						if (pe=="PE" or pe=="MP") and (not (method in ("CLR","CCS"))) :
							with open(_reads.name+"1.fq",'rU') as _homolo_1, open(_reads.name+"2.fq",'rU') as _homolo_2:
								for _read in _homolo_1:
									raw_reads1.write(_read)
								for _read in _homolo_2:
									raw_reads2.write(_read)
						elif (not (method in ("CLR","CCS"))):
							with open(_reads.name+".fq",'rU') as _homolo_1:
								for _read in _homolo_1:
									raw_reads1.write(_read)
						else: # Reads come from PBSIM
							with open(_reads.name+"_0001.fastq",'rU') as _homolo_1:
								for _lnum, _read in enumerate(_homolo_1):
									if _lnum % 4 == 0 or _lnum % 4 == 2: # Make the read id's unique in the final fastq file
										raw_reads1.write(_read[0]+_reads.name+_read.rstrip()[1:]+"\n") # id lines divideable by 4 (before the read seq) or with remainder 2 (before the read quality line) 
									else:	
										raw_reads1.write(_read)
								del _lnum
				_reads, _homolo_1, _homolo_2 = (None, None, None)
				if (pe!="PE" and pe!="MP") or (method in ("CLR","CCS")):
					os.remove(raw_reads2.name)
				if report_aln:
					for _i, _fastqfile in enumerate((raw_reads1.name, raw_reads2.name)):
						if os.path.exists(_fastqfile):
							z.write(_fastqfile,"myreads"+str(_n+1)+"_"+str(_m+1)+"."+str(_i+1)+".fastq")
				for _reads in tmp_reads:			# Clean-up, delete the read files of individual homologues
					if method in ("CLR","CCS"):
						os.remove(_reads.name+'_0001.fastq')
						os.remove(_reads.name+'_0001.maf')
						os.remove(_reads.name+'_0001.ref')
					elif pe=="PE" or pe=="MP":
						os.remove(_reads.name+'1.fq')
						os.remove(_reads.name+'2.fq')
					else:
						os.remove(_reads.name+'.fq')
				del _reads, tmp_reads
				if (not (method in ("CLR","CCS"))) and (not bwa_mem): # In case Illumina is the sequencer and bowtie2 is used for alignment and therefore reference indexes should be made for bowtie
					with NamedTemporaryFile() as bt2_indx, NamedTemporaryFile() as aligner_sam:
						_btindxdir=os.path.realpath(bt2_indx.name)
						_dir='/'.join(_btindxdir.split('/')[:-1])
						tmpfiles.append(aligner_sam.name)
					cmd=Command(bowtiePath+"bowtie2-build -f "+tmp_ref.name+" "+bt2_indx.name+" -q", outfile=subprocess_err.name)
					try:
						exit_msg=cmd.run(timeout=timeout)
					except subprocess.CalledProcessError:
						with open(subprocess_err.name, 'rU') as subprocess_err:
							subprocess_err_msg=[_x if ('*' not in _x and _x!='\n') else '' for _x in subprocess_err.readlines()]
							subprocess_err_msg='\n'+'Original error message:'+'\t'+''.join(''.join(str(_y) for _y in _x) for _x in subprocess_err_msg)+'\n'+"-------------------------------"+'\n'
						raise CustomException("Failed to build bowtie reference! Check the help for bowtie2!"+subprocess_err_msg)
				elif (not (method in ("CLR","CCS"))):# In case Illumina is the sequencer and bwa-mem is used for alignment and therefore reference indexes should be made for bwa-mem
					with NamedTemporaryFile() as aligner_sam:
						tmpfiles.append(aligner_sam.name)
				else:
					tmpfiles.append(raw_reads1.name+".fastq")
					copyfile(raw_reads1.name,raw_reads1.name+".fastq")
					if _m==0:
						copyfile(tmp_ref.name, tmp_ref.name+".fasta")
						copyfile(tmp_ref.name+".fai", tmp_ref.name+".fasta.fai")
					with NamedTemporaryFile() as aligner_sam:
						tmpfiles.append(aligner_sam.name)
				if verbose:
					print("Reads being aligned to the reference...")
				if method in ("CLR", "CCS") and bwa_mem:
					if _m==0:
						cmd=Command(bwaPath+"bwa index " + tmp_ref.name + ".fasta;" +
						bwaPath + "bwa mem -x pacbio " + tmp_ref.name + ".fasta " + raw_reads1.name+".fastq >" + aligner_sam.name, outfile=subprocess_err.name)	
					else:
						cmd=Command(bwaPath + "bwa mem -x pacbio " + tmp_ref.name + ".fasta " + raw_reads1.name+".fastq >" + aligner_sam.name, outfile=subprocess_err.name)
				elif method in ("CLR", "CCS"):
					cmd=Command(blasrPath+"blasr " + raw_reads1.name + ".fastq" + " " + tmp_ref.name+".fasta -bestn 1\
 -minPctIdentity 80 -sam -out " + aligner_sam.name, outfile=subprocess_err.name)
				elif (pe=="PE" or pe=="MP") and bwa_mem:
					if _m==0:
						cmd=Command(bwaPath+"bwa index " + tmp_ref.name + ";" +
						bwaPath + "bwa mem " + tmp_ref.name + " " + raw_reads1.name+
						" "+raw_reads2.name + " > " + aligner_sam.name, outfile=subprocess_err.name)	
					else:
						cmd=Command(bwaPath + "bwa mem " + tmp_ref.name + " " + raw_reads1.name+ 
						" "+raw_reads2.name + " > " + aligner_sam.name, outfile=subprocess_err.name)
				elif (pe=="PE" or pe=="MP"):
					cmd=Command(bowtiePath+"bowtie2 --sensitive -x "+bt2_indx.name+" -1 "+raw_reads1.name+
						" -2 "+raw_reads2.name+" -S "+aligner_sam.name+" -q", outfile=subprocess_err.name)
				elif not bwa_mem:
					cmd=Command(bowtiePath+"bowtie2 --sensitive -x "+bt2_indx.name+" -U "+raw_reads1.name+
						" -S "+aligner_sam.name+" -q", outfile=subprocess_err.name)
				else:
					if _m==0:
						cmd=Command(bwaPath+"bwa index " + tmp_ref.name + ";" +
						bwaPath + "bwa mem " + tmp_ref.name + " " + raw_reads1.name+
						" > " + aligner_sam.name, outfile=subprocess_err.name)	
					else:
						cmd=Command(bwaPath + "bwa mem " + tmp_ref.name + " " + raw_reads1.name+ 
						" > " + aligner_sam.name, outfile=subprocess_err.name)
				try:
					exit_msg=cmd.run(timeout=timeout)
				except subprocess.CalledProcessError:
					with open(subprocess_err.name, 'rU') as subprocess_err:
						subprocess_err_msg=[_x if ('*' not in _x and _x!='\n') else '' for _x in subprocess_err.readlines()]
						subprocess_err_msg='\n'+'Original error message:'+'\t'+''.join(''.join(str(_y) for _y in _x) for _x in subprocess_err_msg)+'\n'+"-------------------------------"+'\n'
					if method in ("CLR","CCS") and bwa_mem:
						raise CustomException("Failed to align PacBio reads in sam format! Check the help for bwa-mem!"+subprocess_err_msg)
					elif method in ("CLR","CCS"):
						raise CustomException("Failed to align PacBio reads in sam format! Check the help for BLASR!"+subprocess_err_msg)
					elif not bwa_mem:
						raise CustomException("Failed to align Illumina reads in sam format! Check the help for bowtie2!"+subprocess_err_msg)
					else:
						raise CustomException("Failed to align Illumina reads in sam format! Check the help for bwa-mem!"+subprocess_err_msg)
				else:
					with NamedTemporaryFile(delete=False) as RGremoved_sam, open(aligner_sam.name,'rU') as RG_sam:
						tmpfiles.append(RGremoved_sam.name)
						for _line in RG_sam:
							if '@HD\t' in _line:# Correct @HD orientation field to be compatible with Picardtools
								RGremoved_sam.write(substitute(r"SO[:].*\t", r"SO:coordinate\t", _line)) 
							elif '@RG\t' in _line:# Remove any RG, as it will be later added by Picardtools
								pass
							else:
								RGremoved_sam.write(_line)
					move(RGremoved_sam.name, aligner_sam.name)
					RGremoved_sam, RGsam, _line = (None, None, None)
					if report_aln:
						z.write(tmp_ref.name,"myref"+str(_n+1)+"_"+str(_m+1)+".fa")
						z.write(tmp_ref.name+".fai","myref"+str(_n+1)+"_"+str(_m+1)+".fa.fai")
						z.write(aligner_sam.name,"mysam"+str(_n+1)+"_"+str(_m+1)+".sam")
				finally:		
					if (not (method in ("CLR","CCS"))) and (not bwa_mem):
						tmpfiles.extend(_dir+"/"+_x for _x in os.listdir(_dir) if _btindxdir.split('/')[-1] in _x) # If bowtie2 is used as aligner, its index files should be scheduled for deletion
						for _file in (_dir+"/"+_x for _x in os.listdir(_dir) if _btindxdir.split('/')[-1] in _x):  # Early Clean-up, remove the bowtie index files as no longer needed
							os.remove(_file)																   	   # In case of the failure of the clean-up, these files will be still deleted when tmp files are finally cleaned
				with NamedTemporaryFile() as aligner_bam, NamedTemporaryFile() as tmp_vcf: 
					tmpfiles.extend([aligner_bam.name+".bam", aligner_bam.name+".bam.bai", tmp_vcf.name, 
									 aligner_bam.name+".RG.bam"])		   # Temporary sorted bam and vcf files added to the clean-up list
				if verbose:
					print("Variants being called from the aligned reads...")
				if method in ("CCS","CLR"):
					platform = "PacBio"
				else:
					platform = "Illumina"
				with NamedTemporaryFile(delete=True) as _tmp_inbam:
					tmpfiles.append(_tmp_inbam.name)
				cmd=[samtoolsPath+"samtools view -bS -o "+_tmp_inbam.name+" "+aligner_sam.name+";"+samtoolsPath+"samtools sort -o " + aligner_bam.name+".bam "+_tmp_inbam.name+";rm "+_tmp_inbam.name,# Convert sam to sorted bam file
					"java -Xmx8g -jar "+PicardtoolsPath+"picard.jar AddOrReplaceReadGroups RGLB="+# Add or replace read group with a single read group in the bam file using Picardtools
					method+"_"+contig_name+"_"+tmp_ref.name+" RGPL="+platform+" RGID=1 RGPU=95 RGSM=L1P1 I="+
					aligner_bam.name+".bam O="+aligner_bam.name+".RG.bam SORT_ORDER=coordinate CREATE_INDEX=FALSE VALIDATION_STRINGENCY=LENIENT;mv "+
					aligner_bam.name+".RG.bam "+aligner_bam.name+".bam",
					samtoolsPath+"samtools index -b "+aligner_bam.name+".bam", # Index the (duplicate removed) bam file
					 FreebayesPath+"freebayes -f "+tmp_ref.name+" -p "+str(ploidy)+" --min-coverage 0 " +
					 " --min-base-quality "+str(3)+" --min-mapping-quality "+str(1)+" --no-indels --no-mnps --no-complex -F "+str(round(1./ploidy*0.4,2))+" "] # Call the variants by freebayes
				errmsg=["Failed to convert sam to bam! Check the help for samtools!",
					"Failed to add RG to the bam file! Check the help for Picardtools!",
					"Failed to index the bam file! Check the help for samtools!",
					"Failed to call the variants! Check freebayes and the bam files!"]
				del platform
				if method in ("CLR"):
					cmd[-1]+="--use-duplicate-reads --theta 0.5 "+ aligner_bam.name+".bam > "+tmp_vcf.name
				else:
					cmd[-1]+=aligner_bam.name+".bam > "+tmp_vcf.name
				if remove_duplicate:
					#tmpfiles.extend([aligner_bam.name+".md.bam", aligner_bam.name+".duplicates"])
					#cmd.insert(2, "java -Xmx8g -jar "+PicardtoolsPath+"picard.jar MarkDuplicates INPUT="+  # Remove duplicate reads from the aligned reads
					#			aligner_bam.name+".bam OUTPUT="+aligner_bam.name+".md.bam METRICS_FILE="+
					#			aligner_bam.name+".duplicates REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT;mv "+
					#			aligner_bam.name+".md.bam "+aligner_bam.name+".bam")
					#errmsg.insert(2, "Failed to remove duplicates from the bam file! Check the help for Picardtools!")
					tmpfiles.extend([aligner_bam.name+".md.bam"])
					cmd.insert(2, samtoolsPath+"samtools rmdup -S "+aligner_bam.name+".bam "+aligner_bam.name+".md.bam;mv "+
								aligner_bam.name+".md.bam "+aligner_bam.name+".bam")
					errmsg.insert(2, "Failed to remove duplicates from the bam file! Check the help for samtools rmdup!") 				
				for _i, _cmd in enumerate(cmd):
					try:
						command=Command(_cmd, outfile=subprocess_err.name)
						exit_msg=command.run(timeout=timeout)
					except subprocess.CalledProcessError:
						if remove_duplicate and _i==2: # If Picardtools fails to remove duplicates, we shall continue without it!
							copyfile(aligner_bam.name+".bam", aligner_bam.name+".md.bam")
							print(errmsg[_i])
							print("Variant calling will be performed without beforehand duplicate removal!")
						else:
							with open(subprocess_err.name, 'rU') as subprocess_err:
								subprocess_err_msg=[_x if ('*' not in _x and _x!='\n') else '' for _x in subprocess_err.readlines()]
								subprocess_err_msg='\n'+'Original error message:'+'\t'+''.join(''.join(str(_y) for _y in _x) for _x in subprocess_err_msg)+'\n'+"-------------------------------"+'\n'
							raise CustomException(errmsg[_i]+subprocess_err_msg)
					else:
						if report_aln:
							if _i==1:
								z.write(aligner_bam.name+".bam", "mybam"+str(_n+1)+"_"+str(_m+1)+".bam")
							elif _i==len(cmd)-2:
								z.write(aligner_bam.name+".bam.bai", "mybam"+str(_n+1)+"_"+str(_m+1)+".bam.bai")
							elif _i==len(cmd)-1:
								z.write(tmp_vcf.name, "myvcf"+str(_n+1)+"_"+str(_m+1)+".vcf")
							else:
								pass
				with NamedTemporaryFile() as hapcompass_vcf, NamedTemporaryFile() as hapcut_vcf, NamedTemporaryFile() as haptree_vcf:
					tmpfiles.extend([hapcompass_vcf.name, hapcut_vcf.name, haptree_vcf.name])  # Adjusted VCF files for the haplotyping software
				with NamedTemporaryFile() as hapcompass_out, NamedTemporaryFile() as sdhap_out: # Temporary haplotype files and folders 
					_compassdir=os.path.realpath(hapcompass_out.name)
					_dir='/'.join(_compassdir.split('/')[:-1])
				haptree_out=tempfile.mkdtemp()
				tmpfiles.extend(sdhap_out.name+_suffix for _suffix in ("", "_Converted"))
				tmpdirs.append(haptree_out)
				with NamedTemporaryFile() as fragment_matrix, NamedTemporaryFile() as fragment_haptree, NamedTemporaryFile() as fragment_sdhap:
					tmpfiles.extend([fragment_matrix.name,fragment_haptree.name,fragment_sdhap.name]) # Temporary fragment matrix for FragmentPoly, fragment file for HapTree and SDhaP 
				cmd=('grep -v "##" '+tmp_vcf.name+' | awk -F" "'+" '"+'BEGIN{OFS="\t";}($4=="REF" \
|| (length($4)==1 && length($5)==1)) {print $1,$2,$3,$4,$5,$6,$7,$8,$9, $10}'+"' >"+hapcompass_vcf.name) # Pre-process the vcf file to throw away any variant other than SNPs!
				if verbose:
					print("VCF file being processed to throw indels and MNPs out...")
				try:
					command=Command(cmd, outfile=subprocess_err.name)
					exit_msg=command.run(timeout=timeout)
				except subprocess.CalledProcessError:      # Exceptions thrown by the haplotyping software are handled so that the simulation pipeline is not interrupted							
					with open(subprocess_err.name, 'rU') as subprocess_err:
						subprocess_err_msg=[_x if ('*' not in _x and _x!='\n') else '' for _x in subprocess_err.readlines()]
						subprocess_err_msg='\n'+'Original error message:'+'\t'+''.join(''.join(str(_y) for _y in _x) for _x in subprocess_err_msg)+'\n'+"-------------------------------"+'\n'
					raise CustomException("Unexpected shell error! Try again later"+subprocess_err_msg)
				with open(hapcompass_vcf.name,'rU') as _vcf1, open(hapcut_vcf.name,'w') as _vcf2, open(haptree_vcf.name,'w') as _vcf3:
					for _line in _vcf1:
						if ('#' not in _line) and (('0/'*ploidy)[:-1] not in _line) and (('1/'*ploidy)[:-1] not in _line):
								_vcf3.write(_line)
								_genotype=_line.rstrip().split()[9]
								_genotype=':'.join(_y if _x>0 else '0/1' for _x, _y in enumerate(_genotype.split(':'))) 
								_vcf2.write('\t'.join(_line.rstrip().split()[:-1])+'\t'+_genotype+'\n')
				_vcf1, _vcf2, _vcf3, _genotype, _line = None, None, None, None, None
				if verbose:
					VCF_called = set()
					Original_var = set()
					with open(haptree_vcf.name,'rU') as _vcf, open(tmp_hap.name+"_varianthaplos.txt",'rU') as _orig, open(errors_in_haplotyping.name,'a') as append_error:
						for _line in _vcf:
							_li = _line.lstrip()
							VCF_called.add(_li.split()[1])
						next(_orig)
						for _line in _orig:
							Original_var.add(_line.split()[2])
						append_error.write("Summary of variant calling (only for heterozygous SNPs) performed by FreeBayes for region {region}, mutation {mut}:\n\tSensitivity = {s:0.3f}\n\
\tPositive predictive value = {ppv:0.3f}\n".format(region=str(_n+1), mut=str(_m+1), s=float(len(Original_var.intersection(VCF_called)))/len(Original_var) if len(Original_var)>0 else float('NaN'), ppv=float(len(Original_var.intersection(VCF_called)))/len(VCF_called) if len(VCF_called)>0 else float('NaN')))
					del Original_var 
					del VCF_called
					_vcf, _orig, append_error = (None, None, None)
				dontrun=[]
				for _algor in ('HapTree', 'SDhaP'):	
					if verbose:
						print("Input being prepared for {0}...".format(_algor))
					if _algor=='HapTree':
						cmd=(HapCutPath+"extractHAIRS --VCF "+hapcut_vcf.name+" --bam "+aligner_bam.name+".bam --maxIS 3000 --out "+
fragment_matrix.name+" --ref "+tmp_ref.name+";FragmentPoly.py -f "+fragment_matrix.name+" -o "+fragment_haptree.name+" -x HapTree")
					else:
						cmd="FragmentPoly.py -f "+fragment_matrix.name+" -o "+fragment_sdhap.name+" -x SDhaP;"#+"cp "+fragment_sdhap.name+" getfragment.frag;" to check the prepared input	
					try:
						command=Command(cmd, outfile=subprocess_err.name)
						exit_msg=command.run(timeout=timeout)
					except subprocess.CalledProcessError:
						dontrun.append('{0}'.format(_algor))
						print >> sys.stderr, "Failed to prepare {0} input...".format(_algor)
						print >> sys.stderr, "Original error message:\n{0}\n{1}\n-------------------------------\n".format(command.std_out, command.std_err)
				cmd=[] # Empty list to later include commands for estimating the haplotypes			
				cmd.append(["java","-Xmx8g","-jar",HapCompassPath+"hapcompass.jar","--ploidy",str(ploidy),"--iterations","19","--bam",
aligner_bam.name+".bam","--vcf",tmp_vcf.name,"-o",hapcompass_out.name]) # Command to estimate the haplotypes by HapCompass
				cmd.append(["HapTree_multiblock.py","-f",fragment_haptree.name,"-v",haptree_vcf.name,"-o",haptree_out,"-t", str(max(1, 0.95*timeout, timeout-50))]) # Command to estimate the haplotypes by HapTree_multiblock
				cmd.append(["hap2","-f",fragment_sdhap.name,"-o",sdhap_out.name,"-p",str(ploidy),"-v",hapcompass_vcf.name]) # Command to estimate the haplotypes by SDhaP
				errmsg=["Failed to run HapCompass! Check HapCompass, VCF and bam files!",
					"Failed to run HapTree! Check HapTree, HapCut, VCF and bam files!",
					"Failed to run SDhaP! Check SDhaP, VCF and bam files!"]
				_haplo_error=False
				out_path=(hapcompass_out.name+"_MWER_solution.txt", haptree_out+"/HapTreeSolution", sdhap_out.name+"_Converted")
				try:
					if timem:
						maxrss=[str(float('NaN'))]*3 	 # list to store the maximum resident set size, i.e. physical memory, used by each haplotyping algorithm 
						algortimes=[str(float('NaN'))]*3 # list to store the total CPU time, i.e. user+system time, used by each haplotyping algorithm
						err_status=['FALSE']*3
					for _i, _cmd in enumerate(cmd):
						try:
							if _i==0:
								if verbose:
									print("HapCompass being run...")
							elif _i==1:
								if 'HapTree' in dontrun:
									raise subprocess.CalledProcessError(-19, 'HapTree', None)
								elif verbose:
									print("HapTree being run...")
							else:
								if 'SDhaP' in dontrun:
									raise subprocess.CalledProcessError(-19, 'SDhaP', None)
								elif verbose:
									print("SDhaP being run...")	
							command=Command(_cmd, shell=False, outfile=subprocess_err.name)
							exit_msg=command.run(timeout=timeout)# command.set_verbose(True) could be used before this line for debugging
							if (not os.path.exists(out_path[_i])) or os.stat(out_path[_i]).st_size==0:
								raise subprocess.CalledProcessError(1, _cmd, None)
						except subprocess.CalledProcessError as e:      # Exceptions thrown by the haplotyping software are handled so that the simulation pipeline is not interrupted
							HapCompass_HapTree_SDhaP_errors[_i]+=1 # Keep track of the number of errors caused by each algorithm
							_haplo_error=True
							if timem:
								err_status[_i]='TRUE' # means the failure of the algorithm
								if e.returncode!=-19:
									maxrss[_i]=('{0:.3f}'.format(round(command.mem.ru_maxrss,3)))
									algortimes[_i]=('{0:.3f}'.format(round(command.mem.ru_utime+command.mem.ru_stime,3)))
							if e.returncode!=-19:
								print >> sys.stderr, errmsg[_i]
								with open(subprocess_err.name, 'rU') as subprocess_err:
									subprocess_err_msg=[_x if ('*' not in _x and _x!='\n') else '' for _x in subprocess_err.readlines()]
									subprocess_err_msg='\n'+'Original error message:'+'\t'+''.join(''.join(str(_y) for _y in _x) for _x in subprocess_err_msg)+'\n'+"-------------------------------"+'\n'
								print >>sys.stderr, subprocess_err_msg
							else:
								print >>sys.stderr, ("ERROR:{0} will not be run as its input was not available in proper format!".format(e.cmd))
							with open(errors_in_haplotyping.name,'a') as append_error:
								append_error.write("\nERROR: an error occured for "+estimate[_i]+
										" at region "+str(_n+1)+", mutant "+str(_m+1)+"!\n")
						else:
							if timem:
								maxrss[_i]=('{0:.3f}'.format(round(command.mem.ru_maxrss,3)))
								algortimes[_i]=('{0:.3f}'.format(round(command.mem.ru_utime+command.mem.ru_stime,3)))
						finally:
							command = None
							del command
				except:
					raise
				finally:
					if timem:
						with open(timem_dat,'a') as _timem_dat:
							for _lst in zip(maxrss, algortimes, ["HapCompass", "HapTree", "SDhaP"], 
									["{0:d}".format(_n+1)]*3, ["{0:d}".format(_m+1)]*3, err_status):
								write_list(_lst, _timem_dat)
					tmpfiles.extend(_dir+"/"+_x for _x in os.listdir(_dir) if _compassdir.split('/')[-1] in _x)
					if (verbose and _haplo_error) or report_aln:
						with open(errors_in_haplotyping.name,'a') as append_error:
							append_error.write("The true haplotypes (generated by haplogenerator), region {0}, mutant {1}:\n".format(_n+1, _m+1))
							if os.path.exists(tmp_hap.name+"_varianthaplos.txt"):
								block_num = 0
								number_of_ones_per_homologue = deque() 
								_nline = 0
								with open(tmp_hap.name+"_varianthaplos.txt",'rU') as _hap_estimate:
									for _line in _hap_estimate:
										append_error.write(_line)
										if not _line.strip():
											if _nline == 0:
												_per_homologue = list(0 for _x in range(0, ploidy))
										elif _line.split()[0]=='Block':
											if _nline == 0:
												_per_homologue = list(0 for _x in range(0, ploidy))
											else:
												_nline=0
											if block_num>0:
												number_of_ones_per_homologue.append(_per_homologue)
												_per_homologue = list(0 for _x in range(0, ploidy))
											block_num+=1
										else:
											_homolo_cols = list(int(_x!='0') for _x in _line.split()[5:])
											_per_homologue = list(_x+_y for _x, _y in zip(_per_homologue, _homolo_cols))
											_nline+=1
									if block_num>0:
										number_of_ones_per_homologue.append(_per_homologue)
										append_error.write('Number of alternative alleles per homolgue in\
 the simulated genome for region {0}, mutation {1}:\n'.format(_n+1, _m+1))
										for _x in range(1, block_num+1):
											append_error.write('\tBlock\t{0}:\t'.format(_x))
											append_error.write(','.join(str(_y) for _y in number_of_ones_per_homologue.popleft())+'\n')	
								_hap_estimate, block_num, number_of_ones_per_homologue, _per_homologue, _nline = None, None, None, None, None
							else:
								append_error.write('\n')
							append_error.write("Haplotypes reported by HapCompass, region {0}, mutant {1}:\n".format(_n+1, _m+1))
							if os.path.exists(hapcompass_out.name+"_MWER_solution.txt"):
								block_num=0
								number_of_ones_per_homologue = deque() 
								_nline = 0
								with open(hapcompass_out.name+"_MWER_solution.txt",'rU') as _hap_estimate:
									for _line in _hap_estimate:
										append_error.write(_line)
										if not _line.strip():
											pass
										elif _line.split()[0]=='BLOCK':
											if _nline == 0:
												_per_homologue = list(0 for _x in range(0, ploidy))
											else:
												_nline=0
											if block_num>0:
												number_of_ones_per_homologue.append(_per_homologue)
												_per_homologue = list(0 for _x in range(0, ploidy))
											block_num+=1
										else:
											_homolo_cols = list(int(_x!='0') for _x in _line.split()[3:])
											_per_homologue = list(_x+_y for _x, _y in zip(_per_homologue, _homolo_cols))
											_nline+=1
									if block_num>0:
										number_of_ones_per_homologue.append(_per_homologue)
										append_error.write('Number of alternative alleles per homolgue in\
 the HapCompass estimate for region {0}, mutation {1}:\n'.format(_n+1, _m+1))
										for _x in range(1, block_num+1):
											append_error.write('\tBlock\t{0}:\t'.format(_x))
											append_error.write(','.join(str(_y) for _y in number_of_ones_per_homologue.popleft())+'\n')	
								_hap_estimate, block_num, number_of_ones_per_homologue, _per_homologue, _nline = None, None, None, None, None
							else:
								append_error.write('\n')
							append_error.write("Haplotypes reported by HapTree, region {0}, mutant {1}:\n".format(_n+1, _m+1))
							if os.path.exists(haptree_out+"/HapTreeSolution"):
								block_num=0
								number_of_ones_per_homologue = deque() 
								_nline = 0
								with open(haptree_out+"/HapTreeSolution",'rU') as _hap_estimate: 
									for _line in _hap_estimate:
										append_error.write(_line)									
										if not _line.strip():
											pass
										elif _line.split()[0]=='BLOCK':
											if _nline == 0:
												_per_homologue = list(0 for _x in range(0, ploidy))
											else:
												_nline=0
											if block_num>0:
												number_of_ones_per_homologue.append(_per_homologue)
												_per_homologue = list(0 for _x in range(0, ploidy))
											block_num+=1
										else:
											_homolo_cols = list(int(_x!='0') for _x in _line.split()[2:])
											_per_homologue = list(_x+_y for _x, _y in zip(_per_homologue, _homolo_cols))
											_nline+=1
									if block_num>0:
										number_of_ones_per_homologue.append(_per_homologue)
										append_error.write('Number of alternative alleles per homolgue in\
 the HapTree estimate for region {0}, mutation {1}:\n'.format(_n+1, _m+1))
										for _x in range(1, block_num+1):
											append_error.write('\tBlock\t{0}:\t'.format(_x))
											append_error.write(','.join(str(_y) for _y in number_of_ones_per_homologue.popleft())+'\n')	
								_hap_estimate, block_num, number_of_ones_per_homologue, _per_homologue, _nline = None, None, None, None, None
							else:
								append_error.write('\n')
							append_error.write("Haplotypes reported by SDhaP:, region {0}, mutant {1}:\n".format(_n+1, _m+1))
							if os.path.exists(sdhap_out.name+"_Converted"):
								block_num=0
								number_of_ones_per_homologue = deque() 
								_nline = 0
								with open(sdhap_out.name+"_Converted",'rU') as _hap_estimate: 
									for _line in _hap_estimate:
										append_error.write(_line)									
										if not _line.strip():
											pass
										elif _line.split()[0]=='Block':
											if _nline == 0:
												_per_homologue = list(0 for _x in range(0, ploidy))
											else:
												_nline=0
											if block_num>0:
												number_of_ones_per_homologue.append(_per_homologue)
												_per_homologue = list(0 for _x in range(0, ploidy))
											block_num+=1
										else:
											_homolo_cols = list(int(_x!='0') for _x in _line.split()[2:])
											_per_homologue = list(_x+_y for _x, _y in zip(_per_homologue, _homolo_cols))
											_nline+=1
									if block_num>0:
										number_of_ones_per_homologue.append(_per_homologue)
										append_error.write('Number of alternative alleles per homolgue in\
 the SDhaP estimate for region {0}, mutation {1}:\n'.format(_n+1, _m+1))
										for _x in range(1, block_num+1):
											append_error.write('\tBlock\t{0}:\t'.format(_x))
											append_error.write(','.join(str(_y) for _y in number_of_ones_per_homologue.popleft())+'\n')	
								_hap_estimate, block_num, number_of_ones_per_homologue, _per_homologue, _nline = None, None, None, None, None
							else:
								append_error.write('\n')
				if _n==0 and _m==0:
					with NamedTemporaryFile(delete=False) as compare_hapcompass, NamedTemporaryFile(delete=False) as compare_HapTree, NamedTemporaryFile(delete=False) as compare_SDhaP:
						method_report=(compare_hapcompass.name, compare_HapTree.name, compare_SDhaP.name) # Files containing haplotyping evaluation results 
					with NamedTemporaryFile() as compare:
						pass
					dummies=[]
					for _i in (0,1,2):
						with NamedTemporaryFile() as dummy:
							dummies.append(dummy.name)
				else:
					pass
				valid_solutions=[]	# The list of valid phasing solutions (maximum 3 equal to the number of estimation algorithms)
				if no_missing or (not _haplo_error): # Only compare the haplotypes if all the three haplotyping algorithms have produced estimates or if no_missing is set!
					_compare_error=False 	# Comparison results only used if all the three haplotype estimates are have been successfully evaluated by hapcomapre or if no_missing is set!
					if verbose:
						print("Estimated haplotypes being compared to the original...")
					for _i, _estim in enumerate(out_path):
						if ploidy>8:
							if verbose:
								print >>sys.stderr, "WARNING: the given ploidy level is greater than 8! HapCompare will use the quick method to calculate Allelic Correlation!"
							cmd = Command("hapcompare.py "+tmp_hap.name+"_varianthaplos.txt "+_estim+" -t -q > "+compare.name, outfile=subprocess_err.name)  # Finally, Comapre the estimated haplotype with the true haplotype! (Ploidy>8) 
						else:
							cmd = Command("hapcompare.py "+tmp_hap.name+"_varianthaplos.txt "+_estim+" -t > "+compare.name, outfile=subprocess_err.name)  # Finally, Comapre the estimated haplotype with the true haplotype! (Ploidy<=8) 
						try:
							exit_msg = cmd.run(timeout=timeout)
						except subprocess.CalledProcessError:
							with open(subprocess_err.name, 'rU') as subprocess_err:
								subprocess_err_msg=[_x if ('*' not in _x and _x!='\n') else '' for _x in subprocess_err.readlines()]
								subprocess_err_msg='\n'+'Original error message:'+'\t'+''.join(''.join(str(_y) for _y in _x) for _x in subprocess_err_msg)+'\n'+"-------------------------------"+'\n'
							print >>sys.stderr, "Failed to run Hapcompare!", subprocess_err_msg
							_compare_error=True
							HapCompass_HapTree_SDhaP_hapcompare_errors[_i]+=1
							with open(errors_in_haplotyping.name,'a') as append_error:
								append_error.write("\nERROR: hapcompare failed to evaluate "+
									estimate[_i]+" estimate for region "+str(_n+1)+" and mutant "+
									str(_m+1)+"!\n") # Reported numbers are extracted from compare file, added column-wise to the report file of each method
								if verbose and (not report_aln): # report_aln = True means that the haplotypes are already written to the log file anyway!
									append_error.write(" The failed estimate:\n")
									if os.path.exists(_estim):
										with open(_estim,'rU') as _hap_estimate:
											for _line in _hap_estimate:
												append_error.write(_line)
										_hap_estimate=None
									else:
										append_error.write('\n')
									append_error.write(" The true haplotypes:\n")
									if os.path.exists(tmp_hap.name+"_varianthaplos.txt"):
										with open(tmp_hap.name+"_varianthaplos.txt",'rU') as _hap_estimate:
											for _line in _hap_estimate:
												append_error.write(_line)
										_hap_estimate=None
									else:
										append_error.write('\n')					
						else:
							with NamedTemporaryFile(delete=True) as compare2:
								tmpfiles.append(compare2.name)
							with open(compare.name,'rU') as _compare, open(compare2.name,'w') as _compare2:
								for _line in _compare:
									_compare2.write(substitute('^.*:','', _line).lstrip()) # Each row will give the values of one of the comparison measures for all of the simulations
								_compare2.write("{0}\n{1}\n".format(_n+1, _m+1)) # Add the id of the reference genomic region, as well as the id of the mutation to quality measures
							move(compare2.name, compare.name)
							cmd=Command("paste -d' ' {0} {1} > {2}".format(compare.name, method_report[_i], dummies[_i])) # Add the processed hapcomapre output to the final report file for each method  
							try:
								exit_msg=cmd.run(timeout=timeout)
							except subprocess.CalledProcessError:
								exit_msg=2
								raise RuntimeError("Unexpected Shell Error!")
							else:
								valid_solutions.append(_i)
					if no_missing or (not _compare_error): # Write the comparisons to file if no error has occured during evaluation or if no_missing is set!
						for _i in valid_solutions:
							move(dummies[_i], method_report[_i])
					else:
						pass
				clean_up_files(tmpfiles)
				clean_up_folders(tmpdirs)
			tmpfiles.extend((tmp_ref.name, tmp_ref.name+'.fai', tmp_ref.name+'.fasta', tmp_ref.name+'.fasta.fai'))
			if method in ("CLR", "CCS") and bwa_mem:
				tmpfiles.extend((tmp_ref.name+".fasta.amb", tmp_ref.name+".fasta.ann", 
					tmp_ref.name+".fasta.bwt", tmp_ref.name+".fasta.pac", tmp_ref.name+".fasta.sa"))
			elif bwa_mem:
				tmpfiles.extend((tmp_ref.name+".amb", tmp_ref.name+".ann", 
					tmp_ref.name+".bwt", tmp_ref.name+".pac", tmp_ref.name+".sa"))
			tmpfiles.append(compare.name)
			tmpfiles.extend(dummies)
			clean_up_files(tmpfiles)
			clean_up_folders(tmpdirs)
			tmpfiles=[]
			tmpdirs=[]
	except CustomException as err:
		print >> sys.stderr, err
		exit_msg=2
	except IOError as e:											
		print >> sys.stderr, "I/O error({0}): {1} {2}".format(e.errno, e.strerror, e.filename)
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
		traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stderr)
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
	else:
		try:
			final_files=list(output+"_"+_suffix for _suffix in ("HapCompass.dat", "HapTree.dat","SDhaP.dat"))
			header=("CORR_ALL_DOSAGES","PSR_ALL_DOSAGES","VER","PSR_CORRECT_DOSAGES",
				"CORR_CORRECT_DOSAGES","MISSING_ESTIMATE","SPURIOUS_ESTIMATE","WRONG_DOSAGES_COMMON_SNPS",
				"PPV_ESTIMATE","NUM_GAPS", "CORRECT_HOMOLOGUE_RATE", "ID_REGION", "ID_MUTATION")
			for _i, _file in enumerate(method_report):
				if os.path.exists(_file):
					with open(_file,'rU') as matrix, open(final_files[_i],'w') as invert: 
						write_list(header, invert) # write the header of the final report for each method
						lis = [line.split() for line in matrix]
						for x in zip(*lis): # transpose the values written to the dummy report of each method, and write it to the final report 
							write_list(x,invert)
				else:
					with open(final_files[_i],'w') as dummy: # Report an empty file
						if verbose:
							print >> sys.stderr, "WARNING: Error in writing {0}.".format(final_files[-i])
						else:
							pass
			with open(errors_in_haplotyping.name,'a') as append_error: # Generate haplosim log file
				end_time_cpu = os.times()
				observed_dosage=list('{0:.3f}'.format(_x*100/float(sum(observed_dosage))) for _x in observed_dosage)
				append_error.write("\nDistribution of simplex, ..., ploidy-plex alleles: "+'\t'.join(str(_x+1)+':'+_y+'%' for _x, _y in enumerate(observed_dosage))+"\n")
				append_error.write("\nTotal CPU time used by haplosim: %s \n" % round(sum(y-x for y, x in zip(end_time_cpu[:-1], start_time_cpu[:-1])), 4))
				append_error.write("\nSummary of haplotyping errors:\n")
				for _i in (0,1,2):
					append_error.write("\t"+str(HapCompass_HapTree_SDhaP_errors[_i])+" exceptions thrown by "+estimate[_i]+".\n")
				if no_missing:
					append_error.write("\nSummary of haplotype evaluation errors (including invalid phasing solutions):\n")
				else:
					append_error.write("\nSummary of haplotype evaluation errors (excluding invalid phasing solutions):\n")
				for _i in (0,1,2):
					append_error.write("\t"+str(HapCompass_HapTree_SDhaP_hapcompare_errors[_i])+
							" exceptions thrown by hapcompare for "+estimate[_i]+" estimates.\n")
			if verbose:
				print >>sys.stdout, "Saving the log report to {0}.log...".format(output)
			move(errors_in_haplotyping.name, output+".log")
			if plotting:
				if verbose:
					print >> sys.stdout, "Saving the variant density histograms to {0}_VariantsHistograms.pdf...".format(output)
				pdfpages.savefig()
				pdfpages.close()
				move(pdffile.name, output+"_VariantsHistograms.pdf")
			exit_msg=0
		except:
			exit_msg=2
			raise
	finally:
		if report_aln:
			if verbose and zip_arc:
				print >>sys.stdout, "Closing the zip archive {0}...".format(zip_arc)
			if z:
				z.close()
		if plotting:
			try:
				pdfpages.close()
			except AttributeError: # Exception only should occur if pdfpages is already closed or has not been created at the first place
				pass
			if pdffile is not None:
				tmpfiles.append(pdffile.name)
		tmpfiles.extend(method_report)
		clean_up_files(tmpfiles)
		if tmp_ref:
			if os.path.exists(tmp_ref.name):
				os.remove(tmp_ref.name)
			if os.path.exists(tmp_ref.name+'.fai'):
				os.remove(tmp_ref.name+'.fai')
			if os.path.exists(tmp_ref.name+'.fasta'):
				os.remove(tmp_ref.name+'.fasta')
			if os.path.exists(tmp_ref.name+'.fasta.fai'):
				os.remove(tmp_ref.name+'.fasta.fai')
			for _indx_bwa_mem_pacbio in (tmp_ref.name+".fasta.amb", tmp_ref.name+".fasta.ann", 
					tmp_ref.name+".fasta.bwt", tmp_ref.name+".fasta.pac", tmp_ref.name+".fasta.sa"):
				if os.path.exists(_indx_bwa_mem_pacbio):
					os.remove(_indx_bwa_mem_pacbio)
			for _indx_bwa_mem in (tmp_ref.name+".amb", tmp_ref.name+".ann", 
								  tmp_ref.name+".bwt", tmp_ref.name+".pac", tmp_ref.name+".sa"):
				if os.path.exists(_indx_bwa_mem):
					os.remove(_indx_bwa_mem)
		if pophap:
			pophap.close()
			os.remove(pophap.filename)
		elif pophapname:
			os.remove(pophapname)
		if errors_in_haplotyping:
			if os.path.exists(errors_in_haplotyping.name):
				os.remove(errors_in_haplotyping.name)
		clean_up_folders(tmpdirs)
		sys.exit(exit_msg)
