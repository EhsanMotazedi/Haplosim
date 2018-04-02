#!/usr/bin/env python
"""A wrapper for HapTree. Takes the input fragment file and splits the fragments into connected sets. 
Each set is given sperately (and in parallel using multiprocessing module) to HapTree to estimate its 
haplotype block, and the blocks are finally written to a single output file. 
Written by Ehsan Motazedi, 10-12-2015, Wageningen UR
Last modified: 03-06-2016 """
import argparse
import ast
import copy
import multiprocessing as mpthread
import networkx as nx
import os 
import os.path
import Queue 
import re
import signal
import subprocess
import sys
sys.path.append('/Users/ehsan/Desktop/Python Scripten')
import tempfile
import threading
import time
import traceback
from functools import partial
from itertools import chain
from shutil import copyfile, rmtree, move

def check_pid(pid):
	""" Checks if a unix process exists or not."""
	try:
		os.kill(pid,0)
	except OSError:
		return False
	else:
		return True

def check_positive_float(value):
	""" Checks if the passed argument is a positive real number."""
	fvalue=None
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

class Command(object):
	""" Handles the spawning of UNIX command line processes."""
	
	def __init__(self, cmd, outfile = None):
		self.cmd = cmd
		self.process = None
		self.outfile = outfile
		self.std_out = None
		self.std_err = None
		self.cpid = None
		self.mem = None 
	
	def run(self):	

		def reader(iostream, outlist, chunk=64):
			"""function to read the I/O byte stream and append it to outlist in blocks of size chunk.
			   The I/O pipe is closed when read."""
			read_block = partial(iostream.read, chunk) #at most 64 bytes are read and stored in each block
			for _block in iter(read_block, b''):
				outlist.append(_block)
			iostream.close()

		try:
			self.std_err, self.std_out = [], []
			self.process = subprocess.Popen(self.cmd, shell=False, close_fds=False,
				preexec_fn = preexec, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			reader1 = threading.Thread(target=reader, args=(self.process.stdout, self.std_out))
			reader2 = threading.Thread(target=reader, args=(self.process.stderr, self.std_err))
			for _reader in (reader1, reader2):
				_reader.start()
			for _reader in (reader1, reader2):
				_reader.join()
			_pid, self.process.returncode, self.mem = os.wait4(self.process.pid, 0) # used instead of the standard (self.std_out, self.std_err) = self.process.communicate(), as wait4() reports resource use
		except:
			if self.process:
				if check_pid(self.process.pid):
					os.kill(self.process.pid, signal.SIGTERM)
				print >> sys.stderr, b''.join(self.std_err).decode('utf-8')	
			raise subprocess.CalledProcessError(1, self.cmd, None)
		finally:
			self.std_err = b''.join(self.std_err).decode('utf-8')
			self.std_out = b''.join(self.std_out).decode('utf-8')
			if self.outfile is not None:
				with open(self.outfile,'w') as fileobject:
					fileobject.write(self.std_out)
					fileobject.write(self.std_err)

		return(self.process.returncode)
				
def component_thread(indx, reads, vcf, outp, error, timeout, lock, sema, qmem, q):
	""" Each component will be phased in parallel by HapTree_v0, using multiple threads. """
	try:
		sema.acquire()
		cmd = Command(cmd=[HapTreePath+'HapTree_v0',reads, vcf, outp], outfile = error)
		t = threading.Thread(target = single_thread, args = (cmd, ))
		t.start()
		t.join(timeout)
		if t.is_alive():
			try:
				os.kill(cmd.process.pid, signal.SIGTERM)
			except OSError as e:
				if e.errno==12:
					raise subprocess.CalledProcessError(1, cmd.cmd, None)
				else:	
					pass # The process finished between the 'is_alive()' and 'terminate()'
			except AttributeError as e:
				pass # The process finished between the 'is_alive()' and 'terminate()'
			else:
				print >> sys.stderr, 'Process #%d was not responding! Timeout termination occured after %f seconds!' %(cmd.process.pid, round(float(timeout), 3))
				raise subprocess.CalledProcessError(1, cmd.cmd , None)
			finally:
				t.join()
		if os.stat(error).st_size > 0:
			lock.acquire()
			print >>sys.stderr, 'WARNING:'
			with open(error, 'rU') as _err:
				print >>sys.stderr,  _err.read() # Print the error message reported for the component by HapTree_v0
			_err = None
			lock.release()
			raise subprocess.CalledProcessError(1, cmd , None)
	except subprocess.CalledProcessError:
		q.put(['False', indx])
	else:
		q.put(['True', -1]) # This means that at least one block has been estimated and should be reported in output
	finally:
		qmem.put(cmd.mem)
		sema.release()
				
class CustomException(Exception):
	def __init__(self, value):
		self.parameter = value
	def __str__(self):
		return repr(self.parameter)

def preexec():
	signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def single_thread(command):
	""" Single thread for a connected SNP graph. """
	command.run()

if __name__ == "__main__":	  # <<The __main__ module starts here!>>
	try:
		qmem = None           # Queue to store memory and CPU time report for each call to HapTree
		processfiles = []     # The list of the name of the files for each individual component of SNP matrix, if multiple components exist 
		exit_msg = 0		
		vcf = None        	  # The VCF file, needed for HapTree
		fragmentfile = None   # The fragment file for HapTree
		output = None         # The output folder name for HapTree
		frag = None           # The input stream of the fragments
		haptree_error_files = [] # List of the error files saving the messages reported by HapTree to sys.stderr
		tmp1, tmp2, tmp3, tmp4 = (None, None, None, None)  # Temporary files for extracting connected fragments
		tmpvcf = None
		_outputs = []			# Temporary folders that store the haplotype block of connected components
		tmpoutput = ''		    # Temporary folder to combine separate HapTree outputs, which belong to different blocks 
		MAXCORES = mpthread.cpu_count()-1 # Max number of available cores
		HapTreePath=os.getenv("HapTreePath")
		NCORES = os.getenv("NCORES") # Desired number of CPU cores to be used. MUST be exported by OS, otherwise is set to one (no parallel execution!) 
		if not NCORES:
			NCORES = 1
			garbage = sys.stderr.write("WARNING: \"NCORES\" not found in the OS-environment! HapTree_multiblock will be run on just one CPU-core!")
			os.environ["NCORES"] = str(NCORES)
		else:
			NCORES = int(NCORES)	
		if NCORES<1:
			raise CustomException("ERROR: Number of cores given to run \"HapTree_multiblock\" must be at least 1!") 
		if not HapTreePath:
			HapTreePath=''
		else:
			HapTreePath=HapTreePath.rstrip('/')+'/'
		Parallel = False # Only gets true if more than one cores are desired as well as available 
		parser = argparse.ArgumentParser(description="A Wrapper for HapTree (Berger et al. 2014). Takes the same inputs as HapTree, and\
 splits the input fragment file into connected sets. Each connected set's hapltype block is separately given to HapTree,\
 and the blocks are finally written in a single output file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
		epilog='In case multiple cores are to be used, their desired number should be set under the name NCORES in the OS environment.\
 Maximum allowed value is the total number of CPU cores minus one. Otherwise, just one core will be used!' )
		parser.add_argument('-f','--fragment', dest='fragmentfile', metavar='Fragment file', type=str, nargs=1,
					   help='The input fragment file in the same format as for HapTree', default=None)
		parser.add_argument('-v','--vcf', dest='vcf', metavar='VCF file', type=str, nargs=1,
					   help='The VCF file related to the fragment file', default=None)
		parser.add_argument('-o','--output', dest='output', metavar='Output folder', type=str, nargs=1,
					   help='The name of output folder to write HapTree output', default=None)
		parser.add_argument('-t','--timeout', dest='timeout', metavar='Time limit', type=check_positive_float, nargs=1,
					   help='The time limit in seconds after which the algorithm is halted by external signal if still running.', default='500')
		try:
			if not len(sys.argv)>1:
				parser.parse_args('-h'.split())
			else:
				args = parser.parse_args() 
		except SystemExit:
			exit_msg=3
			raise
		if not args.fragmentfile:
			raise CustomException("ERROR: no fragment file has been specified!")
		elif not args.vcf:
			raise CustomException("ERROR: no vcf file has been specified!")
		elif not args.output:
			raise CustomException("ERROR: no output name has been specified!")
		else:
			fragmentfile=args.fragmentfile[0]
			vcf=args.vcf[0]
			output=args.output[0]
			timeout=args.timeout if isinstance(args.timeout, float) else args.timeout[0]
		if not fragmentfile:
			raise CustomException("ERROR: no fragment file has been specified!")
		elif not vcf:
			raise CustomException("ERROR: no vcf file has been specified!")
		elif not output:
			raise CustomException("ERROR: no output name has been specified!")
		else:
			pass
		for _file in (fragmentfile, vcf): 
			if not os.path.exists(_file):
				raise IOError(2,"No such file or directory",_file)
		del _file
		Read_Graph=nx.Graph() # The SNP-fragment graph to be grown from the input fragment file
		_fragdict, _edges, _ls = (dict(), [], []) # Dictionary and lists to add each fragment to the SNP-fragment graph
		with open(fragmentfile,'rU') as frag:	# Add the fragment variants to the graph nodes, and an edge between every to variants that lie on the same fragment
			for _frag in frag:
				_fragdict=ast.literal_eval(_frag)
				Read_Graph.add_nodes_from(_fragdict.keys()) # Adding the fragment variants as node to the graph
				_ls=list((x,y) for x in _fragdict.keys() for y in _fragdict.keys()) # Considering an edge between each two variants within the same fragment
				_edges=list(edge for edge in _ls if edge[0]!=edge[1]) # Adding the edges to the graph
				Read_Graph.add_edges_from(_edges)
		_edges, _ls, _fragdict = (None, None, None)
		if len(Read_Graph.nodes())<2:
			try:
				raise IOError()
			except IOError as err:
				err.errno=3
				err.strerror='The input fragment file, {0}, could not be used to build the Read Graph!'.format(fragmentfile)
				err.filename=''
				raise
		else:
			tmpoutput=tempfile.mkdtemp()
		if nx.is_connected(Read_Graph):
			try:
				os.removedirs(tmpoutput)
				included_var=sorted(Read_Graph.nodes()) # Filter the VCF file to include only the variants which are present in the fragment matrix
				with open(vcf,'rU') as _vcf, tempfile.NamedTemporaryFile(delete=False) as tmpvcf:
					_num = -1
					for _var in _vcf:
						if _var.lstrip()[0]=='#': # do not go through the header lines
							pass
						else:
							_num+=1
							if _num in included_var:
								tmpvcf.write(_var)
							else:
								pass
				_vcf, _num, _var = (None, None, None) # Correct the SNP numbers, so that they are mapped to 0,1,2,... in the ascending order
				with open(fragmentfile,'rU') as reads, tempfile.NamedTemporaryFile(delete=False) as tmp1:
					correction_dict=dict()
					for _i, _snp in enumerate(included_var):
						correction_dict[int(_snp)]=_i
					for _frag in reads:
						_fragdict=ast.literal_eval(_frag)
						_fragdict_correct=dict()
						for _item in _fragdict.iteritems():
							_fragdict_correct.update({correction_dict[int(_item[0])]:_item[1]})
						tmp1.write(str(_fragdict_correct)+'\n')
				included_var=None
				with tempfile.NamedTemporaryFile() as haptree_error: # IN SPECIFIC CASES, HapTree_v0 does not throw exceptions. It reports instead an error message which will be saved in this error file.
					haptree_error_files.append(haptree_error.name)   # In case error messages have been reported by HapTree, an exception will be thrown by the spawning process.
				cmd = Command(cmd=[HapTreePath+'HapTree_v0', tmp1.name, tmpvcf.name, tmpoutput], outfile = haptree_error.name)
				t = threading.Thread(target = single_thread, args = (cmd, ))
				t.start()
				time.sleep(0.1)
				t.join(timeout)
				if t.is_alive():
					try:
						os.kill(cmd.process.pid, signal.SIGTERM)
					except OSError as e:
						if e.errno==12:
							raise subprocess.CalledProcessError(1, cmd.cmd, None)
						else:	
							pass # The process finished between the 'is_alive()' and 'terminate()'
					except AttributeError as e:
						pass # The process finished between the 'is_alive()' and 'terminate()'
					else:
						print >> sys.stderr, 'Process #%d was not responding! Timeout termination occured after %f seconds!' %(cmd.process.pid, round(float(timeout), 3))
						raise subprocess.CalledProcessError(1, cmd.cmd , None)
					finally:
						t.join()
				if os.stat(haptree_error.name).st_size > 0:
					raise subprocess.CalledProcessError(1, cmd.cmd , None)
			except subprocess.CalledProcessError:
				exit_msg=2
				with open(haptree_error.name, 'rU') as _err:
					print >> sys.stderr, _err.read()
				_err = None
				raise CustomException("Failed to call HapTree! Check the help for HapTree!")
			finally:
				qmem=Queue.Queue()
				qmem.put(cmd.mem)
		else:
			SNP_frag_components = sorted(nx.connected_components(Read_Graph), key=lambda x: min(x))
			N_components=nx.number_connected_components(Read_Graph)
			_outputs=list(tempfile.mkdtemp() for x in range(1, N_components+1))
			good_outputs=copy.deepcopy(_outputs)
			for _dir in _outputs:
				os.removedirs(_dir)
			with tempfile.NamedTemporaryFile(delete=False) as tmp1, tempfile.NamedTemporaryFile(delete=False) as tmp2, tempfile.NamedTemporaryFile(delete=False) as tmp3, tempfile.NamedTemporaryFile(delete=False) as tmp4:
					pass
			copyfile(frag.name, tmp1.name)
			for _i in range(0, N_components):
				with tempfile.NamedTemporaryFile(delete=False) as _file1, tempfile.NamedTemporaryFile(delete=False) as _file2:
					processfiles.append((_file1.name, _file2.name))
			n_component = 0
			while n_component < N_components:
				connected_reads, remain_reads, vcf_reads, vcf_connected, haptree_error = (None, None, None, None, None) 
				component = set(SNP_frag_components[n_component])
				correction_dict=dict()
				for _i, _snp in enumerate(sorted(component)):
					correction_dict[int(_snp)]=_i
				with open(tmp1.name, 'rU') as reads, open(tmp2.name,'w') as connected_reads, open(tmp3.name,'w') as remain_reads:
					for _frag in reads:
						_fragdict = ast.literal_eval(_frag)
						if set(_fragdict.keys()).issubset(component):
							_fragdict_correct=dict()
							for _item in _fragdict.iteritems():
								_fragdict_correct.update({correction_dict[int(_item[0])]:_item[1]})
							connected_reads.write(str(_fragdict_correct)+'\n')
						else:
							remain_reads.write(_frag)
				_num, _variant= (None, None)
				with open(vcf, 'rU') as vcf_reads, open(tmp4.name, 'w') as vcf_connected:
					_num=-1	
					for _variant in vcf_reads:
						if _variant.lstrip()[0]=='#': # do not go through the header lines
                                                        pass
                                                else:
                                                        _num+=1
						if _num in component:
							vcf_connected.write(_variant)
						else:
							pass	 
				with tempfile.NamedTemporaryFile() as haptree_error:
					haptree_error_files.append(haptree_error.name)
				copyfile(connected_reads.name, processfiles[n_component][0])
				copyfile(vcf_connected.name, processfiles[n_component][1])
				move(tmp3.name, tmp1.name) 
				n_component+=1
			component_mpthreads = []
			partial_solution = []
			remove_lst = []
			if min(NCORES, MAXCORES)>1:
				sema = mpthread.BoundedSemaphore(min(NCORES, MAXCORES))
				lock = mpthread.Lock()
				q = mpthread.Queue()  # It is possible that HapTree fails for some blocks. The other blocks, however must be estimated!	q gather the success info for each block!
				qmem = mpthread.Queue() # The queue to collect memory info for each call to HapTree
				Parallel = True
			else:	# If parallel execution is not possible (as just one core could be used), concurrent execution is performed using threads.
				sema = threading.BoundedSemaphore(NCORES)
				lock = threading.Lock()
				q = Queue.Queue()	# It is possible that HapTree fails for some blocks. The other blocks, however must be estimated! q gather the success info for each block!
				qmem = Queue.Queue() # The queue to collect memory info for each call to HapTree
			for n_component in range(0, N_components):
				if Parallel:
					t = mpthread.Process(target = component_thread,
										 args = (n_component, processfiles[n_component][0], processfiles[n_component][1], _outputs[n_component],
										   haptree_error_files[n_component], timeout, lock, sema, qmem, q))
				else:
					t = threading.Thread(target = component_thread,
										 args = (n_component, processfiles[n_component][0], processfiles[n_component][1], _outputs[n_component],
										   haptree_error_files[n_component], timeout, lock, sema, qmem, q))
	
				component_mpthreads.append(t)
				t.start()
			time.sleep(0.1)
			for x in component_mpthreads: # Wait until all the components have been phased, i.e. their thread is phinished
				x.join()
			while not q.empty():
				_report=q.get()
				remove_lst.append(_report[1])
				partial_solution.append(_report[0])
			del _report
			if 'True' not in partial_solution:
				exit_msg=2
				print >>sys.stderr, "ERROR: HapTree reported errors for all of the SNP matrix components!\n\
No solution could be obtained!\n\
The log report of HapTree:"
				for _file in haptree_error_files:
					with open(_file, 'rU') as _err:
						print >> sys.stderr, _err.read()
				_file, _err = (None, None)
				raise CustomException("Failed to call HapTree! Check the help for HapTree!") 
			else:
				for _i in remove_lst:
					if _i!=-1:
						good_outputs.remove(_outputs[_i])
					else:
						pass
			try: # Now that all components have been tried for phasing, report the available phasings in a single file.
				cmd="cat "+' '.join(x+"/HapTreeSolution" for  x in good_outputs)+" > "+tmpoutput+"/HapTreeSolution"
				exit_msg=subprocess.check_call(cmd, shell=True)
				cmd="cat "+' '.join(x+"/MECScore" for  x in good_outputs)+" > "+tmpoutput+"/MECScore"
				exit_msg=subprocess.check_call(cmd, shell=True)
			except subprocess.CalledProcessError:
				exit_msg=2
				raise RuntimeError('Unexpected Shell Error Occured!')
			with open(tmpoutput+"/HapTreeSolution", 'rU') as _haptree_solution, open(tmp1.name, 'w') as _haptree_solution_header_corrected:
				_size_of_previous_blocks=0
				for _num, _frag in enumerate(_haptree_solution):
					if _frag=="\n":
						_haptree_solution_header_corrected.write(_frag)
					else:
						columns=re.split('\s*|\t',_frag.rstrip())
						if columns[0]=='BLOCK': # Correct the number of the first and last variants in the header of each block
							columns[2]=str(int(columns[2])+_size_of_previous_blocks)
							columns[3]=str(int(columns[3])+_size_of_previous_blocks)
							_size_of_previous_blocks=int(columns[3])
						else:
							pass
						_haptree_solution_header_corrected.write('\t'.join(columns)+'\n')
			os.remove(tmpoutput+"/HapTreeSolution")
			copyfile(tmp1.name, tmpoutput+"/HapTreeSolution")
	except CustomException as err:
		print >> sys.stderr, err
		exit_msg=2
	except IOError as e:											
		print >> sys.stderr,"I/O error({0}): {1} {2}".format(e.errno, e.strerror,e.filename)
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
			if not os.path.exists(output):
				os.makedirs(output)
			else:
				pass
			copyfile(tmpoutput+"/HapTreeSolution", output+"/HapTreeSolution")
			copyfile(tmpoutput+"/MECScore", output+"/MECScore")
			exit_msg=0
		except:
			exit_msg=2
			raise
	finally:
		if qmem and (not qmem.empty()):
			sys.stdout.write("\nRESOURCE USAGE REPORT:\n")
			header='utime,stime,maxrss,ixrss,idrss,isrss,minflt,majflt,\
nswap,inblock,oublock,msgsnd,msgrcv,nsignals,nvcsw,nivcsw'.split(',')
			header=['ru_'+_x for _x in header]
			sys.stdout.write('\t'.join(_x for _x in header)+'\n')
			while not qmem.empty():
				_mem=qmem.get()
				sys.stdout.write('\t'.join('{0:<.9f}'.format(_val) for _val in _mem[:])+'\n')
		if processfiles:
			for _file in chain(*processfiles):
				if os.path.isfile(_file):
					os.remove(_file)
		if haptree_error_files:
			for _file in haptree_error_files:
				if os.path.isfile(_file):
					os.remove(_file)
		for _tmpout in _outputs:
			if os.path.exists(_tmpout):
				rmtree(_tmpout)
		if os.path.exists(tmpoutput):
			rmtree(tmpoutput)
		for _tmp in (tmp1, tmp2, tmp3, tmp4, tmpvcf):
			if _tmp:
				if os.path.isfile(_tmp.name):
					os.remove(_tmp.name)	
		sys.exit(exit_msg)
