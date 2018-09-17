import os,sys
from os import path

#----------------------------- Colors For Print Statements ------------------------------#
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   ORANGE = '\033[38;5;214m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

#------------------------------ Checks the Input Arguments ------------------------------#


if len(sys.argv) != 4:
	print color.BOLD + '\n\nThis script will just loop the actual scripts 1 to 7 for you. Please make sure you have gone through the individual scripts a couple of times\n'
	print color.BOLD + '\tand understand what is going on\n\n'
	print color.BOLD + '\n\nDouble check that you have added all the necessary command-line inputs! (see usage below for an example)\n\n'
	print  color.RED + 'Example Usage:\n\n\t' + color.CYAN + 'katzlab$ python PhyloTOL_postAssembling_LKH.py LKH 200 Table\n\n' + color.END	
	print color.BOLD + 'The current script controls for names, for example, LKH is the shared starting of each file/folder, and 200 the minimum transcript cutoff, a table with code, PhyloToL code and genetic code\n\n'
	sys.exit()
	
### Table: LKH00000 \t Op_me_Jd01 \t Universal	
	
else:
	code = sys.argv[1]
	seqlength = sys.argv[2]
	table = sys.argv[3]
	code10dic={}; gcodedic = {} ; listtaxa= []
	for line in open(table,'r'):
		code10dic[line.split('\t')[0]] = line.split('\t')[1]
		gcodedic[line.split('\t')[1]] = line.split('\t')[2].split('\n')[0]
		listtaxa.append(line.split('\t')[1])
	
def loop_1():
	if not os.listdir('..') == 'Done':
		os.system('mkdir ../Done/')
	for filename in os.listdir(os.curdir+'./'):
		print(filename)
		if filename[0] != "." and filename.endswith('.fasta') and filename.split('-')[1][0:3] == code:	
			os.system('python3 1a_ContigFiltStats.py -in ../' + filename + ' -out ' + filename.split('_transcripts.fasta')[0] + '_rnaSPAdes --minLen ' + seqlength + ' --spades')
	os.system('python3 1b_XSpeciesContamination.py ../XlaneBleeding')
				
def loop_2():
	print('rRNA and BvE')
	os.system('python3 2_Auto_rRNA_BvE.py')
	
def loop_3():
	for foldername in os.listdir(os.curdir+'./'):
		if '-' in foldername and foldername.split('-')[1][0:3] == code:	
			os.system('python3 3_CountOGsUsearch.py -in ../' + foldername + '/' + foldername + '_WTA_NBU.fasta')	

def loop_4():
	for foldername in os.listdir(os.curdir+'./'):
		if '-' in foldername and foldername.split('-')[1][0:3] == code:	
			os.system('python3 4_InFrameStopFreq.py -in ../' + foldername + '/' + foldername + '_WTA_NBU.Renamed.fasta')	

def loop_5():
	for foldername in os.listdir(os.curdir+'./'):
		if '-' in foldername and foldername.split('-')[1][0:3] == code and foldername.split('-')[0] in listtaxa:	
			os.system('python3 5_GCodeTranslate.py -in ../' + foldername + '/' + foldername + '_WTA_NBU.Renamed.fasta -g '+gcodedic[foldername.split('-')[0]])	

def loop_6():
	for foldername in os.listdir(os.curdir+'./'):
		if '-' in foldername and foldername.split('-')[1][0:3] == code and foldername.split('-')[0] in listtaxa:	
			os.system('python3 6_FilterPartials.py -fp ../'+foldername+'/')
			os.system('python3 6b_update_cov_post_removepartials.py '+foldername.split('-')[0])
			

def loop_7():
	for foldername in os.listdir(os.curdir+'./ToRename/'):
		if foldername.endswith('.Final.AA.ORF.fasta'):
			newname = ('_').join(foldername.split('_')[0:3])
			os.system('python3  7_FinalizeName.py -in ../ToRename/' + foldername + ' -n '+ newname)	



def main():
# 	loop_1()
# 	loop_2()
# 	loop_3()
# 	loop_4()
	loop_5()
	loop_6()
	loop_7()

main()