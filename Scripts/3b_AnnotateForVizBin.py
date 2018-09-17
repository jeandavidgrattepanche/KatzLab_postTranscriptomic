import os, sys
from Bio import SeqIO
from Bio.SeqUtils import GC


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

#------------------------------- Main Functions of Script -------------------------------#

def prep_data(fasta_file, tsv_file,MjrC_MnrC):

	inFasta = [i for i in SeqIO.parse(fasta_file,'fasta')]
	inTSV = [i for i in open(tsv_file).read().split('\n') if i != '']
	info_dict = {}
	for i in inTSV:
		if i.split('\t')[1].split('_')[-1] != 'group':
			info_dict.setdefault(i.split('\t')[0],[]).append('_'.join(i.split('\t')[1].split('_')[:2]))
			info_dict[i.split('\t')[0]].append(float(i.split('\t')[2]))
	for i in inFasta:
		if i.description in info_dict.keys():
			info_dict[i.description].append("%.4f" % GC(i.seq))
			info_dict[i.description].append(str(len(i.seq)))
			info_dict[i.description].append(str(i.seq))

	with open('../ForVizBin/'+fasta_file.split('/')[-1].split('.fa')[0]+'.ForVizBin.fasta','w+') as w:
		for k, v in info_dict.items():
			w.write('>'+k+'\n'+v[-1]+'\n')
	with open('../ForVizBin/'+fasta_file.split('/')[-1].split('.fa')[0]+'.AnnotationsForVizBin.txt','w+') as w:
		w.write('label,gc,length\n')
		for k, v in info_dict.items():
			w.write(v[0]+','+','.join(v[2:4])+'\n')
	with open('../ForVizBin/'+fasta_file.split('/')[-1].split('.fa')[0]+'.AnnotationsForVizBin.StrongContam.txt','w+') as w:
		w.write('label,gc,length\n')
		for k, v in info_dict.items():
			if v[1] >= 80 and v[0] != MjrC_MnrC:
				w.write(v[0]+','+','.join(v[2:4])+'\n')
			else:
				w.write('NotStrong,'+','.join(v[2:4])+'\n')

def main():
	if len(sys.argv) != 3:
		print color.BOLD + '\n\nThis script is intended to provide annotations for'\
		+color.RED+' VizBin: '+color.CYAN+'http://claczny.github.io/VizBin/'+color.END
		print color.BOLD + '\n\nThese annotations are based upon the OUTPUT spreadsheet of:'\
		+color.CYAN+' 3_CountOGsUsearch.py'+color.END+color.BOLD+'\n\nYou will need to include'\
		+' the '+color.RED+' MAJOR and MINOR CLADE'+color.END+color.BOLD+' of your taxon for this'\
		+' script to be helpful!\n\n'
		print  color.BOLD + color.RED + 'Two Example Usages Below:\n\n\t' + color.CYAN + 'katzlab$ 3b_AnnotateForVizBin.py'\
		+' ../Op_me_Xxma/Op_me_Xxma_WTA_NBU.Renamed.fasta Op_me\n\n\tkatzlab$ 3b_AnnotateForVizBin.py'\
		+' ../LKH093_Dileptus_mucronatus/LKH093_Dileptus_mucronatus_WTA_NBU.Renamed.fasta Sr_ci\n\n'+color.END
		print color.ORANGE + color.BOLD + '\t\tQuestions/Comments? Email Xyrus (author) at maurerax@gmail.com\n\n'+color.END
		sys.exit()
	elif len(sys.argv) == 3:
		fasta_file = sys.argv[1]
		tsv_file = fasta_file.split('.fa')[0]+'_allOGCleanresults.tsv'
		MjrC_MnrC = sys.argv[2]
		os.system('mkdir ../ForVizBin/')
		print fasta_file
		print tsv_file
		prep_data(fasta_file, tsv_file, MjrC_MnrC)
		
main()