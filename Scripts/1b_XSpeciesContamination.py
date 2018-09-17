#!/usr/bin/python3

__author__ = 'Jean-David Grattepanche'
__version__ = '4, August 21, 2018'
__email__ = 'jeandavid.grattepanche@gmail.com'



import sys
import os
import re
import time
import string
import os.path
from Bio import SeqIO
from sys import argv
listtaxa=[]
toosim = 0.99
seqcoverage = 0.7

def merge_files(folder):
	mergefile = open('../forclustering.fasta','w+')
	print("MERGE following files")
	for taxafile in os.listdir(folder):
		if taxafile[0] != ".":
			taxaname = ('-').join(taxafile.split('-')[0:2])
			print(taxaname)
			listtaxa.append(taxaname)
			mergefile = open('../forclustering.fasta','a')
			for line2 in SeqIO.parse(folder+'/'+taxafile, 'fasta'):
				if int(len(str(line2.seq))) >= 200:
					mergefile.write('>'+taxaname + '_' + line2.description.split('assembled_')[1] + '\n' + str(line2.seq) + '\n')
				else:
					print(line2, " is too short")
			mergefile.close()

	sort_cluster(folder, listtaxa)


def sort_cluster(folder, listtaxa):
	if not os.path.exists('../clusteringresults_vsearch/'):
		os.makedirs('../clusteringresults_vsearch/')
	fastalist = []; fastadict= {}
# 	print("SORT the sequences and create a dictionnary of sequences")
# 	out = open('temp/forclustering_sorted.fasta','w+')
	print("CREATE a dictionnary of sequences")
	for record in SeqIO.parse(open('../forclustering.fasta','r'),'fasta'):
#		IDL  = record.description, int(len(record.seq))
		IDL  = record.description, int(record.description.split('_Cov')[1].replace('\n',''))
		fastalist.append(IDL)
		fastadict[record.description] = record.seq
# 	for seqlength in sorted(fastalist, reverse=True, key=lambda x: x[1]):
# 		out = open('temp/forclustering_sorted.fasta','a')
# 		out.write('>'+seqlength[0] + '\n'+ str(fastadict[seqlength[0]])+'\n')
# 		out.close()	
	print("CLUSTER for too similar sequences that overalp at least 70%")
	print('vsearch --cluster_fast ../forclustering.fasta --strand both --usersort --query_cov '+str(seqcoverage)+' --id '+str(toosim) +' --uc ../clusteringresults_vsearch/results_forclustering.uc' )
	os.system('vsearch --cluster_fast ../forclustering.fasta --strand both --usersort --query_cov '+str(seqcoverage)+' --id '+str(toosim) +' --uc ../clusteringresults_vsearch/results_forclustering.uc' )
	input2 = open('../clusteringresults_vsearch/results_forclustering.uc','r')
	out2 = open('../fastatokeep.fas','w+')
	out3 = open('../fastatoremoved.fas','w+')
	out4 = open('../fastatoremoved.uc','w+')
	print("CREATE a dictionary with clustering results")
	clustdict= {}; clustlist = []; allseq = []; clustline = {}; list= []; i=0; j=0
	for row2 in input2:
		if row2.split('\t')[0] == 'C' and int(row2.split('\t')[2]) < 2: # keep all unique sequences
			out2.write('>'+row2.split('\t')[8] + '\n' + str(fastadict[row2.split('\t')[8]])+ '\n')
		if row2.split('\t')[0] == 'C' and int(row2.split('\t')[2]) > 1: # create another dictionary
#			print("create dico: ", row2.split('\t')[8])
			clustdict.setdefault(row2.split('\t')[8], [row2.split('\t')[8]])
			clustlist.append(row2.split('\t')[8])

	for row3 in open('../clusteringresults_vsearch/results_forclustering.uc','r'):
		if row3.split('\t')[0] == 'H':
#			print("add dico: ", row3.split('\t')[9], row3.split('\t')[8])
			clustdict[row3.split('\t')[9].replace('\n','')].append(row3.split('\t')[8].replace('\n',''))
			clustline[row3.split('\t')[8].replace('\n','')] = row3.replace('\n','')
			clustline[row3.split('\t')[9].replace('\n','')] = row3.replace('\n','')


	print("PARSE the clusters: keep seed sequences (highest coverage) for each cluster")
	for clust in clustlist:
		list = sorted(clustdict[clust], reverse = True, key=lambda x: int(x.split('_Cov')[1]))
		master = list[0]
		Covmaster = int(list[0].split('_Cov')[1])
		master8dig = ('_').join(list[0].split('_')[0:3])[:-2]
		for seq in list:
			clustered =  seq.replace('\n','')
			Covclustered = int(clustered.split('_Cov')[1])
			clustered8dig = ('_').join(clustered.split('_')[0:3])[:-2]
# 			print(master8dig, Covmaster, '//', clustered8dig, Covclustered)
			if float(Covmaster/Covclustered) < 10:
				out2.write('>'+clustered + '\n' + str(fastadict[clustered])+ '\n')
				i +=1
			elif master8dig == clustered8dig:
				out2.write('>'+clustered + '\n' + str(fastadict[clustered])+ '\n')
				i +=1
			elif Covclustered >= 50:				
				out2.write('>'+clustered + '\n' + str(fastadict[clustered])+ '\n')
				i +=1
			else:
				j +=1
				out4 = open('../fastatoremoved.uc','a')
				out3.write('>'+clustered + '\n' + str(fastadict[clustered])+ '\n')
				print(clustline[clustered],'\t' , master )
				out4.write(clustline[clustered]+ '\t' + master + '\n')
				out4.close()
	print('there are ', str(i),' sequences kept for ',str(j),' sequences removed')
	out2.close()
	out3.close()
	splittaxa(listtaxa)
	
def splittaxa(listtaxa):
	for taxa in listtaxa:
		out = open('../'+taxa+'_200bp_Xcont.fasta','w+')
		for kept in SeqIO.parse('../fastatokeep.fas','fasta'):
			if taxa in kept.description:
# 				out = open(taxa+'_200bp_Xcont.fasta','a')
				out.write('>' + kept.description + '\n' + str(kept.seq) + '\n')
# 				os.system('mv ../'+taxa+'_200bp_Xcont.fasta ../'+taxa+'/')
		out.close()
	os.system('mv ../fastatokeep.fas ../clusteringresults_vsearch/')
	os.system('mv ../fastatoremoved.fas ../clusteringresults_vsearch/')
	os.system('mv ../fastatoremoved.uc ../clusteringresults_vsearch/')
	os.system('mv ../forclustering.fasta ../clusteringresults_vsearch/')
	
def main():
	print('\n\n*******************************************************************')
	print('This script:  \n\t Create cluster and remove sequence too similar \n')
	print('\n\npython3 XSpeciesContamination.py folder\n\n')
	print('\nFor example\npython3 XSpeciesContamination.py trimmed200p5 \n\n')
	print('*******************************************************************\n\n')
	script, folder = argv
	merge_files(folder)
# 	sort_cluster(folder)

main()