#python setAnalyser.py --c --wc --hist -gpa '/home/inesmendes/Dropbox/Tese/roary/roary_ines_n61/gene_presence_absence.csv' -uso '/home/inesmendes/Dropbox/Tese/roary/roary_ines_n61/set_difference_unique_set_one' -ust '/home/inesmendes/Dropbox/Tese/roary/roary_ines_n61/set_difference_unique_set_two' -cs '/home/inesmendes/Dropbox/Tese/roary/roary_ines_n61/set_difference_common_set' -cg /home/inesmendes/Dropbox/Tese/roary/roary_ines_n61/gene_presence_absence.Rtab'

import sys
import argparse
import csv
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import wordcloud as WordCloud #https://github.com/amueller/word_cloud
import operator
from collections import Counter
import clustergram as ClusterGram
import pandas as pd


def parse_gene_presence_absence(gene_filename):
	geneID={}
	with open(gene_filename, 'r') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			geneID[row[0]]=row[2]
	return geneID

def parse_set(gene_filename):
	unique_set={}
	setfile=open(gene_filename,'r')
	for line in setfile:
		line=line.replace(': ','\t')
		line=line.replace('\n','')
		line=line.split('\t')
		unique_set[line[0]]=line[1:]
	return unique_set

def control(set_dic, setname):
	species=[]
	for key,value in set_dic.iteritems():
		for item in value:
			item=item.split('_')
			if len(item)>3:
				item=item[1:-1]
				temp=''
				for word in item:
					temp+=str(word)+' '
				item=temp
			else:
				item=item[1]
			if item not in species:
				species.append(item)
	print "Analyzing set: " +str(setname)
	print "Number of isolates in set: " + str(len(species))
	print "Species present:"
	for item in species:
		print str(item)
	print '\n'

def makeHistogram(uniqueDic_one, uniqueDic_two):
	values_one=[]
	for key, value in uniqueDic_one.iteritems():
		values_one.append(len(value))
	values_two=[]
	for key, value in uniqueDic_two.iteritems():
		values_two.append(len(value))

	plt.hist(values_one, bins=range(min(values_one), max(values_one),1),color='blue', alpha=0.5, label='Set One')
	plt.hist(values_two, bins=range(min(values_one), max(values_one),1),color='crimson', alpha=0.5, label='Set Two')
	plt.xticks(np.arange(min(values_one), max(values_one)+1,1))
	plt.title("Gene Distribution Across Samples")
	plt.xlabel("Number of Samples")
	plt.ylabel("Frequency")
	plt.grid(True)
	plt.legend()
	plt.show()

def print_report(genes, set, filename):
	toPrint={}
	c=Counter()
	for key,value in set.iteritems():
		gene_annotation=genes[key]
		gene_annotation=gene_annotation.split('[')
		gene_annotation=gene_annotation[0].strip()
		if key not in toPrint:
			toPrint[key, gene_annotation]=len(value)
		'''else:
			c.update(key)
			key=key+'_'+str(c[key])
			toPrint[key]=[gene_annotation,len(value)]'''
	toPrint_sorted=sorted(toPrint.items(), key=operator.itemgetter(1), reverse=True) #tuple

	report=file(filename+'.txt','w')
	string="Gene Key\tGene Annotation\tNumber of Isolates\n"
	report.write(string)
	for key,value in toPrint_sorted:
		string=str(key[0])+'\t'+str(key[1])+'\t'+str(value)+'\n'
		report.write(string)
	report.close()

	value=toPrint_sorted[0][1]
	report=file('present_in_all_'+filename+'.txt','w')
	for key,values in toPrint_sorted:
		if values == value:
			string=str(key[0])+'\t'+str(key[1])+'\n'
			report.write(string)
	report.close()

def make_wordCloud(genes, uniqueID, filename):
	import matplotlib.pyplot as plt
	from wordcloud import WordCloud

	text=""
	for key,value in uniqueID.iteritems():
		gene_annotation=genes[key]
		if "[" in str(gene_annotation):
			gene_annotation=gene_annotation.split('[')
			gene_annotation=gene_annotation[0]
		gene_annotation=gene_annotation.split(' ')
		for i in range(len(gene_annotation)):
			if gene_annotation[i] == 'hypothetical' or gene_annotation[i] == 'protein':
				break
			else:
				text+= str(gene_annotation[i]) + ' '

	cloudfile=file(filename,"w")
	cloudfile.write(text)
	cloudfile.close()

	#Create wordcloud
	wordcloud = WordCloud().generate(text)
	img=plt.imshow(wordcloud)
	plt.axis("off")
	plt.title(filename)
	plt.show() 
	#or sve as png
	#img.write_png(filename+".png") #TODO grava ao contrario!

def getFileWordCloudExpression(genes, uniqueID,filename):
	expressions={}
	for key,value in uniqueID.iteritems():
		gene_annotation=str(genes[key])
		if "[" in gene_annotation:
			if 'hypothetical' not in gene_annotation:
				gene_annotation=gene_annotation.split('[')
				gene_annotation=gene_annotation[0]
				gene_annotation=gene_annotation.strip()
				if gene_annotation not in expressions:
					expressions[gene_annotation]=len(value)
				else:
					old_value=expressions[gene_annotation]
					expressions[gene_annotation]=int(old_value)+int(len(value))
	fh=open(filename,"w")
	for key,value in expressions.iteritems():
		fh.write(str(key)+':'+str(value)+'\n')
	fh.close()
'''
def getClustergram(gene_presence_absence_file):
	from clustergram import clustergram as ClusterGram
	from heatmapcluster import heatmapcluster
	print 'lala'
	df=pd.read_csv(gene_presence_absence_file, header=0, index_col=0, sep='\t' )
	#print len(df)
	#print list(df)
	header_col=list(df)
	#print len(header_col)
	header_line = df.index.values
	#print header_line
	header_line=header_line.tolist()
	m_data = df.as_matrix()
	#lala=open('matrix.csv','w')
	#lala.write(m_data)
	#lala.close()
	#print m_data[-1,:]
	#h=heatmapcluster(m_data, header_line, header_col)
	#print len(m_data)
	#cg=ClusterGram(data=m_data, row_labels=header_line, col_labels=header_col, figname='clustergram.png')
	#plt.show()

	h = heatmapcluster(m_data, header_line, header_col,
                   num_row_clusters=3, num_col_clusters=0,
                   label_fontsize=6,
                   xlabel_rotation=-75,
                   cmap=plt.cm.coolwarm,
                   show_colorbar=True,
                   top_dendrogram=True)
	plt.show()
'''
def runProgram(args):

	#print args

	genes=parse_gene_presence_absence(args.gpa)
	unique_set_one=parse_set(args.uso)
	unique_set_two=parse_set(args.ust)
	common_set=parse_set(args.cs)

	if args.c:
		control(unique_set_one, 'set one, horses')
		control(unique_set_two, 'set one, humans')
		control(common_set, 'common set')

	set_one=print_report(genes,unique_set_one,"test_horse")
	set_two=print_report(genes,unique_set_two,'test_human')
	set_common=print_report(genes,common_set,'test_common')

	if args.wc:
		make_wordCloud(genes, unique_set_one, "wordcloud_set_one")
		make_wordCloud(genes, unique_set_two, "wordcloud_set_two")

	if args.hist:
		makeHistogram(unique_set_one, unique_set_two)

	if args.awc:
		getFileWordCloudExpression(genes, unique_set_one,"wordCloudExpression.setOne.txt")
		getFileWordCloudExpression(genes, unique_set_two,"wordCloudExpression.setTwo.txt")

def main():

	parser = argparse.ArgumentParser(description="Analyzer for Roary's difference query results.")

	parser.add_argument('-gpa', nargs='?', type=str, help='Gene presence absence csv file, generated by roary.', required=True)
	parser.add_argument('-uso', nargs='?', type=str, help='Set difference unique set one file, generated by roary difference between sets of isolates query.', required=True)
	parser.add_argument('-ust', nargs='?', type=str, help='Set difference unique set one file, generated by roary difference between sets of isolates query.', required=True)
	parser.add_argument('-cs', nargs='?', type=str, help='Set difference common set file, generated by roary difference between sets of isolates query.', required=True)
	parser.add_argument('-cg', nargs='?', type=str, help='Generate Clustergram. gene_presence_absence.Rtab file necessary.', required=False)
	parser.add_argument('--c',help='Print sample control on console.', required=False, action='store_true', default=False)
	parser.add_argument('--wc', help='Create a wordcloud.', required=False, action="store_true", default=False)
	parser.add_argument('--hist', help='Create histogram.', required=False, action="store_true", default=False)
	parser.add_argument('--awc', help='Create expression wordcloud text for wordcloud generation in http://www.wordle.net/', required=False, action="store_true", default=False)

	args = parser.parse_args()

	runProgram(args)

if __name__ == "__main__":
    main()