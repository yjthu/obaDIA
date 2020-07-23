#Web-visualization for pathway enrichment result
#yj 2020-03-02
import sys
import os
import os.path
import urllib
#import socket
import time
import eventlet
eventlet.monkey_patch()
import re
import argparse
from bs4 import BeautifulSoup as bs
#from multiprocessing import Pool, cpu_count
from django.template import Template, Context, loader
from django.conf import settings
settings.configure(DEBUG=True, TEMPLATE_DEBUG=True,TEMPLATE_DIRS=(sys.path[0],))

def file_wash(filename):
	washed=[]
	k=0
	for eachLine in open(filename):
		if (eachLine.strip() != '') and (not eachLine.startswith('#Term')) and (len(eachLine.split('\t')) == 10):
			washed.append(eachLine.strip()+'\t'+str(k))
			k=k+1
	return washed

	
def main_flow(row):
	each=row.strip().split('\t')
	pathway=each[2].strip().lower()
	link=each[9]
	API=link.strip().replace("http","https")

	ko=[one_ko.strip() for one_ko in each[8].strip().split('|')]  #ko must be lower case in KEGG API,but not in title attribute
	gene=[one_gene.strip() for one_gene in each[7].strip().split('|')]
	for col in range(7,9):
		temp=[each_temp.strip() for each_temp in each[col].split('|') if each_temp.strip() != '']
		if list(set(temp)) == [] or list(set(temp)) == ['NA']:
			each[col]='NA'
		else:
			each[col]=', '.join(temp)
	try:
		assert len(ko) == len(gene)
	except:
		print 'assert error'+pathway

	
	term=[API,each[0]]+each[1:-1]

	
	return term


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='KEGG pathway enrichment web visualization')
	parser.add_argument('--table',required=True,help='KOBAS (add KO and color) pathway enrichment input file')
	argv=vars(parser.parse_args())
	filename=argv['table']
	name=os.path.basename(filename).replace("xls","html")
	assert not os.system('cp -r %s .' % (sys.path[0]+'/src'))
	
	if not os.path.exists('src'):
		assert not os.system('mkdir src')	
	parallel_result=[]
	row_pathway=file_wash(filename)
	result=[]
	for eachrow in row_pathway:
		temp = main_flow(eachrow)
		result.append(temp)
	t=loader.get_template('Table_template.html')
	c=Context({'terms':result})
	html=t.render(c)
	open(name,'w').write(html)
	#assert not os.system('rm src/Table_template.html')
	#assert not os.system('rm src/Kegg_map_template.html')	
	
	
