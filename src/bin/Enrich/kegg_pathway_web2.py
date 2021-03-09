#Web-visualization for pathway enrichment result
#yj 2019-03-23
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

def updown(file):
	up={}
	down={}
	deg_file=open(file).readlines()
	
	for eachLine in deg_file:
		if eachLine.strip() != '':
			temp=eachLine.split()
			if temp[1] == "down":
				down[temp[0].strip()]=temp[0].strip()
			else:
				up[temp[0].strip()]=temp[0].strip()
	return (up,down)
	
def main_flow(row):
	each=row.strip().split('\t')
	pathway=each[2].strip().lower()
	API=each[9].strip().replace("http","https")

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
	ko_gene={}
	ko_red=[]
	ko_yellow=[]
	ko_green=[]
	#download the webpage and motify the links
	for i,each_ko in enumerate(ko):
		if each_ko != '':
			if each_ko not in ko_gene:
				ko_gene[each_ko]={}
				ko_gene[each_ko]['up']=[]
				ko_gene[each_ko]['down']=[]
			if gene[i] in up:
				ko_gene[each_ko]['up'].append(gene[i])
			if gene[i] in down:
				ko_gene[each_ko]['down'].append(gene[i])
	
	term=[API,each[0]]+each[1:-1]

	with eventlet.Timeout(300,False):
		try:
			content=urllib.urlopen(API)
			soup=bs(content,"html.parser")
			img=soup.find_all('img')
			for each_img in img:
				if each_img.has_attr('id') and each_img['id'] == 'pathwayimage':
					pic_url='https://www.genome.jp'+ each_img['src']
					break
			urllib.urlretrieve(pic_url,filename='src/%s.png' % pathway)

			map=soup.map
			for each_area in map.find_all('area'):
				if not (each_area.has_attr('id') and each_area.has_attr('href') and each_area.has_attr('class') and each_area['class'] == 'original'):
					continue
				ko_set=[each_ko.strip() for each_ko in re.search(r'\?(.*)',each_area['href']).group(1).split('+')]
				each_area['href']='https://www.genome.jp'+each_area['href']
				inner_html='<ul>'
				for each_ko in ko_gene:
					if each_ko in ko_set:
						inner_html+='<li>%s</li>' % each_ko
						if ko_gene[each_ko]['up'] != []:
							inner_html+='<ul><li><font color=\\"red\\">Up regulated</font></li><ul><font color=\\"red\\">%s</font></ul></ul>' % ' '.join([each_gene for each_gene in ko_gene[each_ko]['up']])
						if ko_gene[each_ko]['down'] != []:
							inner_html+='<ul><li><font color=\\"blue\\">Down regulated</font></li><ul><font color=\\"blue\\">%s</font></ul></ul>' % ' '.join([each_gene for each_gene in ko_gene[each_ko]['down']])
				inner_html+='</ul>'
				if inner_html != '<ul></ul>':
					each_area['onmouseover']='javascript: showInfo("'+inner_html+'");'
					del each_area['onmouseout']
			t=loader.get_template('src/Kegg_map_template.html')
			c=Context({'title':pathway,'map_content':str(map.prettify(formatter=None)),'image':pathway+'.png'})
			html=t.render(c)
			open('src/'+pathway+'.html','w').write(html)
			link='src/'+pathway+'.html'
			term=[link,each[0]]+each[1:-1]
		except:
			print pathway+" fail..."
	return term


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='KEGG pathway enrichment web visualization')
	parser.add_argument('--table',required=True,help='KOBAS (add KO and color) pathway enrichment input file')
	parser.add_argument('--diff',required=True,help="gene UpDown file")
	argv=vars(parser.parse_args())
	filename=argv['table']
	DEG=argv['diff']
	name=os.path.basename(filename).replace("xls","html")
	assert not os.system('cp -r %s .' % (sys.path[0]+'/src'))
	
	if not os.path.exists('src'):
		assert not os.system('mkdir src')	
	parallel_result=[]
	row_pathway=file_wash(filename)
	up,down=updown(DEG)
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
	
	
