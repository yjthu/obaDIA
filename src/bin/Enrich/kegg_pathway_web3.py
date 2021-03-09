#Web-visualization for pathway enrichment result
#guoyang 2012-11-06
#yj 20209-03-08
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

def file_fc(file):
	fc={}
	fdr={}
	deg_file=open(file).readlines()
	for eachLine in deg_file:
		if (eachLine.strip() != '') and (not eachLine.startswith('#')):
			temp=eachLine.split('\t')
			fc[temp[0].strip()]=temp[3].strip()
			fdr[temp[0].strip()]=temp[5].strip()
	return (fc,fdr)

def anno(file):
	annomap={}
	genemap={}
	deg_file=open(file).readlines()
	for eachLine in deg_file:
		if (eachLine.strip() != '') and (not eachLine.startswith('#')):
			temp=eachLine.split('\t')
			genemap[temp[0].strip()] = temp[1].strip()
			if(temp[1].strip() in annomap.keys()):
				annomap[temp[1].strip()] = annomap[temp[1].strip()]+','+temp[0].strip()
			else:
				annomap[temp[1].strip()]=temp[0].strip()
	return (genemap,annomap)
	
def main_flow(row):
	each=row.strip().split('\t')
	pathway=each[2].strip().lower()
	link=each[9]
	API=link.strip().replace("http","https")
	
	ko1=[one_ko.strip() for one_ko in each[8].strip().split('|')]  #ko must be lower case in KEGG API,but not in title attribute
	gene1=[one_gene.strip() for one_gene in each[7].strip().split('|')]
	
	bg=link.strip().replace('http://www.genome.jp/kegg-bin/show_pathway?','').replace("%09cyan",'').replace("%09red",'').replace("%09blue",'').replace("%09yellow",'').replace("%09green",'').replace("%09pink",'')
	kob=bg.strip().split('/')[1:-1]
	ko2 = list(set(kob) ^ set(ko1))
	
	gene2 = [annomap[e].strip().split(',') for e in ko2 if e in annomap.keys()]
	gene3 = [i for item in gene2 for i in item]
	ko3 = [genemap[e].strip() for e in gene3]
	
	ko = ko1 + ko3
	gene = gene1 + gene3
	
	
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

	#download the webpage and motify the links
	for i,each_ko in enumerate(ko):
		if each_ko != '':
			if each_ko not in ko_gene:
				ko_gene[each_ko]={}
				ko_gene[each_ko]['up']=[]
				ko_gene[each_ko]['down']=[]
				ko_gene[each_ko]['mid']=[]	
			if gene[i] in up:
				ko_gene[each_ko]['up'].append(gene[i])
			elif gene[i] in down:
				ko_gene[each_ko]['down'].append(gene[i])
			else:
				ko_gene[each_ko]['mid'].append(gene[i])
	
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
							inner_html+='<ul><li><font color=\\"red\\">Up regulated</font></li>'
							for each_gene in ko_gene[each_ko]['up']:
								inner_html+='<ul><font color=\\"red\\"><a href=\\"https://www.uniprot.org/uniprot/%s\\">%s</a>: FC=%s, FDR=%s</font></ul>' % (each_gene, each_gene, fc[each_gene], fdr[each_gene])
							inner_html+='</ul>'
						if ko_gene[each_ko]['down'] != []:
							inner_html+='<ul><li><font color=\\"blue\\">Down regulated</font></li>'
							for each_gene in ko_gene[each_ko]['down']:
								inner_html+='<ul><font color=\\"blue\\"><a href=\\"https://www.uniprot.org/uniprot/%s\\">%s</a>: FC=%s, FDR=%s</font></ul>' % (each_gene, each_gene, fc[each_gene], fdr[each_gene])
							inner_html+='</ul>'
						if ko_gene[each_ko]['mid'] != []:
							inner_html+='<ul><li><font color=\\"black\\">Not significant</font></li>'
							for each_gene in ko_gene[each_ko]['mid']:
								inner_html+='<ul><font color=\\"black\\"><a href=\\"https://www.uniprot.org/uniprot/%s\\">%s</a>: FC=%s, FDR=%s</font></ul>' % (each_gene, each_gene, fc[each_gene], fdr[each_gene])
							inner_html+='</ul>'
							
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
	parser.add_argument('--diff2',required=True,help="diff full file")
	parser.add_argument('--anno',required=True,help="ko anno map file")
	argv=vars(parser.parse_args())
	filename=argv['table']
	DEG=argv['diff']
	Diff=argv['diff2']
	annof=argv['anno']
	name=os.path.basename(filename).replace("xls","html")
	assert not os.system('cp -r %s .' % (sys.path[0]+'/src'))
	
	if not os.path.exists('src'):
		assert not os.system('mkdir src')	
	parallel_result=[]
	row_pathway=file_wash(filename)
	up,down=updown(DEG)
	fc,fdr=file_fc(Diff)
	genemap,annomap=anno(annof)
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
	
	
