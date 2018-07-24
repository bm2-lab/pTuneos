#!/usr/bin/python
# -*- coding: UTF-8 -*-
######annotate neoantigens related gene using GEPIA,THPA,
import sys,getopt,os,subprocess
import pandas as pd 
from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
import time,os
import urllib
from pyvirtualdisplay import Display
from selenium import webdriver
from pyper import *
binary = FirefoxBinary('/usr/bin/firefox')

opts,args=getopt.getopt(sys.argv[1:],"hi:o:s:",["input_neo_file","out_dir","sample_id"])
input_neo_file =""
out_dir=""
sample_id=""
USAGE='''
	This script convert deletion VCF derived VEP result to fasta format file for netMHC
	usage: python deletion2fasta.py -i <input_neo_file> -o <outdir> -s <sample_id>
		required argument:
			-i | --input_neo_file : input final neoantigen file
			-o | --out_dir : output directory
			-s | --sample_id : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_neo_file"):
		input_neo_file=value
	elif opt in ("-o","--out_dir"):
		out_dir =value
	elif opt in ("-s","--sample_id"):
		sample_id =value  
	
#print coverage
if (input_neo_file =="" or out_dir =="" or sample_id==""):
	print USAGE
	sys.exit(2)


data_final=pd.read_table(input_neo_file,header=0,sep='\t')
gene_list=data_final.Gene.drop_duplicates()

def expdiy_url_get(genename):
    display = Display(visible=0, size=(800, 600))
    display.start()
    browser = webdriver.Firefox(firefox_binary=binary)
    try:
        browser.get("http://gepia.cancer-pku.cn/detail.php?gene={name}&tag=expdiy".format(name=genename))
        time.sleep(2)
        expdiy_pic_url=browser.find_element_by_xpath('//*[@id="pre_download"]/a').get_attribute('href')
        browser.quit()
        display.stop()
        return expdiy_pic_url
    except:
        browser.quit()
        display.stop()
        return "null_url"

def barplot_url_get(genename):
    display = Display(visible=0, size=(800, 600))
    display.start()
    browser = webdriver.Firefox(firefox_binary=binary)
    try:
        browser.get("http://gepia.cancer-pku.cn/detail.php?gene={name}&tag=expdiy".format(name=genename))
        time.sleep(2)
        barplot_pic_url=browser.find_element_by_xpath('//*[@id="general_load_barplot"]/iframe').get_attribute('src')
        browser.quit()
        display.stop()
        return barplot_pic_url
    except:
        browser.quit()
        display.stop()
        return "null_url"
        

def survival_url_get(genename):
    display = Display(visible=0, size=(800, 600))
    display.start()
    browser = webdriver.Firefox(firefox_binary=binary)
    try:
        browser.get("http://gepia.cancer-pku.cn/detail.php?gene={name}&tag=survival".format(name=genename))
        time.sleep(2)
        survival_pic_url=browser.find_element_by_xpath('//*[@id="iframe"]').get_attribute('src')
        browser.quit()
        display.stop()
        return survival_pic_url
    except:
        browser.quit()
        display.stop()
        return "null_url"
        

def bodymap_url_get(genename):
    display = Display(visible=0, size=(800, 600))
    display.start()
    browser = webdriver.Firefox(firefox_binary=binary)
    try:
        browser.get("http://gepia.cancer-pku.cn/detail.php?gene={name}&tag=bodymap".format(name=genename))
        time.sleep(2)
        bodymap_tumor_pic_url = browser.find_element_by_xpath('//*[@id="general_load_bodymap_tumor"]/iframe').get_attribute('src')
        bodymap_normal_pic_url = browser.find_element_by_xpath('//*[@id="general_load_bodymap_normal"]/iframe').get_attribute('src')
        browser.quit()
        display.stop()
        return bodymap_tumor_pic_url,bodymap_normal_pic_url
    except:
        browser.quit()
        display.stop()
        return "null_url","null_url"
        


def Download_pic(url,savePath,file_pre):
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    try:
        urlopen = urllib.URLopener()
        fp = urlopen.open(url)
        data = fp.read()
        fp.close()
        file=open(savePath + file_pre,'w+b')
        file.write(data)
        file.close()
    except IOError, error:
        print "DOWNLOAD %s ERROR!==>>%s" % (url, error)
    except Exception, e:
        print "Exception==>>" + e

for gene in gene_list:
    expdiy_url = expdiy_url_get(gene)
    barplot_url = barplot_url_get(gene)
    survival_url = survival_url_get(gene)
    bodymap_tumor_url,bodymap_normal_url = bodymap_url_get(gene)
    Download_pic(bodymap_tumor_url,out_dir+'/'+gene+'/',gene+'_bodymap_tumor')
    Download_pic(bodymap_normal_url,out_dir+'/'+gene+'/',gene+'_bodymap_normal')
    Download_pic(expdiy_url,out_dir+'/'+gene+'/',gene+'.svg')
    Download_pic(barplot_url,out_dir+'/'+gene+'/',gene+'_expression_barplot')
    Download_pic(survival_url,out_dir+'/'+gene+'/',gene+'.pdf')    




def Plot_Protein_Expression(pro_exp_database,gene_list_input,out_dir):
    r=R()
    str_pro='''
        pro_exp_db=\"%s\"
        gene_list=\"%s\"
        outdir=\"%s\"
        library(ggplot2)
        data_tumor<-read.table(,header=TRUE,sep',')
        gene_list=
        for (gene in gene_list){
            new_data_tumor<-subset(data_tumor,data_tumor$Gene.name==gene,select=c('Tumor','Level','Count.patients'))
            p<-ggplot(new_data_tumor,aes(x=Tumor,y=Count.patients,fill=Level))+geom_bar(stat="identity")+xlab("Tumor Type")+ylab("Patient Count")+ggtitle(paste(paste("Protein expression of",gene_select),'in different tumor'))+theme_bw()+theme(axis.line=element_line(colour="black",size=1))+theme(axis.text.x=element_text(face="italic",colour="black",size=10,angle=40,hjust=1))
            file=paste(gene_select,'_Protein_Exp.png')
            ggsave(p,filename=paste(ourdir,file,sep='/'))
}
'''%(pro_exp_database,gene_list_input,out_dir)
path_pro_exp_database="${DATABASE}/Annotation/protein/cancer.csv"#Path to the cancer.csv file
#Plot_Protein_Expression(path_pro_exp_database,gene_list,out_dir)








	

