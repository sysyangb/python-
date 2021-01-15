#-*- coding:utf-8 -*-
import re,sys,os
from collections import OrderedDict
import rpy2.robjects as robjects
import requests
from bs4 import BeautifulSoup
from retrying import retry

@retry(stop_max_attempt_number=5)
def getSummary(geneId):

    res = requests.get('https://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&dopt=full_report&list_uids={}'.format(geneId),timeout=10)

    res.encoding = 'utf-8'

    soup = BeautifulSoup(res.text,'html.parser')

    summarylist = []

    for key in soup.find('dl',id='summaryDl'):
        if key.string!="\n":
            newstr = re.sub("<.*?>"," ",str(key).replace("\n",""))
            temp = []
            for key in newstr.split(' '):
                if key != "":
                    temp.append(key)
            summarylist.append(" ".join(temp))

    summary = ""

    summaryindex = "NA"

    if "Summary" in summarylist:
        summaryindex = summarylist.index("Summary")+1

    if summaryindex!="NA":
        summary = summarylist[summaryindex]
    else:
        summary = "No Summary"

    genetype = summarylist[summarylist.index("Gene type")+1]
    genename = summarylist[summarylist.index("Official Symbol")+1]
    fullname = summarylist[summarylist.index("Official Full Name")+1]
    if re.search("provided by",genename):
        genename = genename.split("provided by")[0]
    if re.search("provided by",fullname):
        fullname = fullname.split("provided by")[0]

    return '\t'.join([genename,fullname,genetype,summary])

_input = sys.argv[1]
_num = int(sys.argv[2])
ensgList = []
for line in open(_input,"r"):
    if re.search("^ensgId", line):
        continue
    ensgList.append(line.rstrip("\n").split('\t')[_num])

out = open("genes.txt", "a")
headers = ["species", "ensgId", "entezeId", "hgncId", "strand", "length", "position", "summary"]
out.write('\t'.join(headers) + '\n')
for key in ensgList:
    temp = ["Homo_sapiens", key, "-", "-", "-", "-", "-", "-"]
    out.write("\t".join(temp) + '\n')
out.close()

rscript = '''
library(XML)
library(RCurl)
library(httr)
library(xml2)

allgene <- read.table("genes.txt",header = T,sep = "\t",stringsAsFactors = F)
for(i in 1:nrow(allgene)){
  if(allgene[i,3]=="-"){
    tryCatch({
      url<-paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",allgene[i,2],"&keywords=",allgene[i,1])
      print(paste0(i," --- scrapy genecards information of ",allgene[i,2]))
      geneinformation <- GET(url,content_type_xml())
      html <- htmlParse(geneinformation,encoding = 'utf-8')
      #entrze id 
      allgene[i,3] <- as.character(xpathSApply(doc=html,path = "//div[@class='gc-subsection']/div[@class='gc-subsection-inner-wrap']/ul/li/a/text()",fun = xmlValue)[[2]])
      #hgnc id 
      allgene[i,4] <- as.character(xpathSApply(doc=html,path = "//div[@class='gc-subsection']/div[@class='gc-subsection-inner-wrap']/ul/li/a/text()",fun = xmlValue)[[1]])
      #strand
      allgene[i,5] <- as.character(xpathSApply(doc = html,path = "//div[@class='gc-subsection']/div[@class='gc-subsection-inner-wrap']/div[@class='row']/div/dl/dd/text()",fun=xmlValue)[[3]])
      #length
      allgene[i,6] <- as.character(xpathSApply(doc = html,path = "//div[@class='gc-subsection']/div[@class='gc-subsection-inner-wrap']/div[@class='row']/div/dl/dd/text()",fun=xmlValue)[[2]])
      #position
      allgene[i,7] <- as.character(xpathSApply(doc = html,path = "//div[@class='gc-subsection']/div[@class='gc-subsection-inner-wrap']/div[@class='row']/div/dl/dt/text()",fun=xmlValue)[[1]])
      #summary
      allgene[i,8] <- as.character(gsub("\r\n","",xpathSApply(doc = html,path = "//*[@id='summaries']/div[1]",fun = xmlValue)))
    }, error = function(e){print(paste0("Error ",i," ",allgene[i,2]))}
    ) 
    Sys.sleep(runif(1,min=1,max=2))}
}
write.table(allgene[,c(1:7)],sep = "\t",file = "genecards.txt")
'''
###开始爬去genecards信息
robjects.r(rscript)

###开始爬去NCBI信息
ensgIDdict = OrderedDict()
for line in open("genecards.txt","r"):
    if re.search("species",line):
        continue
    info = line.rstrip('\n').split('\t')
    ensgid = info[2].split('"')[1]
    ensgIDdict[ensgid] = OrderedDict()
    ensgIDdict[ensgid]['geneid'] = info[3].split('"')[1]
    ensgIDdict[ensgid]['genecards'] = '\t'.join(info[3:])

out = open(sys.argv[3],"a")
out.write("ensgId"+'\t'+"geneSymbol"+'\t'+"full name"+'\t'+"gene type"+'\t'+"summary"+'\t'+"entrzeId"+'\t'+"HGNCId"+'\t'+"strand"+'\t'+"length"+"\t"+"position"+'\n')
n = 1
for key in ensgIDdict:
    ncbisummary = ""
    print("{} --- scrapy NCBI information of {}".format(n,key))
    n += 1
    try:
        ncbisummary = getSummary(ensgIDdict[key]['geneid'])
        out.write(key+'\t'+ncbisummary+'\t'+ensgIDdict[key]['genecards']+'\n')
    except:
        out.write(key+'\t'+"\t".join([".",".",".","."])+'\t'+ensgIDdict[key]['genecards']+'\n')
out.close()

os.system("rm genes.txt genecards.txt")