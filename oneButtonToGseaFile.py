#-*- coding:utf-8 -*-
import re,sys,os,argparse

def prepareFile(file,annotationFile,outputPrefix):
    filename = file.split(".txt")[0]
    if re.search("/",outputPrefix):                         ##如果用户输入带路径的输出前缀, 则需要检查文件夹是否存在, 如果不存在则创建.
        path = "/".join(outputPrefix.split("/")[:-1])
        if not os.path.exists(path):
            os.system("mkdir -p {}".format(path))

    ## 读取注释的文件
    annotDict = {}
    for line in open(annotationFile,"r"):
        if re.search("#",line):
            continue
        info = line.rstrip("\n").split("\t")
        annotDict[info[0]] = info[1]

    ## output
    gmt = open("{}.gct".format(outputPrefix),"a")
    cls = open("{}.cls".format(outputPrefix),"a")
    chip = open("{}.chip".format(outputPrefix),"a")
    #genenum = 0
    geneDict ={}
    ensgidToGenesymbol = {}
    classNum = 0
    classList =[]                                           ## 写在cls的表型
    classList2 = []                                         ## 写在gct文件的样品名
    totalgenesymbol = []
    for line in open(file,"r"):
        if re.search("^#",line):     ## 第一行用#开头注释起来
            for index ,element in enumerate(line.rstrip("\n").split("\t")[1:]):
                classList.append(annotDict[element])
                classList2.append(element)
                #classList.append(''.join(re.findall(r'[A-Za-z]', element)))
            continue
        info = line.rstrip("\n").split('\t')
        genesymbol = info[0].split("(")[0]
        if genesymbol not in totalgenesymbol:                ##只能有一个基因symbol值
            totalgenesymbol.append(genesymbol)
            ensgid = info[0].split("(")[1].split(")")[0]
            if re.search(".",ensgid):                        ##防止ensgId附带小数点
                ensgid = ensgid.split(".")[0]
            dataall = []
            for key in info[1:]:
                dataall.append(float(key))
            if set(dataall)=={0}:
                continue
            #genenum += 1
            geneDict[ensgid]= ensgid+"\t"+'\t'.join(info[1:])
            ensgidToGenesymbol[ensgid] = genesymbol

    ## write result
    chip.write("Probe Set ID"+'\t'+"Gene Symbol"+"\t"+"Gene Title"+'\n')
    for key in ensgidToGenesymbol:
        chip.write(key+"\t"+ensgidToGenesymbol[key]+'\t'+"NA"+'\n')
    chip.close()
    lableUnique = list(set(classList))                      ##保持cls文件的顺序
    lableUnique.sort(key=classList.index)
    cls.write("{} {} 1".format(len(classList),len(lableUnique))+'\n')
    cls.write("# {}".format(" ".join(lableUnique)) + '\n')
    cls.write(" ".join(classList) + '\n')
    cls.close()
    gmt.write("#1.2"+'\n')
    gmt.write(str(len(geneDict))+'\t'+str(len(classList))+ '\n')
    gmt.write("NAME"+"\t"+"Description"+'\t'+"\t".join(classList2)+ '\n')
    for key in geneDict:
        gmt.write(key+'\t'+geneDict[key]+'\n')

def main():
    parser = argparse.ArgumentParser(usage='\n\npython oneButtonToGseaFile.py --input [] --annotation [] --outputPrefix []', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--input', '-I', required=True, help='整理好的输入文件,第一行以#注释')
    parser.add_argument('--annotation', '-A', required=True, help='分组注释文件')
    parser.add_argument('--outputPrefix', '-O', required=True, help='输出文件的前缀,可以附带目录')
    args = parser.parse_args()
    inputfile = args.input
    outputPrefix = args.outputPrefix
    annotionFile = args.annotation
    prepareFile(file=inputfile,annotationFile=annotionFile,outputPrefix=outputPrefix)

if __name__ == '__main__':
    main()
