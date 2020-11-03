**01: oneButtonToGseaFile.py 一键处理表达谱生成可以用来做GSEA分析的脚本**
- 要注意的是 input file 和 annotation file 的第一行都需要用 <span style="color:blue">#</span> 注释, 详见实例文件.
- 由于编者的习惯, 将表达谱第一列的格式规定为 genesymbol(ensgId)(gene_type), 这个在gtf文件中可以提取到的.

```
usage:

python oneButtonToGseaFile.py --input [] --annotation [] --outputPrefix []

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -I INPUT
                        整理好的输入文件,第一行以#注释
  --annotation ANNOTATION, -A ANNOTATION
                        分组注释文件
  --outputPrefix OUTPUTPREFIX, -O OUTPUTPREFIX
                        输出文件的前缀,可以附带目录
```

    
 
