# This script reproduces bioinformatics analyses 
# for the purpose of comparing patient-derived pancreatic 
# cancers that proliferate in the liver early or late 
# after splenic injection

# Install necessary packages (only once)
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("hgu133plus2.db")
# biocLite("genefilter")
# biocLite("annotate")
# biocLite("GEOquery")
# biocLite("dplyr")
# biocLite("reshape2")
# biocLite("GSVA")

library(limma)
library(hgu133plus2.db)
library(genefilter)
library(annotate)
library(GEOquery)
library(dplyr)
library(reshape2)
library(GSVA)

# Specify a target folder for downloaded gene expression data
geo.directory = ""

# Download gene expression data from:
# Gene Expression Omnibus Accession GSE46385
# This dataset was originally referenced in a PLoS ONE article 
# with the PubMed ID number: PMID24204737
eset.dl = getGEO('GSE46385',destdir = geo.directory)
eset = eset.dl[[1]]
annotation(eset) = "hgu133plus2"

# Use boxplots as a quality control to verify that the data is normalized
boxplot(exprs(eset))

# Simplify sample information related 
# to patient-derived tumor lines
pData(eset)$Line = gsub("patient/line: ","",pData(eset)$`characteristics_ch1.3`,fixed=T)
pData(eset)$SampleType = gsub("( ).*","",pData(eset)$title)
pData(eset)$Line[grep("mus",pData(eset)$characteristics_ch1.3)] = "mPanc96"

# Select the subset of patient-derived 
# tumor lines analyzed in this study:
# Fast (early proliferators in the liver): T608, T215, mPanc96
# Slow (late proliferators in the liver): T366, T738, T654
eset.splenic = eset[,pData(eset)$SampleType == "Tumor" & pData(eset)$Line %in% c(
  "UVA366", "UVA738", "UVA654","UVA608", "UVA215", "mPanc96")]
pData(eset.splenic)$LiverRate = ifelse(pData(eset.splenic)$Line %in% c(
  "UVA608","mPanc96","UVA215"),"fast","slow")

# Apply feature filtering to gene expression data
# so that there is only one probeset per gene.
# Default parameters choose probesets with the 
# highest inter-quartile range across all samples.
eset.filtered = featureFilter(eset.splenic)

# Annotate probesets with feature information about 
# unique ENTREZ genes and GENE SYMBOLS
eset.annotation = annotation(eset.filtered)
eset.feature.names = featureNames(eset.filtered)
eset.fdata = data.frame(
  ID = eset.feature.names,
  ENTREZ = unlist(lookUp(eset.feature.names, eset.annotation, c("ENTREZID"),load=T)),
  SYMBOL = unlist(lookUp(eset.feature.names, eset.annotation, c("SYMBOL"),load=T)),
  GENE = unlist(lookUp(eset.feature.names, eset.annotation, c("GENENAME"),load=T)),
  stringsAsFactors=F)
eset.fdata$URL = paste0("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",eset.fdata$SYMBOL)
rownames(eset.fdata) = as.character(eset.fdata$ENTREZ)
featureNames(eset.filtered) = as.character(eset.fdata$ENTREZ)
# Replace gene annotation in the filtered dataset
fData(eset.filtered) = eset.fdata

# Specify patient-derived tumors as independent variables
eset.formula = (~0+Line)

# Construct a design matrix associating each sample 
# with one of the 6 patient-derived tumor lines.
eset.design = model.matrix(eset.formula, data = pData(eset.filtered))
colnames(eset.design) = gsub("Line","",colnames(eset.design))

# Formulate contrast between early and late 
# proliferating patient-derived tumor lines
eset.contrast = c(LMets = "(UVA608+ mPanc96 + UVA215) / 3 - (UVA366 + UVA738 + UVA654) / 3")

eset.make.contrast = makeContrasts(contrasts = eset.contrast, levels = eset.design)
colnames(eset.make.contrast) = names(eset.design)

# Apply robust regression using linear models for microarrays (limma)
# to create a table of differentially expressed genes.
efit = lmFit(eset.filtered, eset.design, method="robust", maxit=1000) %>% 
  contrasts.fit(eset.make.contrast) %>% eBayes()

etable = topTable(efit,coef=1,number=Inf,sort.by= "none") %>% 
  mutate(FC = 2^logFC, pctFC = 100*round(FC-1, 3),
         FDR = adj.P.Val, DIR = sign(logFC),
         SIG = -log10(FDR), MAG = SIG*DIR)

# (optional) Create a table for all differentially expressed probesets
efit.all = lmFit(eset.splenic, eset.design, method="robust", maxit=1000) %>% 
  contrasts.fit(eset.make.contrast) %>% eBayes()
etable.all = topTable(efit.all,coef=1, number=Inf, sort.by= "none") %>% 
  mutate(FC = 2^logFC, pctFC = 100*round(FC-1, 3), 
         FDR = adj.P.Val, DIR = sign(logFC),
         SIG = -log10(FDR), MAG = SIG*DIR)

# Load gene set enrichment database for KEGG pathways
# Downloaded from the Broad Institute's Molecular Signatures Database
msig.filename = "KEGG.gmt"
msig.data.read = readLines(msig.filename)
msig.data.list = strsplit(msig.data.read, split="\t")
msig.data.id = sapply(msig.data.list,function(x) x[[1]])
msig.data.sets = lapply(setNames(msig.data.list,msig.data.id),function(x) x[c(-1,-2)])
msig.data.annotation = data.frame(
  row.names = msig.data.id,
  ID = msig.data.id,
  SET = msig.data.id,
  PATHWAY = tolower(gsub("_"," ",gsub("KEGG_","",msig.data.id))),
  URL = sapply(msig.data.list,function(x) x[[2]]),
  ENTREZLIST = sapply(msig.data.sets,function(x) paste(x,collapse=";")),
  ENTREZN = sapply(msig.data.sets,length), stringsAsFactors = F)
msig.data = list(sets = msig.data.sets,annotation = msig.data.annotation)

# Filter out genes from gene sets that do not have expression data
kegg.sets = lapply(msig.data$sets,intersect,featureNames(eset.filtered))

# Merge gene expression changes with gene sets
kegg.melt = melt(kegg.sets) %>% filter(!is.na(value)) %>%
  transmute(ENTREZ = as.character(value), SET = as.character(L1)) %>% 
  left_join(etable %>% select(ENTREZ, SYMBOL, logFC, FDR, SIG, MAG), by="ENTREZ")

# Summarize expression changes across gene sets
kegg.annotation = kegg.melt %>% 
  arrange(FDR,SET,ENTREZ) %>% group_by(SET) %>% summarize(
    MeanFC = mean(logFC),MedianFC = median(logFC),
    NGenesSig = sum(FDR<=0.05),
    NGenesUp = sum(FDR<=0.05 & logFC>0),
    NGenesDn = sum(FDR<=0.05 & logFC<0),
    NGenes = n_distinct(ENTREZ),
    Genes = paste(ENTREZ, collapse=";"),
    GenesSig = paste(ENTREZ[FDR <= 0.05], collapse=";"),
    GenesUp = paste(ENTREZ[FDR <= 0.05 & logFC > 0], collapse=";"),
    GenesDn = paste(ENTREZ[FDR <= 0.05 & logFC < 0], collapse=";"),
    GeneSymbols = paste(sort(SYMBOL), collapse=";"),
    GeneSymbolsSig = paste(SYMBOL[FDR <= 0.05], collapse=";"),
    GeneSymbolsUp = paste(SYMBOL[FDR <= 0.05 & logFC > 0], collapse=";"),
    GeneSymbolsDn = paste(SYMBOL[FDR <= 0.05 & logFC < 0], collapse=";")) %>% 
  ungroup %>% mutate(
    PctSig = 100 * round(NGenesSig / NGenes, 2),
    PctUp = 100 * round(NGenesUp / NGenes, 2),
    PctDn = 100 * round(NGenesDn / NGenes, 2)) %>%
  left_join(msig.data$annotation, by="SET") %>% select(
    Set = SET,Pathway = PATHWAY,Link = URL, MeanFC, MedianFC,
    NGenes, NGenesSig, NGenesDn, NGenesUp,
    PctSig, PctUp, PctDn,
    Genes, GenesSig, GenesDn, GenesUp,
    GeneSymbols, GeneSymbolsSig, GeneSymbolsDn, GeneSymbolsUp) %>%
  data.frame(stringsAsFactors = F, check.names = F, check.rows = F)
rownames(kegg.annotation) = kegg.annotation$Set

gene.kegg.annotation = kegg.melt %>% 
  left_join(kegg.annotation %>% select(SET = Set,Pathway), by="SET") %>%
  arrange(ENTREZ,Pathway) %>% group_by(ENTREZ) %>% 
  summarize(Pathways = paste(Pathway,collapse=";")) %>% ungroup

# Use similar experimental setup for gene set enrichment analysis
gset.design = eset.design
gset.make.contrast = eset.make.contrast
gset.eset = eset.filtered

# Apply the "single sample gene set enrichment analysis" method
gset.result = gsva(gset.eset, kegg.sets, method="ssgsea")
gset = gset.result$es.obs
fData(gset) = kegg.annotation[featureNames(gset),]
gfit = lmFit(gset, gset.design, method="robust", maxit=1000) %>% 
  contrasts.fit(gset.make.contrast) %>% eBayes()
gtable = gfit %>% topTable(coef = 1, number=Inf, sort.by = "none") %>%
  mutate(logFC, FC = 2 ^ logFC, pctFC = 100 * round(FC - 1, 3), FDR = adj.P.Val, 
         DIR = sign(logFC), SIG = -log10(FDR), MAG = SIG*DIR) #%>%

# Some interesting genes to discuss
discussion.genes = c("GMDS","PMM1","PMM2","NID1","PTK2","CDH17","ERBB3","EGFR")
discussion.table = etable %>% 
  filter(SYMBOL %in% discussion.genes) %>% 
  left_join(gene.kegg.annotation, by="ENTREZ") %>% 
  arrange(FDR, desc(abs(logFC))) %>%
  mutate(Pathways = ifelse(is.na(Pathways), "", Pathways))

expression.table.all = etable.all %>% 
  mutate(ENTREZ = ENTREZ_GENE_ID) %>% 
  left_join(gene.kegg.annotation, by="ENTREZ") %>% 
  arrange(FDR, desc(abs(logFC))) %>%
  mutate(Pathways = ifelse(is.na(Pathways), "", Pathways))

expression.table = etable %>% 
  left_join(gene.kegg.annotation, by="ENTREZ") %>% 
  arrange(FDR, desc(abs(logFC))) %>%
  mutate(Pathways = ifelse(is.na(Pathways), "", Pathways))

enrichment.table = gtable %>% arrange(FDR, desc(abs(logFC)))

# Export tables to text files which can be easily imported into Excel.
# Tables 2 and 3 are subsets of splenicInjectionPaperExprsTable.txt
write.table(enrichment.table,file="splenicInjectionPaperKeggTable.txt",
            sep="\t", row.names = F, quote = F)
write.table(expression.table,file="splenicInjectionPaperExprsTable.txt",
            sep="\t", row.names = F, quote = F)
write.table(expression.table.all,file="splenicInjectionPaperAllExprsTable.txt",
            sep="\t", row.names = F, quote = F)
write.table(discussion.table,file="splenicInjectionPaperDiscussionTable.txt",
            sep="\t", row.names = F, quote = F)

