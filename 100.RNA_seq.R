rm(list=ls())
gc()

library(stringr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(reshape2)
library(limma)

source('/data/panzhong/genome_sequence/000Seurat.functions.panzhong.version5.R')
source('/data/jianzhongxiang/genome_sequence/000.Function.jzx.R')

config<-read.delim('config.tsv')

#读取参数
GSE_number = config[1,2]
species = config[2,2]
projectid = config[3,2]
plotid = as.numeric(config[4,2])
ID_type = config[5,2]
g_test = config[6,2]
g_control = config[7,2]
pvalue_column = config[8,2]
pvalue_cutoff = as.numeric(config[9,2])
fc_cutoff = as.numeric(config[10,2])
go_analysis = ifelse(config[11,2] == 'T',T,F)
analysis_type = config[12,2]


df_sample <- config[-1:-14,-3]
colnames(df_sample) <- c('sample','group')
rownames(df_sample) <- df_sample[,1]

if(grepl("GSE", GSE_number)){
  df_count <- get_array_data(GSE_number = GSE_number)
}else{
  df_count <- tryCatch({
    read.delim(file = paste0(projectid, '.txt'), row.names = 1)
  }, error = function(e) {
    read.delim(file = paste0(projectid, '.tsv'), row.names = 1)
  })
  df_count <- ID_convert(df_count,fromType = ID_type, toType="SYMBOL", species=species)
}

df_count<-na.omit(df_count)
df_count <- rawdata_colname_match(data = df_count,column_order = rownames(df_sample))

if(analysis_type == 'DESeq2'){
  
  diff_genes <- rna_seq_count_deseq2(
    df_count,
    df_sample,
    gcontrol = g_control,
    gtest = g_test,
    plotid = plotid,
    fc_cutoff = fc_cutoff,
    fc_column = 'log2FoldChange',
    pvalue_cutoff = pvalue_cutoff,
    pvalue_column = pvalue_column,
    go_analysis = go_analysis,
    species = species
  )
  
  #df_count <- read.delim(file = paste0('fig',plotid,'c.',projectid,".DESeq2.ncount_vst.tsv"),row.names = 1)
}else if(analysis_type == 'limma'){
  
  pvalue_column <- switch(pvalue_column,
                          pvalue = "P.Value",
                          "adj.P.Val")
  
  diff_genes <- rna_seq_tpm_limma_test(
    df_tpm = df_count,          # TPM数据框
    df_sample = df_sample,      # 样本信息数据框
    plotid = plotid,            # 绘图ID
    gcontrol = g_control,       
    gtest = g_test,             
    fc_cutoff = fc_cutoff,      
    fc_column = 'logFC',        
    pvalue_cutoff = pvalue_cutoff,  
    pvalue_column = pvalue_column,  
    go_analysis = go_analysis,  # 是否进行GOKEGG分析
    species = species,          # 物种名称
    recutoff = F                
  )
}else if(analysis_type == 'Wilcoxon'){
  
  pvalue_column <- switch(pvalue_column,
                          pvalue = "Pvalue",
                          "padj")
  
  diff_genes <- rna_seq_tpm_wilcox_test(
    df_tpm = df_count,
    df_sample,
    plotid = plotid,
    gcontrol = g_control,
    gtest = g_test,
    fc_cutoff = fc_cutoff,
    fc_column = 'log2FoldChange',
    pvalue_cutoff = pvalue_cutoff,
    pvalue_column = pvalue_column,
    go_analysis = go_analysis,
    species = species,
    recutoff = F
  )
}





