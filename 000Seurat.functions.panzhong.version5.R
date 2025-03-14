### source('g:/BaiduSyncdisk/100sop/sop_perl/000Seurat.functions.panzhong.version5.R')
### source('/data/panzhong/genome_sequence/000Seurat.functions.panzhong.version5.R')


rds_gene_set<-'/data/panzhong/genome_sequence/scrna/list_gene_set.YB20250208.rds'

output_table_with_rowname<-function(df_data,name='gene',output='output')
{
	df_data<-as.data.frame(df_data)
  df_data[[name]]<-rownames(df_data)
  nc<-ncol(df_data)
  df_data<-df_data[,c(nc,1:(nc-1))]
  write.table(df_data,paste0(output,'.tsv'),quote=F,sep='\t',col.names = T,row.names = F)
}

file2genes<-function(filename)
{
  df_string<-read.delim(filename,header=F)
  head(df_string)
  colnames(df_string)[1]<-'gene'
  string_genes<-unique(df_string$gene)
  return(string_genes)
}


## mycolors<-get_colors(2)
get_colors<-function(n,style='defaut')
{
  if(n<=9)
  {
    library(ggsci)
    library("scales")
    library(RColorBrewer)
    
    #mycolors= pal_lancet('lanonc')(9)
    mycolors<-brewer.pal(9,"Set1")
    mycolors<-mycolors[1:n]
    #show_col(mycolors)
				if(n==2)
				{
					if(style=='defaut')
					{
						mycolors<-brewer.pal(9,"Set1")
						mycolors<-mycolors[1:n]
					}else if(style=='style1')
					{
						mycolors<-c('blue','red')
					}else if(style=='style2')
					{
						mycolors<-c('orange','red')
					}
				}
  }else if(n<=20)
  {
    library(ggsci)
    library("scales")
    mycolors= pal_d3('category20')(20)
    # show_col(mycolors)
    mycolors<-mycolors[1:n]
  }else if(n<=51)
  {
    library(ggsci)
    library("scales")
    mycolors= pal_igv()(51)
    #show_col(mycolors)
    mycolors<-mycolors[1:n]
  }else
  {
    library(viridis)
    mycolors<-viridis(n, option = "D")
    #show_col(mycolors)
  }
  return(mycolors)
}




#  ids=idmap('GPL6947',type = 'soft')
#  head(ids)
## expr<-matrix_probeid2symbol(exp,ids)
matrix_probeid2symbol<-function(mt_data,ids){
  colnames(ids)<-c('probeID','symbol')
  head(ids)
  ids_used=ids[ids$symbol != '',]
  ids_used=ids_used[ids_used$probeID %in%  rownames(mt_data),]
  
  mt_data=mt_data[ids_used$probeID,] 
  ids_used$median=apply(mt_data,1,median)
  ids_used=ids_used[order(ids_used$symbol,ids_used$median,decreasing = T),]
  ids_used=ids_used[!duplicated(ids_used$symbol),]
  mt_data=mt_data[ids_used$probeID,]
  rownames(mt_data)=ids_used$symbol
  return(mt_data)
}

matrix_fpkm2tpm<-function(df_fpkm)
{
  fpkm2tpm <- function(fpkm)
  {
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  df_tpm <- apply(df_fpkm,2,fpkm2tpm)
  return(df_tpm)
}



rna_seq_qc<-function(df_data,phenotype,group='group',sampl_cutoff=30,plotid='10')
{
  # df_data<-data_ref
  #  phenotype<-phenotype_ref
  if(1)   #########  绘制boxplot和logCount的boxplot  #########  绘制样本相关性矩阵
  {
    nsample<-dim(df_data)[2]
    if(nsample>sampl_cutoff)
    {
      samples<-sample(colnames(df_data),sampl_cutoff)
      df_data_boxplot<-as.matrix(df_data[,samples])
      class(df_data_boxplot)
    }else{
      df_data_boxplot<-as.matrix(df_data)
    }
    df_data_boxplot[1:3,1:3]
    dim(df_data_boxplot)
    df_data_boxplot<-as.data.frame(df_data_boxplot)
    phenotype_boxplot<-phenotype[colnames(df_data_boxplot),,drop=F]
    ggplot2_boxplot_rna_seq(df_data_boxplot,phenotype_boxplot,group='group',samplename = T,output=paste0('fig',plotid,'a','.boxplot'))
    ggcorrplot_1matrix(df_data_boxplot,output=paste0('fig',plotid,'b',".sample_correlation"),save.data=T)
  }
  
  if(1)  #########  绘制PCA图
  {
    class(df_data)
    sds<-apply(df_data,1,sd)
    df_data_pca<-df_data[sds>0,]
    ggplot2_pca(df_data_pca,phenotype,output=paste0('fig',plotid,'c',".PCA"))
    ggplot2_pca_3d(df_data_pca,phenotype,output=paste0('fig',plotid,'d',".PCA-3d"))
  }
}

rna_seq_qc_combat<-function(df_data,phenotype,group='project',sampl_cutoff=30,plotid='10')
{
  # df_data<-data_ref
  #  phenotype<-phenotype_ref
  if(1)   #########  绘制boxplot和logCount的boxplot  #########  绘制样本相关性矩阵
  {
    nsample<-dim(df_data)[2]
    if(nsample>sampl_cutoff)
    {
      samples<-sample(colnames(df_data),sampl_cutoff)
      df_data_boxplot<-as.matrix(df_data[,samples])
      class(df_data_boxplot)
    }else{
      df_data_boxplot<-as.matrix(df_data)
    }
    df_data_boxplot[1:3,1:3]
    dim(df_data_boxplot)
    df_data_boxplot<-as.data.frame(df_data_boxplot)
    phenotype_boxplot<-phenotype[colnames(df_data_boxplot),,drop=F]
    ggplot2_boxplot_rna_seq(df_data_boxplot,phenotype_boxplot,group=group,samplename = T,output=paste0('fig',plotid,'a','.boxplot'))
    ggcorrplot_1matrix(df_data_boxplot,output=paste0('fig',plotid,'b',".sample_correlation"),save.data=T)
  }
  
  if(1)  #########  绘制PCA图
  {
    class(df_data)
    sds<-apply(df_data,1,sd)
    df_data_pca<-df_data[sds>0,]
    ggplot2_pca(df_data_pca,phenotype,output=paste0('fig',plotid,'c',".PCA"))
    ggplot2_pca_3d(df_data_pca,phenotype,output=paste0('fig',plotid,'d',".PCA-3d"))
    ggplot2_pca_2factor(df_data_pca,phenotype,output=paste0('fig',plotid,'e',".PCA-2factor"))
  }
}

combat_with_reference<-function(test_project='',test_data='',test_phenotype='',ref_project='',ref_data='',ref_phenotype='',plotid=10)
{
  if(identical(rownames(test_phenotype),colnames(test_data)))
  {
    print('test data good!')
  }
  if(identical(rownames(ref_phenotype),colnames(ref_data)))
  {
    print('ref data good!')
  }  
  
  test_phenotype$project=test_project
  ref_phenotype$project=ref_project
  combi_phenotype<-rbind(test_phenotype,ref_phenotype)
  output_table_with_rowname(combi_phenotype,name='sample',output=paste0('fig',plotid,'.',ref_project,'-',test_project,'-combat.phenotype'))
  venn_genes<-intersect(rownames(test_data),rownames(ref_data))
  combi_data<-cbind(test_data[venn_genes,],ref_data[venn_genes,])
  rna_seq_qc_combat(combi_data,combi_phenotype,group='project',sampl_cutoff=200,plotid=paste0(plotid,'a.before.'))
  
  library(sva)
  model <- model.matrix(~group, data = combi_phenotype)
  df_expr_combat <- ComBat(dat = combi_data, batch = combi_phenotype$project, mod = model,ref.batch=ref_project)
  rna_seq_qc_combat(df_expr_combat,combi_phenotype,group='project',sampl_cutoff=200,plotid=paste0(plotid,'b.after.'))
  df_expr_combat<-as.data.frame(df_expr_combat)
  output_table_with_rowname(df_expr_combat,name='gene',output=paste0('fig',plotid,'.',ref_project,'-',test_project,'-combat.data'))
  return(df_expr_combat)
}

### input主要有两个矩阵：df_count和df_sample。
### 要求:1,df_count是计数矩阵，通常行名为基因symbol，列名为sampleid
### 要求:2,df_sample是样本分组矩阵，行名为sampleid。其中一列为group。矩阵中可以包含其它列。
### 要求:3,df_count的列名与df_sample的行名一一对应，且顺序一致（此条经优化后，可以不严格遵循）。
###  deseq2的fc列pvalue列fdr列：log2FoldChange、P.Value、padj
rna_seq_count_deseq2<-function(df_count,df_sample,plotid=4,gcontrol='Control',gtest='ALS',
fc_cutoff=2,fc_column='log2FoldChange',pvalue_cutoff=0.05,pvalue_column='padj',
go_analysis=F,species='human',recutoff=F,countfilter=1)
{
library(stringr)
library(DESeq2)
library(BiocParallel)

  df_sample_used<-df_sample[df_sample$group %in% c(gcontrol,gtest),,drop=F]
  df_sample_used$group<-factor(df_sample_used$group,levels=c(gcontrol,gtest))
  df_sample_used<-df_sample_used[order(df_sample_used$group),,drop=F]
  
  samples_df_count<-colnames(df_count)
  samples_df_sample<-rownames(df_sample_used)
  samples_venn<-samples_df_sample[samples_df_sample %in% samples_df_count]
  df_sample_used<-df_sample_used[samples_venn,,drop=F]
  group_list<-df_sample_used$group
  df_group<-table(group_list)
  df_group<-as.data.frame(df_group)
  colnames(df_group)<-c('group','n_sample')
  write.table(df_group,file = paste0('fig',plotid,'a.',projectid,".sample_number.tsv"),row.names = F,sep='\t',quote = FALSE)
  
  if(df_group[1,2]>0 & df_group[2,2]>0)
  {
    print(df_group)
  }else{
    print('sample number good!')
    return()
  }
    
  df_count_used<-df_count[,samples_venn]
  logFC_cutoff<-log2(fc_cutoff)
  comparison<-paste0(gtest,'_vs_',gcontrol)
	if(!recutoff)
	{
		  if(1)   #####  使用DESeq2进行差异分析，并且导出标准化count
		  {
			# 去除表达量过低的基因
			df_count_used <- df_count_used[rowMeans(df_count_used)>countfilter,]
			dim(df_count)
			# df_count<-round(df_count)
			#DEseq2均一化
			colData <- data.frame(row.names=colnames(df_count_used), group_list)
			head(colData)
			dds <- DESeqDataSetFromMatrix(df_count_used, DataFrame(group_list), design= ~ group_list)
			
			dds_norm <- DESeq(dds,parallel = T) 
			sizeFactors(dds_norm)
			head(dds_norm)
			
			normalized_count<-as.data.frame(counts(dds_norm,normalized=TRUE))
			normalized_count[1:3,1:3]
			sum_col<-apply(normalized_count,2,sum)
			head(sum_col,30)
			write.table(sum_col,file = paste0('fig',plotid,'b.',projectid,".count_colsum.tsv"),row.names = T,sep='\t',quote = FALSE)
			summary(sum_col)
			
			rld<-vst(dds_norm,blind=F)
			log_normalized_count<-assay(rld)

			#log_normalized_count<-log2(normalized_count+1)
			normalized_count_out<-cbind(rownames(normalized_count),normalized_count)
			colnames(normalized_count_out)[1]<-'gene'
			normalized_count_out[1:3,1:3]
			write.table(normalized_count_out,file = paste0('fig',plotid,'c.',projectid,".DESeq2.ncount.tsv"),row.names = FALSE,sep='\t',quote = FALSE)
			write.table(log_normalized_count,file = paste0('fig',plotid,'c.',projectid,".DESeq2.ncount_vst.tsv"),row.names = T,sep='\t',quote = FALSE)
			
			group_mean=apply(normalized_count,1,function(x) aggregate(x,by=list(group_list),mean)$x)
			group_mean<-t(group_mean)
			head(group_mean)
			colnames(group_mean)<-levels(group_list)
			write.table(group_mean,file = paste0('fig',plotid,'c.',projectid,".DESeq2.ncount.mean.tsv"),row.names = T,sep='\t',quote = FALSE)
		  }
		  
		  if(1)   #########  绘制boxplot和logCount的boxplot  #########  绘制样本相关性矩阵
		  {
			nsample<-dim(df_count_used)[2]
			if(nsample>30)
			{
				samples<-sample(colnames(normalized_count),30)
				log_normalized_count_boxplot<-log_normalized_count[,samples]
				class(log_normalized_count_boxplot)
			}else{
			  log_normalized_count_boxplot<-log_normalized_count
			}
			log_normalized_count_boxplot<-as.data.frame(log_normalized_count_boxplot)
			ggplot2_boxplot_matrix(log_normalized_count_boxplot,group='sample',ytitle='log2(Normalized Count)',output=paste0('fig',plotid,'d.',projectid,'.boxplot'))
			ggcorrplot_1matrix(log_normalized_count_boxplot,output=paste0('fig',plotid,'e.',projectid,".sample_correlation"),save.data=T)
		  }
		  
		  if(1)  #########  绘制PCA图
		  {
			phenotype<-df_sample_used[,'group',drop=F]			
			sds<-apply(log_normalized_count,1,sd)
			log_normalized_count_pca<-log_normalized_count[sds>0,]
			ggplot2_pca(log_normalized_count_pca,phenotype,output=paste0('fig',plotid,'f.',projectid,".PCA"))
			ggplot2_pca_3d(log_normalized_count_pca,phenotype,output=paste0('fig',plotid,'f.',projectid,".PCA-3d"))
		  }
		  
		  res <- results(dds_norm,contrast=c("group_list",gtest,gcontrol))
		  res[1:3,]
		  padj_max<-max(res$padj,na.rm=TRUE)
		  res$padj[is.na(res$padj)]<-padj_max
		  if(1)   ######  输出差异表达基因表格
		  {
			dim(res)
			res$change = as.factor(
			  ifelse(res[[pvalue_column]] < pvalue_cutoff & abs(res[[fc_column]]) > logFC_cutoff,
					 ifelse(res[[fc_column]] > logFC_cutoff ,'UP','DOWN'),'NOT'))
			DEG_full<-merge(res,group_mean[,c(gcontrol,gtest)],by='row.names')       
			colnames(DEG_full)[1]<-'gene'
			write.table(DEG_full,file=paste0('fig',plotid,'g.',projectid,'.count.DESeq2.',comparison,'.deg.all.tsv'),row.names = F,sep='\t',quote = FALSE)
			DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]
			write.table(DEG_sig,file=paste0('fig',plotid,'g.',projectid,'.count.DESeq2.',comparison,'.deg.sig.tsv'),row.names = F,sep='\t',quote = FALSE)
			dim(DEG_sig)
			
				df_change<-table(DEG_sig$change)
			df_change<-as.data.frame(df_change)
			colnames(df_change)<-c('change','n_gene')
			write.table(df_change,file = paste0('fig',plotid,'g.',projectid,".DEG_gene_number.tsv"),row.names = F,sep='\t',quote = FALSE)
		  }
	}else{
		
			DEG_full <- read.table(file = paste0('fig',plotid,'g.',projectid,'.count.DESeq2.',comparison,'.deg.all.tsv'),header=T,sep='\t',check.names=FALSE)  
			DEG_full$change = as.factor(
			  ifelse(DEG_full[[pvalue_column]] < pvalue_cutoff & abs(DEG_full[[fc_column]]) > logFC_cutoff,
					 ifelse(DEG_full[[fc_column]] > logFC_cutoff ,'UP','DOWN'),'NOT'))
			write.table(DEG_full,file=paste0('fig',plotid,'g.',projectid,'.count.DESeq2.',comparison,'.deg.all.tsv'),row.names = F,sep='\t',quote = FALSE)
			DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]    
			write.table(DEG_sig,file=paste0('fig',plotid,'g.',projectid,'.count.DESeq2.',comparison,'.deg.sig.tsv'),row.names = F,sep='\t',quote = FALSE)
			dim(DEG_sig)
			df_change<-table(DEG_sig$change)
			df_change<-as.data.frame(df_change)
			colnames(df_change)<-c('change','n_gene')
			print(df_change)
			write.table(df_change,file = paste0('fig',plotid,'g.',projectid,".DEG_gene_number.tsv"),row.names = F,sep='\t',quote = FALSE)			
			log_normalized_count<-read.delim(paste0('fig',plotid,'c.',projectid,'.DESeq2.ncount_vst.tsv'),check.names=F,row.names=1)
	}
  if(nrow(DEG_sig)>1)   ######  绘制差异表达基因的热图（不区分PCG和LNC）
  {
    # log_normalized_count_heatmap<-read.delim('fig5.GSE153960.Heatmap.ALS_vs_Control.tsv',row.names=1)
    # colnames(log_normalized_count_heatmap)<-gsub('[.]','-',colnames(log_normalized_count_heatmap))
    # df_sample<-read.delim(paste0('fig0.',projectid,".sample_sheet_used.tsv"))
    # df_sample$group<-factor(df_sample$group,levels=c('Control','ALS'))
    phenotype<-df_sample_used[,'group',drop=F]
    head(phenotype)
    log_normalized_count_heatmap<-log_normalized_count[rownames(log_normalized_count) %in% DEG_sig$gene,]
    log_normalized_count_heatmap[1:3,1:3]
    
    ggheatmap_yingbio_2groups(log_normalized_count_heatmap,phenotype,geneid='gene',pheight=10,pw=10,output=paste0('fig',plotid,'h.',projectid,".Heatmap.",comparison),save.data=T)
  }
  
  if(1) ########## 火山图
  {
    head(DEG_full)
    ggplot2_volcano(DEG_full[,c(pvalue_column,fc_column)],comparison=comparison,output=paste0('fig',plotid,'v.',projectid,".volcano.",comparison),
    fc_threshold=2^logFC_cutoff,pvalue_threshold=pvalue_cutoff,ylab='padj',ymax=300)
  }
  
  if(nrow(DEG_sig)>1)  ######  进行GO和KEGG pathway分析
  {
    library(stats)
    if(0)
    {
      DEG_full <- read.table(file = paste0('fig4.',projectid,'.count.DESeq2.',comparison,'.deg.all.tsv'),header=T,sep='\t',check.names=FALSE)  
      DEG_sig <- read.table(file = paste0('fig4.',projectid,'.count.DESeq2.',comparison,'.deg.sig.tsv'),header=T,sep='\t',check.names=FALSE)  
    }
    #go_analysis=T
    head(DEG_sig)
    plotid<-'17'
    if(go_analysis)
    {
      aaa<-DEG_sig$gene[DEG_sig[[fc_column]]>0]
      aaa
      gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UP'),paste0('kegg_',comparison,'_UP')),species = species,go_pvalueCutoff=1,kegg_pvalueCutoff=1)
      
      aaa<-DEG_sig$gene[DEG_sig[[fc_column]]<0]
      gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_DOWN'),paste0('kegg_',comparison,'_DOWN')),species = species,go_pvalueCutoff=1,kegg_pvalueCutoff=1)
    }
  }
  
}


### input主要有两个矩阵：df_count和df_sample。
### 要求:1,df_tpm是log转化后的tpm，通常为log(TPM+1)，通常行名为基因symbol，列名为sampleid
### 要求:2,df_sample是样本分组矩阵，行名为sampleid。其中一列为group。矩阵中可以包含其它列。
### 要求:3,df_count的列名与df_sample的行名一一对应，且顺序一致（此条经优化后，可以不严格遵循）。
### input主要有两个矩阵：df_count和df_sample。
### 要求:1,df_tpm是log转化后的tpm，通常为log(TPM+1)，通常行名为基因symbol，列名为sampleid
### 要求:2,df_sample是样本分组矩阵，行名为sampleid。其中一列为group。矩阵中可以包含其它列。
### 要求:3,df_count的列名与df_sample的行名一一对应，且顺序一致（此条经优化后，可以不严格遵循）。
###  wilcox_test的fc列pvalue列fdr列：logFC、P.Value、adj.P.Val
rna_seq_tpm_wilcox_test<-function(df_tpm,df_sample,plotid=2,gcontrol='Normal',gtest='Tumor',
fc_cutoff=2,fc_column='log2FoldChange',pvalue_cutoff=0.05,pvalue_column='padj',
go_analysis=F,species='human',recutoff=F)
{
  
  df_sample_used<-df_sample[df_sample$group %in% c(gcontrol,gtest),,drop=F]
  class(df_sample_used)
  df_sample_used$group<-factor(df_sample_used$group,levels=c(gcontrol,gtest))
  df_sample_used<-df_sample_used[order(df_sample_used$group),,drop=F]
  
  samples_df_count<-colnames(df_tpm)
  samples_df_sample<-rownames(df_sample_used)
  samples_venn<-samples_df_sample[samples_df_sample %in% samples_df_count]
  df_sample_used<-df_sample_used[samples_venn,,drop=F]
  group_list<-df_sample_used$group
  df_group<-table(group_list)
  df_group<-as.data.frame(df_group)
  colnames(df_group)<-c('group','n_sample')
  write.table(df_group,file = paste0('fig',plotid,'a.',projectid,".sample_number.tsv"),row.names = F,sep='\t',quote = FALSE)
  
  if(df_group[1,2]>0 & df_group[2,2]>0)
  {
    print(df_group)
  }else{
    print('sample number good!')
    return()
  }
  
  df_tpm_used<-df_tpm[,samples_venn]
  logFC_cutoff<-log2(fc_cutoff)
  comparison<-paste0(gtest,'_vs_',gcontrol)
  if(!recutoff)
  {
    if(1)   #####  弱信号筛选和平均TPM计算
    {
      # 去除表达量过低的基因
      df_tpm_used <- df_tpm_used[rowMeans(df_tpm_used)>0.1,]
      dim(df_tpm_used)
      # df_count<-round(df_count)
      #DEseq2均一化
      
      group_mean=apply(2^df_tpm_used-1,1,function(x) aggregate(x,by=list(group_list),mean)$x)
      group_mean<-t(group_mean)
      head(group_mean)
      colnames(group_mean)<-levels(group_list)
      write.table(group_mean,file = paste0('fig',plotid,'b.',projectid,".tpm.mean.tsv"),row.names = T,sep='\t',quote = FALSE)
    }
    
    if(1)   #########  绘制boxplot和logCount的boxplot  #########  绘制样本相关性矩阵
    {
      nsample<-dim(df_tpm_used)[2]
      if(nsample>30)
      {
        samples<-sample(colnames(df_tpm_used),30)
        df_tpm_used_boxplot<-as.matrix(df_tpm_used[,samples])
        class(df_tpm_used_boxplot)
      }else{
        df_tpm_used_boxplot<-as.matrix(df_tpm_used)
      }
      df_tpm_used_boxplot[1:3,1:3]
      dim(df_tpm_used_boxplot)
      df_tpm_used_boxplot<-as.data.frame(df_tpm_used_boxplot)
      ggplot2_boxplot_matrix(df_tpm_used_boxplot,group='sample',output=paste0('fig',plotid,'c.',projectid,'.boxplot'))
      ggcorrplot_1matrix(df_tpm_used_boxplot,output=paste0('fig',plotid,'d.',projectid,".sample_correlation"),save.data=T)
    }
    
    if(1)  #########  绘制PCA图
    {
		class(df_tpm_used)
		phenotype<-df_sample_used[,'group',drop=F]
		sds<-apply(df_tpm_used,1,sd)
		df_tpm_pca<-df_tpm_used[sds>0,]
		ggplot2_pca(df_tpm_pca,phenotype,output=paste0('fig',plotid,'e.',projectid,".PCA"))
		ggplot2_pca_3d(df_tpm_pca,phenotype,output=paste0('fig',plotid,'e.',projectid,".PCA-3d"))
    }
    
    if(1)   ######  输出差异表达基因表格
    {
      
      # Run the Wilcoxon rank-sum test for each gene
      pvalues <- sapply(1:nrow(df_tpm_used),function(i){
        data<-cbind.data.frame(gene=as.numeric(t(df_tpm_used[i,])),group_list)
        p=wilcox.test(gene~group_list, data)$p.value
        return(p)
      })
      fdr=p.adjust(pvalues,method = "fdr")
      
      # Calculate fold-change for each gene
      dataCon1=df_tpm_used[,c(which(group_list==gcontrol))]
      dataCon2=df_tpm_used[,c(which(group_list==gtest))]
      foldChanges=rowMeans(dataCon2)-rowMeans(dataCon1)
      # Output results based on FDR threshold
      res<-data.frame(log2FoldChange=foldChanges, Pvalue=pvalues, padj=fdr)
      rownames(res)=rownames(df_tpm_used)
      head(res)
      res$change = as.factor(
        ifelse(res[[pvalue_column]] < pvalue_cutoff & abs(res[[fc_column]]) > logFC_cutoff,
               ifelse(res[[fc_column]] > logFC_cutoff ,'UP','DOWN'),'NOT'))
      DEG_full<-merge(res,group_mean,by='row.names')
      colnames(DEG_full)[1]<-'gene'
      write.table(DEG_full,file=paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.all.tsv'),row.names = F,sep='\t',quote = FALSE)
      DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]
      write.table(DEG_sig,file=paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.sig.tsv'),row.names = F,sep='\t',quote = FALSE)
      dim(DEG_sig)
      
      df_change<-table(DEG_sig$change)
      df_change<-as.data.frame(df_change)
      colnames(df_change)<-c('change','n_gene')
      print(df_change)
      write.table(df_change,file = paste0('fig',plotid,'f.',projectid,".DEG_gene_number.tsv"),row.names = F,sep='\t',quote = FALSE)
    }
  }else{
    DEG_full <- read.table(file = paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.all.tsv'),header=T,sep='\t',check.names=FALSE)  
    DEG_full$change = as.factor(
      ifelse(DEG_full[[pvalue_column]] < pvalue_cutoff & abs(DEG_full[[fc_column]]) > logFC_cutoff,
             ifelse(DEG_full[[fc_column]] > logFC_cutoff ,'UP','DOWN'),'NOT'))
    write.table(DEG_full,file=paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.all.tsv'),row.names = F,sep='\t',quote = FALSE)
    DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]
    write.table(DEG_sig,file=paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.sig.tsv'),row.names = F,sep='\t',quote = FALSE)
    dim(DEG_sig)
    df_change<-table(DEG_sig$change)
    df_change<-as.data.frame(df_change)
    colnames(df_change)<-c('change','n_gene')
    print(df_change)
    write.table(df_change,file = paste0('fig',plotid,'f.',projectid,".DEG_gene_number.tsv"),row.names = F,sep='\t',quote = FALSE)
  }
  
  if(nrow(DEG_sig)>1)   ######  绘制差异表达基因的热图（不区分PCG和LNC）
  {
    # log_normalized_count_heatmap<-read.delim('fig5.GSE153960.Heatmap.ALS_vs_Control.tsv',row.names=1)
    # colnames(log_normalized_count_heatmap)<-gsub('[.]','-',colnames(log_normalized_count_heatmap))
    # df_sample<-read.delim(paste0('fig0.',projectid,".sample_sheet_used.tsv"))
    # df_sample$group<-factor(df_sample$group,levels=c('Control','ALS'))
    phenotype<-df_sample_used[,'group',drop=F]
    head(phenotype)
    df_tpm_used_heatmap<-df_tpm_used[rownames(df_tpm_used) %in% DEG_sig$gene,]
    df_tpm_used_heatmap[1:3,1:3]
    
    ggheatmap_yingbio_2groups(df_tpm_used_heatmap,phenotype,geneid='gene',pheight=10,pw=10,output=paste0('fig',plotid,'h.',projectid,".Heatmap.",comparison),save.data=T)
  }
  
  if(1) ########## 火山图
  {
    head(DEG_full)
    ggplot2_volcano(DEG_full[,c(pvalue_column,fc_column)],comparison=comparison,output=paste0('fig',plotid,'g.',projectid,".volcano.",comparison),
    fc_threshold=2^logFC_cutoff,pvalue_threshold=pvalue_cutoff,ylab=pvalue_column,ymax=300)
  }
  
  if(nrow(DEG_sig)>1)  ######  进行GO和KEGG pathway分析
  {
    library(stats)
    if(0)
    {
      DEG_full <- read.table(file = paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.all.tsv'),header=T,sep='\t',check.names=FALSE)  
      DEG_sig <- read.table(file = paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.sig.tsv'),header=T,sep='\t',check.names=FALSE)  
    }
    #go_analysis=T
    head(DEG_sig)
    plotid<-'17'
    if(go_analysis)
    {
      aaa<-DEG_sig$gene[DEG_sig[[fc_column]]>0]
      aaa
      gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UP'),paste0('kegg_',comparison,'_UP')),species = species,go_pvalueCutoff=1,kegg_pvalueCutoff=1)
      
      aaa<-DEG_sig$gene[DEG_sig[[fc_column]]<0]
      gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_DOWN'),paste0('kegg_',comparison,'_DOWN')),species = species,go_pvalueCutoff=1,kegg_pvalueCutoff=1)
    }
  }
  
}


wilcox_test_matrix<-function(df_data,group_list)
{
	  if(class(group_list) =='factor')
	  {
		print('good')
	  }else{
		group_list<-factor(group_list,levels=unique(group_list))
	  }
		head(df_data)
	  group_mean=apply(df_data,1,function(x) aggregate(x,by=list(group_list),mean)$x)
	  group_mean<-t(group_mean)
	  head(group_mean)
	  colnames(group_mean)<-levels(group_list)
	  head(group_mean)
	  Pvalue <- sapply(1:nrow(df_data),function(i){
		data<-cbind.data.frame(gene=as.numeric(t(df_data[i,])),group_list)
		p=wilcox.test(gene~group_list, data)$p.value
		return(p)
	  })
	  FDR=p.adjust(pvalues,method = "fdr")
	  res<-cbind(group_mean,Pvalue,FDR)
	  head(res)
	  
	  return(res)
}
kruskal_test_matrix<-function(df_data,group_list)
{
	if(class(group_list) =='factor')
	{
	  print('good')
	}else{
	  group_list<-factor(group_list,levels=unique(group_list))
	}
	head(df_data)
	group_mean=apply(df_data,1,function(x) aggregate(x,by=list(group_list),mean)$x)
	group_mean<-t(group_mean)
	head(group_mean)
	colnames(group_mean)<-levels(group_list)
	head(group_mean)
	Pvalue <- sapply(1:nrow(df_data),function(i){
	  data<-cbind.data.frame(gene=as.numeric(t(df_data[i,])),group_list)
	  p=kruskal.test(gene~group_list, data)$p.value
	  return(p)
	})
	FDR=p.adjust(pvalues,method = "fdr")
	res<-cbind(group_mean,Pvalue,FDR)
	head(res)

	return(res)
}

### input主要有两个矩阵：df_tpm和df_sample。
### 要求:1,df_tpm是"未log转化的tpm"，通常行名为基因symbol，列名为sampleid
### 要求:2,df_sample是样本分组矩阵，行名为sampleid。其中一列为group。矩阵中可以包含其它列。
### 要求:3,df_tpm的列名与df_sample的行名一一对应，且顺序一致（此条经优化后，可以不严格遵循）。
###  limma的fc列pvalue列fdr列：logFC、P.Value、adj.P.Val
rna_seq_tpm_limma<-function(df_tpm,df_sample,plotid=2,gcontrol='Normal',gtest='Tumor',
                            fc_cutoff=2,fc_column='logFC',pvalue_cutoff=0.05,pvalue_column='adj.P.Val',
                            go_analysis=F,species='human',recutoff=F)
{
  # df_tpm<-df_count
  
  library(limma)
  df_sample_used<-df_sample[df_sample$group %in% c(gcontrol,gtest),]
  df_sample_used$group<-factor(df_sample_used$group,levels=c(gcontrol,gtest))
  df_sample_used<-df_sample_used[order(df_sample_used$group),]
  
  samples_df_count<-colnames(df_tpm)
  samples_df_sample<-rownames(df_sample_used)
  samples_venn<-samples_df_sample[samples_df_sample %in% samples_df_count]
  df_sample_used<-df_sample_used[samples_venn,]
  group_list<-df_sample_used$group
  df_group<-table(group_list)
  df_group<-as.data.frame(df_group)
  colnames(df_group)<-c('group','n_sample')
  write.table(df_group,file = paste0('fig',plotid,'a.',projectid,".sample_number.tsv"),row.names = F,sep='\t',quote = FALSE)
  
  if(df_group[1,2]>0 & df_group[2,2]>0)
  {
    print(df_group)
  }else{
    print('sample number good!')
    return()
  }
  
  df_tpm_used<-df_tpm[,samples_venn]
  logFC_cutoff<-log2(fc_cutoff)
  comparison<-paste0(gtest,'_vs_',gcontrol)
  if(!recutoff)
  {
    if(1)   #####  弱信号筛选和平均TPM计算
    {
      # 去除表达量过低的基因
      df_tpm_used <- df_tpm_used[rowMeans(df_tpm_used)>1,]
      dim(df_tpm_used)
      # df_count<-round(df_count)
      #DEseq2均一化
    }
    
    
    if(1)  ## voom标准化
    {
      design <- model.matrix(~0+group_list)
      colnames(design)=levels(group_list)
      rownames(design)=colnames(df_tpm_used)
      
      png(paste0('fig',plotid,'c.',projectid,'.limma.voom.png'),height=4096,width=4096,units="px",res=300)
      v <- voom(df_tpm_used,design, normalize="quantile",plot=T)
      dev.off()
      
      #voomE
      voom_matrix<-v$E
      head(voom_matrix)
      voom_matrix_out<-cbind(rownames(voom_matrix),voom_matrix)
      colnames(voom_matrix_out)[1]<-"gene_name"
      write.table(voom_matrix_out,file=paste0('fig',plotid,'c.',projectid,'.limma.voomE.tsv'),row.names = F,sep='\t',quote = FALSE)
      
      if(1)  ## 计算各组的平均强度，并构建线性模型
      {
        group_mean=apply(2^voom_matrix-1,1,function(x) aggregate(x,by=list(group_list),mean)$x)
        group_mean<-t(group_mean)
        colnames(group_mean)<-levels(group_list)
        write.table(group_mean,file = paste0('fig',plotid,'c.',projectid,".limma.voomE.mean.tsv"),row.names = T,sep='\t',quote = FALSE)
        head(group_mean)
      }
      
    }
    
    
    if(1)   #########  绘制boxplot和logCount的boxplot  #########  绘制样本相关性矩阵
    {
      nsample<-dim(df_tpm_used)[2]
      if(nsample>30)
      {
        samples<-sample(colnames(voom_matrix),30)
        df_tpm_used_boxplot<-voom_matrix[,samples]
        class(df_tpm_used_boxplot)
      }else{
        df_tpm_used_boxplot<-voom_matrix
      }
      df_tpm_used_boxplot<-as.data.frame(df_tpm_used_boxplot)
      ggplot2_boxplot_matrix(df_tpm_used_boxplot,group='sample',output=paste0('fig',plotid,'c.',projectid,'.boxplot'))
      ggcorrplot_1matrix(df_tpm_used_boxplot,output=paste0('fig',plotid,'d.',projectid,".sample_correlation"),save.data=T)
    }
    
    if(1)  #########  绘制PCA图
    {
      phenotype<-df_sample_used[,'group',drop=F]  
      sds<-apply(voom_matrix,1,sd)
      voom_matrix_used<-voom_matrix[sds>0,]
      ggplot2_pca(voom_matrix_used,phenotype,output=paste0('fig',plotid,'e.',projectid,".PCA"))
      ggplot2_pca_3d(voom_matrix_used,phenotype,output=paste0('fig',plotid,'e.',projectid,".PCA-3d"))
    }
    
    
    if(1)   ######  输出差异表达基因表格
    {
      
      if(1)  ## 差异分析
      {
        #比较矩阵，进行两两之间的比较
        fit <- lmFit(v, design)
        contrast <- makeContrasts(contrasts=paste0(gtest,"-",gcontrol),levels=design)
        fit2 <- contrasts.fit(fit, contrast)
        fit2=eBayes(fit2)
      }
      
      if(1)  ## 筛选出上调和下调基因
      {
        # topTable
        DEG = topTable(fit2, coef=1,number=Inf)
        # 去掉那些NA值
        DEG = na.omit(DEG)
        head(DEG,6)
        
        #筛选出上调和下调基因
        DEG$change = as.factor(
          ifelse(DEG[[pvalue_column]] < pvalue_cutoff & abs(DEG[[fc_column]]) > logFC_cutoff,
                 ifelse(DEG[[fc_column]] > logFC_cutoff ,'UP','DOWN'),'NOT'))
        
        group_mean<-group_mean[rownames(DEG),]
        df_tpm_used<-df_tpm_used[rownames(DEG),]
        colnames(df_tpm_used)<-paste0(colnames(df_tpm_used),'_tpm')
        DEG_full<-cbind(DEG,group_mean)
        DEG_full<-cbind(DEG_full,df_tpm_used)            
        output_table_with_rowname(DEG_full,name='gene',output=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.all'))
        DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]
        print(dim(DEG_sig))
        output_table_with_rowname(DEG_sig,name='gene',output=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.sig'))
      }          
      df_change<-table(DEG_sig$change)
      df_change<-as.data.frame(df_change)
      colnames(df_change)<-c('change','n_gene')
      print(df_change)
      write.table(df_change,file = paste0('fig',plotid,'f.',projectid,".DEG_gene_number.tsv"),row.names = F,sep='\t',quote = FALSE)
    }
  }else{
    DEG_full <- read.delim(file = paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.all.tsv'),row.names=1,header=T,sep='\t',check.names=FALSE)  
    head(DEG_full)
    DEG_full$change = as.factor(
      ifelse(DEG_full[[pvalue_column]] < pvalue_cutoff & abs(DEG_full[[fc_column]]) > logFC_cutoff,
             ifelse(DEG_full[[fc_column]] > logFC_cutoff ,'UP','DOWN'),'NOT'))
    output_table_with_rowname(DEG_full,name='gene',output=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.all'))
    DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]    
    output_table_with_rowname(DEG_sig,name='gene',output=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.sig'))
    dim(DEG_sig)
    head(DEG_sig)
    df_change<-table(DEG_sig$change)
    df_change<-as.data.frame(df_change)
    colnames(df_change)<-c('change','n_gene')
    print(df_change)
    write.table(df_change,file = paste0('fig',plotid,'f.',projectid,".DEG_gene_number.tsv"),row.names = F,sep='\t',quote = FALSE)        
    voom_matrix<-read.delim(paste0('fig',plotid,'c.',projectid,'.limma.voomE.tsv'),check.names=F,row.names=1)
    
  }
  
  if(nrow(DEG_sig)>1)   ######  绘制差异表达基因的热图（不区分PCG和LNC）
  {
    phenotype<-df_sample_used[,'group',drop=F]
    voom_matrix_heatmap<-voom_matrix[rownames(voom_matrix) %in% rownames(DEG_sig),]
    head(voom_matrix_heatmap)
    dim(voom_matrix_heatmap)
    nr<-nrow(voom_matrix_heatmap)
    ns<-ncol(voom_matrix_heatmap)
    show_r<-F
    show_c<-F
    if(ns<=20)
    {show_c=T}
    if(nr<=100)
    {show_r=T}
    
    ggheatmap_yingbio_2groups(voom_matrix_heatmap,phenotype,geneid='gene',show_colname = show_c,show_rowname = show_r,pheight=10,pw=10,output=paste0('fig',plotid,'h.',projectid,".Heatmap.",comparison),save.data=T)
  }
  
  if(1) ########## 火山图
  {
    head(DEG_full)
    ggplot2_volcano(DEG_full[,c(pvalue_column,fc_column)],comparison=comparison,output=paste0('fig',plotid,'g.',projectid,".volcano.",comparison),
    fc_threshold=2^logFC_cutoff,pvalue_threshold=pvalue_cutoff,ylab=pvalue_column,ymax=300)
  }
  
  if(nrow(DEG_sig)>1)  ######  进行GO和KEGG pathway分析
  {
    library(stats)
    if(0)
    {
      DEG_full <- read.table(file = paste0('fig',plotid,'.f.',projectid,'.tpm.wilcox.',comparison,'.deg.all.tsv'),header=T,sep='\t',check.names=FALSE)  
      DEG_sig <- read.table(file = paste0('fig',plotid,'.f.',projectid,'.tpm.wilcox.',comparison,'.deg.sig.tsv'),header=T,sep='\t',check.names=FALSE)  
    }
    #go_analysis=T
    head(DEG_sig)
    plotid<-'17'
    if(go_analysis)
    {
      aaa<-rownames(DEG_sig[DEG_sig[[fc_column]]>0,])
      aaa
      gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UP'),paste0('kegg_',comparison,'_UP')),species = species,go_pvalueCutoff=1,kegg_pvalueCutoff=1)
      
      aaa<-rownames(DEG_sig[DEG_sig[[fc_column]]<0,])
      gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_DOWN'),paste0('kegg_',comparison,'_DOWN')),species = species,go_pvalueCutoff=1,kegg_pvalueCutoff=1)
    }
  }
  
}

rna_seq_gsva_custom<- function(expr=expr,genesets=gsc,tf_plot=T,phenotype=phenotype,group.by='group',title='Value',plotid='14'){
    library(GSVA)
    library(msigdbr)
    library(pheatmap)
    library(ggplot2)
    library(reshape2)
    
    expr<-as.matrix(expr)      
    gsvaPar <- ssgseaParam(expr, genesets,minSize=2)
    gsva.res <- gsva(gsvaPar, verbose=FALSE)
    # gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
    gsva.df <- data.frame(gsva.res, check.names = F)
    write.table(gsva.df, paste0("fig",plotid,".gsva.tsv"),sep='\t',row.names = T,quote = FALSE)
    ngene<-nrow(gsva.df)
    if(tf_plot)
    {
		
        ggheatmap_yingbio_2groups_for_gsva(gsva.df,phenotype,geneid='genesets',pheight=6,
                                           output=paste0('fig',plotid,'a.gsva.heatmap'),save.data=T)
        df_gsva<-as.matrix(gsva.df)
        df_plot<-melt(df_gsva)
        head(df_plot)
        colnames(df_plot)<-c('celltype','sample','score')
        df_plot_used<-merge(df_plot,phenotype,by.x='sample',by.y=0)
        head(df_plot_used)
        if(ngene>1)
		{
			ggplot2_boxplot_gsva(df_plot_used,groupby='celltype',score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'b.',projectid,".gsva.boxpot")) 
			ggplot2_boxplot_paired(df_plot_used,groupby=c('group','celltype'),score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'c.',projectid,".gsva.boxplot_paired"))
		}else{
			ymax<-max(df_plot_used[,3])
			ggboxplot_2groups_comparison(df_plot_used,group='group',value='score',title=title,pw=4,
                                       output=paste0('fig',plotid,'b.boxplot.comparison.',title),label_y=ymax+1)
		}
    }
    return(gsva.df)
}

### df_score<-rna_seq_CRDscore(expr=df_logtpm,genes=hubgenes,name='RloopScore',plotid=plotid,tf_plot=F)
### 结果是1行的矩阵，列名为样本名
rna_seq_CRDscore<- function(expr=expr,genes=hubgenes,name='CRDscore',tf_plot=T,phenotype=phenotype,group.by='group',title='Value',plotid='14'){
library(pheatmap)
library(ggplot2)
library(reshape2)
library(CRDscore)
expr <- as.data.frame(expr)
score <- cal_CRDscore(expr = expr, n.bins = 50, circadians = genes, study.type = "bulk_RNAseq")
score = as.data.frame(t(score))
rownames(score)<-name
dim(score)
write.table(score, paste0("fig",plotid,".CRDscore.tsv"),sep='\t',row.names = T,quote = FALSE)     
		if(tf_plot)
		{
		  score<-as.matrix(score)
		  df_plot<-melt(score)
		  head(df_plot)
		  colnames(df_plot)<-c('celltype','sample','score')
		  df_plot_used<-merge(df_plot,phenotype,by.x='sample',by.y=0)
		  head(df_plot_used)
		  ggplot2_boxplot_gsva(df_plot_used,groupby='celltype',score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'b.',projectid,".gsva.boxpot")) 
		  ggplot2_boxplot_paired(df_plot_used,groupby=c('group','celltype'),score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'c.',projectid,".gsva.boxplot_paired"))
		}
return(score)
}


rna_seq_groupby_score<-function(df_data,groupby='score',cutoff='median',levels=c('low','high'))
{
	df_group<-data.frame(df_data[,groupby,drop=F])
	colnames(df_group)<-groupby
  if(cutoff=='median')
  {
    cutoff<-median(df_group[[groupby]])
  }
  
  df_group$group<-levels[1]
  df_group[df_group[[groupby]]>cutoff,'group']<-levels[2]
  return(df_group)
}
 
 
### 芯片数据使用limma进行差异分析的函数
### input主要有两个矩阵：df_tpm和df_sample。
### 要求:1,df_tpm是标准化后的Intensity，通常为log(Raw+1)再进行芯片间的quantile，通常行名为基因symbol，列名为sampleid
### 要求:2,df_sample是样本分组矩阵，行名为sampleid。其中一列为group。矩阵中可以包含其它列。
### 要求:3,df_count的列名与df_sample的行名一一对应，且顺序一致（此条经优化后，可以不严格遵循）。
###  limma的fc列pvalue列fdr列：logFC、P.Value、adj.P.Val
microarray_limma_deg<-function(df_tpm,df_sample,plotid=2,gcontrol='Normal',gtest='Tumor',
                                       fc_cutoff=2,fc_column='logFC',pvalue_cutoff=0.05,pvalue_column='adj.P.Val',
                                       go_analysis=F,species='human',recutoff=F)
{
  df_sample_used<-df_sample[df_sample$group %in% c(gcontrol,gtest),]
  df_sample_used$group<-factor(df_sample_used$group,levels=c(gcontrol,gtest))
  df_sample_used<-df_sample_used[order(df_sample_used$group),]
  samples_df_count<-colnames(df_tpm)
  samples_df_sample<-rownames(df_sample_used)
  samples_venn<-samples_df_sample[samples_df_sample %in% samples_df_count]
  df_sample_used<-df_sample_used[samples_venn,]
  group_list<-df_sample_used$group
  df_group<-table(group_list)
  df_group<-as.data.frame(df_group)
  colnames(df_group)<-c('group','n_sample')
  write.table(df_group,file = paste0('fig',plotid,'a.',projectid,".sample_number.tsv"),row.names = F,sep='\t',quote = FALSE)
  
  if(df_group[1,2]>0 & df_group[2,2]>0)
  {
	print(df_group)
  }else{
	print('sample number good!')
	return()
  }
  
  df_tpm_used<-df_tpm[,samples_venn]
  logFC_cutoff<-log2(fc_cutoff)
  comparison<-paste0(gtest,'_vs_',gcontrol)
  if(!recutoff)
  {
	if(1)   #####  弱信号筛选和平均TPM计算
	{
	  # 去除表达量过低的基因
	  #df_tpm_used <- df_tpm_used[rowMeans(df_tpm_used)>0.1,]
	  dim(df_tpm_used)
	  # df_count<-round(df_count)
	  #DEseq2均一化
	  
	  group_mean=apply(2^df_tpm_used-1,1,function(x) aggregate(x,by=list(group_list),mean)$x)
	  group_mean<-t(group_mean)
	  head(group_mean)
	  colnames(group_mean)<-levels(group_list)
	  write.table(group_mean,file = paste0('fig',plotid,'b.',projectid,".Intensity.mean.tsv"),row.names = T,sep='\t',quote = FALSE)
	}
	
	if(1)   #########  绘制boxplot和logCount的boxplot  #########  绘制样本相关性矩阵
	{
	  nsample<-dim(df_tpm_used)[2]
	  if(nsample>30)
	  {
		samples<-sample(colnames(df_tpm_used),30)
		df_tpm_used_boxplot<-df_tpm_used[,samples]
		class(df_tpm_used_boxplot)
	  }else{
		df_tpm_used_boxplot<-df_tpm_used
	  }
	  df_tpm_used_boxplot<-as.data.frame(df_tpm_used_boxplot)
	  ggplot2_boxplot_matrix(df_tpm_used_boxplot,group='sample',output=paste0('fig',plotid,'c.',projectid,'.boxplot'))
	  ggcorrplot_1matrix(df_tpm_used_boxplot,output=paste0('fig',plotid,'d.',projectid,".sample_correlation"),save.data=T)
	}
	
	if(1)  #########  绘制PCA图
	{
	  phenotype<-df_sample_used[,'group',drop=F] 
      sds<-apply(df_tpm_used,1,sd)
      df_tpm_used_pca<-df_tpm_used[sds>0,]
	  ggplot2_pca(df_tpm_used_pca,phenotype,output=paste0('fig',plotid,'e.',projectid,".PCA"))
	  ggplot2_pca_3d(df_tpm_used_pca,phenotype,output=paste0('fig',plotid,'e.',projectid,".PCA-3d"))
	}
	
	
	if(1)   ######  输出差异表达基因表格
	{
	  
	  design <- model.matrix(~0+group_list)
	  colnames(design) <- levels(group_list)
	  fit <- lmFit(df_tpm,design)
	  
	  png(paste0('fig',plotid,'fa.',projectid,"plotSA.fit.png"),height=4096,width=4096,units="px",res=300)
	  plotSA(fit, main="Final model: Mean-variance trend")
	  dev.off()
	  
	  cont.matrix<-makeContrasts(contrasts=paste0(gtest,"-",gcontrol),levels=design)
	  fit2 <- contrasts.fit(fit, cont.matrix)
	  fit2 <- eBayes(fit2)
	  png(paste0('fig',plotid,'fb.',projectid,"plotSA.fit2.png"),height=4096,width=4096,units="px",res=300)
	  plotSA(fit2, main="Final model: Mean-variance trend")
	  dev.off()
	  DEG=topTable(fit2,adjust='fdr',number=Inf)
	  print(head(DEG))
	  DEG$change = as.factor(
		ifelse(DEG[[pvalue_column]] < pvalue_cutoff & abs(DEG[[fc_column]]) > logFC_cutoff,
			   ifelse(DEG[[fc_column]] > logFC_cutoff ,'UP','DOWN'),'NOT'))
	  print(table(DEG$change))
	  
	  group_mean<-group_mean[rownames(DEG),]
	  df_tpm_add<-df_tpm_used[rownames(DEG),]        
	  colnames(df_tpm_add)<-paste0(colnames(df_tpm_add),'_Normalized')
	  DEG_full<-cbind(DEG,group_mean)
	  DEG_full<-cbind(DEG_full,df_tpm_add)            
	  output_table_with_rowname(DEG_full,name='gene',output=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.all'))
	  DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]
	  print(dim(DEG_sig))
	  output_table_with_rowname(DEG_sig,name='gene',output=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.sig'))
	  
	  df_change<-table(DEG_sig$change)
	  df_change<-as.data.frame(df_change)
	  colnames(df_change)<-c('change','n_gene')
	  write.table(df_change,file = paste0('fig',plotid,'f.',projectid,".DEG_gene_number.tsv"),row.names = F,sep='\t',quote = FALSE)
	}
  }else
  {
	DEG_full <- read.delim(file = paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.all.tsv'),row.names=1,header=T,sep='\t',check.names=FALSE)  
	head(DEG_full)
	DEG_full$change = as.factor(
	  ifelse(DEG_full[[pvalue_column]] < pvalue_cutoff & abs(DEG_full[[fc_column]]) > logFC_cutoff,
			 ifelse(DEG_full[[fc_column]] > logFC_cutoff ,'UP','DOWN'),'NOT'))
	output_table_with_rowname(DEG_full,name='gene',output=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.all'))
	DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]    
	output_table_with_rowname(DEG_sig,name='gene',output=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.sig'))
	dim(DEG_sig)
	head(DEG_sig)
	df_change<-table(DEG_sig$change)
	df_change<-as.data.frame(df_change)
	colnames(df_change)<-c('change','n_gene')
	print(df_change)
	write.table(df_change,file = paste0('fig',plotid,'f.',projectid,".DEG_gene_number.tsv"),row.names = F,sep='\t',quote = FALSE)        
  }
  
  
  
  
  if(nrow(DEG_sig)>1)   ######  绘制差异表达基因的热图（不区分PCG和LNC）
  {
	phenotype<-df_sample_used[,'group',drop=F]
	head(phenotype)
	df_tpm_used_heatmap<-df_tpm_used[rownames(df_tpm_used) %in% rownames(DEG_sig),]
	df_tpm_used_heatmap[1:3,1:3]
	
	nr<-nrow(df_tpm_used_heatmap)
	ns<-ncol(df_tpm_used_heatmap)  
	show_r<-F
	show_c<-F
	if(ns<=20)
	{show_c=T}
	if(nr<=100)
	{show_r=T}
	
	ggheatmap_yingbio_2groups(df_tpm_used_heatmap,phenotype,geneid='gene',show_colname = show_c,show_rowname = show_r,pheight=10,pw=10,output=paste0('fig',plotid,'g.',projectid,".Heatmap.",comparison),save.data=T)
  }
  
  if(1) ########## 火山图
  {
	head(DEG_full)
	ggplot2_volcano(DEG_full[,c(pvalue_column,fc_column)],comparison=comparison,output=paste0('fig',plotid,'g.',projectid,".volcano.",comparison),
	fc_threshold=2^logFC_cutoff,pvalue_threshold=pvalue_cutoff,ylab=pvalue_column,ymax=300)
  }
  
  if(nrow(DEG_sig)>1)  ######  进行GO和KEGG pathway分析
  {
	library(stats)
	if(0)
	{
	  DEG_full <- read.table(file = paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.all.tsv'),header=T,sep='\t',check.names=FALSE)  
	  DEG_sig <- read.table(file = paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.sig.tsv'),header=T,sep='\t',check.names=FALSE)  
	}
	#go_analysis=T
	head(DEG_sig)
	plotid<-'17'
	if(go_analysis)
	{
	  aaa<-rownames(DEG_sig[DEG_sig[[fc_column]]>0,])
	  aaa
	  gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UP'),paste0('kegg_',comparison,'_UP')),species = species,go_pvalueCutoff=1,kegg_pvalueCutoff=1)
	  
	  aaa<-rownames(DEG_sig[DEG_sig[[fc_column]]<0,])
	  gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_DOWN'),paste0('kegg_',comparison,'_DOWN')),species = species,go_pvalueCutoff=1,kegg_pvalueCutoff=1)
	}
  }
  
}



#####  行是样本，列是细胞类型
cibersort_plot<-function(df_infiltration,phenotype,groupby='group',plotid=10)
{
	library(ggplot2)
library(pheatmap)
    head(df_infiltration)
    if(1)#堆积柱状图
    {
      data=t(df_infiltration)
      col=rainbow(nrow(data),s=1,v=1)
      #绘制柱状图
      pdf(paste0('fig',plotid,"a.cibersort",".barplot2.pdf"),height=10,width=22)
      par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
      a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8,axisnames=F)
      a2=axis(2,tick=F,labels=F)
      axis(2,a2,paste0(a2*100,"%"))
      #axis(1,a1,labels=F)
      par(srt=60,xpd=T);
      #text(a1,-0.02,colnames(data),adj=1,cex=0.6);
      par(srt=0)
      ytick2 = cumsum(data[,ncol(data)])
      ytick1 = c(0,ytick2[-length(ytick2)])
      legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
      dev.off()
      
      tiff(paste0('fig',plotid,"a.cibersort",".barplot2.tiff"),height=10*300,width=22*300,units="px",res=300,compression='lzw')
      par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
      a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8,axisnames=F)
      a2=axis(2,tick=F,labels=F)
      axis(2,a2,paste0(a2*100,"%"))
      #axis(1,a1,labels=F)
      par(srt=60,xpd=T);
      #text(a1,-0.02,colnames(data),adj=1,cex=0.6);
      par(srt=0)
      ytick2 = cumsum(data[,ncol(data)])
      ytick1 = c(0,ytick2[-length(ytick2)])
      legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
      dev.off()
    }
    if(1)   ##### 使用生信自学网的脚本绘制相关性图形
    {
      library(corrplot)  
      #绘制相关性图形
      pdf(paste0('fig',plotid,"a.cibersort",".corrplot.pdf"),height=13,width=13)              #保存图片的文件名称
      par(oma=c(0.5,1,1,1.2))
      data=df_infiltration
      head(data)
      dim(data)
      data=data[,colMeans(data)>0]
      M=cor(data)
      corrplot(M,
               order="hclust",
               method = "color",
               addCoef.col = "black",
               diag = TRUE,
               tl.col="black",
               #col=colorRampPalette(c("blue", "white", "red"))(50)
      )
      dev.off()
      
      tiff(paste0('fig',plotid,"a.cibersort",".corrplot.tiff"),height=13*300,width=13*300,units="px",res=300,compression='lzw')              #保存图片的文件名称
      par(oma=c(0.5,1,1,1.2))
      data=df_infiltration
      head(data)
      dim(data)
      data=data[,colMeans(data)>0]
      M=cor(data)
      corrplot(M,
               order="hclust",
               method = "color",
               addCoef.col = "black",
               diag = TRUE,
               tl.col="black",
               #col=colorRampPalette(c("blue", "white", "red"))(50)
      )
      dev.off()
    }
    if(1)  ##### 使用生信自学网的脚本绘制热图
    {
      library(stringr)
      data=t(df_infiltration)
      phenotype_used<-phenotype[,groupby,drop=F]
      groups<-levels(phenotype_used[,groupby])
      
      annotation_col=phenotype_used
      head(annotation_col)
      
      ann_colors = list(groupby= c('a' ="#ce6700",'b'="#1b8c19") )
      class(ann_colors$groupby)
      names(ann_colors$groupby)<-groups
      
      pdf(paste0('fig',plotid,"c.cibersort",".heatmap.pdf"),height=6,width=12)
      pheatmap::pheatmap(data, 
                         annotation=annotation_col,
                         annotation_colors = ann_colors,
                         color = colorRampPalette(c('black',"red", "yellow"))(50),
                         show_colnames=F,
                         cluster_cols =F,
                         fontsize = 8,
                         fontsize_row=8,
                         fontsize_col=5)
      dev.off()
      tiff(paste0('fig',plotid,"c.cibersort",".heatmap.tiff"),height=6*300,width=12*300,units="px",res=300,compression='lzw')
      pheatmap::pheatmap(data, 
                         annotation=annotation_col,
                         annotation_colors = ann_colors,
                         color = colorRampPalette(c('black',"red", "yellow"))(50),
                         show_colnames=F,
                         cluster_cols =F,
                         fontsize = 8,
                         fontsize_row=8,
                         fontsize_col=5)
      dev.off()
    }
    if(1)   #### boxplot 像tcga一样带tumor和normal两种配对的样本。
    {
      data<-t(df_infiltration)
      phenotype_used<-phenotype[,groupby,drop=F]
      #####   data的行为celltype
      ggplot2_boxplot_paired_cibersort(data, phenotype_used,groupby=groupby,plotid=paste0('fig',plotid,'e'),ymax=0.9,label_y=0.75)
    }
}


#### 每种细胞类型分开作boxpot，并且比较两组间的差异。   
cibersort_celltype_boxplot<-function(df_infiltration,phenotype,groupby='group',plotid=10)
{
	library(ggplot2)
	library(pheatmap)
 
  data=df_infiltration
  head(data)
  dim(data)
  class(data)
  head(phenotype)
  dim(phenotype)
  for(i in 1:ncol(data))
  {
     #i<-6
     mydata<-cbind(data[,i,drop=F],phenotype[,groupby,drop=F])
     head(mydata)
     class(mydata)
     dim(mydata)
     
     groups<-levels(phenotype[,groupby])
     
     mydata_low<-mydata[mydata[,groupby]==groups[1],]
     dim(mydata_low)
     mydata_high<-mydata[mydata[,groupby]==groups[1],]
     meana<-mean(mydata_low[,1])
     meanb<-mean(mydata_high[,1])
     celltype<-colnames(data)[i]
     #celltype<-gsub('[.]',' ',celltype)
     #celltype<-gsub('CD8  ','CD8+ ',celltype)
     #celltype<-gsub('CD4  ','CD4+ ',celltype)
     ymax<-max(mydata[,1])
     ggboxplot_2groups_comparison(mydata,group=groupby,value=colnames(data)[i],title=celltype,output=paste0('fig',plotid,'.boxplot.',i),label_y=ymax+0.02,pw=4)
  }

}





seurat_qc_plot <- function(pbmc,plotid='01')
{
library(reshape2)
Idents(pbmc) <- "orig.ident" 
#绘制基因特征的小提琴图
violin<-VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","percent.rb"),pt.size=0, ncol = 3)
ggsave(paste0("fig",plotid,"a1.Violin-percentmt.tiff"), plot = violin, width = 12, height = 12,compression='lzw') 
ggsave(paste0("fig",plotid,"a1.Violin-percentmt.pdf"), plot = violin, width = 12, height = 12)


#测序深度的相关性图
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
plot3 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.HB",pt.size=1.5)
plot4 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.rb",pt.size=1.5)
pearplot <-CombinePlots(plots = list(plot1, plot2,plot3, plot4))
ggsave(paste0("fig",plotid,"a2.lm-Ngene_vs_Ncount.tiff"), plot = pearplot, width = 12, height = 6,compression='lzw') 
ggsave(paste0("fig",plotid,"a2.lm-Ngene_vs_Ncount.pdf"), plot = pearplot, width = 12, height = 6)

if(1)
{
library(cowplot)
violin<-plot_grid(VlnPlot(pbmc, features=c('nFeature_RNA'),pt.size=0, ncol = 1)+theme(axis.text.x = element_text(angle = 70, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE),
          VlnPlot(pbmc, features=c('nCount_RNA'),pt.size=0, ncol = 1)+theme(axis.text.x = element_text(angle = 70, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE),
          VlnPlot(pbmc, features=c('percent.mt'),pt.size=0, ncol = 1)+theme(axis.text.x = element_text(angle = 70, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE), ncol = 1)
ggsave(paste0("fig",plotid,"a1b.Violin-percentmt.angle.tiff"), plot = violin, width = 12, height = 18,compression='lzw') 
ggsave(paste0("fig",plotid,"a1b.Violin-percentmt.angle.pdf"), plot = violin, width = 12, height = 18)          

plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
pearplot <-CombinePlots(plots = list(plot1, plot2),ncol=1)
ggsave(paste0("fig",plotid,"a2b.lm-Ngene_vs_Ncount.tiff"), plot = pearplot, width = 12, height = 12,compression='lzw') 
ggsave(paste0("fig",plotid,"a2b.lm-Ngene_vs_Ncount.pdf"), plot = pearplot, width = 12, height = 12)
}


pbmc$log10GenesPerUMI <- log10(pbmc$nFeature_RNA) / log10(pbmc$nCount_RNA)
metadata <- pbmc@meta.data
df_sample_count<-as.data.frame(table(metadata$orig.ident))
colnames(df_sample_count)<-c('sample','Ncells')
write.table(df_sample_count, file=paste0("fig",plotid,".sample.Ncells.tsv"), quote=F, sep="\t", row.names=F)

if('group' %in% colnames(metadata))
{
df_sample_count<-as.data.frame(table(metadata$group))
colnames(df_sample_count)<-c('sample','Ncells')
write.table(df_sample_count, file=paste0("fig",plotid,".group.Ncells.tsv"), quote=F, sep="\t", row.names=F)

plota<-ggplot(metadata,aes(x=group, fill=group)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("Ncells")
ggsave(paste0("fig",plotid,"b1.group.Ncells.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"b1.group.Ncells.pdf"), plot = plota, width = 8, height = 8)
}

plota<-ggplot(metadata,aes(x=orig.ident, fill=orig.ident)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("Ncells")
ggsave(paste0("fig",plotid,"b1.Ncells.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"b1.Ncells.pdf"), plot = plota, width = 8, height = 8)

plotb<-ggplot(metadata,aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Ncount density") +
  	geom_vline(xintercept = 300)
ggsave(paste0("fig",plotid,"b2.Ncount_density.tiff"), plot = plotb, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"b2.Ncount_density.pdf"), plot = plotb, width = 8, height = 8)

plotc<-ggplot(metadata,aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 50)
ggsave(paste0("fig",plotid,"b3.Ngene_density.tiff"), plot = plotc, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"b3.Ngene_density.pdf"), plot = plotc, width = 8, height = 8)


plotd<-ggplot(metadata,aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("nFeature_RNA")
ggsave(paste0("fig",plotid,"b4.boxplot-Ngene.tiff"), plot = plotd, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"b4.boxplot-Ngene.pdf"), plot = plotd, width = 8, height = 8)


plote<-ggplot(metadata,aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~orig.ident)
ggsave(paste0("fig",plotid,"b5.lm-Ngene.tiff"), plot = plote, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"b5.lm-Ngene.pdf"), plot = plote, width = 8, height = 8)

plotf<-ggplot(metadata,aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.01)
ggsave(paste0("fig",plotid,"b6.mt_density.tiff"), plot = plotf, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"b6.mt_density.pdf"), plot = plotf, width = 8, height = 8)

plotg<-ggplot(metadata,aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)
ggsave(paste0("fig",plotid,"b7.GenesPerUMI.tiff"), plot = plotg, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"b7.GenesPerUMI.pdf"), plot = plotg, width = 8, height = 8)
}

#seurat_plot_deg_multicelltypes(pbmc,plotid='04',groupby='celltype',splitby='group',
#                               group_test='CD',group_control='control',
#                               deg.calculte=T,logFCfilter=1,adjPvalFilter=0.05)
seurat_plot_deg_multicelltypes<-function(pbmc,deg.calculte=T,groupby='celltype',splitby='group',group_test='CD',group_control='control',plotid='04',logFCfilter=1,adjPvalFilter=0.05)
{
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(reshape2)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)  
  library(magrittr)
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  
  library(viridis)
  library(future)
  plan("multicore", workers = 6)
  options(future.globals.maxSize = 2000 * 1024^2)
  
  
  Idents(pbmc) <- groupby
  type=unique(pbmc@meta.data[[groupby]])
  ntype<-length(type)
  if(deg.calculte)
  {
        r.deg=data.frame()
        for (i in 1:ntype) {
          Idents(pbmc) <- groupby
          deg=FindMarkers(pbmc,ident.1 = group_test,ident.2 = group_control,
                          group.by = splitby,subset.ident =type[i])
          write.table(deg, file =paste0("fig",plotid,".",type[i],'deg.tsv'),sep='\t')
          
          deg[[groupby]]=type[i]
          deg$unm=i-1
          deg$gene=rownames(deg)
          r.deg=rbind(deg,r.deg)
        }
        head(r.deg)
        table(r.deg$unm)
  }else{
          r.deg=read.delim(paste0("fig",plotid,".deg.multicelltypes.tsv"),header = T)
  }

  #############################################################################
  # 根据自己计算的marker基因数量确定log2FC的阈值，这里先定为1.5
  r.deg <- subset(r.deg, p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter)
  r.deg$threshold <- as.factor(ifelse(r.deg$avg_log2FC > 0 , 'Up', 'Down'))
  dim(r.deg)
  r.deg$adj_p_signi <- as.factor(ifelse(r.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
  r.deg$thr_signi <- paste0(r.deg$threshold, "_", r.deg$adj_p_signi)
  r.deg$unm %<>% as.vector(.) %>% as.numeric(.)
  ##自定义显示想要展示的基因名
  ##这里挑选log2FC为top5的基因进行展示
  top_up_label <- r.deg %>% 
    subset(., threshold%in%"Up") %>% 
    group_by(unm) %>% 
    top_n(n = 5, wt = avg_log2FC) %>% 
    as.data.frame()
  top_down_label <- r.deg %>% 
    subset(., threshold %in% "Down") %>% 
    group_by(unm) %>% 
    top_n(n = -5, wt = avg_log2FC) %>% 
    as.data.frame()
  top_label <- rbind(top_up_label,top_down_label)
  top_label$thr_signi %<>% 
    factor(., levels = c("Up_Highly","Down_Highly","Up_Lowly","Down_Lowly"))
  # 保存到文件，便于小白套用格式
  write.table(r.deg, file =paste0("fig",plotid,".deg.multicelltypes.tsv"),sep='\t')
  ##也可以基于output_pbmc.markers.csv文件，手动挑选出想要标注名字的基因，
  #例如标注参与某一通路的基因，然后将文件命名为easy_input_label.csv
  colnames(r.deg)
  ### 准备绘制暗灰色背景所需数据 
  background_position <- r.deg %>%
    dplyr::group_by(unm) %>%
    dplyr::summarise(Min = min(avg_log2FC) - 0.2, Max = max(avg_log2FC) + 0.2) %>%
    as.data.frame()
  ## `summarise()` ungrouping output (override with `.groups` argument)
  background_position$unm %<>% as.vector(.) %>% as.numeric(.)
  background_position$start <- background_position$unm - 0.4
  background_position$end <- background_position$unm + 0.4
  
  ### 准备绘制中间区域cluster彩色bar所需数据
  cluster_bar_position <- background_position
  cluster_bar_position$start <- cluster_bar_position$unm - 0.5
  cluster_bar_position$end <- cluster_bar_position$unm + 0.5
  cluster_bar_position$unm %<>% 
    factor(., levels = c(0:max(as.vector(.))))
  
  ## 设置填充颜色
  cols_thr_signi <- c("Up_Highly" = "#d7301f",
                      "Down_Highly" = "#225ea8",
                      "Up_Lowly" = "black",
                      "Down_Lowly" = "black")
  if(ntype<=7){
    cols_cluster <- c("0" = "#35978f",
                      "1" = "#8dd3c7",
                      "2" = "#ffffb3",
                      "3" = "#bebada",
                      "4" = "#fb8072",
                      "5" = "#80b1d3",
                      "6" = "#fdb462",
                      "7" = "#b3de69")
  }else if(ntype<=32){
    cols_cluster <- c("#cf4b35","#4ca8bf","#19937d","#394f7d","#e39279","#7e8aaa","#7a5f47",
                               "#8bc3b6","#ae9a84","#41529b","#bb3e31","#709356","#e2da84","#44657d",
                               "#af5f39","#76215f","#71b666","#77181c","#d4d4cc","#edab4c","#a7ab7b",
                               "#58869c","#c98d62","#a8706b","#878a78","#6f5461","#e26b21","#b11f23",
                               "#7b3d87","#568d90","#e68b9c","#040000")
  }else{
    cols_cluster <- turbo(ntype)
  }
  plota= ggplot() +
    geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,
                                              ymax = Max),
              fill = "#525252", alpha = 0.1) + ###添加灰色背景色
    geom_jitter(data = r.deg, aes(x =unm, y = avg_log2FC, colour = thr_signi),
                size = 1,position = position_jitter(seed = 1)) +
    scale_color_manual(values = cols_thr_signi) +
    scale_x_continuous(limits = c(-0.5, max(r.deg$unm) + 0.5),
                       breaks = seq(0, max(r.deg$unm), 1),
                       label = seq(0, max(r.deg$unm),1)) + #修改坐标轴显示刻度
    
    # 根据top_label标注基因名
    geom_text_repel(data = top_label, aes(x =unm, y = avg_log2FC, label = gene),
                    position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                    box.padding = unit(0, "lines")) +
    
    geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.4,
                                               ymax = 0.4, fill = unm), color = "black", alpha = 1, show.legend = F) +
    scale_fill_manual(values = cols_cluster) +
    labs(x = "Cluster", y = "average log2FC") +
    theme_bw()
  
    plotb <- plota + theme(panel.grid.minor = element_blank(), ##去除网格线
                           panel.grid.major = element_blank(),
                           axis.text.y = element_text(colour = 'black', size = 14),
                           axis.text.x = element_text(colour = 'black', size = 14, vjust = 80), #调整x轴坐标,vjust的值按照最终结果稍加调整
                           panel.border = element_blank(), ## 去掉坐标轴
                           axis.ticks.x = element_blank(), ## 去掉的坐标刻度线
                           axis.line.y = element_line(colour = "black")) #添加y轴坐标轴
    ggsave(filename = paste0("fig",plotid,".deg_pointplot.tiff"), plot = plotb, width = 9, height = 6,compression='lzw')
    ggsave(filename = paste0("fig",plotid,".deg_pointplot.pdf"), plot = plotb, width = 9, height = 6)
}



## seurat_plot_ncell_barplot_with_errorbar(metadata,sample='orig.ident',group='group',celltype='celltype_sub',legend_pos=c(0.8,0.95),plotid='08')
seurat_plot_ncell_barplot_with_errorbar <- function(metadata,sample='orig.ident',group='group',celltype='celltype_sub',legend_pos=c(0.8,0.95),plotid='08')
{
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(grid)
  
  head(metadata)
  df_ncell<-as.data.frame(table(metadata[c(sample,celltype)]))
  head(df_ncell)
  df_group<-unique(metadata[,c(sample,group)])
  df_plot<-merge(df_ncell,df_group,by=sample)
  head(df_plot)
  
  formula <- as.formula(paste0(sample,'~', celltype))
  df_ncell_matrix <- dcast(df_ncell, formula)
  colnames(df_ncell_matrix)[1]<-x
  write.table(df_ncell_matrix,paste0("fig",plotid,".",sample,"-",celltype,".Ncell.tsv"),sep="\t",row.names=F,quote=F)
  
  ncelltype<-length(unique(metadata[[celltype]]))
  ngroup<-length(unique(metadata[[group]]))
  mycolors<-get_colors(ngroup)
  ggplot2_barplot_paired_with_errorbar(df_plot,value='Freq',x=group,y=celltype,y_title='Cells',output=paste0("fig",plotid,".barplot.",group,"-",celltype,".Ncell"))
}


seurat_markers_plot_5<- function(pbmc,markers,group.by='celltype',plotid='04',assay='RNA')
{
	DefaultAssay(pbmc) <- assay
	Idents(pbmc) <- group.by
	ntype<-length(levels(pbmc))
	markers<-unique(markers)
	ngene<-length(markers)

	plota<-NULL
	plota<-FeaturePlot(object = pbmc, features = markers, cols = c("green", "red"),ncol=ngene,raster=F)
	ggsave(paste0("fig",plotid,".FeaturePlot.markers.tiff"), plot = plota, width = ngene*6, height = 5,compression='lzw',limitsize=F)
	ggsave(paste0("fig",plotid,".FeaturePlot.markers.pdf"), plot = plota, width = ngene*6, height = 5,limitsize=F)

	plota<-NULL
	plota<-VlnPlot(object = pbmc, features = markers,pt.size = 0,ncol=ngene)
	ggsave(paste0("fig",plotid,".Violin.markers.tiff"), plot = plota, width = ngene*6, height = 5,compression='lzw',limitsize=F)
	ggsave(paste0("fig",plotid,".Violin.markers.pdf"), plot = plota, width = ngene*6, height = 5,limitsize=F)

	plota<-NULL
	plota<-RidgePlot(pbmc, features = markers, ncol = ngene)
	ggsave(paste0("fig",plotid,".RidgePlot.markers.tiff"), plot = plota, width = ngene*6, height = 5,compression='lzw',limitsize=F)
	ggsave(paste0("fig",plotid,".RidgePlot.markers.pdf"), plot = plota, width = ngene*6, height = 5,limitsize=F)
}

seurat_markers_plot_25<- function(pbmc,markers,group.by='celltype',plotid='04',assay='RNA')
{
	DefaultAssay(pbmc) <- assay
	Idents(pbmc) <- group.by
	ntype<-length(levels(pbmc))
	markers<-unique(markers)
	ngene<-length(markers)
	nrow<-ceiling(ngene/5)

	plota<-NULL
	plota<-FeaturePlot(object = pbmc, features = markers, cols = c("green", "red"),ncol=5,raster=F)
	ggsave(paste0("fig",plotid,".FeaturePlot.markers.tiff"), plot = plota, width = 30, height = nrow*5,compression='lzw',limitsize=F)
	ggsave(paste0("fig",plotid,".FeaturePlot.markers.pdf"), plot = plota, width = 30, height = nrow*5,limitsize=F)

	plota<-NULL
	plota<-VlnPlot(object = pbmc, features = markers,pt.size = 0,ncol=5)
	ggsave(paste0("fig",plotid,".Violin.markers.tiff"), plot = plota, width = 30, height = nrow*5,compression='lzw',limitsize=F)
	ggsave(paste0("fig",plotid,".Violin.markers.pdf"), plot = plota, width = 30, height = nrow*5,limitsize=F)

	plota<-NULL
	plota<-RidgePlot(pbmc, features = markers, ncol = 5)
	ggsave(paste0("fig",plotid,".RidgePlot.markers.tiff"), plot = plota, width = 30, height = nrow*5,compression='lzw',limitsize=F)
	ggsave(paste0("fig",plotid,".RidgePlot.markers.pdf"), plot = plota, width = 30, height = nrow*5,limitsize=F)
}
# seurat_markers_plot(pbmc,markers,group.by='celltype',plotid='04',dotplot.only=F,plot.heatmap=F)
seurat_markers_plot<- function(pbmc,markers,assay='RNA',group.by='celltype',ngene_cutoff=20000,plotid='04',dotplot.only=F,plot.heatmap=F,pw=NA,ph=NA){	
	#Idents(pbmc) <- "celltype"
	DefaultAssay(pbmc) <- assay
	Idents(pbmc) <- group.by
	ntype<-length(levels(pbmc))
	all_genes <- rownames(pbmc)[rowSums(pbmc) > 0]
	markers<-unique(markers)
	markers_used<-markers[markers %in% all_genes]
	ngene<-length(markers_used)	

	plota<-NULL
	if(is.na(ph))
	{
		ph<-max(ntype*0.5,5)
	}
	if(is.na(pw))
	{
		pw<-max(ngene*0.4+3,5)
	}
	if(! dotplot.only)
	{
		if(plot.heatmap)
		{
			seurat_markers_plot_heatmap(pbmc,markers_used,group.by=group.by,output=paste0('fig',plotid,'.Heatmap.markers'))
		}
	}
	if(ngene < ngene_cutoff)
	{
		plota<-DotPlot(object = pbmc, features = markers_used)+ theme(axis.text.x =element_text(angle=45,hjust=1))
		ggsave(paste0("fig",plotid,".Bubble.markers.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
		ggsave(paste0("fig",plotid,".Bubble.markers.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
		if(ngene>1)
		{
			plota<-NULL
			plota<-VlnPlot(object = pbmc, features = markers_used,stack=T,flip=T,pt.size = 0)+ theme(legend.position = "none",axis.text.x =element_text(angle=45,hjust=1))
			ggsave(paste0("fig",plotid,".Violin.markers",".tiff"), plot = plota, width = ph, height = pw,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,".Violin.markers",".pdf"), plot = plota, width = ph, height = pw,limitsize=F)
		}
			if(! dotplot.only)
			{
				ncut<-ceiling(ngene/25)
				i<-0
				while(i<ncut)
				{
					markersi<-markers_used[(i*25+1):(i*25+25)]
					markersi<-na.omit(markersi)
					ngenei<-length(markersi)
					if(ngenei<6)
					{
						seurat_markers_plot_5(pbmc,markersi,group.by=group.by,plotid=paste0(plotid,i),assay=assay)
					}else{
						seurat_markers_plot_25(pbmc,markersi,group.by=group.by,plotid=paste0(plotid,i),assay=assay)
					}
					i=i+1
				}
			}
	}
}

seurat_marker_dotplot_paired<- function(pbmc,markers,group.by='celltype',paired_group='group',plotid='04',pw=NA,ph=NA,plot_bygene=T){
		library(tidyverse)
	#Idents(pbmc) <- "celltype"
	DefaultAssay(pbmc) <- "RNA"
	Idents(pbmc) <- group.by
	ntype<-length(levels(pbmc))
	groups<-unique(pbmc@meta.data[,paired_group])
	ngroup<-length(groups)
	ndot<-ntype*ngroup
	
	ngene<-length(markers)
	
	plota<-NULL
	if(is.na(ph))
	{
		ph<-max(ntype*0.5,8)
	}
	if(is.na(pw))
	{
		pw<-max(ngene*0.4+3,6)
	}
	plota<-DotPlot(object = pbmc, features = markers,split.by=paired_group)+ theme(axis.text.x =element_text(angle=45,hjust=1))
	ggsave(paste0("fig",plotid,".dotplot.a.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
	ggsave(paste0("fig",plotid,".dotplot.a.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
	
	plotb<-NULL
	plotb<-FeaturePlot(object = pbmc, features = markers,split.by=paired_group, cols = c("green", "red"),ncol=2,raster=F)
	ggsave(paste0("fig",plotid,".FeaturePlot.markers.split.tiff"), plot = plotb, width = 2*6, height = ngene*5,compression='lzw',limitsize=F)
	ggsave(paste0("fig",plotid,".FeaturePlot.markers.split.pdf"), plot = plotb, width = 2*6, height = ngene*5,limitsize=F)
	
	if(ngene>1)
	{
		plotc<-NULL
		plotc<-VlnPlot(object = pbmc, features = markers,stack=T,flip=T,pt.size = 0)+ theme(legend.position = "none",axis.text.x =element_text(angle=45,hjust=1))
		ggsave(paste0("fig",plotid,".Violin.markers",".tiff"), plot = plotc, width = 4, height = ngene*3,compression='lzw',limitsize=F)
		ggsave(paste0("fig",plotid,".Violin.markers",".pdf"), plot = plotc, width = 4, height = ngene*3,limitsize=F)
	}

	
	df_plot<-plota$data
	
		df_plot$group<-'unknown'
		for(i in 1:ngroup)
		{
		  #i=1
		  group<-groups[i]
		  group_<-paste0('_',group)
		  
		  df_plot[grepl(paste0(group_,'$'),df_plot$id),'group']<-group
		}

		for(i in 1:nrow(df_plot))
		{
		  #i=1
		  group<-df_plot[i,'group']
		  group_<-paste0('_',group)
		  df_plot[i,'celltype']<-sub(paste0(group_,'$'),'',df_plot[i,'id'])
		}
	write.table(df_plot, paste0("fig",plotid,".dotplot.a.tsv"), sep='\t', quote=F,row.names=T)
	if(plot_bygene)
	{
		ph2<-max(ngroup*0.5,4)
		pw2<-max(ntype*0.4+3,6)

		for(i in 1:length(markers))
		{
			genei<-markers[i]
			df_ploti<-df_plot[df_plot$features.plot==genei,]
			plota<-ggplot(data = df_ploti, mapping = aes(x=celltype,y=group)) +
			  geom_point(aes(size=pct.exp,col=avg.exp))+ 
			  #scale_color_viridis(name='Average Expression')+
			  scale_color_gradient(name='Average Expression',low = "grey", high = "blue")+
			  guides(size = guide_legend(title = 'Percent Expressed'))+
			theme_classic()+
			theme(text=element_text(size=15,color='black'),axis.title.x =element_blank(),axis.text.x =element_text(angle=30,hjust=1),
					axis.title.y =element_blank())+
			  theme(panel.border=element_blank())+
			  theme(panel.background = element_blank(),plot.title = element_text(hjust = 0.5,size=30),
					legend.title = element_text(size = 10, face = 2),
					legend.text = element_text( size = 10,face = 'bold'),legend.background=element_blank())
			ggsave(paste0("fig",plotid,".dotplot.b.",genei,".tiff"), plot = plota, width = pw2, height = ph2,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,".dotplot.b.",genei,".pdf"), plot = plota, width = pw2, height = ph2,limitsize=F)
		}
	}

}


seurat_markers_plot_too_much<- function(pbmc=pbmc,markers,group.by='celltype',plotid='04',output_dir='output'){	
#Idents(pbmc) <- "celltype"
Idents(pbmc) <- group.by
ntype<-length(levels(pbmc))
markers<-unique(markers)
ngene<-length(markers)
nrow<-ceiling(ngene/5)
features<-rownames(x=pbmc)
		        if(ngene < 200)
		        {
		        plota<-NULL
				ph<-max(ntype*0.5,5)
				plota<-DotPlot(object = pbmc, features = markers)+ theme(axis.text.x =element_text(angle=45,hjust=1))
				ggsave(paste0("fig",plotid,".Bubble.markers.tiff"), plot = plota, width = ngene*0.4+3, height = ph,compression='lzw',limitsize=F)
				ggsave(paste0("fig",plotid,".Bubble.markers.pdf"), plot = plota, width = ngene*0.4+3, height = ph,limitsize=F)
			}
				
		#gene='CXCL8'
		for(gene in markers)
		{
			if(gene %in% features)
			{
				plota<-NULL                                         
				
				                                         ##  cols = c("green", "red")
				                                         #   scale_colour_gradient2(low = "gray", mid = "gray", high = "red",midpoint = 0.5)
				                                         #   & scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values = c(1.0,0.8,0.6,0.4,0.2,0))
				                                         ##  & scale_colour_gradientn(colours = c("red",'red',"gray"),values = c(1.0,0.2,0.1,0))
				                                         
				plota<-FeaturePlot(object = pbmc, features = gene)
				ggsave(paste0(output_dir,"/fig",plotid,".FeaturePlot.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
				ggsave(paste0(output_dir,"/fig",plotid,".FeaturePlot.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)

				plota<-NULL
				plota<-RidgePlot(pbmc, features = gene)
				ggsave(paste0(output_dir,"/fig",plotid,".RidgePlot.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
				ggsave(paste0(output_dir,"/fig",plotid,".RidgePlot.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)
				
				plota<-NULL
				plota<-VlnPlot(object = pbmc, features = gene,pt.size = 0)
				ggsave(paste0(output_dir,"/fig",plotid,".Violin.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
				ggsave(paste0(output_dir,"/fig",plotid,".Violin.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)
			}
		}
}

#seurat_cell_density_plot(pbmc.x,group.by='group',plotid='11')
seurat_cell_density_plot<- function(pbmc,group.by='group',plotid='04'){	
Idents(pbmc) <- group.by
ntype<-length(levels(pbmc))
		if(ntype<5)
		{
				 library(viridis)
					metadata <- data.frame(pbmc@meta.data)
					metadata$UMAP1 <- Embeddings(pbmc, "umap")[,1]
					metadata$UMAP2 <- Embeddings(pbmc, "umap")[,2]
                    
                    colnames(metadata)[which(colnames(metadata)==group.by)]<-'xxxxxxxxxxxxxxxx'
                    
					#metadata_bg <- metadata[,-(which(colnames(metadata)=="group"))]
					metadata$xxxxxxxxxxxxxxxx <- factor(metadata$xxxxxxxxxxxxxxxx )
					density_plot<-NULL
					density_plot <- ggplot(metadata, aes(x=UMAP1, y=UMAP2)) +
					  stat_density_2d(geom="raster", aes(fill=stat(ndensity)), contour=F) + 
					  #ndensity calculates the normalized density for each sample--otherwise density would be affected by the number of cells for each sample, which is variable
					  #geom_point(data=metadata_bg, shape=16, size=0.1, alpha=0.2, color="white") +
						#scale_fill_gradientn(colours = c("white",'#ffff80','#a2137b',"Navy"),values = c(0,0.2,0.55,1),name="Density")+
						scale_fill_gradientn(colours = c("white",'#ffff80','#CE0000',"#7d132a"),values = c(0,0.25,0.75,1),name="Density")+
						
					  #scale_x_continuous(expand=c(0,0)) +
					  #scale_y_continuous(expand=c(0,0)) +
					  facet_wrap(~xxxxxxxxxxxxxxxx, ncol=8) +
					  theme_classic() +
					  theme(strip.background = element_blank(),plot.title = element_text(size=25),
							strip.text = element_text(size=25, color="black"),
							axis.text=element_blank(),
							axis.title=element_blank(),
							axis.ticks=element_blank(),
							axis.line = element_blank(),
							plot.background = element_rect(fill = "transparent", color = NA),
							legend.text=element_text(size=15, color="black"),
							legend.title=element_text(size=15, color="black"))
					ggsave(paste0("fig",plotid,".stat_density_2d.raster.tiff"), plot = density_plot, width = ntype*8, height = 6,compression='lzw',limitsize=F)
					ggsave(paste0("fig",plotid,".stat_density_2d.raster.pdf"), plot = density_plot, width = ntype*8, height = 6,limitsize=F)
if(0)
{
					density_plot<-NULL
					density_plot <- ggplot(metadata, aes(x=UMAP1, y=UMAP2)) +
					  stat_density_2d(geom="raster", aes(fill=stat(ndensity)), contour=F) + 
					  scale_fill_viridis(option="plasma",name="Density")+
					  scale_x_continuous(expand=c(0,0)) +
					  scale_y_continuous(expand=c(0,0)) +
					  facet_wrap(~xxxxxxxxxxxxxxxx, ncol=8) +
					  theme_classic() +
					  theme(strip.background = element_blank(),plot.title = element_text(size=20),
							strip.text = element_text(size=12, color="black"),
							axis.text=element_blank(),
							axis.title=element_blank(),
							axis.ticks=element_blank(),
							axis.line = element_blank(),
							plot.background = element_rect(fill = "transparent", color = NA),
							legend.text=element_text(size=15, color="black"),
							legend.title=element_text(size=20, color="black"))
					ggsave(paste0("fig",plotid,".stat_density_2d.raster.viridis.tiff"), plot = density_plot, width = ntype*8, height = 6,compression='lzw',limitsize=F)
					ggsave(paste0("fig",plotid,".stat_density_2d.raster.viridis.pdf"), plot = density_plot, width = ntype*8, height = 6,limitsize=F)

					density_plot<-NULL
					density_plot <- ggplot(metadata, aes(x=UMAP1, y=UMAP2)) +
					  stat_density_2d(geom="raster", aes(fill=stat(ndensity)), contour=F) + 
					  scale_fill_viridis(option="plasma",name="Density",direction=-1)+
					  scale_x_continuous(expand=c(0,0)) +
					  scale_y_continuous(expand=c(0,0)) +
					  facet_wrap(~xxxxxxxxxxxxxxxx, ncol=8) +
					  theme_classic() +
					  theme(strip.background = element_blank(),plot.title = element_text(size=20),
							strip.text = element_text(size=12, color="black"),
							axis.text=element_blank(),
							axis.title=element_blank(),
							axis.ticks=element_blank(),
							axis.line = element_blank(),
							plot.background = element_rect(fill = "transparent", color = NA),
							legend.text=element_text(size=15, color="black"),
							legend.title=element_text(size=20, color="black"))
					ggsave(paste0("fig",plotid,".stat_density_2d.raster.viridis2.tiff"), plot = density_plot, width = ntype*8, height = 6,compression='lzw',limitsize=F)
					ggsave(paste0("fig",plotid,".stat_density_2d.raster.viridis2.pdf"), plot = density_plot, width = ntype*8, height = 6,limitsize=F)
}
		}
}


#seurat_markers_plot_heatmap(pbmc=pbmc.x,transcription_factor_gene,group.by='sub_celltype',plotid='08')
seurat_markers_plot_heatmap<- function(pbmc=pbmc,markers,group.by='celltype',output='Heatmap',pw=NULL,ph=NULL){	
library(circlize)
library(pheatmap)
Idents(pbmc) <- group.by

ntype<-length(levels(pbmc))
markers<-unique(markers)
features<-rownames(x=pbmc)
markers_used<-markers[markers %in% features]
ngene<-length(markers_used)
if(is.null(pw))
{
	pw<-max(5,ntype)
}
if(is.null(ph))
{
	ph<-max(6,ngene*0.2)
}
				if(1)
				{
					pdf(file=paste0(output,".col1.pdf"),width=pw,height=ph)
					print(DoHeatmap(object = pbmc, features = markers,label=F))
					dev.off()
					png(file=paste0(output,".col1.png"),width=pw*300,height=ph*300,res=300)
					#print(DoHeatmap(object = pbmc, features = top5.head) + NoLegend())
					print(DoHeatmap(object = pbmc, features = markers,label=F))
					dev.off()
				}
				if(1)
				{
					pdf(file=paste0(output,".col2.pdf"),width=pw,height=ph)
					print(DoHeatmap(object = pbmc, features = markers,label=F)+scale_fill_gradientn(colors = c("#406aa8", "white", "#d91216")))
					#+theme(text=element_text(size = 20))
					dev.off()
					png(file=paste0(output,".col2.png"),width=pw*300,height=ph*300,res=300)
					print(DoHeatmap(object = pbmc, features = markers,label=F)+scale_fill_gradientn(colors = c("#406aa8", "white", "#d91216")))
					dev.off()
				}
}


#seurat_markers_plot_heatmap_groupmean(pbmc=pbmc.x,transcription_factor_gene,group.by='sub_celltype',plotid='08')
seurat_markers_plot_heatmap_groupmean<- function(pbmc=pbmc,markers,group.by='celltype',plotid='04',output_dir='output'){	

#pbmc=pbmc.x
#markers=list_genesets[[1]]
#group.by='celltype_sub'
#plotid='08'

library(circlize)
library(ComplexHeatmap)
#Idents(pbmc) <- "celltype"
Idents(pbmc) <- group.by

ntype<-length(levels(pbmc))
markers<-unique(markers)
features<-rownames(x=pbmc)
markers_used<-markers[markers %in% features]
ngene<-length(markers_used)

df_data <- AverageExpression(pbmc, assays = "RNA", slot = "data")[[1]]
write.table(df_data, paste0("fig",plotid,".expr_mean.by-",group.by,".tsv"), sep='\t', quote=F,row.names=T)
df_data<-df_data[rownames(df_data) %in% markers_used,]
df_data_output<-cbind(rownames(df_data),df_data)
colnames(df_data_output)[1]<-'gene'
write.table(df_data_output, paste0("fig",plotid,".Heatmap.",group.by,".tsv"), sep='\t', quote=F,row.names=F)

		if(ngene<100)
		{
				cols = colorRamp2(c(min(df_data), (min(df_data)+max(df_data))/2, max(df_data)), c("#377EB8", "white", "#E41A1C"))
				plot.pheatmap<-NULL
				plot.pheatmap<-ComplexHeatmap::pheatmap(df_data,fontsize_row=10,col = cols)
				#dev.off()
				pdf(file=paste0("fig",plotid,".Heatmap.",group.by,".nonscale.pdf"), height=ngene*25/100, width=ntype*0.8+2)
				print(plot.pheatmap)
				dev.off()
				tiff(file=paste0("fig",plotid,".Heatmap.",group.by,".nonscale.tiff"), height=ngene*25/100*300, width=(ntype*0.8+2)*300,compression='lzw')
				print(plot.pheatmap)
				dev.off()
				
				df_data<-as.data.frame(df_data)
				sds<-NULL
				for(i in 1:nrow(df_data))
				{
					sds[i] <- sd(df_data[i,])
				}
				df_data<-df_data[sds>0,]
				df_data<-as.matrix(df_data)
				plot.pheatmap<-NULL
				plot.pheatmap <-pheatmap(df_data,scale="row",fontsize_row=10)
				#dev.off()
				pdf(file=paste0("fig",plotid,".Heatmap.",group.by,".scale_row.pdf"), height=ngene*25/100, width=ntype*0.8+2)
				print(plot.pheatmap)
				dev.off()
				tiff(file=paste0("fig",plotid,".Heatmap.",group.by,".scale_row.tiff"), height=ngene*25/100*300, width=(ntype*0.8+2)*300,compression='lzw',res=300)
				print(plot.pheatmap)
				dev.off()
		}else{
			print("too many genes")
			min(df_data)
			cols = colorRamp2(c(min(df_data), (min(df_data)+max(df_data))/2, max(df_data)), c("#377EB8", "white", "#E41A1C"))
			plot.pheatmap<-NULL
			plot.pheatmap<- pheatmap(df_data,fontsize_row=10,col = cols,show_rownames = F)
			#dev.off()
			pdf(file=paste0("fig",plotid,".Heatmap.",group.by,".nonscale.pdf"), height=ngene*8/100, width=ntype*0.8+2)
			print(plot.pheatmap)
			dev.off()
			#tiff(file=paste0("fig",plotid,".Heatmap.",group.by,".nonscale.tiff"), height=ngene*8/100*300, width=(ntype*0.8+2)*300,compression='lzw',res=300)
			#print(plot.pheatmap)
			#dev.off()
			
			df_data<-as.data.frame(df_data)
			sds<-NULL
			for(i in 1:nrow(df_data))
			{
				sds[i] <- sd(df_data[i,])
			}
			df_data<-df_data[sds>0,]
			df_data<-as.matrix(df_data)
			plot.pheatmap<-NULL
			plot.pheatmap <-pheatmap(df_data,scale="row",fontsize_row=10,show_rownames = F)
			#dev.off()
			pdf(file=paste0("fig",plotid,".Heatmap.",group.by,".scale_row.pdf"), height=ngene*8/100, width=ntype*0.8+2)
			print(plot.pheatmap)
			dev.off()
			tiff(file=paste0("fig",plotid,".Heatmap.",group.by,".scale_row.tiff"), height=ngene*8/100*300, width=(ntype*0.8+2)*300,compression='lzw',res=300)
			print(plot.pheatmap)
			dev.off()
		}
}


seurat_markers_plot2<- function(pbmc,markers,group.by='celltype',plotid='04'){	
#Idents(pbmc) <- "celltype"
Idents(pbmc) <- group.by
ntype<-length(levels(pbmc))
markers<-unique(markers)
ngene<-length(markers)
nrow<-ceiling(ngene/5)
		if(ngene<5)
		{
				                                       plota<-NULL                                         
				
				                                         ##  cols = c("green", "red")
				                                         #   scale_colour_gradient2(low = "gray", mid = "gray", high = "red",midpoint = 0.5)
				                                         #   & scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values = c(1.0,0.8,0.6,0.4,0.2,0))
				                                         ##  & scale_colour_gradientn(colours = c("red",'red',"gray"),values = c(1.0,0.2,0.1,0))
				                                         
				plota<-FeaturePlot(object = pbmc, features = markers,ncol=ngene) & scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values = c(1.0,0.8,0.6,0.4,0.2,0))
				ggsave(paste0("fig",plotid,".FeaturePlot.markers.tiff"), plot = plota, width = ngene*6, height = 5,compression='lzw',limitsize=F)
				ggsave(paste0("fig",plotid,".FeaturePlot.markers.pdf"), plot = plota, width = ngene*6, height = 5,limitsize=F)

				plota<-NULL
				plota<-VlnPlot(object = pbmc, features = markers,pt.size = 0,ncol=ngene)
				ggsave(paste0("fig",plotid,".Violin.markers.tiff"), plot = plota, width = ngene*6, height = 5,compression='lzw',limitsize=F)
				ggsave(paste0("fig",plotid,".Violin.markers.pdf"), plot = plota, width = ngene*6, height = 5,limitsize=F)
				
				plota<-NULL
				ph<-min(ntype*0.5,5)
				plota<-DotPlot(object = pbmc, features = markers)+ theme(axis.text.x =element_text(angle=45,hjust=1))
				ggsave(paste0("fig",plotid,".Bubble.markers.tiff"), plot = plota, width = ngene*0.5+4, height = ph,compression='lzw')
				ggsave(paste0("fig",plotid,".Bubble.markers.pdf"), plot = plota, width = ngene*0.5+4, height = ph)
		}else{
				plota<-NULL
				plota<-FeaturePlot(object = pbmc, features = markers, ncol=5) & scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values = c(1.0,0.8,0.6,0.4,0.2,0))
				ggsave(paste0("fig",plotid,".FeaturePlot.markers.tiff"), plot = plota, width = 30, height = nrow*5,compression='lzw',limitsize=F)
				ggsave(paste0("fig",plotid,".FeaturePlot.markers.pdf"), plot = plota, width = 30, height = nrow*5,limitsize=F)

				plota<-NULL
				plota<-VlnPlot(object = pbmc, features = markers,pt.size = 0,ncol=5)
				ggsave(paste0("fig",plotid,".Violin.markers.tiff"), plot = plota, width = 30, height = nrow*5,compression='lzw',limitsize=F)
				ggsave(paste0("fig",plotid,".Violin.markers.pdf"), plot = plota, width = 30, height = nrow*5,limitsize=F)
				
				plota<-NULL
				ph<-max(ntype*0.5,5)
				plota<-DotPlot(object = pbmc, features = markers)+ theme(axis.text.x =element_text(angle=45,hjust=1))
				ggsave(paste0("fig",plotid,".Bubble.markers.tiff"), plot = plota, width = ngene*0.4+3, height = ph,compression='lzw',limitsize=F)
				ggsave(paste0("fig",plotid,".Bubble.markers.pdf"), plot = plota, width = ngene*0.4+3, height = ph,limitsize=F)
		}
}

#seurat_markers_plot_seperate(pbmc,markers,group.by='celltype',plotid='04')
seurat_markers_plot_seperate<- function(pbmc,markers,group.by='celltype',plotid='04'){	
#Idents(pbmc) <- "celltype"
Idents(pbmc) <- group.by
ntype<-length(levels(pbmc))
markers<-unique(markers)
ngene<-length(markers)
nrow<-ceiling(ngene/5)
		
		for(gene in markers)
		{
				plota<-NULL                                         
				
				##  cols = c("green", "red")
				#   scale_colour_gradient2(low = "gray", mid = "gray", high = "red",midpoint = 0.5)
				#   & scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values = c(1.0,0.8,0.6,0.4,0.2,0))
				##  & scale_colour_gradientn(colours = c("red",'red',"gray"),values = c(1.0,0.2,0.1,0))
				
				plota<-FeaturePlot(object = pbmc, features = gene)
				ggsave(paste0("fig",plotid,".FeaturePlot.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
				ggsave(paste0("fig",plotid,".FeaturePlot.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)

				plota<-NULL
				plota<-RidgePlot(pbmc, features = gene)+theme(legend.position = 'none')
				ggsave(paste0("fig",plotid,".RidgePlot.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
				ggsave(paste0("fig",plotid,".RidgePlot.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)
				
				plota<-NULL
				plota<-VlnPlot(object = pbmc, features = gene,pt.size = 0)+theme(legend.position = 'none')
				ggsave(paste0("fig",plotid,".Violin.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
				ggsave(paste0("fig",plotid,".Violin.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)
		}
}

#seurat_markers_4p(pbmc,markers,group.by='celltype',plotid='04')
seurat_markers_4p<- function(pbmc,markers,group.by='celltype',plotid='04'){	
	#Idents(pbmc) <- "celltype"
	library(ggplot2)
	library(magrittr)
	DefaultAssay(pbmc) <- "RNA"
	Idents(pbmc) <- group.by
	ntype<-length(levels(pbmc))
	markers<-unique(markers)
	markers_used<-markers[markers %in% rownames(pbmc)]
	ngene<-length(markers_used)
		i<-1
		while(i<=ngene)
		{
				genei<-markers_used[i]
			  gene_FeaturePlot<-FeaturePlot(object = pbmc,features = genei,order  = TRUE,min.cutoff = 'q10', label = TRUE,repel = TRUE)
			  gene_VlnPlot<-VlnPlot(object = pbmc, features = genei,pt.size = 0)+theme(legend.position = 'none')
			  gene_RidgePlot<-RidgePlot(pbmc, features = genei)+theme(legend.position = 'none')
			  gene_DotPlot<-DotPlot(pbmc, features = genei) + RotatedAxis()+theme(plot.background = element_rect(fill='white'))
			  gene_plot<-cowplot::plot_grid(gene_FeaturePlot, gene_VlnPlot, gene_RidgePlot, gene_DotPlot, nrow = 2, labels = LETTERS[1:4])
			  ggplot2::ggsave(plot=gene_plot,paste0('fig',plotid,'.4pictures.',genei,'.pdf'),device = 'pdf',width = 32,height = 32*0.8,units = 'cm')
			  ggplot2::ggsave(plot=gene_plot,paste0('fig',plotid,'.4pictures.',genei,'.tiff'),device = 'tiff',width = 32,height = 32*0.8,units = 'cm',compress='lzw')
			i=i+1
		}
}

#  seurat_cell_annotation_first(pbmc,plotid=20,species='human')
seurat_cell_annotation_first<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F,assay='RNA'){
####  此函数中的marker的来自潘重2023年6月之前的单细胞数据分析的总结,后续再不停更新
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  Epithelial=c('WFDC2','PAX8','EPCAM','GPC3','MMP7','KRT19','PROM1','ALDH1A1','CD24','KRT14','KRT5','JUP'),
  Fibroblasts=c('PDPN','FN1','VIM','DCN','GSN','COL1A1','COL1A2','COL4A1','COL4A2','PDGFRA'),
  SMC =  c('DES','ACTA2','TAGLN','MYH11','MYLK','ACTG2'),
  Pericytes=c('BCAM','PDGFRB','KCNJ8','RGS5','ABCC9','MCAM','NOTCH3'),
  Endothelial =  c('MME','PECAM1','THBD','CDH5','ENG','VWF','CCL21'),
  #mesothelial =  c('LY6G','MSIN'),                                                 #### 间皮细胞
  #mesenchymal =  c('PDGFRA','TCF21','CDH11'),                                      ###  间充质细胞
  Adipocytes =  c('GPAM','FASN','LEP'),
  Neuronal=c('NRXN1'),
  Immune = c("PTPRC"),
  Lymphoid = c("CD69",'MKI67'),
  Tcells = c('CD3D','CD3E','CD3G','CD247','CD4','CD8A','CD8B'),
  NK=c('NCAM1','CD19','GNLY','KLRD1','NKG7','FCGR3A','FGFBP2','CX3CR1','KLRF1'),           ####20
  ILCs=c('AREG','TNFRSF18','TNFRSF25'),
  Bcells=c( 'CD79A','MS4A1'),
  Plasma=c( 'MZB1','SDC1','IGHG1'),
  Myeloid =  c("CD14",'LYZ','CD300E'),  
  Macrophages=c('CD68','CD86','CD163','FCER1G','MARCO'),
  Monocytes=c('FCN1','APOBEC3A','THBS1'),
  DC = c('HLA-DQB2','HLA-DPB1','BIRC3'),                                                 
  Mast =  c('MS4A2','TPSB2',"TPSAB1",'CST3','KIT'),                                       ##24
  #Neutrophils =  c("FCGR3B","CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2",'CSF3R'),           ###15
  #####     CD11b-->ITGAM   CD66b-->CEACAM8
  Neutrophils =  c("ITGAM","CEACAM8",'LY6G','Retnlg'),           ###15
  #https://mp.weixin.qq.com/s/lUZeQEPUSaoDumMnpspMQg
  #https://pubmed.ncbi.nlm.nih.gov/32719519/
  Basophils =  c('CD34','CD200R3','FCER1A','CD9',"CPA3",'ITGB2','IL3RA'),
  Eosinophils =  c('SIGLEC5', 'IL5RA', 'CCR3', 'EPX')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
	#markers_first[[21]][1]<-'Fcgr4'
	#markers_first[[21]][1]<-'Siglecf'
}
p <- DotPlot(pbmc, features = markers_first,assay=assay ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.first.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,ngene_cutoff=20000,dotplot.only=F,plot.heatmap=T)
	}

}






gene_mouse2human_df<-function(aaa)
{
library(clusterProfiler)
library(org.Hs.eg.db)
library(homologene)
library(stringr)
#aaa<-URGs_mouse
df_output<-data.frame('mouse'=aaa)
df_output$human<-toupper(aaa)

AAA<-toupper(aaa)
bbb<-bitr(AAA, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Hs.eg.db)
AAA_used<-bbb[,1]
head(df_output)
AAA_lost<-AAA[! AAA %in% AAA_used]
df_output[df_output$human %in% AAA_lost,'human']<-NA

AAA_lost<-str_to_title(AAA_lost)
ccc<-homologene(AAA_lost, inTax = 10090, outTax = 9606)
for(i in 1:nrow(ccc))
{
df_output[df_output$mouse==ccc[i,1],'human']<-ccc[i,2]
}
return(df_output)
}

gene_human2mouse_df<-function(aaa)
{
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(homologene)
    library(stringr)
    #aaa<-URGs_mouse
    df_output<-data.frame('human'=toupper(aaa))
    df_output$mouse<-str_to_title(aaa)
    
    AAA<-str_to_title(aaa)
    bbb<-bitr(AAA, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Mm.eg.db)
    AAA_used<-bbb[,1]
    head(df_output)
    AAA_lost<-AAA[! AAA %in% AAA_used]
    df_output[df_output$mouse %in% AAA_lost,'mouse']<-NA
    
    AAA_lost<-toupper(AAA_lost)
    ccc<-homologene(AAA_lost, inTax = 9606, outTax = 10090)
    for(i in 1:nrow(ccc))
    {
      df_output[df_output$human==ccc[i,1],'mouse']<-ccc[i,2]
    }
    return(df_output)
  }
  
gene_mouse2human_vector<-function(aaa)
{
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(homologene)
	library(stringr)
	df_output<-data.frame('mouse'=aaa)
	df_output$human<-toupper(aaa) 
	aaa_used<-unique(df_output$human)
	bbb<-bitr(aaa_used, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Hs.eg.db)
	AAA_used<-bbb[,1]
	head(df_output)
	AAA_lost<-aaa_used[! aaa_used %in% AAA_used]

	df_output[df_output$human %in% AAA_lost,'human']<-NA

	AAA_lost<-toupper(AAA_lost)
	ccc<-homologene(AAA_lost, inTax = 10090, outTax = 9606)
	if(nrow(ccc)>0)
	{
	  for(i in 1:nrow(ccc))
	  {
		df_output[df_output$mouse==ccc[i,1],'human']<-ccc[i,2]
	  }
	}
	output<-unique(df_output[,2])
	output<-output[!is.na(output)]
	return(output)
}

gene_human2mouse_vector<-function(aaa)
{
library(clusterProfiler)
library(org.Mm.eg.db)
library(homologene)
library(stringr)

DBkeyset='org.Mm.eg.db'
  
df_output<-data.frame('human'=aaa)
df_output$mouse<-str_to_title(aaa) 
aaa_used<-unique(df_output$mouse)

bbb = try(bitr(aaa_used, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=DBkeyset),silent = T)
if(class(bbb) =='try-error')
{
	return('')
	}else{
		AAA_used<-bbb[,1]
		head(df_output)
		AAA_lost<-aaa_used[! aaa_used %in% AAA_used]

		df_output[df_output$mouse %in% AAA_lost,'mouse']<-NA

		AAA_lost<-toupper(AAA_lost)
		ccc<-homologene(AAA_lost, inTax = 9606, outTax = 10090)
		if(nrow(ccc)>0)
		{
			for(i in 1:nrow(ccc))
			{
			   df_output[df_output$human==ccc[i,1],'mouse']<-ccc[i,2]
			}
		}
		output<-unique(df_output[,2])
		output<-output[!is.na(output)]
		return(output)
		}
}


gene_human2rat_vector<-function(aaa)
{
library(clusterProfiler)
library(org.Mm.eg.db)
library(homologene)
library(stringr)

DBkeyset='org.Rn.eg.db'
  
df_output<-data.frame('human'=aaa)
df_output$mouse<-str_to_title(aaa) 
aaa_used<-unique(df_output$mouse)

bbb = try(bitr(aaa_used, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=DBkeyset),silent = T)
if(class(bbb) =='try-error')
{
	return('')
	}else{
		AAA_used<-bbb[,1]
		head(df_output)
		AAA_lost<-aaa_used[! aaa_used %in% AAA_used]

		df_output[df_output$mouse %in% AAA_lost,'rat']<-NA

		AAA_lost<-toupper(AAA_lost)
		ccc<-homologene(AAA_lost, inTax = 9606, outTax = 10116)
		if(nrow(ccc)>0)
		{
			for(i in 1:nrow(ccc))
			{
			   df_output[df_output$human==ccc[i,1],'rat']<-ccc[i,2]
			}
		}
		output<-unique(df_output[,2])
		output<-output[!is.na(output)]
		return(output)
		}
}



seurat_plot_umap<-function(pbmc,group.by='celltype',plotid='04',pt.size=0.5,pw=8,ph=6)
{
	levels(pbmc)
	ncell<-dim(pbmc)[2]
	if(ncell>10000)
	{
		pt.size=0.2
	}
					plota<-NULL
					plota<-DimPlot(pbmc, reduction = "umap", group.by = group.by,pt.size = pt.size,label=T,raster=F)
					ggsave(paste0("fig",plotid,"a.umap.",group.by,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
					ggsave(paste0("fig",plotid,"a.umap.",group.by,".pdf"), plot = plota, width = pw, height = ph)
					plota<-NULL
					plota<-DimPlot(pbmc, reduction = "umap", group.by = group.by,pt.size = pt.size,label=F,raster=F)
					ggsave(paste0("fig",plotid,"b.umap.",group.by,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
					ggsave(paste0("fig",plotid,"b.umap.",group.by,".pdf"), plot = plota, width = pw, height = ph)
					
	library(viridis)
	ntype<-length(unique(pbmc@meta.data[[group.by]]))
	mycolors<-get_colors(ntype)
	#cols_cluster <- rev(turbo(ntype+2))
	cols_cluster <- mycolors
					plota<-NULL
					plota<-DimPlot(pbmc, reduction = "umap", group.by = group.by,pt.size = pt.size,label=F,cols=cols_cluster,raster=F)
					ggsave(paste0("fig",plotid,"c.umap.",group.by,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
					ggsave(paste0("fig",plotid,"c.umap.",group.by,".pdf"), plot = plota, width = pw, height = ph)
					plota<-NULL
					plota<-DimPlot(pbmc, reduction = "umap", group.by = group.by,pt.size = pt.size,label=T,cols=cols_cluster,raster=F)
					ggsave(paste0("fig",plotid,"d.umap.",group.by,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
					ggsave(paste0("fig",plotid,"d.umap.",group.by,".pdf"), plot = plota, width = pw, height = ph)
}


seurat_cell_annotation_pbmc<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F,assay='RNA'){
#####  https://www.nature.com/articles/s41467-023-38356-1
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  Immune = c("PTPRC"),
  Lymphoid = c("CD69",'MKI67'),
  Tcells = c('CD3D','CD3E','CD3G','CD247','CD4','CD8A','CD8B','CD8B1'),
  NK=c('NCAM1','CD19','GNLY','KLRD1','NKG7','FGFBP2','CX3CR1','KLRF1'),
  ILCs=c('AREG','TNFRSF18','TNFRSF25'),
  Bcells=c( 'CD79A','MS4A1'),
  Plasma=c( 'MZB1','SDC1','IGHG1'),  
  Megakaryocytes=c('PPBP'), 
  Platelet =  c("GP9",'PF4'), 
  myeloid =  c("CD14",'VCAN','LYZ'),  
  Monocytes=c('FCGR3A','MS4A7','FCN1','APOBEC3A','THBS1','S100A8'),
  Macrophages=c('CD68','CD86','CD163','FCER1G','MARCO'),             ###39
  DC = c('HLA-DQA1','HLA-DQB1','CLEC10A','HLA-DQB2','HLA-DPB1','BIRC3','FCER1A','CST3'), 
  Mast =  c('MS4A2','TPSB2',"TPSAB1",'KIT'),
  Neutrophils =  c('ITGAM','CEACAM8','LY6G',"FCGR3B","CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2",'CSF3R','HCAR3'),
  Basophils =  c('CD34','CD200R3','CD9',"CPA3",'ITGB2','PECAM1','IL3RA'),   #c('CD34','CD200R3','FCER1A','CD9',"CPA3",'ITGB2','PECAM1','IL3RA'),
  Eosinophils =  c('SIGLEC5', 'IL5RA', 'CCR3', 'EPX'),
  Erythrocytes=c('HBB', 'HBA2'),
  HSC=c('CYTL1', 'GATA2','CRHBP','HLF','DUSP1','ADGRG6','PCDH9')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
	markers_first[[14]][4]<-'Fcgr4'
	markers_first[[16]][1]<-'Siglecf'
}
p <- DotPlot(pbmc, features = markers_first,assay=assay ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.pbmc.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}


seurat_cell_annotation_granulocyte<-function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
	Idents(pbmc) <- group.by
	nc<-length(levels(pbmc))
	ph<-max(nc/2,6)
	th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust=0.5))
	markers_first = list(
	Immune = c("PTPRC"),
	myeloid =  c('LYZ','LYZ1','LYZ2'),
	granulocyte=c('CCR3','CD33','IL5RA'),
    ### CD11C-->ITGAX
	Mast =  c('CMA1','CPE',"CXCR4",'IL10RA','PTGER2','MS4A2','TPSB2',"TPSAB1",'CST3','KIT'),
	Neutrophils =  c("CD14","ELANE","FCGR3A","MME","MPO",'ITGAM','CEACAM8','LY6G','Retnlg','S100A8','S100A9'),           ###16
	Basophils =  c('CD9','CD22','CD36','CD40LG',"FCGR2B"),
	####      CD31-->PECAM1   CD203c-->ENPP3
	Bas_m_Mast =  c("ITGB2",'PECAM1','IL3RA'),
	Bas_m_Neu =  c("ENPP3"),
	Eosinophils =  c('CCR1', 'FCAR', 'IFNAR1', 'ITGA4','SIGLEC8')
	)
if(species =='mouse')
{
	library(stringr)
	for(i in 1:length(markers_first))
	{
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
	#markers_first[[3]][5]<-'Fcnb'
}
  plota<-NULL
  plota <- DotPlot(pbmc,features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=plota,filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 20,height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,ngene_cutoff=20000,dotplot.only=F,plot.heatmap=T)
	}
}

seurat_cell_annotation_Neutrophils<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
#####  
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  ##  CD45-->PTPRC   CD16-->FCGR3A;FCGR3B  CD15-->FUT4  CD62L-->SELL
  Immune = c("PTPRC"),
  myeloid =  c("CD14",'VCAN','LYZ'),           
  ## CD11b-->ITGAM          CD66b-->CEACAM8   CD11c-->ITGAX   CD18-->ITGB2
  Neu =  c('FCGR3A',"FCGR3B","CXCR2",'ITGAM','SELL'),
  Neu_other=c('LY6G',"SLC25A37","G0S2","ADGRG3","PROK2","STEAP4","CMTM2",'CSF3R'),
  Neu_activate =  c('ITGB2','CD177','PRTN3'),
  Neu_old=  c('CXCR4','CD47','ITGAX','CD24','ICAM1','TLR4'),                ### 衰老中性粒
  Nbh=c('FUT4'),                                        ### B细胞协助中性粒
  Neu_back=c("CXCR1",'ELANE'),                          ### 回流中性粒
  ####  CD10-->MME   CD31-->PECAM1   CD49d-->ITGA4  CD182-->CXCR2  CD43-->SPN
  LDN=c("ARG1",'CD33','CEACAM8','MME','PECAM1'),                                 ### 低密度中性粒
  #PMN-MDSC=c("ARG1",'CD33','FCGR3A','ITGAM')                                 ### 多形核MDSC;粒细胞MDSC，抑制T细胞反应
  PAN=c('FLT1','ITGA4'),                  # c("CXCR4",'FLT1','ITGA4')              ### 促血管生成中性粒细胞
  PMN=c("CD63",'SPN','IL17RA'),              #  c("CD63",'CD49','CXCR2','IL17RA')                   ### 多形核MDSC;粒细胞MDSC，抑制T细胞反应
  TAN=c("CCR7",'CD86')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
	markers_first[[3]][1]<-'Fcgr3'
	markers_first[[3]][2]<-'Fcgr4'
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,ngene_cutoff=20000,dotplot.only=F,plot.heatmap=T)
	}

}

seurat_cell_annotation_Neutrophils_PMID38447573<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
#####  https://doi.org/10.1016/j.cell.2024.02.005  IF:64.5
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  ##  H3F3A-->H3-3A
  ARG1 = c("ARG1","CD163","H3F3A"),
  NFKBIZ_HIF1A =  c("NFKBIZ","HIF1A","MAPK6"),          
  MMP9 =  c("MMP9","CD177","FCN1"),
  IFIT1_ISG15 = c("IFIT1","ISG15","IFIT2"),
  CXCR2 =  c("CXCR2","CXCR1","P2RY13"),           
  CXCL8_IL1B =  c("CXCL8","IL1B","FRMD4B"),
  TXNIP = c("TXNIP","MALAT1","IFITM2"),
  VEGFA_SPP1 =  c("VEGFA","SPP1","CCL3"),           
  HLA_DR_CD74 =  c("HLA-DRA","CD74","HLA-DQB1"),
  S100A12 = c("S100A12","S100A6","S100A4")
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,ngene_cutoff=20000,dotplot.only=F,plot.heatmap=T)
	}

}



seurat_cell_annotation_liver<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
####  https://www.sciencedirect.com/science/article/pii/S2589004221012013
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
Hepatocytes=c('ALB','SERPINA3K'),                        ####  肝实质细胞
HPS_Cho=c('EPCAM','TSPAN8'),                             ####  肝祖细胞、胆管细胞
HSC=c('RELN','RGS5'),                                    ####  肝星状细胞
  Fibroblasts=c('PDPN','FN1','VIM','DCN','GSN','COL1A1','COL1A2','MSLN','MYL7'),   ##15
  SMC =  c('ACTA2','TAGLN','MYH11','MYLK','ACTG2','DES'),                           ##6
  Pericyte=c('PDGFRB','KCNJ8','ABCC9'),                                      ##
  Endothelial =  c('MME','PECAM1','THBD','CDH5','ENG','VWF','CCL21','CLEC4G','FABP4'),
  Immune = c("PTPRC"),
  Lymphoid = c("CD69",'MKI67'),
  Tcells = c('CD3D','CD3E','CD3G','CD247','CD4','CD8A','CD8B'),
  NK=c('NCAM1','CD19','GNLY','KLRD1','NKG7','FCGR3A','FGFBP2','CX3CR1','KLRF1','XCL1','NCR1'),
  ILCs=c('AREG','TNFRSF18','TNFRSF25'),
  Bcells=c( 'CD79A','MS4A1'),
  Plasma=c( 'MZB1'),
  myeloid =  c("CD14",'LYZ'),  
  Mac=c('CD68','CD86','CD163','FCER1G','MARCO'),
  KC=c('CLEC4F','VSIG4','LYZ1'),
  mono=c('FCN1','APOBEC3A','THBS1','CCR2'),
  DC = c('HLA-DQB2','HLA-DPB1','BIRC3','FLT3','H2-OA'), 
  Mast =  c('MS4A2','TPSB2',"TPSAB1",'CST3','KIT'),
  Neutrophils =  c("FCGR3B","CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2",'CSF3R','LY6G','RETNLG','CCL4')
  #https://mp.weixin.qq.com/s/lUZeQEPUSaoDumMnpspMQg
  #https://pubmed.ncbi.nlm.nih.gov/32719519/
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
	markers_first[[21]][1]<-'Fcgr4'
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}


seurat_cell_annotation_colon2<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
	### https://www.cell.com/iscience/fulltext/S2589-0042(23)00908-2
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  Colonocyte=c('CA1','SELENBP1','SLC26A2','KRT20','SLC26A3'),
  CT=c('AQP8','GUCA2A','CEACAM7','GUCA2B'),                 ###    crypt top colonocytes       隐窝顶部结肠细胞
  GOB=c('SPINK4','MUC2','ITLN1','TFF3','ZG16'),                                ###    goblet cells                杯状细胞
  EE =  c('CA7','BEST4','PCSK1N','OTOP2','SPIB')                                     ###    enteroendocfine cells       内分泌细胞  
)
if(species =='mouse')
{
    library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}

seurat_cell_annotation_colon<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
	### https://mp.weixin.qq.com/s/gLaV-biRq4LMFgISdLeC0A
	### https://mp.weixin.qq.com/s/0ukdzOzePAFUvYGTsrzfmA
	### https://www.cmghjournal.org/article/S2352-345X(22)00032-7/fulltext
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  Epi=c('EPCAM','KRT19'),
  colonocytes=c('ALPI','CD24','SLC26A2','CA1','BEST4','HES1','SPIB','HES4','CEACAM7','CA4','CEACAM6'),
  CT=c("OTOP2","MEIS1",'FAM162A','HSD17B2'),                 ###    crypt top colonocytes       隐窝顶部结肠细胞
  ABS=c("KRT20","GUCA2A","ALDOB"),                                         ###    absorptive cell             吸收细胞
  SSC=c("MSLN","MUC5AC","AQP5","TACSTD2","FSCN1","TFF2","ANXA1","ANXA10","MUC17",
        "S100P","GSDMB","GSDMD","L18","RELB","MDK","RARA","RXRA","AHR","AGRN","PDX1"),
  ASC=c("CLDN2","CD44","AXIN2","RNF43","TGFBI","EPHB2","TEAD2","CDX2"),
  STM =  c('LGR5','OLFM4',"ASCL2"),                                                    ###    stem cells                  干细胞
  TAC=c('PCNA','MKI67','CCNB1','CENPA','PTTG1','MCM7','CDK2'),                                                        ###    transit amplifing cells     过渡放大细胞
  GOB=c('ATOH1','MUC2','TFF3','RETNLB','SPDEF'),                                ###    goblet cells                杯状细胞
  EE =  c('CHGA','NEUROD1','GCG','HES6'),                                     ###    enteroendocfine cells       内分泌细胞  
  TUF = c("POU2F3","SOX9","DCLK1"),                                                       ###    tuft cells                  肠道簇细胞
  DCS=c("REG4",'KIT')                                                    #####  colon deep crypt secretory cells   结肠深隐窝分泌细胞
)
if(species =='mouse')
{
    library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}

seurat_cell_annotation_prostate_cancer<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
#####  https://pubmed.ncbi.nlm.nih.gov/36750562/

#####  https://www.jianshu.com/p/20385a420733

Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
	Fibroblasts=c('DCN','LUM','PTN','IGF1','APOD','COL1A2','FBLN1','MEG3','CXCL12'),
	Epithelial=c('EPCAM','PCA3','KRT18','KRT8','MA'),
	Epi_Club=c('SCGB1A3','WFDC2','LCN2','MMP7','KRT4','TACSTD2'),                    #### 俱乐部细胞，是圆顶状细胞具有短绒毛，分布在肺部细支气管，分泌糖胺聚糖。
	Epi_Hillock=c('KRT13','S100A16','S100A14','KRT19'),                                 #### 小丘细胞, 是免疫反应性气道上皮细胞。
	Epi_Basal=c('TP63','KRT14','KRT5'),                                                  #### 基底上皮细胞，与基膜直接接触的上皮细胞。
	Epi_Luminal=c('KLK4','KLK2','KLK3','ACPP','AR'),                                     #### 管腔细胞。
	Endothelial=c('RAMP2','TM4SF1','RNASE1','EGFL7','RAMP3','PLVAP','AQP1','ECSCR','FKBP1A','AC011526.1','EMP1','DARC','VWF','EMCN')
)

if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}

seurat_cell_annotation_prostate_cancer_epithelial<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
#####  https://pubmed.ncbi.nlm.nih.gov/36750562/
#####  https://www.jianshu.com/p/20385a420733
#####  https://www.atsjournals.org/doi/pdf/10.1164/ajrccm-conference.2023.207.1_MeetingAbstracts.A4479
#####  https://mp.weixin.qq.com/s/ZUkw189-p4OMqsTqu91I9w

Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
	Epithelial=c('EPCAM','PCA3','KRT18','KRT8','MA'),
	Epi_Club=c('WFDC2','LCN2','MMP7'),                    #### 俱乐部细胞，是圆顶状细胞具有短绒毛，分布在肺部细支气管，分泌糖胺聚糖。
	Epi_Hillock=c('KRT4','S100A16','TACSTD2','KRT13'),                                 #### 小丘细胞, 是免疫反应性气道上皮细胞。
	Epi_Basal=c('TP63','KRT14','KRT5'),                                                  #### 基底上皮细胞，与基膜直接接触的上皮细胞。
	Epi_Luminal=c('KLK4','KLK2','KLK3','ACPP','AR')                                     #### 管腔细胞。
)

if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}
seurat_cell_annotation_breast_epithelial<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){

Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
	Epithelial=c('EPCAM','PCA3','MA'),
	Myo=c("KRT17", "KRT14", "KRT5", "ACTA2", "MYL9", "MYLK", "MYH11"),
	Lum=c("KRT19", "KRT18", "KRT8"),
	Hs=c("PRLR", "CITED1", "PGR", "PROM1", "ESR1"),
	AV=c("MFGE8", "TRF", "CSN3", "WFDC18", "ELF5", "LTF"),
	Lp=c("KIT", "ALDHLA3", "CD14")
)

if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}

seurat_cell_annotation_liver_epithelial<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){

Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
	Epithelial=c('EPCAM','PCA3','KRT18','KRT8','MA'),
	cholangiocytes=c("FYXD2", "TM4SF4", "ANXA4"),
	hepatocytes =c("APOC3", "APOC1", "FABP1"),
	STEM=c('KRT19','PROM1','ALDH1A1','CD24','ANPEP','CD44','ICAM1','CD47','SOX9')
)

if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}

seurat_cell_annotation_lung_epithelial<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){

Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))

markers_first = list(
	Epithelial=c('EPCAM','PCA3','KRT8','MA'),
	AT1=c("AGER"),
	AT2=c("SFTPA1"),
	Club=c("SCGB1A1"), 
	Basal=c("KRT17"),
	Ciliated=c("TPPP3")
)

if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,"a.dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

markers_first = list(
	Epithelial=c('EPCAM','PCA3','KRT8','MA'),
	squamous= c("SPRR3","GDPD3","SPRR1A","SPRR2A"),
lonocyte=c("RARRES2","TMPRSS11E","ASCL3","CFTR","FOXI2","FOXI1"),
cil_acute_phase=c("1SG20",'ISG15',"SAA4","SAA2","SAA1"),
Ciliated=c("EFHC1","CCDC153","CCDC113","FOXJ1"),
FOXN4=c("CDC20B","MYCL","FOXN4","CCNO"),
Goblet=c("PIGR","BP1","MUC5A","VMO1"),
Club=c("SCGB3A1","CYP2A13","CYP2B6","SCGB1A1"),
Basal=c("BCAM","KRT15","KRT5","TP63")
)

if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,"b.dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}
	
	
	
	
	
	
	

}

seurat_cell_annotation_gastric_epithelial<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))

markers_first = list(
	Epithelial=c('EPCAM','PCA3','KRT8','MA','MUC1','KRT18'),
	goblet=c("MUC2", "ITLN1"),
	enterocytes=c("FABP1", "APOA1"),
	GMCs=c("MUC6", "TFF2"),
	PMCs=c("MUC5AC", "TFF1"),
	chief=c("PGA4", "PGA3"),
	enteroendocrine=c("CHGA", "CHGB")
)

if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}

seurat_cell_annotation_kidney_epithelial<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){

Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))

markers_first = list(
	Epithelial=c('EPCAM','PCA3','KRT8','MA'),
	PT_VCAM1=c('HAVCR1'),
	PC=c('AQP2'),
	PCT=c('CUBN','LRP2','SLC34A1','SLC5A12','SLC5A2','ALDOB'),
	CFH=c('CFH'),
	LOH=c('SLC12A1'),
	DCT=c('SLC12A3','SLC12A2'),
	DCT_CT=c('SLC8A1'),	
	CD_PC=c('AQP3'),
	CD_ICA=c('AQP6','KIT','SLC26A7'),
	CD_ICB=c('SLC26A4','ATP6V0D2'),
	PODO=c('NPHS1','NPHS2')
)

if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}




#df_genes<-list2df(markers_first)
#write.table(df_genes,file=paste0("Cell_markers_first20231116.tsv"),quote=F,sep="\t", row.names=F,na = "")
list2df<- function(markers_first){
  nlist<-length(markers_first)
  nmax<-0
  for(i in 1:nlist)
  {
    #i=1
    genes<-markers_first[[i]]
    name<-names(markers_first)[i]
    if(nmax<length(genes))
    {
      nmax<-length(genes)
    }
  }
  df_genes<-data.frame(tmp=rep('',nmax))
  for(i in 1:nlist)
  {
    #i=1
    genes<-markers_first[[i]]
    name<-names(markers_first)[i]
    df_genes[,name]<-rep('',nmax)
    for(j in 1:length(genes))
    {
      df_genes[j,name]<-genes[j]
    }
  }
  df_genes$tmp<-NULL
  for(i in 1:nlist)
  {
    #i=1
    for(j in 1:nmax)
    {
      if(df_genes[j,i]=='')
      {
         df_genes[j,i]<-NA
      }
    }
  }
  return(df_genes)
}


seurat_cell_annotation_uterus<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
    #参考来源：https://www.nature.com/articles/s41588-021-00972-2
	#桂婷--子宫腺肌病项目--20230927
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
	Epi=c('EPCAM','MMP7'),
	Epi_SOX9=c('SOX9'),
	Lumenal=c('LGR5','PTGS1'),
	Glandular=c('SCGB2A2','SLC18A2','PAEP'),
	Ciliated=c('FOXJ1','PIFO'),
  SMC =  c('ACTA2','LEFTY2','ACTG2','RGS5','NTRK2','FHL5','GUCY1A2'),
  PV_MYH11 =  c('MYH11'),
  PV_STEAP4 =  c('STEAP4'),
  Fibroblasts=c('DCN','GSN','COL1A1','COL1A2','PDGFRA'),
  Fibro_C7 =  c('C7','OGN'),
  eS=c('IGF1','PCOLCE','MMP11','SFRP1'),
  dS=c('DKK1','FOXO1'),  
  Endothelial =  c('MME','PECAM1','THBD','CDH5','ENG','VWF','CCL21'),
  Endothelial_artery =  c('CD34','SEMA3G','GJA5'),
  Endothelial_vein =  c('ACKR1','PLVAP'),
  Myeloid =  c("CD14",'CSF1R'),  
  Lymphoid =  c('CD8A','IL7R','CD40LG','CD3G','NCAM1'),  
  Immune = c("PTPRC")
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}



seurat_cell_annotation_heart<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
    #参考来源：https://pubmed.ncbi.nlm.nih.gov/32971526/
	#杨瑶琳--心力衰竭项目--20230919
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 1, hjust=1))
markers_first = list(
  Fibroblasts=c('DCN','GSN','COL1A1','COL1A2','PDGFRA'),
  SMC =  c('DES','TAGLN','MYH11','MYLK','ACTA2'),
  Pericyte=c('PDGFRB','KCNJ8','RGS5','ABCC9'),
  Endothelial =  c('PECAM1','CDH5','VWF'),
  Mesothelial =  c('MSLN','WT1','BNC1'),
  Adipocytes =  c('GPAM','FASN','LEP'),
  Neuronal =  c('PLP1','NRXN1','NRXN3'),
  Cardiomyocyte=c('MYH7','MYL2','FHL2','NPPA','MYL4'),
  Myeloid =  c("CD14",'CD68','C1QA'),  
  Lymphoid =  c('CD8A','IL7R','CD40LG'),  
  Immune = c("PTPRC")
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}


seurat_cell_annotation_retina<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
	#参考来源：https://pubmed.ncbi.nlm.nih.gov/34364890/
	#王军--糖尿病视网膜病变项目--20230919
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 1, hjust=1))
markers_first = list(
  Cone=c('OPN1MW'),
  Rod=c('NRL'),
  Bipolar =  c('VSX2'),
  RGC = c("THY1"),
  Amacrine = c('PAX6'),
  Muller=c( 'RLBP1'),
  Microglia=c( 'CX3CR1'),
  Astrocyte =  c("GFAP"),
  Endothelial =  c("PECAM1"),
  Pericyte =  c('KCNJ8'),
  RPE=c('RPE65')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}



#  seurat_cell_annotation_first(pbmc,plotid=21)
seurat_cell_annotation_neuronal<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
### https://www.cellsignal.com/science-resources/hallmarks-of-ndg/neuronal-markers
### https://www.abcam.com/neuroscience/microglia-markers
### https://www.cellsignal.jp/science-resources/hallmarks-of-ndg/microglial-markers
### https://mp.weixin.qq.com/s/yLEk3ZtqWjk1cuk1NMQUwQ
### 神经元标志物：
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 1, hjust=1))
markers_first = list(
  Neurons=c('RBFOX3','TUBB3','MAP2','ENO2','NEFL','MAPT'),                                                       ####  Neurons  神经元细胞     TUJ1-->TUBB3  NeuN-->RBFOX3  tau-->MAPT
  Astrocyte=c('GFAP','ALDH1L1',"ADGRV1", "GPC5", "RYR3",'AQP4'),                                                           ####  Astrocyte  星形胶质细胞
  Oligodendrocytes =  c('NKX2-2','NKX6-2','PDGFA','SOX10','OLIG1','OLIG2',"MBP","PLP1","CNP",'ST18'),                               ####  Oligodendrocytes  少突胶质细胞
  Microglia =  c("TMEM119",'P2RY12','CD68','AIF1','ITGAM','HCLS1','PYCARD','CX3CR1','SALL1','C3','LRMDA','DOCK8')      ####  Microglia  小胶质细胞(巨噬细胞)
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}
}


seurat_cell_annotation_Neuronal_sub_old<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
  #参考来源：PMID：30385464 https://zhuanlan.zhihu.com/p/46201084
  Idents(pbmc) <- group.by
  nc<-length(levels(pbmc))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,
                                      vjust = 1, hjust=1))
  markers_first = list(
     Oligodendrocytes = c("OLIG1", "OLIG2", "OLIG3","SOX10"),
     Astrocytes = c("GFAP", "ALDH1L1", "EAAT1", "EAAT2",'S100B'),
     Microglia = c("AIF1", "CD40LG", "CD68", "PTPRC", "SLC2A5"),
     Mature_Neurons = c("MAP2"),
     Glu_Neurons = c("VGLUT1", "VGLUT2","FOLH1", "GOT1","HDLBP", "SLC17A8", "SLC1A1", "SLC14A2", "SLC1A6", "SLC17A6", "SLC17A7"),
     GABA_Neurons = c("ABAT", "GAD1", "GAD2", "PPP1R1B", "SLC6A1", "SLC32A1", "SLC6A13"),
     DA_Neurons = c("DBH", "FOXA2", "KCNJ6", "LMX1B", "NR4A2", "PRKN", "SLC6A2", "SLC6A3", "TH", "TOR1A",'SOX6'),
     Serotonergic_Neurons = c("FEV","SLC6A4","TPH2"),
     Cholinergic_Neurons = c("ACHE", "CHAT", "SLC5A7", "SLC18A3")
   )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-str_to_title(markers_first[[i]])
    }
  }
  p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}

seurat_cell_annotation_Neuronal_sub<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
  #参考来源：PMID：30385464 https://zhuanlan.zhihu.com/p/46201084 PMID：38981007
  #https://www.cellsignal.cn/pathways/neuronal-and-glial-cell-markers
  #http://xteam.xbio.top/CellMarker/
  Idents(pbmc) <- group.by
  nc<-length(levels(pbmc))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,
                                      vjust = 1, hjust=1))
  markers_first = list(
     Neu_Stem=c('SOX2','MSI1','ID1','ID3','PROM1'),
     Neuroepithelial=c('HES1','SOX10','NES'),
     Radial_glial=c('PAX6','PCNA','HOPX'),
     Olig = c("OLIG1", "OLIG2",'CNP'),
     Astrocytes = c("GFAP", "ALDH1L1", "SLC1A2", "SLC1A3",'S100B'),
     Microglia = c("AIF1", "CD40LG", "CD68","SLC2A5"),
     Neurons=c('DCX'),
     Glu_Neurons = c('VGLUT1','VGLUT2','NEUROD2','NEUROD6','TBR1',"SLC1A6", "SLC17A6", "SLC17A7"),
     GABA_Neurons = c("GAD1", "GAD2", "DLX1", "DLX5","SLC6A1","SLC6A13"),
     DA_Neurons = c("DBH","FOXA2","LMX1B", "NR4A2", "PRKN", "SLC6A2", "SLC6A3", "TH"),
     Serotonergic_Neurons = c("FEV","SLC6A4","TPH2"),
     Cholinergic_Neurons = c("ACHE", "CHAT", "SLC5A7", "SLC18A3")
   )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-str_to_title(markers_first[[i]])
    }
  }
  p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}


seurat_cell_annotation_brain_organoids<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
#####  https://stemcellres.biomedcentral.com/articles/10.1186/s13287-023-03302-x
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
	Epithelial=c('WFDC2','PAX8','EPCAM'),
	Fibroblasts=c('PDPN','FN1','VIM','DCN','GSN','COL1A1','COL1A2','ACTA2'),
	Endothelial =  c('MME','PECAM1','THBD','CDH5','ENG','VWF','CCL21'),
	Neuronal=c('NRXN1'),
	Immune = c("PTPRC"),
	Precursors = c('PAX6','SOX1','SOX2','NES','GFAP','HOPX','MKI67','FABP7','PTPRZ1','ARL13B','FAM107A','TUB','ASCL1','EOMES'),
	
	Neurons = c('NEFH','NEUROD1','DCX','ELAVL3','ELAVL4','TUBB3','TBR1','RELN','RBFOX3','MAP2','BCL11B','SATB2',
'POU3F2','CUX1','CUX2','FOXP2','SLC17A7','TH','GPHN','GAD2','GAD1','SLC32A1','CALB1','PVALB','CALB2','SST'),
	AS = c('S100B'),
	OLIGO = c('OLIG2','SOX10','MBP','FOXO4'),
	Other = c('TTR','AIF1')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
	markers_first[[14]][1]<-'Fcgr4'
	markers_first[[17]][1]<-'Siglecf'
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}


seurat_cell_annotation_brain_organoids2<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
#####  https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(19)30337-6
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  GABAerg=c('GAD2','DLX1','DLX5'),
  Glutamat=c('NEUROD2','MAPT','SNAP25','GRIA2'),
  Glial =  c('TTYH1','MT2A','SLC1A3',"NR2F2"),
  Int.Prog = c('GADD45G','EOMES','TAC3','NHLH1'),
  Progenitor = c('ID3','ID1','PAX6','CA2','OTX2'),
  MC = c('MKI67')  
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}


seurat_cell_annotation_brain_organoids_PMID37889749<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
####  此函数中的marker的来自潘重2023年6月之前的单细胞数据分析的总结,后续再不停更新
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  genes=c("VIM", "HES5", "SOX2", "FABP7", "HES1", "NES", "MKI67", 
  "ASPM", "UBE2C", "OLIG1", "OLIG2", "STMN2", "RELN", "NEUROD1", 
  "GRIA2", "SLC17A6", "BCL11B", "FEZF2", "PCP4", "FOXP2", "GAD1", "FLT1", "DDIT3")
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
	#markers_first[[21]][1]<-'Fcgr4'
	#markers_first[[21]][1]<-'Siglecf'
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.first.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}

seurat_cell_annotation_brain_organoids_PMID37605285<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
####  此函数中的marker的来自潘重2023年6月之前的单细胞数据分析的总结,后续再不停更新
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
genes=c("TOP2A", "MKI67", "HES1", "LGALS3", "PAX6", "SOX2", "HOPX", "TNC", 
"S100B", "GFAP", "AQP4", "EOMES", "PPP1R17", "NHLH1", "STMN2", "GAP43", "DCX",
 "SLC17A7", "TBR1", "BCL11B", "SATB2", "POU3F2", "SLC32A1", "GAD1", "GAD2")
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
	#markers_first[[21]][1]<-'Fcgr4'
	#markers_first[[21]][1]<-'Siglecf'
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.first.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}



##  seurat_cell_annotation_myeloids(pbmc,plotid=20)
seurat_cell_annotation_myeloids_sxjns<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  此函数中的marker的来源为生信技能树的公众号文章：https://mp.weixin.qq.com/s/lUZeQEPUSaoDumMnpspMQg
	Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 0.5, hjust=0.5))
markers_first = list(
  Mac=c("C1QA","C1QB","C1QC","SELENOP","RNASE1","DAB2","LGMN","PLTP","MAF","SLCO2B1"),
  mono=c("VCAN","FCN1","CD300E","S100A12","EREG","APOBEC3A","STXBP2","ASGR1","CCR2","NRG1"),
  neutrophils =  c("FCGR3B","CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2" ),
  pDC = c("GZMB","SCT","CLIC3","LRRC26","LILRA4","PACSIN1","CLEC4C","MAP1A","PTCRA","C12orf75"),
  DC1 = c("CLEC9A","XCR1","CLNK","CADM1","ENPP1","SNX22","NCALD","DBN1","HLA-DOB","PPY"),
  DC2=c( "CD1C","FCER1A","CD1E","AL138899.1","CD2","GPAT3","CCND2","ENHO","PKIB","CD1B"),
  DC3 =  c("HMSD","ANKRD33B","LAD1","CCR7","LAMP3","CCL19","CCL22","INSM1","TNNT2","TUBB2B")
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
plota<-NULL
plota <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=plota, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 20, height =ph)
}

##  seurat_cell_annotation_myeloids(pbmc,plotid=20)
seurat_cell_annotation_myeloids_PMC9829493<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  此函数中的marker的来源为论文“Heterogeneity of tumor-infiltrating myeloid cells in era of single-cell genomics”
####  期刊为Chin J Cancer Res. 2022 Dec 30。影响因子为4分
####  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9829493/
	Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,
                                      vjust = 0.5,hjust=0.5))
  myeloids = list(
    Mast_KIT=c('TPSAB1','CPA3','KIT','MS4A2'),
    pDC_LILRA4=c('MZB1','LILRA4','IRF7','IL3RA'),
    cDC1_CLEC9A=c('XCR1','CLEC9A','BATF3','IRF8'),
    cDC2_CD1C=c("CLEC10A",'FCER1A','CD1C','CD1E'),
    DC_LAMP3 =  c('LAMP3','BIRC3','CCR7','FSCN1' ),
    Mono_CD14 = c('CD14','S100A8','VCAN'),
    Mono_CD16 = c('FCGR3A','LILRB2','LST1'),
    Macro_FCN1=c("FCN1","CD163",'CD68'),
    Macro_LYVE1=c("LYVE1",'PLTP','SEPP1'),
    Macro_NLRP3 =  c("NLRP3","IL1B",'CXCL2','EREG' ),
    Macro_PPARG = c("PPARG","MARCO",'MRC1','MSR1'),
    Macro_ISG15 = c("ISG15","CXCL10",'IFITM3','GBP1'),
    Macro_C1QC = c( "C1QA","C1QB",'C1QC','APOE'),
    Macro_SPP1 =  c("SPP1","VEGFA",'GPNMB','FN1')
  )
  if(species =='mouse')
{
library(stringr)
   for(i in 1:length(myeloids))
   {
	   #i=1
	   myeloids[[i]]<-str_to_title(myeloids[[i]])
	}
	myeloids[[8]][1]<-'Fcnb'
}
  plota<-NULL
  plota <- DotPlot(pbmc,features = myeloids,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=plota,filename=paste0("fig",plotid,".PMC9829493.tiff"),compression='lzw',limitsize=F,width = 20,height =ph)
}

##  seurat_cell_annotation_myeloids(pbmc,plotid=20)
seurat_cell_annotation_myeloids_PMID33545035_first<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  此函数中的marker的来源为论文“A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells”
####  期刊为《Cell》. 2021_2_4。影响因子为64.5分
####  https://pubmed.ncbi.nlm.nih.gov/33545035/
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 0.5, hjust=0.5))
myeloids = list(
  Mast	=c("TPSB2","TPSAB1","KIT","CLU","VWA5A"),
  cDC	=c("CLEC10A","CD1C","CD74","FCER1A","HLA-DQB1"),
  pDC 	=c("TCF4","IRF8","IRF7","GZMB","LILRA4" ),
  Mo_Mq =c("FCN1","CD68","FCGR3A","S100A9","CD14")
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(myeloids))
   {
	   #i=1
	   myeloids[[i]]<-str_to_title(myeloids[[i]])
	}
}
plota<-NULL
plota <- DotPlot(pbmc, features = myeloids,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=plota, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 20, height =ph)
}

##  seurat_cell_annotation_myeloids(pbmc,plotid=20)
seurat_cell_annotation_myeloids_PMID33545035<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  此函数中的marker的来源为论文“A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells”
####  期刊为《Cell》. 2021_2_4。影响因子为64.5分
####  https://pubmed.ncbi.nlm.nih.gov/33545035/
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 0.5, hjust=0.5))
myeloids = list(
  Mast_KIT		=c("KIT","TPSAB1","CPA3"),
  pDC_LILRA4	=c("LILRA4","GZMB","IL3RA","FAM129C","IGJ","PTPRCAP","IGKC","IRF4"),
  cDC1_CLEC9A	=c("CLEC9A","FLT3","IDO1","XCR1","CADM1","C1orf54","SLAMF8"),
  cDC2_CD1C		=c("CD1C","FCER1A","HLA-DQA1","CD1E","CLEC10A"),
  cDC3_LAMP3	=c("LAMP3","CCR7","FSCN1","CST7","IL4I1"),
  Mono_CD14		=c("CD14","FCN1","S100A9","S100A8"),
  Mono_CD16		=c("FCGR3A","LST1","LILRB2","HK3","PILRA"),
  Macro_C1QC	=c("C1QC","C1QA","APOE","TREM2","SLCO2B1","APOC1","RNASE1"),
  Macro_CX3CR1=c("CX3CR1"),
  Macro_FN1		=c("FN1"),
  Macro_GPNMB	=c("GPNMB"),
  Macro_IL1B	=c("IL1B"),
  Macro_INHBA	=c("INHBA","IL1RN","CCL4"),
  Macro_ISG15	=c("ISG15","CCL2","CXCL10"),
  Macro_LYVE1	=c("LYVE1","PLTP","SEPP1"),
  Macro_NLRP3	=c("NLRP3","EREG"),
  Macro_PPARG	=c("PPARG"),
  Macro_SPP1	=c("SPP1"),
  Macro_VCAN	=c("Vcan","VEGFA","OLR1","TREM1"),
  Macro_rtms	=c("THBS1")
  ####  https://mp.weixin.qq.com/s/-w0e9IhNcYyWUMqajJ9xqQ
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(myeloids))
   {
	   #i=1
	   myeloids[[i]]<-str_to_title(myeloids[[i]])
	}
}
plota<-NULL
plota <- DotPlot(pbmc, features = myeloids,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=plota, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 30, height =ph)
}


seurat_cell_annotation_myeloids_panzhong<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
	Idents(pbmc) <- group.by
	nc<-length(levels(pbmc))
	ph<-max(nc/2,6)
	th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust=0.5))
	markers_first = list(
	Immune = c("PTPRC"),
	myeloid =  c('LYZ','LYZ1','LYZ2'),
	Mono_C1=c("CD14",'LY6C'),         
	Mono_C2=c('FCGR3A','MS4A7','FCGR3B','Fcgr4','CX3CR1'),         
	Mono_C3=c('HLA-DRA','ITGAX','Cd209','H2'),         ##13
    ### CD11C-->ITGAX
	Macro_m_Mono=c('Adgre1','CD68','FCER1G','MARCO'),
	Mac=c('CD86','CD163','TREM2'),
	preDC=c('FLT3','CSF1R','BIRC3'),            #c('FLT3','CSF1R','ITGAX','BIRC3')
	cDC1 = c('CD8A','CLEC9A','ITGAE','THBD','XCR1'),                                        ##8
	cDC2 = c('CD1C','CD207','ITGAM','NOTCH2','SIRPA'), 
	pDC = c('CLEC4C','LILRB4','NRP1','CCR7','B220','SiglecH','KDR','FCRLA'), 
	moDC = c('MRC1'),                   #  c('MRC1','CD209'),                            ##13
	LC = c('CD1A','ID2'),                                                  
	Mast =  c('MS4A2','TPSB2',"TPSAB1",'CST3','KIT'),
	Neutrophils =  c("CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2",'CSF3R'),           ###16
	Basophils =  c('CD34','CD200R3','FCER1A','CD9',"CPA3",'ITGB2','IL3RA'),
	Eosinophils =  c('SIGLEC5', 'IL5RA', 'CCR3', 'EPX')
	)
if(species =='mouse')
{
	library(stringr)
	for(i in 1:length(markers_first))
	{
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
	#markers_first[[3]][5]<-'Fcnb'
}
  plota<-NULL
  plota <- DotPlot(pbmc,features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=plota,filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 20,height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}
}

seurat_cell_annotation_myeloids_panzhong_dup<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
	Idents(pbmc) <- group.by
	nc<-length(levels(pbmc))
	ph<-max(nc/2,6)
	th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust=0.5))
	markers_first = list(
	Immune = c("PTPRC"),
	myeloid =  c('LYZ','LYZ1','LYZ2'),
	Mono_Classical=c("CD14",'LY6C1','CCR2','CCR5','CD62L'),         
	Mono_NC=c('FCGR3A','FCGR3B','Fcgr4','CX3CR1','HLA-DR','CD43','Treml4'),         
	Mono_Med=c("CD14",'FCGR3A','FCGR3B','Fcgr4','HLA-DR','CD68','ITGAX','Cd209','H2'),         ##13
    ### CD11C-->ITGAX
	Macro_m_Mono=c('Adgre1','CD68','FCER1G','MARCO','CD80'),
	Macrophages=c('CD68','CD86','CD163','FCER1G','MARCO','RETNLB','Retnla','CLEC7A','Mgl2'),     
	preDC=c('FLT3','CSF1R','ITGAX','BIRC3'),
	cDC1 = c('CD8A','CLEC9A','ITGAE','THBD','XCR1'),                                        ##8
	cDC2 = c('CD1C','CD207','ITGAM','NOTCH2','SIRPA'), 
	pDC = c('CLEC4C','LILRB4','NRP1','CCR7','B220','SiglecH'), 
	moDC = c('MRC1'),                   #  c('MRC1','CD209'),                            ##13
	LC = c('CD1A','ID2'),                                                  
	Mast =  c('MS4A2','TPSB2',"TPSAB1",'CST3','KIT'),
	Neutrophils =  c("CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2",'CSF3R'),           ###16
	#https://mp.weixin.qq.com/s/lUZeQEPUSaoDumMnpspMQg
	#https://pubmed.ncbi.nlm.nih.gov/32719519/
	Basophils =  c('CD34','CD200R3','FCER1A','CD9',"CPA3",'ITGB2','IL3RA'),
	####      CD31-->PECAM1
	Bas_m_Mast =  c("ITGB2",'PECAM1','IL3RA'),
	####        CD203c-->ENPP3
	Bas_m_Neu =  c("ENPP3"),
	Eosinophils =  c('SIGLEC5', 'IL5RA', 'CCR3', 'EPX')
	)
if(species =='mouse')
{
	library(stringr)
	for(i in 1:length(markers_first))
	{
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
	#markers_first[[3]][5]<-'Fcnb'
}
  plota<-NULL
  plota <- DotPlot(pbmc,features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=plota,filename=paste0("fig",plotid,".markers.tiff"),compression='lzw',limitsize=F,width = 20,height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}
}



seurat_cell_annotation_Macrophages_M1M2<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
#####  https://stemcellres.biomedcentral.com/articles/10.1186/s13287-023-03302-x
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))

df_genes<-get_M1M2_marker()
m1_genes<-df_genes[df_genes$geneset=='M1','gene']
m2_genes<-df_genes[df_genes$geneset=='M2','gene']

markers_first = list(
	M1=m1_genes,
	M2=m2_genes
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.M1M2.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}




seurat_cell_annotation_T_cells_panzhong<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  此函数中的marker的来自潘重2023年6月之前的单细胞数据分析中T细胞的总结
####  其中大部分来自Nature review
####  印象笔记的链接为：https://app.yinxiang.com/fx/ef58d501-3d50-4486-a29a-e07cc1bfd265
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  general=c('PTPRC','CD3D','CD3E','CD4','IL7R','CD8A','CD8B'),
  γδT=c('TRGV9','TRDV2','TRDC'),
  Naive=c('CCR7','SELL','CD5'),  
  Tm=c('CD44','S100A4','GPR183'),
  Te=c('FASLG','FAS'),
  Th1=c('IFNG','IL2','LTA'),
  Th2=c('IL4','IL5','IL13'),
  Th17=c('IL17A','RORA','RORC'),
  Th9=c('IL9'),
  Th22=c('IL22'),
  Tfh=c('CXCR5','ICOS','IL21','PDCD1','BCL6'),
  Treg=c('FOXP3','IL2RA'),
  NKT=c('CD1D','KLRF1','NKG7','GNLY'),
  #https://www.nature.com/articles/s41598-020-76659-1
  MAIT=c('SLC4A10','KLRB1'),
  Tex=c('CXCL13'),
  Tgd=c('TRA','TRB','TRG','TRD'),
  ILCs=c('AREG','TNFRSF25')  
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
	markers_first[[5]][1]<-'Fasl'
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}


seurat_cell_annotation_T_cells_panzhong_simple<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
####  https://www.jianshu.com/p/0127c9b380c9
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  general=c('PTPRC','CD3D','CD3E','CD4','CD8A','CD8B','CD8B1'),
  γδT=c('TRGV9','TRDV2','TRDC'),
  Naive=c('TCF7','CCR7','SELL','CD5'),  
  Tcm=c('IL2'),
  Tem=c('IFNG','IL4','CX3CR1','IL7R'),
  Temra=c('KLRG1'),
  Tm=c('CD44','S100A4','GPR183'),
  Trm=c('CD69','ITGAE','ITGA1'),
  Te=c('FASLG','FAS'),
  Treg=c('FOXP3','IL2RA'),
  Tfh=c('IL21'),
  Th1=c('LTA','TNFA'),
  Th2=c('IL5','IL13'),
  Th17=c('IL17A','IL17F','RORA','RORC'),
  Th9=c('IL9','IL10'),
  Th22=c('IL22'),
  NKT=c('CD1D','KLRF1','NKG7','GNLY'),
  #https://www.nature.com/articles/s41598-020-76659-1
  MAIT=c('SLC4A10','KLRB1'),
  Tex=c('CXCL13','CTLA4'),
  ILCs=c('AREG','TNFRSF25')  
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
	markers_first[[9]][1]<-'Fasl'
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}
}


seurat_cell_annotation_T_cells_PMID33958794<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
  ####  https://pubmed.ncbi.nlm.nih.gov/33958794/
  Idents(pbmc) <- group.by
  nc<-length(levels(pbmc))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
  markers_first = list(
    CD4_Tn=c("LEF1","TCF7","SELL"),
    CD4_Tem=c("IL7R","CD40LG","ANXA1","FOS","JUN"),
    CD4_Treg=c("FOXP3","SAT1","IL2RA"),
    CD4_Tex=c("CTLA4","PDCD1","CXCL13","CD200","TNFRSF18"),
    CD8_Tn=c("CCR7","NELL2","CD55","KLF2"),
    CD8_Trm=c("TOB1","ZNF683","CCL5"),
    CD8_Tem=c("GZMK","EOMES","ITM2C"),
    CD8_Temra=c("CX3CR1","GNLY","GZMH"),
    CD8_Tex=c("GZMB","LAG3","CCL4L2"),
    NK_cyto=c("FCGR3A","FGFBP2","TYROBP"),
    NK_rest=c("AREG","XCL1","KLRC1"),
    Vγ9_Vδ2_Tγδ=c("TRDV2","TRGV9","MTRNR2L8","KLRD1"),
    Tγδ=c("TRDV1","KLRC3","CTSW","CD7"),
    Proliferating=c("MKI67","STMN1","TUBA1B","HIST1H4C")
  )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-str_to_title(markers_first[[i]])
    }
    markers_first[[5]][1]<-'Fasl'
  }
  p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}

seurat_cell_annotation_T_cells_zhangzhemin<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  此函数中的marker的来自张泽民文章：Lineage tracking reveals dynamic relationships of T cells in colorectal cancer
####  Nature. 2018 Dec   IF：69.5
####  https://pubmed.ncbi.nlm.nih.gov/30479382/
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 1, hjust=1))
markers_first<-list(
general=c('PTPRC','CD3D','CD3E','CD4','CD8A','CD8B'),
#7
cd8.Tn=c('CCR7','SELL','LEF1','TCF7'),
#11
cd8.Tcm=c('IL7R','CD27','CD28','PRF1','GZMA','CCL5','GPR183','S1PR1'),
#19
cd8.Temra=c('KLRG1','CX3CR1','FCGR3A','FGFBP2','GZMH','TBX21','EOMES','S1PR5'),
#27
cd8.Tem=c('GZMK','CXCR4','CXCR3','CD44'),
#31
cd8.Trm=c('CD6','XCL1','XCL2','MYADM','CAPG','RORA','NR4A1','NR4A2','NR4A3','CD69','ITGAE'),
#42
cd8.IEL=c('CD160','KIR2DL4','TMIGD2','KLRC1','KLRC2','KLRC3','IKZF2','ENTPD1'),
#50
cd8.Tex=c('HAVCR2','CXCL13','PDCD1','LAYN','TOX','IFNG','GZMB','MIR155HG','TNFRSF9'),
#59
cd8.MAIT=c('SLC4A10','KLRB1','ZBTB16','NCR3','RORC')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,"a.dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)

markers_first<-list(
general=c('PTPRC','CD3D','CD3E','CD4','CD8A','CD8B'),
#7
cd4.Tn=c('CCR7','SELL','LEF1','TCF7'),
#11
cd4.BTcm=c('PTGER2','ICAM2','ANXA1','ANXA2','S1PR1'),
#16
cd4.Temra=c('KLRG1','CX3CR1','NKG7','PRF1','GNLY','GZMH','TBX21','CTSW','S1PR5'),
#25
cd4.NTcm=c('RGS1','CD69'),
#27
cd4.Trm=c('KLRB1','PTGER4','IL7R','CXCR6','NR4A1','NR4A2','NR4A3','MYADM'),
#35
cd4.Tfh=c('CXCR5','BCL6','ICA1','TOX','TOX2','IL6ST','MAGEH1','BTLA','ICOS','PDCD1','CD200'),
#46
cd4.Tem=c('GZMK','GZMA','CCL5','IFNG','RUNX3','EOMES','CXCR3','CXCR4','CD44'),
#55
cd4.Th17=c('IL23R','RORC','IL17A','FURIN','CTSH','CCR6','CAPG','ITGAE'),
#63
cd4.Th1=c('CXCL13','BHLHE40','GZMB','HAVCR2','IGFLR1'),
#68
cd4.Treg=c('FOXP3','IL2RA','IL10RA','IKZF2','RTKN2','CDC25B','S1PR4'),
#75
cd4.Tfr=c('IL10','CCR4'),
#77
cd4.Ttreg=c('CCR8','TNFRSF18','LAYN','TNFRSF9','CTLA4','BATF','IL21R')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,"b.dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)

}


seurat_cell_annotation_T_cells_PMID32385277<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
####  https://pubmed.ncbi.nlm.nih.gov/32385277/
####  https://pubmed.ncbi.nlm.nih.gov/38200551/
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  general=c('PTPRC','CD3D','CD3E','CD4','IL7R','CD8A','CD8B','CD8B1'),
  NK=c('CD1D','KLRF1','NKG7','GNLY'),
  γδT=c('TRDC','TRGC1','TRGC2'),  
  Naive=c("TCF7", "SELL", "LEF1", "CCR7"),  
  Exhausted=c("LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4",'LAYN'),     ## 23
  Cytotoxic=c("IL2", "GZMA", "PRF1", "GZMB", "GZMK"),
  Treg=c("IL2RA", "FOXP3", "IKZF2", "TGFB1", "TGFB3", "TGFBI", "TGFBR1"),    #36
  Tfh=c("MAF", "CXCL13", "CXCR5"),
  Th17=c('IL17A','IL17F','RORA','RORC',"IRF4", "CREM", "NR4A2","BHLHE40"),   ##46
  Th1=c("STAT4", "IL12RB2", "IFNG"),
  Th2=c("GATA3", "STAT6", "IL4")  
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}
}


seurat_cell_annotation_NK<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,
                                    vjust = 1, hjust=1))
markers_first = list(
  Immune = c("PTPRC"),
  Lymphoid = c("CD69",'MKI67'),
  Tcells = c('CD3D','CD3E','CD3G','CD247','CD4','CD8A','CD8B'),
  NK=c('CD14','CD19','NCAM1','FCGR3A','GNLY','KLRD1','NKG7','FGFBP2','CX3CR1','KLRF1','NCR3'),
  #https://www.biocompare.com/Editorial-Articles/576304-A-Guide-to-NK-Cell-Markers/
  ILCs=c('AREG','TNFRSF18','TNFRSF25'),
  Bcells=c( 'CD79A','MS4A1')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}



seurat_cell_annotation_fibroblasts_panzhong<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  https://mp.weixin.qq.com/s/2byK9WIdiUNDPpNlutdyJw
####  https://mp.weixin.qq.com/s/QED03-6_y9SZGUXYXoF96Q
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  general=c('S100A4','FGF7','ACTA2','COL1A1','COL1A2','TAGLN'),
  vascular_CAF=c('GJA4','RGS5'),
  matrix_CAF=c('LUM','DCN','VCAN'),  
  inflammatory_CAF=c('C3','C7'),
  antigen_presenting_CAF=c('CD74','HLA-DRA','HLA-DRB1'),
  EMT_like_CAF=c('KRT19','KRT8'),
  lip_ofibroblast=c('APOA2','FABP1','FABP4','FRZB')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = 'seurat_clusters' )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}

seurat_cell_annotation_fibroblasts_panzhong2<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
####  此函数中的marker的来自潘重2023年6月之前的单细胞数据分析的总结,后续再不停更新
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  Epithelial=c('WFDC2','PAX8','EPCAM','GPC3','MMP7'),
  Fibroblasts=c('PDPN','FN1','VIM','DCN','GSN','COL1A1','COL1A2','COL4A1','COL4A2'),
  SMC =  c('DES','ACTA2','TAGLN','MYH11','MYLK','ACTG2'),
  Pericyte=c('BCAM','PDGFRB','KCNJ8','RGS5','ABCC9'),
  Endothelial =  c('MME','PECAM1','THBD','CDH5','ENG','VWF','CCL21'),
  mesothelial =  c('LY6G','MSIN'),
  mesenchymal =  c('PDGFRA','TCF21','CDH11'),                   ###  32
  Neuronal=c('NRXN1'),
  Immune = c("PTPRC")
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
	#markers_first[[21]][1]<-'Fcgr4'
	#markers_first[[21]][1]<-'Siglecf'
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}

seurat_cell_annotation_fibroblasts_mhnb<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  https://www.jianshu.com/p/5efd228b9b76
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  general=c('S100A4','FGF7','COL1A1','COL1A2'),
  myofibroblast=c('ACTA2','FAP','PDPN'),
  iCAFs=c('C3','C7','CXCL12'),
  dPVAL=c('TAGLN','MYH11','MYLK'),
  imPVL=c('PDGFRB','CD36','RGS5')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = 'seurat_clusters' )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}

seurat_cell_annotation_fibroblasts_nature_cancer<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  https://mp.weixin.qq.com/s/R3Ui3Kg5WAFcd4PEY3NxCg
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  general=c('S100A4','FGF7','COL1A1','COL1A2'),
  myofibroblast=c('ACTA2','MMP3','IL4'),
  iCAFs=c('CXCL1','CXCL2','CXCL12','IL6','CFD','C1QA'),
  aCAFs=c('CD74','HLA-DRA','HLA-DRB'),
  vCAF=c('MCAM'),
  meCAF=c('PLA2G2A'),
  cd63CAF=c('CD63'),
  cox2CAF=c('PTGS2')
)
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = 'seurat_clusters' )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}


seurat_cell_annotation_fibroblasts_PMID33981032<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  https://www.nature.com/articles/s41586-021-03549-5
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  Fib_PI16=c('PI16','CD34','HAS1','PLIN2'),
  Fib_NPNT=c('NPNT','CES1','ADH1B','FGFR4'),
  Fib_LRRC15=c('LRRC15','COL11A1','CTHRC1'),
  Fib_ADAMDEC1=c('ADAMDEC1','CXCL14'),
  Fib_CCL19=c('CCL19','GREM1','TNFSF13B'),
  Fib_COL3A1=c('COL3A1','POSTN')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = 'seurat_clusters' )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}


seurat_cell_annotation_fibroblasts_PMID34857954<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  https://www.nature.com/articles/s41588-021-00972-2
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  fibroblasts=c('IGF1','PCOLCE','C7','OGN','MMP11','SFRP1','DKK1','FOXO1')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-str_to_title(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = 'seurat_clusters' )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}



seurat_cell_annotation_MSC<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
####  MSC Mesenchymal Stem Cell
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  Epithelial=c('WFDC2','PAX8','EPCAM','GPC3','MMP7'),
  Fibroblasts=c('PDPN','FN1','VIM','DCN','GSN','COL1A1','COL1A2','COL4A1','COL4A2'),
  SMC =  c('DES','ACTA2','TAGLN','MYH11','MYLK','ACTG2'),
  Pericyte=c('BCAM','PDGFRB','KCNJ8','RGS5','ABCC9'),  
  Endothelial =  c('MME','PECAM1','THBD','CDH5','VWF','CCL21'),        ## 31
  ## CD73 -->NT5E  CD90 --> THY1   CD105 --> ENG
  MSC =  c('ENG','NT5E','THY1','BSG'),                   ###  35
  MSC_neg =  c('PTPRC','CD34','CD14','CD19','HLA-DRA'),                   ###  40
  mesothelial =  c('LY6G','MSIN'),
  mesenchymal =  c('PDGFRA','TCF21','CDH11'),                   ###  45
  ###  CD29 --> ITGB1   CD73 -->NT5E  CD90 --> THY1   CD105 --> ENG
  progenitor =  c('ITGB1','CD44')                   ###  50
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
	if(plot.full)
	{
		markers<-as.character(unlist(markers_first))
		seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
	}

}



seurat_cell_annotation_Endothelial_panzhong1<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  https://zhuanlan.zhihu.com/p/614994990
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  general=c('CDH5','PECAM1'),
  aCap_EC=c('AW112021','EDNRB','CA4'),
  gCap_EC=c('SEMA3C','KIT','GPIHBP1'),  
  artery_EC=c('MGP','BMX'),
  vein_EC=c('PRSS23','SLC6A2'),
  lymphatic_ECs=c('PROX1','PDPN','MMRN1')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = 'seurat_clusters' )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}


seurat_cell_annotation_Endothelial_panzhong2<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  https://zhuanlan.zhihu.com/p/614994990
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
  general=c('CDH5','PECAM1'),
  lymphatic_ECs=c('CCL21','PROX1'),
  arteries_EC=c('HEY1', 'IGFBP3'), 
  capillaries_EC=c('CD36', 'CA4'),
  veins_EC=c('ACKR1')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = 'seurat_clusters' )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}


seurat_cell_annotation_Endothelial_PMID36127427<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
####  https://pubmed.ncbi.nlm.nih.gov/36127427/
####  Single cell atlas identifies lipid-processing and immunomodulatory endothelial cells in healthy and malignant breast
####  见附件“41467_2022_33052_MOESM5_ESM.xlsx”
Idents(pbmc) <- group.by
nc<-length(levels(pbmc))
ph<-max(nc/2,6)
th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
markers_first = list(
	general=c('MME','VWF','PECAM1','THBD','CDH5','ENG'),
EC1_LEC=c('CCL21','PROX1','LYVE1', 'PDPN', 'TFF3'),
EC2_artery=c('HEY1', 'CLU','ELN','IGFBP3'), 
EC3_artery=c('CXCL12', 'AMD1', 'FN1'),
EC4_angiogenic=c('KDR', 'RGCC', 'EDNRB','VWA1', 'ESM1', 'INSR','COL4A2', 'COL4A1', 'MCAM', 'NOTCH4'),
EC5_transitioning=c('PLCG2', 'MTRNR2L1'), 
EC6_vein=c('ACKR1','POSTN', 'VCAN', 'NR2F2', 'C7', 'CYP1B1', 'APLNR'),
EC7_vein=c('LAMA5'),
EC8_vein=c('HLA-DQA2', 'HLA-DQA1', 'VCAM1', 'PLAT'),
EC9_vein=c('EDN1', 'CYTL1'),
EC10_capillary=c('CD36','ICAM1', 'CCL2', 'SELE'),
EC11_capillary=c('ID2', 'SPP1'),
EC12_capillary=c('CA4', 'FABP4', 'LPL', 'CLDN5', 'FABP5')
)
if(species =='mouse')
{
library(stringr)
   for(i in 1:length(markers_first))
   {
	   #i=1
	   markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
	}
}
p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = 'seurat_clusters' )  +th
ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}





## pbmc.x<-seurat_cell_annotation_manual(pbmc.x,filename='celltype_manual.txt',group.by = 'seurat_clusters',celltype_column='celltype_sub')
## table(pbmc.x@meta.data$celltype_sub)
seurat_cell_annotation_manual<- function(pbmc,filename='celltype_manual.txt',group.by = 'seurat_clusters',celltype_column='celltype',plotid='04'){
	
df_celltype_manual <- read.delim(file = filename,sep='\t',header=F)
pbmc@meta.data[[celltype_column]]<- 'unknown'
pbmc@meta.data[[celltype_column]]<-as.character(pbmc@meta.data[[celltype_column]])
pbmc@meta.data[[celltype_column]]<- 'unknown'
for(i in 1:nrow(df_celltype_manual))
{
	pbmc@meta.data[which(pbmc@meta.data[[group.by]] ==df_celltype_manual[i,1]),celltype_column] <- df_celltype_manual[i,2]
}
plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = celltype_column,pt.size = 0.3,label=T,raster=F)
ggsave(paste0("fig",plotid,'d.umap.',celltype_column,'.tiff'), plot = plota, width = 8, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,'d.umap.',celltype_column,'.pdf'), plot = plota, width = 8, height = 6)

return(pbmc)
}


# pbmc.x<-seurat_cell_annotation_bycluster(pbmc.x,group.by = 'seurat_clusters',celltype_column='celltype_sub',plotid='04')
seurat_cell_annotation_bycluster<- function(pbmc,group.by = 'seurat_clusters',celltype_column='celltype',plotid='04'){
class(pbmc@meta.data$seurat_clusters)
	clusters<-levels(pbmc.x@meta.data$seurat_clusters)
	class(clusters)
	clusters_num<-as.numeric(clusters)
	clusters_level<-paste0('C', clusters_num+1)
	
pbmc@meta.data$seurat_clusters<-as.character(pbmc@meta.data$seurat_clusters)
pbmc@meta.data$seurat_clusters<-as.numeric(pbmc@meta.data$seurat_clusters)
pbmc@meta.data[[celltype_column]]<-paste0('C', pbmc@meta.data$seurat_clusters+1)
pbmc@meta.data$seurat_clusters<-factor(pbmc@meta.data$seurat_clusters)
pbmc@meta.data[[celltype_column]]<-factor(pbmc@meta.data[[celltype_column]],levels=clusters_level)

plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = celltype_column,pt.size = 0.3,label=T)
ggsave(paste0("fig",plotid,'d.umap.',celltype_column,'.tiff'), plot = plota, width = 8, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,'d.umap.',celltype_column,'.pdf'), plot = plota, width = 8, height = 6)
return(pbmc)
}


#seurat_celltype_check(pbmc,plotid=20,group.by = 'group',celltype='celltype_manual',celltypes=c('Myeloid cells','Endothelial cells'))
#seurat_celltype_check(pbmc,plotid=22,group.by = 'group',celltype='celltype_manual',celltypes=c('T cells','Fibroblasts'))

seurat_celltype_check<- function(pbmc,plotid=20,group.by = 'group',celltype='celltype_manual',celltypes=NULL,species='human'){
		print(table(pbmc@meta.data[,celltype]))
		list_sc<-list()
		if(is.null(celltypes))
		{
			celltypes<-unique(pbmc@meta.data[,celltype])
		}
		print(celltypes)
		ncelltype=length(celltypes)
		metadata<-pbmc@meta.data
		for(i in 1:ncelltype)
		{
			#i<-1
			plotid_new=plotid+i
			celltype_name<-gsub(' ','_',celltypes[i])
			print(celltype_name)
			cells_used<-rownames(metadata[metadata[[celltype]]==celltypes[i],])
			list_sc[[i]]<-subset(x = pbmc, cells = cells_used)
			ncells<-dim(list_sc[[i]])[2]
			if(ncells>500)
			{
				list_sc[[i]]<-seurat_integrate_harmony(list_sc[[i]],group.by=group.by,res=resolution,plotid=plotid_new,species=species)
				seurat_cell_annotation_first(list_sc[[i]],plotid=paste0(plotid_new,'g'),species=species)
				seurat_cell_annotation_pbmc(list_sc[[i]],plotid=paste0(plotid_new,'h'),species=species)
				metadatai<-list_sc[[i]]@meta.data
				write.table(metadatai,file=paste0("fig",plotid_new,'.',celltype_name,".metadata.txt"),quote=F,sep="\t", row.names=T)
				output_rds=paste0("fig",plotid_new,'.',celltype_name,".rds")
				saveRDS(list_sc[[i]],file=output_rds)
			}else if(ncells>100){
				list_sc[[i]]<-seurat_integrate_null(list_sc[[i]],group.by=group.by,res=resolution,plotid=plotid_new,species=species)
				seurat_cell_annotation_first(list_sc[[i]],plotid=paste0(plotid_new,'g'),species=species)
				seurat_cell_annotation_pbmc(list_sc[[i]],plotid=paste0(plotid_new,'h'),species=species)
				metadatai<-list_sc[[i]]@meta.data
				write.table(metadatai,file=paste0("fig",plotid_new,'.',celltype_name,".metadata.txt"),quote=F,sep="\t", row.names=T)
			}
			else{
				seurat_cell_annotation_first(pbmc,plotid=paste0(plotid_new,'g'),group.by = 'seurat_clusters',species=species)
				seurat_cell_annotation_first(pbmc,plotid=paste0(plotid_new,'h'),group.by = celltype,species=species)
			}
		}
}

seurat_splitby_celltype<- function(pbmc,projectid='projectid',celltype='celltype_manual',celltypes=c('Myeloid cells'),plotid=10){
		list_sc<-list()
		ncelltype=length(celltypes)
		metadata<-pbmc@meta.data
		for(i in 1:ncelltype)
		{
			#i<-1
			plotid_new=plotid+i
			celltype_name<-gsub(' ','_',celltypes[i])
			cells_used<-rownames(metadata[metadata[[celltype]]==celltypes[i],])
			list_sc[[i]]<-subset(x = pbmc, cells = cells_used)
			counts <- GetAssayData(object = list_sc[[i]], layer = "counts")
			nonzero <- counts > 0
			selected_gene <- rownames(nonzero)[Matrix::rowSums(nonzero) >= min.cells]
			list_sc[[i]] <- subset(list_sc[[i]], features = selected_gene)
			rm(counts,nonzero,selected_gene)
			gc()
			dir.create(celltype_name)
			output_rds=paste0(celltype_name,"/Seurat.",projectid,'.',celltype_name,".rds")
			saveRDS(list_sc[[i]],file=output_rds)
			file.copy('seurat_parameters.tsv',paste0(celltype_name,'/','seurat_parameters.tsv'))
		}
	}




#filename='fig21.Myeloid_cells.metadata.txt'
#pbmc<-seurat_celltype_update(pbmc,filename=filename,group= c(9),celltype_column='celltype_manual',celltype='T cells')

seurat_celltype_update<- function(pbmc,filename='',group= c(9),celltype_column='celltype_manual',celltype='Myeloid cells'){
			metadata<-read.delim(filename,sep='\t',header=T,row.names=1)
			cells_used<-rownames(metadata[metadata$seurat_clusters %in% group,])
			pbmc@meta.data[which(rownames(pbmc@meta.data) %in% cells_used),celltype_column] <- celltype
			return(pbmc)
}



# pbmc@meta.data$celltype_nmf<-pbmc@meta.data$celltype_manual
# pbmc@meta.data$celltype_nmf<-as.character(pbmc@meta.data$celltype_nmf)

# filename='fig04.GSE216651.Macrophages.NMF.metadata.tsv'
# pbmc<-seurat_celltype_update_sub(pbmc,filename=filename,celltype_column='celltype_nmf',celltype_sub='celltype_nmf',prefix='Macro_')

seurat_celltype_update_sub<- function(pbmc,filename='',celltype_column='celltype_sub',celltype_sub='Myeloid cells',prefix='Macro_')
{
metadata<-read.delim(filename,sep='\t',header=T,row.names=1)
cells_used<-rownames(metadata)
for(i in 1:nrow(metadata))
{
	# i=1
	pbmc@meta.data[which(rownames(pbmc@meta.data) ==rownames(metadata)[i]),celltype_column] <- paste0(prefix,metadata[i,celltype_sub])
}
return(pbmc)
}



#pbmc<-seurat_celltype_callback(pbmc=pbmc,pbmc.x=pbmc.y,celltypea='celltype4',celltypeb='celltype4',plotid='14e')
seurat_celltype_callback<- function(pbmc,pbmc.x,celltypea='celltype',celltypeb='celltype',plotumap=T,plotid='06',prefix='Macro_'){	
		if(! celltypea %in% colnames(pbmc@meta.data))
		{
			print(paste0('not find ',celltypea))
			pbmc@meta.data[[celltypea]]<-'unknown'
		}		
		pbmc.x@meta.data[[celltypeb]]<-as.character(pbmc.x@meta.data[[celltypeb]])		
		for(i in 1:nrow(pbmc.x@meta.data))
		{
			pbmc@meta.data[rownames(pbmc@meta.data)==rownames(pbmc.x@meta.data)[i],celltypea]<-paste0(prefix,pbmc.x@meta.data[i,celltypeb])
		}
		if(plotumap)
		{
        plota<-NULL
		plota<-DimPlot(pbmc, reduction = "umap", group.by = celltypea,pt.size = 0.5,label=F)
		ggsave(paste0("fig",plotid,"d.umap.",celltypea,".tiff"), plot = plota, width = 8, height = 8,compression='lzw')
		ggsave(paste0("fig",plotid,"d.umap.",celltypea,".pdf"), plot = plota, width = 8, height = 8)
	}
		return(pbmc)
}
	
	
seurat_integrate_harmony <- function(pbmc,group.by='orig.ident',species = "human",
                         mt.pattern="^MT-",mt.list=NULL,dim.use=20,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=2000,
                         res=0.5,plotid='02') {
library(harmony)
Idents(pbmc) <- group.by
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nfeatures)
features_variable<-VariableFeatures(object = pbmc)
features_all <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = features_variable)
pbmc <- RunPCA(pbmc, features = features_variable)

plota<-NULL
plota<-VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
ggsave(paste0("fig",plotid,"b.pcaGene.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"b.pcaGene.pdf"), plot = plota, width = 10, height = 10)

plota<-NULL
plota<-ElbowPlot(pbmc, ndims = 50)
ggsave(paste0("fig",plotid,"c.pca50.tiff"), plot = plota, width = 10, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"c.pca50.pdf"), plot = plota, width = 10, height = 6)

#绘制主成分分析图形
plota<-NULL
plota<-DimPlot(object = pbmc, reduction = "pca")
ggsave(paste0("fig",plotid,"d.PCA.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"d.PCA.pdf"), plot = plota, width = 10, height = 10)

pbmc.integrated <- RunHarmony(pbmc, group.by , plot_convergence = F,dims.use = 1:dim.use)
pbmc.integrated = RunUMAP(pbmc.integrated, reduction = "harmony", dims = 1:dim.use)

pbmc.integrated = FindNeighbors(pbmc.integrated, reduction = "harmony",dims = 1:dim.use)
pbmc.integrated = FindClusters(pbmc.integrated, resolution = res)

plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 6)

meta<-pbmc.integrated@meta.data
if('group' %in% colnames(meta))
{
plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e1.umap.group.pdf"), plot = plota, width = 8, height = 8)

plota<-NULL
plota <- DimPlot(pbmc.integrated, reduction = "umap", split.by = "group")
ggsave(paste0("fig",plotid,"e2.umap.group.split.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e2.umap.group.split.pdf"), plot = plota, width = 12, height = 6)
}

plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e3.umap.sample.pdf"), plot = plota, width = 8, height = 8)

nsample<-length(unique(pbmc.integrated@meta.data$orig.ident))
nh<-ceiling(nsample/4)
plotc<-NULL
plotc <- DimPlot(pbmc.integrated, reduction = "umap", split.by = "orig.ident",ncol=4)
ggsave(paste0("fig",plotid,"e4.umap.sample.split.tiff"), plot = plotc, width = 12, height = nh*3,compression='lzw')
ggsave(paste0("fig",plotid,"e4.umap.sample.split.pdf"), plot = plotc, width = 12, height = nh*3)

plota<-NULL
plota<-FeaturePlot(object = pbmc.integrated, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
ggsave(paste0("fig",plotid,"e5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)

DefaultAssay(pbmc.integrated) <- "RNA"
pbmc.integrated<-NormalizeData(pbmc.integrated)
all.genes<-rownames(pbmc.integrated)
pbmc.integrated<-ScaleData(pbmc.integrated,features=all.genes)

if(species=='human')
{
pbmc.integrated <- CellCycleScoring(
  object = pbmc.integrated,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
}else{
	library(stringr)
	pbmc.integrated <- CellCycleScoring(
  object = pbmc.integrated,
  g2m.features = gene_human2mouse_vector(cc.genes$g2m.genes),
  s.features = gene_human2mouse_vector(cc.genes$s.genes)
)
	}
plota<-NULL
plota<-VlnPlot(pbmc.integrated, features = c("S.Score","G2M.Score"),group.by="orig.ident")
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
#head(pbmc@meta.data)
#table(pbmc$Phase)
plota<-NULL
plota<-DimPlot(pbmc.integrated,group.by='Phase',split.by = "Phase")
ggsave(paste0("fig",plotid,"f2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"f2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)
DefaultAssay(pbmc.integrated) <- "RNA"
return(pbmc.integrated)
}



seurat_integrate_harmony_v5_notuse <- function(pbmc,group.by='orig.ident',species = "human",
                         mt.pattern="^MT-",mt.list=NULL,dim.use=20,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=2000,
                         res=0.5,plotid='02') {
library(harmony)
Idents(pbmc) <- group.by
DefaultAssay(pbmc) <- "RNA"
pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc[[group.by]])

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nfeatures)
features_variable<-VariableFeatures(object = pbmc)
features_all <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = features_variable)
pbmc <- RunPCA(pbmc, features = features_variable)

plota<-NULL
plota<-VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
ggsave(paste0("fig",plotid,"b.pcaGene.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"b.pcaGene.pdf"), plot = plota, width = 10, height = 10)

plota<-NULL
plota<-ElbowPlot(pbmc, ndims = 50)
ggsave(paste0("fig",plotid,"c.pca50.tiff"), plot = plota, width = 10, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"c.pca50.pdf"), plot = plota, width = 10, height = 6)

#绘制主成分分析图形
plota<-NULL
plota<-DimPlot(object = pbmc, reduction = "pca")
ggsave(paste0("fig",plotid,"d.PCA.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"d.PCA.pdf"), plot = plota, width = 10, height = 10)


pbmc.integrated <- IntegrateLayers(object = pbmc, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",  verbose = FALSE)
#pbmc.integrated <- RunHarmony(pbmc, group.by , plot_convergence = F,dims.use = 1:dim.use)
pbmc.integrated = RunUMAP(pbmc.integrated, reduction = "harmony", dims = 1:dim.use)
pbmc.integrated = FindNeighbors(pbmc.integrated, reduction = "harmony",dims = 1:dim.use)
pbmc.integrated = FindClusters(pbmc.integrated, resolution = res)


plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 6)

meta<-pbmc.integrated@meta.data
if('group' %in% colnames(meta))
{
plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e1.umap.group.pdf"), plot = plota, width = 8, height = 8)

plota<-NULL
plota <- DimPlot(pbmc.integrated, reduction = "umap", split.by = "group")
ggsave(paste0("fig",plotid,"e2.umap.group.split.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e2.umap.group.split.pdf"), plot = plota, width = 12, height = 6)
}

plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e3.umap.sample.pdf"), plot = plota, width = 8, height = 8)

nsample<-length(unique(pbmc.integrated@meta.data$orig.ident))
nh<-ceiling(nsample/4)
plotc<-NULL
plotc <- DimPlot(pbmc.integrated, reduction = "umap", split.by = "orig.ident",ncol=4)
ggsave(paste0("fig",plotid,"e4.umap.sample.split.tiff"), plot = plotc, width = 12, height = nh*3,compression='lzw')
ggsave(paste0("fig",plotid,"e4.umap.sample.split.pdf"), plot = plotc, width = 12, height = nh*3)

plota<-NULL
plota<-FeaturePlot(object = pbmc.integrated, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
ggsave(paste0("fig",plotid,"e5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)

DefaultAssay(pbmc.integrated) <- "RNA"
pbmc.integrated<-NormalizeData(pbmc.integrated)
all.genes<-rownames(pbmc.integrated)
pbmc.integrated<-ScaleData(pbmc.integrated,features=all.genes)

if(species=='human')
{
pbmc.integrated <- CellCycleScoring(
  object = pbmc.integrated,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
}else{
	library(stringr)
	pbmc.integrated <- CellCycleScoring(
  object = pbmc.integrated,
  g2m.features = str_to_title(cc.genes$g2m.genes),
  s.features = str_to_title(cc.genes$s.genes)
)
	}
plota<-NULL
plota<-VlnPlot(pbmc.integrated, features = c("S.Score","G2M.Score"),group.by="orig.ident")
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
#head(pbmc@meta.data)
#table(pbmc$Phase)
plota<-NULL
plota<-DimPlot(pbmc.integrated,group.by='Phase',split.by = "Phase")
ggsave(paste0("fig",plotid,"f2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"f2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)
DefaultAssay(pbmc.integrated) <- "RNA"
return(pbmc.integrated)
}


seurat_integrate_harmony_cellcycle <- function(pbmc,group.by='orig.ident',species = "human",
                         mt.pattern="^MT-",mt.list=NULL,dim.use=20,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=2000,
                         res=0.5,plotid='02') {
library(harmony)
Idents(pbmc) <- group.by
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nfeatures)
features_variable<-VariableFeatures(object = pbmc)
features_all <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = features_variable)
pbmc <- RunPCA(pbmc, features = features_variable)
pbmc<-ScaleData(pbmc,features=features_all)

s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes

if(species=='human')
{
		pbmc <- CellCycleScoring(
		  object = pbmc,
		  g2m.features = g2m.genes,
		  s.features = s.genes
		)
}else{
		library(stringr)
		g2m.genes = gene_human2mouse_vector(g2m.genes)
		s.genes = gene_human2mouse_vector(s.genes)
		pbmc <- CellCycleScoring(
		  object = pbmc,
		  g2m.features = g2m.genes,
		  s.features = s.genes
		)
}
## https://satijalab.org/seurat/articles/cell_cycle_vignette.html
pbmc <- RunPCA(pbmc, features = c(s.genes, g2m.genes))
pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), nfeatures.print = 10)

plota<-NULL
plota<-VlnPlot(pbmc, features = c("S.Score","G2M.Score"),group.by="orig.ident")
ggsave(paste0("fig",plotid,"a1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"a1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
plota<-NULL
plota<-DimPlot(pbmc,group.by='Phase',split.by = "Phase")
ggsave(paste0("fig",plotid,"a2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"a2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)

plota<-NULL
plota<-DimPlot(pbmc,group.by='Phase')
ggsave(paste0("fig",plotid,"a3.umap.Phase.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"a3.umap.Phase.pdf"), plot = plota, width = 8, height = 8)


plota<-NULL
plota<-VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
ggsave(paste0("fig",plotid,"b1.pcaGene.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"b1.pcaGene.pdf"), plot = plota, width = 10, height = 10)

plota<-NULL
plota<-ElbowPlot(pbmc, ndims = 50)
ggsave(paste0("fig",plotid,"c2.pca50.tiff"), plot = plota, width = 10, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"c2.pca50.pdf"), plot = plota, width = 10, height = 6)

#绘制主成分分析图形
plota<-NULL
plota<-DimPlot(object = pbmc, reduction = "pca")
ggsave(paste0("fig",plotid,"d.PCA.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"d.PCA.pdf"), plot = plota, width = 10, height = 10)

pbmc <- RunHarmony(pbmc, group.by , plot_convergence = F,dims.use = 1:dim.use)
pbmc = RunUMAP(pbmc, reduction = "harmony", dims = 1:dim.use)

pbmc = FindNeighbors(pbmc, reduction = "harmony",dims = 1:dim.use)
pbmc = FindClusters(pbmc, resolution = res)

plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 6)

meta<-pbmc@meta.data
if('group' %in% colnames(meta))
{
plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"f1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"f1.umap.group.pdf"), plot = plota, width = 8, height = 8)

plota<-NULL
plota <- DimPlot(pbmc, reduction = "umap", split.by = "group")
ggsave(paste0("fig",plotid,"f2.umap.group.split.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"f2.umap.group.split.pdf"), plot = plota, width = 12, height = 6)
}

plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"f3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"f3.umap.sample.pdf"), plot = plota, width = 8, height = 8)

nsample<-length(unique(pbmc@meta.data$orig.ident))
nh<-ceiling(nsample/4)
plotc<-NULL
plotc <- DimPlot(pbmc, reduction = "umap", split.by = "orig.ident",ncol=4)
ggsave(paste0("fig",plotid,"f4.umap.sample.split.tiff"), plot = plotc, width = 12, height = nh*3,compression='lzw')
ggsave(paste0("fig",plotid,"f4.umap.sample.split.pdf"), plot = plotc, width = 12, height = nh*3)

plota<-NULL
plota<-FeaturePlot(object = pbmc, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
ggsave(paste0("fig",plotid,"f5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"f5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)


plota<-NULL
plota<-VlnPlot(pbmc, features = c("S.Score","G2M.Score"),group.by="orig.ident")
ggsave(paste0("fig",plotid,"g1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"g1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
plota<-NULL
plota<-DimPlot(pbmc,group.by='Phase',split.by = "Phase")
ggsave(paste0("fig",plotid,"g2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"g2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)

plota<-NULL
plota<-DimPlot(pbmc,group.by='Phase')
ggsave(paste0("fig",plotid,"g3.umap.Phase.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"g3.umap.Phase.pdf"), plot = plota, width = 8, height = 8)


DefaultAssay(pbmc) <- "RNA"
return(pbmc)
}


seurat_integrate_CCA <- function(pbmc,group.by='orig.ident',species='human',
                         mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=2000,
                         res=0.5,plotid='02') {
Idents(pbmc) <- group.by
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.list <- SplitObject(pbmc, split.by = group.by)
for (i in 1:length(pbmc.list)){
  pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]], verbose = FALSE)
  pbmc.list[[i]] <- FindVariableFeatures(pbmc.list[[i]], selection.method = "vst",nfeatures = nfeatures, verbose = FALSE)}
features <- SelectIntegrationFeatures(object.list = pbmc.list)
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)
infer<-table(pbmc@meta.data$orig.ident)
cell_num_min<-min(infer)
kweight=100
if(cell_num_min<100)
{
	kweight<-cell_num_min
}
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors,k.weight=kweight)
DefaultAssay(pbmc.integrated) <- "integrated"

pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunPCA(pbmc.integrated, npcs = dim.use, verbose = FALSE)

plota<-NULL
plota<-VizDimLoadings(object = pbmc.integrated, dims = 1:4, reduction = "pca",nfeatures = 20)
ggsave(paste0("fig",plotid,"b.pcaGene.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"b.pcaGene.pdf"), plot = plota, width = 10, height = 10)
plota<-NULL
plota<-ElbowPlot(pbmc.integrated, ndims = 50)
ggsave(paste0("fig",plotid,"c.pca50.tiff"), plot = plota, width = 10, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"c.pca50.pdf"), plot = plota, width = 10, height = 6)

#绘制主成分分析图形
plota<-NULL
plota<-DimPlot(object = pbmc.integrated, reduction = "pca")
ggsave(paste0("fig",plotid,"d.PCA.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"d.PCA.pdf"), plot = plota, width = 10, height = 10)

pbmc.integrated <- RunUMAP(pbmc.integrated, reduction = "pca", dims = 1:dim.use)
pbmc.integrated <- FindNeighbors(pbmc.integrated, reduction = "pca", dims = 1:dim.use)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = res)

plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 8)

meta<-pbmc.integrated@meta.data
if('group' %in% colnames(meta))
{
plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e1.umap.group.pdf"), plot = plota, width = 8, height = 8)

plotc<-NULL
plotc <- DimPlot(pbmc.integrated, reduction = "umap", split.by = "group")
ggsave(paste0("fig",plotid,"e2.umap.group.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e2.umap.group.split.pdf"), plot = plotc, width = 12, height = 6)
}
plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e3.umap.sample.pdf"), plot = plota, width = 8, height = 8)
plotc<-NULL
plotc <- DimPlot(pbmc.integrated, reduction = "umap", split.by = "orig.ident",ncol=4)
ggsave(paste0("fig",plotid,"e4.umap.sample.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e4.umap.sample.split.pdf"), plot = plotc, width = 12, height = 6)
plota<-NULL
plota<-FeaturePlot(object = pbmc.integrated, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
ggsave(paste0("fig",plotid,"e5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)

DefaultAssay(pbmc.integrated) <- "RNA"
pbmc.integrated<-NormalizeData(pbmc.integrated)
all.genes<-rownames(pbmc.integrated)
pbmc.integrated<-ScaleData(pbmc.integrated,features=all.genes)

if(species=='human')
{
pbmc.integrated <- CellCycleScoring(
  object = pbmc.integrated,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
}else{
	library(stringr)
	pbmc.integrated <- CellCycleScoring(
  object = pbmc.integrated,
  g2m.features = gene_human2mouse_vector(cc.genes$g2m.genes),
  s.features = gene_human2mouse_vector(cc.genes$s.genes)
)
}

plota<-NULL
plota<-VlnPlot(pbmc.integrated, features = c("S.Score","G2M.Score"),group.by="orig.ident")
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
#head(pbmc@meta.data)
#table(pbmc$Phase)
plota<-NULL
plota<-DimPlot(pbmc.integrated,group.by='Phase',split.by = "Phase")
ggsave(paste0("fig",plotid,"f2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"f2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)
DefaultAssay(pbmc.integrated) <- "RNA"
return(pbmc.integrated)
}






seurat_integrate_SCTransform <- function(pbmc,group.by='orig.ident',species='human',
                         mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=2000,
                         res=1.5,plotid='02') {
Idents(pbmc) <- group.by
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.list <- SplitObject(pbmc, split.by = group.by)

if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
pbmc.list <- lapply(X = pbmc.list, FUN = SCTransform, vars.to.regress = "percent.mt", verbose = FALSE)
}else{
pbmc.list <- lapply(X = pbmc.list, FUN = SCTransform, vars.to.regress = "percent.mt", verbose = FALSE, method = "glmGamPoi")
}
features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = nfeatures)
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = features)
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT",anchor.features = features)
infer<-table(pbmc@meta.data$orig.ident)
cell_num_min<-min(infer)
kweight=100
if(cell_num_min<100)
{
	kweight<-cell_num_min
}
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT",k.weight=kweight)
pbmc.integrated <- RunPCA(pbmc.integrated, verbose = FALSE)
plota<-NULL
plota<-VizDimLoadings(object = pbmc.integrated, dims = 1:4, reduction = "pca",nfeatures = 20)
ggsave(paste0("fig",plotid,"b.pcaGene.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"b.pcaGene.pdf"), plot = plota, width = 10, height = 10)
plota<-NULL
plota<-ElbowPlot(pbmc.integrated, ndims = 50)
ggsave(paste0("fig",plotid,"c.pca50.tiff"), plot = plota, width = 10, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"c.pca50.pdf"), plot = plota, width = 10, height = 6)

#绘制主成分分析图形
plota<-NULL
plota<-DimPlot(object = pbmc.integrated, reduction = "pca")
ggsave(paste0("fig",plotid,"d.PCA.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"d.PCA.pdf"), plot = plota, width = 10, height = 10)
pbmc.integrated <- RunUMAP(pbmc.integrated, reduction = "pca", dims = 1:dim.use)

pbmc.integrated <- FindNeighbors(pbmc.integrated, reduction = "pca", dims = 1:dim.use)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = res)
pbmc.integrated <- RunUMAP(pbmc.integrated, reduction = "pca", dims = 1:dim.use)
plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 8)

meta<-pbmc.integrated@meta.data
if('group' %in% colnames(meta))
{
plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e1.umap.group.pdf"), plot = plota, width = 8, height = 8)

plotc<-NULL
plotc <- DimPlot(pbmc.integrated, reduction = "umap", split.by = "group")
ggsave(paste0("fig",plotid,"e2.umap.group.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e2.umap.group.split.pdf"), plot = plotc, width = 12, height = 6)
}
plota<-NULL
plota<-DimPlot(pbmc.integrated, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e3.umap.sample.pdf"), plot = plota, width = 8, height = 8)
plotc<-NULL
plotc <- DimPlot(pbmc.integrated, reduction = "umap", split.by = "orig.ident",ncol=4)
ggsave(paste0("fig",plotid,"e4.umap.sample.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e4.umap.sample.split.pdf"), plot = plotc, width = 12, height = 6)
plota<-NULL
plota<-FeaturePlot(object = pbmc.integrated, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
ggsave(paste0("fig",plotid,"e5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)

DefaultAssay(pbmc.integrated) <- "RNA"
pbmc.integrated<-NormalizeData(pbmc.integrated)
all.genes<-rownames(pbmc.integrated)
pbmc.integrated<-ScaleData(pbmc.integrated,features=all.genes)


if(species=='human')
{
pbmc.integrated <- CellCycleScoring(
  object = pbmc.integrated,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
}else{
	library(stringr)
	pbmc.integrated <- CellCycleScoring(
  object = pbmc.integrated,
  g2m.features = gene_human2mouse_vector(cc.genes$g2m.genes),
  s.features = gene_human2mouse_vector(cc.genes$s.genes)
)
}
plota<-NULL
plota<-VlnPlot(pbmc.integrated, features = c("S.Score","G2M.Score"),group.by="orig.ident")
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
#head(pbmc@meta.data)
#table(pbmc$Phase)
plota<-NULL
plota<-DimPlot(pbmc.integrated,group.by='Phase',split.by = "Phase")
ggsave(paste0("fig",plotid,"f2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"f2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)
DefaultAssay(pbmc.integrated) <- "RNA"
return(pbmc.integrated)
}




seurat_integrate_SCTransform_v5 <- function(pbmc,group.by='orig.ident',species='human',
											mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=5,
											nf.low=500,nf.high=6000,nfeatures=2000,
											res=1.5,plotid='02') {
  Idents(pbmc) <- group.by
  DefaultAssay(pbmc) <- "RNA"
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc.list <- SplitObject(pbmc, split.by = group.by)
  Idents(pbmc) <- group.by
  pbmc[['RNA']] <- split(pbmc[['RNA']],f=pbmc@meta.data[[group.by]])
  
  options(future.globals.maxSize = 3e+09)
  pbmc <- SCTransform(pbmc)
  pbmc <- RunPCA(pbmc, npcs = dim.use, verbose = F)
  pbmc <- IntegrateLayers(
	object = pbmc,
	method = CCAIntegration,
	normalization.method = "SCT",
	verbose = F
  )
  pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
  pbmc <- FindNeighbors(pbmc, dims = 1:dim.use, reduction = "integrated.dr")
  pbmc <- FindClusters(pbmc, resolution = res)
  pbmc <- RunUMAP(pbmc, dims = 1:dim.use, reduction = "integrated.dr")
  
  plota<-NULL
  plota<-DimPlot(pbmc, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
  ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 8)
  
  meta<-pbmc@meta.data
  if('group' %in% colnames(meta))
  {
	plota<-NULL
	plota<-DimPlot(pbmc, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
	ggsave(paste0("fig",plotid,"e1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
	ggsave(paste0("fig",plotid,"e1.umap.group.pdf"), plot = plota, width = 8, height = 8)
	
	plotc<-NULL
	plotc <- DimPlot(pbmc, reduction = "umap", split.by = "group")
	ggsave(paste0("fig",plotid,"e2.umap.group.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
	ggsave(paste0("fig",plotid,"e2.umap.group.split.pdf"), plot = plotc, width = 12, height = 6)
  }
  plota<-NULL
  plota<-DimPlot(pbmc, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
  ggsave(paste0("fig",plotid,"e3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,"e3.umap.sample.pdf"), plot = plota, width = 8, height = 8)
  plotc<-NULL
  plotc <- DimPlot(pbmc, reduction = "umap", split.by = "orig.ident",ncol=4)
  ggsave(paste0("fig",plotid,"e4.umap.sample.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"e4.umap.sample.split.pdf"), plot = plotc, width = 12, height = 6)
  plota<-NULL
  plota<-FeaturePlot(object = pbmc, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
  ggsave(paste0("fig",plotid,"e5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"e5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)
  
  DefaultAssay(pbmc) <- "RNA"
  all.genes<-rownames(pbmc)
  pbmc<-ScaleData(pbmc,features=all.genes)
  
  if(species=='human')
  {
	pbmc <- CellCycleScoring(
	  object = pbmc,
	  g2m.features = cc.genes$g2m.genes,
	  s.features = cc.genes$s.genes
	)
	
  }else{
	library(stringr)
	pbmc <- CellCycleScoring(
	  object = pbmc,
	  g2m.features = gene_human2mouse_vector(cc.genes$g2m.genes),
	  s.features = gene_human2mouse_vector(cc.genes$s.genes)
	)
  }
  plota<-NULL
  plota<-VlnPlot(pbmc, features = c("S.Score","G2M.Score"),group.by="orig.ident")
  ggsave(paste0("fig",plotid,"f1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,"f1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
  #head(pbmc@meta.data)
  #table(pbmc$Phase)
  plota<-NULL
  plota<-DimPlot(pbmc,group.by='Phase',split.by = "Phase")
  ggsave(paste0("fig",plotid,"f2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"f2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)
  DefaultAssay(pbmc) <- "RNA"
  return(pbmc)
}


seurat_integrate_SCTransform_v5_small <- function(pbmc,group.by='orig.ident',species='human',
											mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=5,
											nf.low=500,nf.high=6000,nfeatures=2000,
											res=1.5,plotid='02') {
  Idents(pbmc) <- group.by
  DefaultAssay(pbmc) <- "RNA"
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc.list <- SplitObject(pbmc, split.by = group.by)
  Idents(pbmc) <- group.by
  pbmc[['RNA']] <- split(pbmc[['RNA']],f=pbmc@meta.data[[group.by]])
  
  options(future.globals.maxSize = 3e+09)
  pbmc <- SCTransform(pbmc)
  pbmc <- RunPCA(pbmc, npcs = dim.use, verbose = F)
  pbmc <- IntegrateLayers(
	object = pbmc,
	method = CCAIntegration,
	k.weight = 50,
	normalization.method = "SCT",
	verbose = F
  )
  pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
  pbmc <- FindNeighbors(pbmc, dims = 1:dim.use, reduction = "integrated.dr")
  pbmc <- FindClusters(pbmc, resolution = res)
  pbmc <- RunUMAP(pbmc, dims = 1:dim.use, reduction = "integrated.dr")
  
  plota<-NULL
  plota<-DimPlot(pbmc, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
  ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 8)
  
  meta<-pbmc@meta.data
  if('group' %in% colnames(meta))
  {
	plota<-NULL
	plota<-DimPlot(pbmc, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
	ggsave(paste0("fig",plotid,"e1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
	ggsave(paste0("fig",plotid,"e1.umap.group.pdf"), plot = plota, width = 8, height = 8)
	
	plotc<-NULL
	plotc <- DimPlot(pbmc, reduction = "umap", split.by = "group")
	ggsave(paste0("fig",plotid,"e2.umap.group.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
	ggsave(paste0("fig",plotid,"e2.umap.group.split.pdf"), plot = plotc, width = 12, height = 6)
  }
  plota<-NULL
  plota<-DimPlot(pbmc, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
  ggsave(paste0("fig",plotid,"e3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,"e3.umap.sample.pdf"), plot = plota, width = 8, height = 8)
  plotc<-NULL
  plotc <- DimPlot(pbmc, reduction = "umap", split.by = "orig.ident",ncol=4)
  ggsave(paste0("fig",plotid,"e4.umap.sample.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"e4.umap.sample.split.pdf"), plot = plotc, width = 12, height = 6)
  plota<-NULL
  plota<-FeaturePlot(object = pbmc, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
  ggsave(paste0("fig",plotid,"e5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"e5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)
  
  DefaultAssay(pbmc) <- "RNA"
  all.genes<-rownames(pbmc)
  pbmc<-ScaleData(pbmc,features=all.genes)
  
  if(species=='human')
  {
	pbmc <- CellCycleScoring(
	  object = pbmc,
	  g2m.features = cc.genes$g2m.genes,
	  s.features = cc.genes$s.genes
	)
	
  }else{
	library(stringr)
	pbmc <- CellCycleScoring(
	  object = pbmc,
	  g2m.features = gene_human2mouse_vector(cc.genes$g2m.genes),
	  s.features = gene_human2mouse_vector(cc.genes$s.genes)
	)
  }
  plota<-NULL
  plota<-VlnPlot(pbmc, features = c("S.Score","G2M.Score"),group.by="orig.ident")
  ggsave(paste0("fig",plotid,"f1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,"f1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
  #head(pbmc@meta.data)
  #table(pbmc$Phase)
  plota<-NULL
  plota<-DimPlot(pbmc,group.by='Phase',split.by = "Phase")
  ggsave(paste0("fig",plotid,"f2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"f2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)
  DefaultAssay(pbmc) <- "RNA"
  return(pbmc)
}



seurat_integrate_null <- function(pbmc,group.by='orig.ident',
                         mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=2000,
                         res=0.5,plotid='02',species = 'human') {
Idents(pbmc) <- group.by
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = nfeatures)
#输出特征方差图
top10 <- head(x = VariableFeatures(object = pbmc), 10)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plota<-CombinePlots(plots = list(plot1, plot2))
ggsave(paste0("fig",plotid,"a.VariableFeatures.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"a.VariableFeatures.pdf"), plot = plota, width = 12, height = 6)
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc,npcs = 50,pc.genes=VariableFeatures(object = pbmc))

plota<-NULL
plota<-VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
ggsave(paste0("fig",plotid,"b.pcaGene.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"b.pcaGene.pdf"), plot = plota, width = 10, height = 10)
plota<-NULL
plota<-ElbowPlot(pbmc, ndims = 50)
ggsave(paste0("fig",plotid,"c.pca50.tiff"), plot = plota, width = 10, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"c.pca50.pdf"), plot = plota, width = 10, height = 6)

#绘制主成分分析图形
plota<-NULL
plota<-DimPlot(object = pbmc, reduction = "pca")
ggsave(paste0("fig",plotid,"d.PCA.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"d.PCA.pdf"), plot = plota, width = 10, height = 10)

pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:dim.use)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:dim.use)
pbmc <- FindClusters(pbmc, resolution = res)

plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 8)

meta<-pbmc@meta.data
if('group' %in% colnames(meta))
{
plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e1.umap.group.pdf"), plot = plota, width = 8, height = 8)

ngroup<-length(unique(meta$group))

plotc<-NULL
plotc <- DimPlot(pbmc, reduction = "umap", split.by = "group")
ggsave(paste0("fig",plotid,"e2.umap.group.split.tiff"), plot = plotc, width = ngroup*6, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e2.umap.group.split.pdf"), plot = plotc, width = ngroup*6, height = 6)
}
plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e3.umap.sample.pdf"), plot = plota, width = 8, height = 8)

nsample<-length(unique(meta$orig.ident))
nh<-ceiling(nsample/4)

plotc<-NULL
plotc <- DimPlot(pbmc, reduction = "umap", split.by = "orig.ident",ncol=4)
ggsave(paste0("fig",plotid,"e4.umap.sample.split.tiff"), plot = plotc, width = 12, height = nh*3,compression='lzw')
ggsave(paste0("fig",plotid,"e4.umap.sample.split.pdf"), plot = plotc, width = 12, height = nh*3)
plota<-NULL
plota<-FeaturePlot(object = pbmc, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
ggsave(paste0("fig",plotid,"e5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)

if(0)
{
#细胞周期回归：上一步找到的高变基因，常常会包含一些细胞周期相关基因。
#它们会导致细胞聚类发生一定的偏移，即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
?CaseMatch
cc.genes
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA1))
#细胞周期评分
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA1))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA1))
scRNA1 <- CellCycleScoring(object=scRNA1,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAb <- ScaleData(scRNA1, vars.to.regress = c("S.Score", "G2M.Score"))
#     然后重新PCA、UMAP等计算。
}


DefaultAssay(pbmc) <- "RNA"
all.genes<-rownames(pbmc)
pbmc<-ScaleData(pbmc,features=all.genes)

if(species=='human')
{
pbmc <- CellCycleScoring(
  object = pbmc,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
}else{
	library(stringr)
	pbmc <- CellCycleScoring(
  object = pbmc,
  g2m.features = str_to_title(cc.genes$g2m.genes),
  s.features = str_to_title(cc.genes$s.genes)
)
	}
	
plota<-NULL
plota<-VlnPlot(pbmc, features = c("S.Score","G2M.Score"),group.by="orig.ident")
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
#head(pbmc@meta.data)
#table(pbmc$Phase)
plota<-NULL
plota<-DimPlot(pbmc,group.by='Phase',split.by = "Phase")
ggsave(paste0("fig",plotid,"f2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"f2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)
DefaultAssay(pbmc) <- "RNA"
return(pbmc)
}



seurat_integrate_null2 <- function(pbmc,group.by='orig.ident',
                         mt.pattern="^MT-",mt.list=NULL,dim.use=30,mt.cutoff=5,
                         nf.low=500,nf.high=6000,nfeatures=2000,
                         res=0.5,plotid='02',species = 'human') {
Idents(pbmc) <- group.by
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = nfeatures)
#输出特征方差图
top10 <- head(x = VariableFeatures(object = pbmc), 10)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plota<-CombinePlots(plots = list(plot1, plot2))
ggsave(paste0("fig",plotid,"a.VariableFeatures.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"a.VariableFeatures.pdf"), plot = plota, width = 12, height = 6)
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc,npcs = 50,pc.genes=VariableFeatures(object = pbmc))

plota<-NULL
plota<-VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
ggsave(paste0("fig",plotid,"b.pcaGene.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"b.pcaGene.pdf"), plot = plota, width = 10, height = 10)
plota<-NULL
plota<-ElbowPlot(pbmc, ndims = 50)
ggsave(paste0("fig",plotid,"c.pca50.tiff"), plot = plota, width = 10, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"c.pca50.pdf"), plot = plota, width = 10, height = 6)

#绘制主成分分析图形
plota<-NULL
plota<-DimPlot(object = pbmc, reduction = "pca")
ggsave(paste0("fig",plotid,"d.PCA.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
ggsave(paste0("fig",plotid,"d.PCA.pdf"), plot = plota, width = 10, height = 10)

pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:dim.use)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:dim.use)
pbmc <- FindClusters(pbmc, resolution = res)

plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 8)

meta<-pbmc@meta.data
if('group' %in% colnames(meta))
{
plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e1.umap.group.pdf"), plot = plota, width = 8, height = 8)

ngroup<-length(unique(meta$group))

plotc<-NULL
plotc <- DimPlot(pbmc, reduction = "umap", split.by = "group")
ggsave(paste0("fig",plotid,"e2.umap.group.split.tiff"), plot = plotc, width = ngroup*6, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e2.umap.group.split.pdf"), plot = plotc, width = ngroup*6, height = 6)
}
plota<-NULL
plota<-DimPlot(pbmc, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
ggsave(paste0("fig",plotid,"e3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"e3.umap.sample.pdf"), plot = plota, width = 8, height = 8)

nsample<-length(unique(meta$orig.ident))
nh<-ceiling(nsample/4)

plotc<-NULL
plotc <- DimPlot(pbmc, reduction = "umap", split.by = "orig.ident",ncol=4)
ggsave(paste0("fig",plotid,"e4.umap.sample.split.tiff"), plot = plotc, width = 12, height = nh*3,compression='lzw')
ggsave(paste0("fig",plotid,"e4.umap.sample.split.pdf"), plot = plotc, width = 12, height = nh*3)
plota<-NULL
plota<-FeaturePlot(object = pbmc, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
ggsave(paste0("fig",plotid,"e5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"e5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)

if(0)
{
#细胞周期回归：上一步找到的高变基因，常常会包含一些细胞周期相关基因。
#它们会导致细胞聚类发生一定的偏移，即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
?CaseMatch
cc.genes
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA1))
#细胞周期评分
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA1))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA1))
scRNA1 <- CellCycleScoring(object=scRNA1,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAb <- ScaleData(scRNA1, vars.to.regress = c("S.Score", "G2M.Score"))
#     然后重新PCA、UMAP等计算。
}


DefaultAssay(pbmc) <- "RNA"
all.genes<-rownames(pbmc)
pbmc<-ScaleData(pbmc,features=all.genes)

if(species=='human')
{
pbmc <- CellCycleScoring(
  object = pbmc,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
}else{
	library(stringr)
	pbmc <- CellCycleScoring(
  object = pbmc,
  g2m.features = str_to_title(cc.genes$g2m.genes),
  s.features = str_to_title(cc.genes$s.genes)
)
	}
	
plota<-NULL
plota<-VlnPlot(pbmc, features = c("S.Score","G2M.Score"),group.by="orig.ident")
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
ggsave(paste0("fig",plotid,"f1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
#head(pbmc@meta.data)
#table(pbmc$Phase)
plota<-NULL
plota<-DimPlot(pbmc,group.by='Phase',split.by = "Phase")
ggsave(paste0("fig",plotid,"f2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
ggsave(paste0("fig",plotid,"f2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)
DefaultAssay(pbmc) <- "RNA"
return(pbmc)
}


seurat_singler<- function(pbmc,reference='BlueprintEncodeData',plotid='04',label='label.main',species='human',plot.only=F){
		
		library(celldex)
		library(SingleR)

		library(reshape2)
		library(ggpubr)		
		if(!plot.only)
		{
		counts<-LayerData(pbmc, assay = "RNA", layer = "counts")
		clusters<-pbmc@meta.data$seurat_clusters
		ann=pbmc@meta.data$orig.ident

		if(species=='human')
		{
			print(species)
			#ref=celldex::BlueprintEncodeData()
			if(reference == 'BlueprintEncodeData')
			{
			   ref=celldex::BlueprintEncodeData()
			}else if(reference == 'MonacoImmuneData')
			{
			   ref=celldex::MonacoImmuneData()
			}else if(reference == 'HumanPrimaryCellAtlasData')
			{
			   ref=celldex::HumanPrimaryCellAtlasData()
			}else if(reference == 'DatabaseImmuneCellExpressionData')
			{
			   ref=celldex::DatabaseImmuneCellExpressionData()
			}else if(reference == 'NovershternHematopoieticData')
			{
			   ref=celldex::NovershternHematopoieticData()
			}else{
				print("unknown singler reference")
				ref=celldex::BlueprintEncodeData()		
			}
			#ref=celldex::MonacoImmuneData()
			#ref=celldex::HumanPrimaryCellAtlasData()	
			#ref=celldex::DatabaseImmuneCellExpressionData()	
			#ref=celldex::NovershternHematopoieticData()
		}else if(species=='mouse'){
			
			#ref=celldex::MouseRNAseqData()
			if(reference == 'MouseRNAseqData')
			{
			   ref=celldex::MouseRNAseqData()
			}else if(reference == 'ImmGenData')
			{
			   ref=celldex::ImmGenData()
			}else{
				print("unknown singler reference")
				ref=celldex::MouseRNAseqData()		
			}
		}

		singler=SingleR(test=counts, ref =ref,labels=ref[[label]], clusters = clusters)

		tiff(file=paste0("fig",plotid,"a.SingleR.celltype.heatmap.tiff"),width=3000,heigh=2400,units='px',res=300,compression='lzw')
		plotScoreHeatmap(singler)
		dev.off()
		pdf(file=paste0("fig",plotid,"a.SingleR.celltype.heatmap.pdf"),width=10,height=6)
		plotScoreHeatmap(singler)
		dev.off()

		tiff(file=paste0("fig",plotid,"b.SingleR.DeltaDistribution.tiff"),width=3000,heigh=2400,units='px',res=300,compression='lzw')
		plotDeltaDistribution(singler,ncol=3)
		dev.off()
		pdf(file=paste0("fig",plotid,"b.SingleR.DeltaDistribution.pdf"),width=10,height=6)
		plotDeltaDistribution(singler,ncol=3)
		dev.off()


		clusterAnn=as.data.frame(singler)
		clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
		clusterAnn=clusterAnn[,c("id", "labels")]
		write.table(clusterAnn,file=paste0("fig",plotid,".clusterAnn.txt"),quote=F,sep="\t", row.names=F)
		if(0)
		{
		singler2=SingleR(test=counts, ref =ref,labels=ref$label.main)
		tiff(file=paste0("fig",plotid,"b.UMAP.annotation.heatmap2.tiff"),width=3000,heigh=2400,units='px',res=300,compression='lzw')
		plotScoreHeatmap(singler2)
		dev.off()
		pdf(file=paste0("fig",plotid,"b.UMAP.annotation.heatmap2.pdf"),width=10,height=6)
		plotScoreHeatmap(singler2)
		dev.off()
		cellAnn=as.data.frame(singler2)
		cellAnn=cbind(id=row.names(cellAnn), cellAnn)
		cellAnn=cellAnn[,c("id", "labels")]
		write.table(cellAnn, file=paste0("fig",plotid,".cellAnn.txt"), quote=F, sep="\t", row.names=F)
		}

		#cluster注释后的可视化
		#newLabels=singler$labels
		#names(newLabels)=levels(pbmc)
		#pbmc=RenameIdents(pbmc, newLabels)

		pbmc@meta.data$celltype_singleR<- 'unknown'
		for(i in 1:nrow(clusterAnn)){
		  pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == clusterAnn[i,1]),'celltype_singleR'] <- clusterAnn[i,2]
		  }
		  
		}
		head(pbmc@meta.data)
		plota<-NULL
		plota<-DimPlot(pbmc, reduction = "umap", group.by = 'celltype_singleR',pt.size = 0.5,label=T)
		ggsave(paste0("fig",plotid,"c.umap.celltype.tiff"), plot = plota, width = 8, height = 6,compression='lzw')
		ggsave(paste0("fig",plotid,"c.umap.celltype.pdf"), plot = plota, width = 8, height = 6)


		metadata <- pbmc@meta.data
		df_celltype_count<-as.data.frame(table(metadata$celltype_singleR))
		colnames(df_celltype_count)<-c('celltype_singleR','Ncells')
		write.table(df_celltype_count, file=paste0("fig",plotid,"c.celltype.Ncells.tsv"), quote=F, sep="\t", row.names=F)
		plot_cellcount_barplot(metadata,plotid=plotid,x='celltype_singleR')
		#plot_celltype_proportion(metadata,plotid=plotid,x='group',y='celltype')
		return(pbmc)
}

seurat_hdWGCNA<-function(pbmc,celltype_column='cell_type',sample_column='Sample',celltype_target='INH',plotid=20)
{
  if(1)  #### 载入hdWGCNA包
  {
    # single-cell analysis package
    library(Seurat)
    
    # plotting and data science packages
    library(tidyverse)
    library(cowplot)
    library(patchwork)
    
    # co-expression network analysis packages:
    library(WGCNA)
    library(hdWGCNA)
    package.version('hdWGCNA')
    
    # using the cowplot theme for ggplot
    theme_set(theme_cowplot())
    
    # set random seed for reproducibility
    set.seed(12345)
    
    # optionally enable multithreading
    enableWGCNAThreads(nThreads = 8)
    
    # load the Zhou et al snRNA-seq dataset
    # pbmc <- readRDS('Zhou_2020.rds')
  }
  
pbmc_wgcna <- SetupForWGCNA(
  pbmc,
  gene_select = "fraction", 
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

# pbmc_wgcna<-MetacellsByGroups(seurat_obj = pbmc_wgcna,group.by = c(celltype_column, sample_column),k = 25,max_shared = 10,  ident.group = celltype_column)
pbmc_wgcna<-MetacellsByGroups(seurat_obj = pbmc_wgcna,group.by = c(celltype_column, sample_column),reduction = 'harmony',k = 25,max_shared = 10,  ident.group = celltype_column)
# normalize metacell expression matrix:
pbmc_wgcna <- NormalizeMetacells(pbmc_wgcna)
### metacell_obj <- GetMetacellObject(pbmc_wgcna)
pbmc_wgcna <- SetDatExpr(
  pbmc_wgcna,
  group_name = celltype_target, # the name of the group of interest in the group.by column
  group.by=celltype_column, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
pbmc_wgcna <- TestSoftPowers(
  pbmc_wgcna,
  networkType = 'signed' # you can also use "signed" or "signed hybrid"
)

plot_list <- PlotSoftPowers(pbmc_wgcna)
plota<-NULL
plota<-wrap_plots(plot_list, ncol=2)
ggsave(paste0("fig",plotid,"a.PlotSoftPowers",".tiff"), plot = plota, width = 10, height = 8,compression='lzw',limitsize=F)
ggsave(paste0("fig",plotid,"a.PlotSoftPowers",".pdf"), plot = plota, width = 10, height = 8,limitsize=F)

power_table <- GetPowerTable(pbmc_wgcna)
head(power_table)

# 如果没有指定软阈值，construcNetwork会自动指定软阈值
# construct co-expression network:
pbmc_wgcna <- ConstructNetwork(
  pbmc_wgcna,
  tom_name = celltype_target, # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE # 允许覆盖已存在的同名文件
)

pdf(file=paste0("fig",plotid,"b.hdWGCNA_Dendrogram.pdf"),width=10,height=8,onefile=F)
print(PlotDendrogram(pbmc_wgcna, main='hdWGCNA Dendrogram'))
dev.off()
tiff(file=paste0("fig",plotid,"b.hdWGCNA_Dendrogram.tiff"),width=10*300,height=8*300,compression='lzw',res=300)
print(PlotDendrogram(pbmc_wgcna, main='hdWGCNA Dendrogram'))
dev.off()

pbmc_wgcna@misc$tutorial$wgcna_modules %>% head


# 可选：检查topoligcal重叠矩阵(TOM)
# TOM <- GetTOM(pbmc_wgcna)
# TOM

# need to run ScaleData first or else harmony throws an error:
# pbmc_wgcna <- ScaleData(pbmc_wgcna, features=VariableFeatures(pbmc_wgcna))

# compute all MEs in the full single-cell dataset
pbmc_wgcna <- ModuleEigengenes(
  pbmc_wgcna,
  group.by.vars=sample_column
)

unique(pbmc_wgcna@meta.data$orig.ident)
unique(pbmc_wgcna@meta.data$Sample)

# harmonized module eigengenes:
# allow the user to apply Harmony batch correction to the MEs, yielding harmonized module eigengenes (hMEs)
hMEs <- GetMEs(pbmc_wgcna)
n_me<-dim(hMEs)[2]-1
print(paste0('get MEs: ',n_me))
# module eigengenes:
#MEs <- GetMEs(seurat_obj, harmonized=FALSE)
# compute eigengene-based connectivity (kME):
# focus on the “hub genes”
pbmc_wgcna <- ModuleConnectivity(
  pbmc_wgcna,
  group.by = celltype_column, 
  group_name = celltype_target
)

# rename the modules
pbmc_wgcna <- ResetModuleNames(
  pbmc_wgcna,
  new_name = paste0(celltype_target,"-M")
)
pbmc_wgcna@misc$tutorial$wgcna_modules[,c(-1,-2,-3)] %>% head

PlotKMEs<-function (seurat_obj, n_hubs = 10, text_size = 2, ncol = 5, plot_widths = c(3, 2), wgcna_name = NULL) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 
                                                             "grey")
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  mod_colors <- modules %>% subset(module %in% mods) %>% dplyr::select(c(module, 
                                                                         color)) %>% dplyr::distinct()
  hub_df <- GetHubGenes(seurat_obj, n_hubs=n_hubs, wgcna_name=wgcna_name)
  plot_list <- lapply(mods, function(x) {
    print(x)
    cur_color <- subset(mod_colors, module == x) %>% .$color
    cur_df <- subset(hub_df, module == x)
    top_genes <- cur_df %>% dplyr::top_n(n_hubs, wt = kME) %>% 
      .$gene_name
    p <- cur_df %>% ggplot(aes(x = reorder(gene_name, kME), 
                               y = kME)) + geom_bar(stat = "identity", width = 1, 
                                                    color = cur_color, fill = cur_color) + ggtitle(x) + 
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
            plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), 
            axis.line.x = element_blank())
    p_anno <- ggplot() + annotate("label", x = 0, y = 0, 
                                  label = paste0(top_genes, collapse = "\n"), size = text_size, 
                                  fontface = "italic", label.size = 0) + theme_void()
    patch <- p + p_anno + plot_layout(widths = plot_widths)
    patch
  })
  wrap_plots(plot_list, ncol = ncol)
}

# plot genes ranked by kME for each module
nr<-ceiling(n_me/5)

plota<-NULL
plota <- PlotKMEs(pbmc_wgcna, n_hubs = 10,ncol=5,wgcna_name = "tutorial")
ggsave(paste0("fig",plotid,"c.PlotKMEs",".tiff"), plot = plota, width = 25, height = nr*4,compression='lzw',limitsize=F)
ggsave(paste0("fig",plotid,"c.PlotKMEs",".pdf"), plot = plota, width = 25, height = nr*4,limitsize=F)



# get the module assignment table:
modules <- GetModules(pbmc_wgcna) %>% 
  subset(module != 'grey')

# show the first 6 columns:
head(modules)
dim(modules)
#          gene_name            module   color     kME_grey kME_CD4+ T-cells_NEW1 kME_CD4+ T-cells_NEW2
# HES4          HES4 CD4+ T-cells_NEW1     red 0.0008704183             0.1680883            0.05642900
# ISG15        ISG15 CD4+ T-cells_NEW1     red 0.1094323715             0.2194204            0.20275178
# TNFRSF18  TNFRSF18 CD4+ T-cells_NEW2   brown 0.2613211844             0.1958609            0.31299310
# TNFRSF4    TNFRSF4 CD4+ T-cells_NEW3 magenta 0.2009701385             0.1884660            0.27907966
# ACAP3        ACAP3 CD4+ T-cells_NEW4    blue 0.0086625773            -0.0524401            0.03869565
# AURKAIP1  AURKAIP1 CD4+ T-cells_NEW2   brown 0.2577977323             0.1609058            0.32569947

# get hub genes
hub_df <- GetHubGenes(pbmc_wgcna, n_hubs = 10)
head(hub_df)
#   gene_name            module       kME
# 1      SRGN CD4+ T-cells_NEW1 0.6773476
# 2      CREM CD4+ T-cells_NEW1 0.6354395
# 3     HSPD1 CD4+ T-cells_NEW1 0.5762528
# 4       UBC CD4+ T-cells_NEW1 0.5738136
# 5  HSP90AA1 CD4+ T-cells_NEW1 0.5646901
# 6     HSPH1 CD4+ T-cells_NEW1 0.5626037

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
pbmc_wgcna <- ModuleExprScore(
  pbmc_wgcna,
  n_genes = 25,
  method='UCell' # Seurat方法(AddModuleScore)
)


# 制作每个模块的hMEs特征图
plot_list <- ModuleFeaturePlot(
  pbmc_wgcna,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
plota<-NULL
plota<-wrap_plots(plot_list, ncol=5)
ggsave(paste0("fig",plotid,"d.hMEs",".tiff"), plot = plota, width = 25, height = nr*4,compression='lzw',limitsize=F)
ggsave(paste0("fig",plotid,"d.hMEs",".pdf"), plot = plota, width = 25, height = nr*4,limitsize=F)

# 制作每个模块的hub scores特征图
plot_list <- ModuleFeaturePlot(
  pbmc_wgcna,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
# stitch together with patchwork
plota<-NULL
plota<-wrap_plots(plot_list, ncol=5)
ggsave(paste0("fig",plotid,"e.hubscorePlot",".tiff"), plot = plota, width = 25, height = nr*4,compression='lzw',limitsize=F)
ggsave(paste0("fig",plotid,"e.hubscorePlot",".pdf"), plot = plota, width = 25, height = nr*4,limitsize=F)

if(1)
{
  # devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
  library(ggradar)
  
    # 每个模块在不同细胞亚群中的情况
    head(pbmc_wgcna@meta.data)
    plota<-ModuleRadarPlot(
      pbmc_wgcna,
      group.by = sample_column,
      barcodes = rownames(pbmc_wgcna@meta.data[pbmc_wgcna@meta.data[[celltype_column]]==celltype_target,]),
      axis.label.size=4,
      grid.label.size=4
    )
    ggsave(paste0("fig",plotid,"e2.ModuleRadarPlot",".tiff"), plot = plota, width = 25, height = nr*4,compression='lzw',limitsize=F)
    ggsave(paste0("fig",plotid,"e2.ModuleRadarPlot",".pdf"), plot = plota, width = 25, height = nr*4,limitsize=F)

}

# 查看模块相关图
pdf(file=paste0("fig",plotid,"f.ModuleCorrelogram.pdf"),width=10,height=10,onefile=F)
print(ModuleCorrelogram(pbmc_wgcna))
dev.off()
tiff(file=paste0("fig",plotid,"f.ModuleCorrelogram.tiff"),width=10*300,height=10*300,compression='lzw',res=300)
print(ModuleCorrelogram(pbmc_wgcna))
dev.off()





# get hMEs from seurat object
MEs <- GetMEs(pbmc_wgcna, harmonized=TRUE)
dim(MEs)
modules <- GetModules(pbmc_wgcna)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
pbmc_wgcna@meta.data <- cbind(pbmc_wgcna@meta.data, MEs)
# plot with Seurat's DotPlot function
ncelltype<-length(unique(pbmc_wgcna@meta.data[[celltype_column]]))
ph<-max(3,ncelltype)
pw<-max(6,n_me)
plota<-NULL
plota <- DotPlot(pbmc_wgcna, features=mods, group.by = celltype_column)
plota <- plota +RotatedAxis() +scale_color_gradient2(high='red', mid='grey95', low='blue')
ggsave(paste0("fig",plotid,"g.mods_DotPlot",".tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
ggsave(paste0("fig",plotid,"g.mods_DotPlot",".pdf"), plot = plota, width = pw, height = ph,limitsize=F)

plota<-NULL
plota<-VlnPlot(object = pbmc_wgcna, features = mods,stack=T,flip=T,pt.size = 0, group.by = celltype_column)+ theme(legend.position = "none",axis.text.x =element_text(angle=45,hjust=1))
ggsave(paste0("fig",plotid,"g2.mods_Violin",".tiff"), plot = plota, width = ph, height = pw,compression='lzw',limitsize=F)
ggsave(paste0("fig",plotid,"g2.mods_Violin",".pdf"), plot = plota, width = ph, height = pw,limitsize=F)



# 使用ModuleNetworkPlot可视化每个模块前50(数值可自定)的hub gene
ModuleNetworkPlot(
  pbmc_wgcna, 
  outdir='ModuleNetworks', # new folder name
  n_inner = 20, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)



# hubgene network(基因数可自定)


pdf(file=paste0("fig",plotid,"h.HubGeneNetworkPlot.pdf"),width=6,height=6,onefile=F)
print(HubGeneNetworkPlot(pbmc_wgcna,n_hubs = 2,n_other=2,edge_prop = 0.75,mods = 'all'))
dev.off()
tiff(file=paste0("fig",plotid,"h.HubGeneNetworkPlot.tiff"),width=6*300,height=6*300,compression='lzw',res=300)
print(HubGeneNetworkPlot(pbmc_wgcna,n_hubs = 2,n_other=2,edge_prop = 0.75,mods = 'all'))
dev.off()




pbmc_wgcna <- RunModuleUMAP(
  pbmc_wgcna,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(pbmc_wgcna)

# plot with ggplot
plota<-NULL
plota<-ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
ggsave(paste0("fig",plotid,"i.ModuleUMAP",".tiff"), plot = plota, width = 10, height = 10,compression='lzw',limitsize=F)
ggsave(paste0("fig",plotid,"i.ModuleUMAP",".pdf"), plot = plota, width = 10, height = 10,limitsize=F)

pdf(file=paste0("fig",plotid,"j.ModuleUMAPPlot.pdf"),width=10,height=10,onefile=F)
print(ModuleUMAPPlot(pbmc_wgcna,edge.alpha=0.25,sample_edges=TRUE,edge_prop=0.1,label_hubs=2 ,keep_grey_edges=FALSE))
dev.off()
tiff(file=paste0("fig",plotid,"j.ModuleUMAPPlot.tiff"),width=10*300,height=10*300,compression='lzw',res=300)
print(ModuleUMAPPlot(pbmc_wgcna,edge.alpha=0.25,sample_edges=TRUE,edge_prop=0.1,label_hubs=2 ,keep_grey_edges=FALSE))
dev.off()

# 保存数据
library(qs)
qsave(pbmc_wgcna, 'hdWGCNA_Zhou_2020_control.qs')


}

### infercnv_obj<-seurat_infercnv_step1(pbmc,group.by='celltype_manual',ref_cells=c('Myeloid_cells','T_cells','B_cells','Plasma'),output='infercnv')
seurat_infercnv_step1_cnv<- function(pbmc_infercnv,group.by='seurat_clusters',ref_cells=c('immnue cells'),file_gene="GRCh38_ensembl_v107_gene_coor.txt",output=paste0('infercnv_',projectid,'_20240314')){
		library(infercnv)
		if(!dir.exists(output))
		{
			dir.create(output)
		}		
		counts<-as.matrix(GetAssayData(object = pbmc, layer = "counts"))
		if(1)
		{
			metadata<-pbmc_infercnv@meta.data
			metadata_out<-cbind(rownames(metadata),metadata[,group.by])
			colnames(metadata_out)<-c('cell_id','celltype')
			write.table(metadata_out,file=paste0(output,"/cell_annotation_infercnv.tsv"),sep="\t",row.names=F,quote=F,col.names=F)
		}
		table(metadata$celltype_manual)
		infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
								annotations_file=paste0(output,"/cell_annotation_infercnv.tsv"),
								delim="\t",
								gene_order_file=file_gene,
								ref_group_names=ref_cells,
								)
		infercnv_obj = infercnv::run(infercnv_obj,
							 cutoff=0.1, 
							 out_dir=output,
							 cluster_by_groups=T, 
							 analysis_mode="subclusters",
							 hclust_method="ward.D2",
							 denoise=T,
							 HMM=F,
							 num_threads=10,
							 plot_steps=F,
							 png_res=150
							 )
		return(infercnv_obj)
}


#infercnv_obj = readRDS("./infercnv_try2/run.final.infercnv_obj")
#CNVscore<-seurat_infercnv_step2_CNVscore(infercnv_obj,infercnv_dir='infercnv_try2',plotid=plotid)
seurat_infercnv_step2_CNVscore<-function(infercnv_obj,infercnv_dir='infercnv_try2',file_gene="GRCh38_ensembl_v107_gene_coor.txt",plotid=1)
{
  expr <- infercnv_obj@expr.data
  dim(expr)
  
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc_vector <- as.vector(unlist(normal_loc))
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  test_loc_vector <- as.vector(unlist(test_loc))
  
  anno_normal<-vector()
  for(i in 1:length(normal_loc))
  {
    namei<-names(normal_loc)[i]
    print(length(normal_loc[[i]]))
    anno_normal<-c(anno_normal,rep(namei,length(normal_loc[[i]])))
  }
  head(anno_normal)
  
  anno_test<-vector()
  for(i in 1:length(test_loc))
  {
    namei<-names(test_loc)[i]
    print(length(test_loc[[i]]))
    anno_test<-c(anno_test,rep(namei,length(test_loc[[i]])))
  }
  head(anno_test)
  
  anno.df=data.frame(
    CB=c(colnames(expr)[normal_loc_vector],colnames(expr)[test_loc_vector]),
    class=c(anno_normal,anno_test)
  )
  head(anno.df)
  
  gn <- rownames(expr)
  geneFile <- read.table(file_gene,header = F,sep = "\t",stringsAsFactors = F)
  rownames(geneFile)=geneFile$V1
  sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
  expr=expr[intersect(gn,geneFile$V1),]
  head(sub_geneFile,4)
  expr[1:4,1:4]
  dim(expr)
  
  expr2=expr-1
  expr2=expr2 ^ 2
  CNV_score=as.data.frame(colMeans(expr2))
  colnames(CNV_score)="CNV_score"
  head(CNV_score)
  write.table(CNV_score,file=paste0(infercnv_dir,"/CNV_score.tsv"),quote=F,sep="\t", row.names=T)
  
  df_annotation<-read.delim(paste0(infercnv_dir,"/cell_annotation_infercnv.tsv"),header=F)
  colnames(df_annotation)<-c('cell_id','celltype')
  head(df_annotation)
  dim(df_annotation)
  
  df_subcluster<-read.delim(paste0(infercnv_dir,'/infercnv_subclusters.observation_groupings.txt'),sep=' ')
  tail(df_subcluster)
  dim(df_subcluster)
  
  df_plot<-merge(CNV_score,df_annotation,by.x='row.names',by.y='cell_id')
  head(df_plot)
  dim(df_plot)
  rownames(df_plot)<-df_plot$Row.names
  df_plot$Row.names<-NULL
  ggplot2_violin(df_plot,score='CNV_score',group='celltype',output=paste0("fig",plotid,'a.violin.celltype'))
  
  df_plot2<-merge(df_plot,df_subcluster[,'Annotation.Group',drop=F],by=0,all.x=T)
  head(df_plot2)
  dim(df_plot2)
  df_plot2[is.na(df_plot2$Annotation.Group),'Annotation.Group']<-0
  
  ggplot2_violin(df_plot2,score='CNV_score',group='Annotation.Group',output=paste0("fig",plotid,'b.violin.subcluster'))
  return(CNV_score)
}

seurat_infercnv_step2_CNVscore_kmeans<-function(infercnv_obj,k=10,infercnv_dir='infercnv_try2',file_gene="GRCh38_ensembl_v107_gene_coor.txt",plotid=1)
{
	library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library("RColorBrewer")

	expr <- infercnv_obj@expr.data
	dim(expr)

	normal_loc <- infercnv_obj@reference_grouped_cell_indices
	normal_loc_vector <- as.vector(unlist(normal_loc))
	test_loc <- infercnv_obj@observation_grouped_cell_indices
	test_loc_vector <- as.vector(unlist(test_loc))

	anno_normal<-vector()
	for(i in 1:length(normal_loc))
	{
	namei<-names(normal_loc)[i]
	print(length(normal_loc[[i]]))
	anno_normal<-c(anno_normal,rep(namei,length(normal_loc[[i]])))
	}
	head(anno_normal)

	anno_test<-vector()
	for(i in 1:length(test_loc))
	{
	namei<-names(test_loc)[i]
	print(length(test_loc[[i]]))
	anno_test<-c(anno_test,rep(namei,length(test_loc[[i]])))
	}
	head(anno_test)

	anno.df=data.frame(
	CB=c(colnames(expr)[normal_loc_vector],colnames(expr)[test_loc_vector]),
	class=c(anno_normal,anno_test)
	)
	head(anno.df)

	gn <- rownames(expr)
	geneFile <- read.table(file_gene,header = F,sep = "\t",stringsAsFactors = F)
	rownames(geneFile)=geneFile$V1
	sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
	expr=expr[intersect(gn,geneFile$V1),]
	head(sub_geneFile,4)
	expr[1:4,1:4]
	dim(expr)
	chroms<-c(1:22,'MT','X','Y')
	chroms<-factor(chroms,levels=chroms)
	chroms_used<-chroms[chroms %in% geneFile$V2]


	kmeans.result <- kmeans(t(expr), k)
	kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
	kmeans_df$CB=rownames(kmeans_df)
	kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
	kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
	rownames(kmeans_df_s)=kmeans_df_s$CB
	kmeans_df_s$CB=NULL
	kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
	head(kmeans_df_s)
	classes<-unique(kmeans_df_s$class)
	col_class<-rainbow(length(classes))
	names(col_class)<-classes
			if(0)
			{
				#定义热图的注释，及配色
				top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = chroms_used,labels_gp = gpar(cex = 1.5)))
				#color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:k] #类别数
				color_v <- turbo(k)
				names(color_v)=as.character(1:k)
				left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=col_class,kmeans_class=color_v))
				
				
				#left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("red","blue"),kmeans_class=color_v))

				#下面是绘图
				pdf(paste0(infercnv_dir,'/fig0.infercnv_kmeans.pdf'),width = 15,height = 10)
				ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
								col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
								cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
								column_split = factor(sub_geneFile$V2, paste("chr",chroms_used,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
								column_gap = unit(2, "mm"),

								heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),

								top_annotation = top_anno,left_annotation = left_anno, #添加注释
								row_title = NULL,column_title = NULL)
				draw(ht, heatmap_legend_side = "right")
				dev.off()
			}
	#每一类对应的CB保存在kmeans_df_s数据框中
	write.table(kmeans_df_s, file = paste0(infercnv_dir,'/fig0.kmeans_df_s.txt'), quote = FALSE, sep = '\t', row.names = T, col.names = T)

	expr2=expr-1
	expr2=expr2 ^ 2
	CNV_score=as.data.frame(colMeans(expr2))
	colnames(CNV_score)="CNV_score"
	head(CNV_score)
	write.table(CNV_score,file=paste0(infercnv_dir,"/CNV_score.tsv"),quote=F,sep="\t", row.names=T)

	df_annotation<-read.delim(paste0(infercnv_dir,"/cell_annotation_infercnv.tsv"),header=F)
	colnames(df_annotation)<-c('cell_id','celltype')
	head(df_annotation)
	dim(df_annotation)


	df_plot<-merge(CNV_score,df_annotation,by.x='row.names',by.y='cell_id')
	head(df_plot)
	dim(df_plot)
	rownames(df_plot)<-df_plot$Row.names
	df_plot$Row.names<-NULL
	ggplot2_violin(df_plot,score='CNV_score',group='celltype',output=paste0("fig",plotid,'a.violin.celltype'))
	write.table(df_plot,file=paste0("fig",plotid,'a.violin.celltype.tsv'),quote=F,sep="\t", row.names=T)


	df_plot2<-merge(df_plot,kmeans_df_s[,'kmeans_class',drop=F],by=0,all.x=T)
	head(df_plot2)
	dim(df_plot2)
	df_plot2[is.na(df_plot2$kmeans_class),'kmeans_class']<-0

	#df_plot2_mean<-aggregate(df_plot2[,'CNV_score'], by = list(df_plot2[,'kmeans_class']), FUN = mean)
	#head(df_plot2_mean)
	#colnames(df_plot2_mean)<-c('kmeans_class','CNV_score_mean')


	plota<-ggplot2_violin(df_plot2,score='CNV_score',group='kmeans_class',output=paste0("fig",plotid,'b.violin.kmeans_class'))
	write.table(df_plot2,file=paste0("fig",plotid,'b.violin.kmeans_class.tsv'),quote=F,sep="\t", row.names=T)
	return(CNV_score)
}


#pbmc<-seurat_infercnv_step3_malignant(pbmc,malignant_clusters=c(1:6,8:9,11:18),method='subclusters',infercnv_dir='infercnv_try2',plotid=1)
seurat_infercnv_step3_malignant<-function(pbmc,malignant_clusters=c(1),method='subclusters',infercnv_dir='infercnv_try2',plotid=1)
{
			if(method == 'subclusters')
			{
				print(method)
				df_subcluster<-read.delim(paste0(infercnv_dir,'/infercnv_subclusters.observation_groupings.txt'),sep=' ')
				tail(df_subcluster)
				dim(df_subcluster)
				malignant_cells<-rownames(df_subcluster[df_subcluster$Annotation.Group %in% malignant_clusters,])
			}else if(method == 'kmeans')
			{
				print(method)
				kmeans_df_s<-read.delim(paste0(infercnv_dir,'/fig0.kmeans_df_s.txt'))
				tail(kmeans_df_s)
				table(kmeans_df_s$kmeans_class)
				dim(kmeans_df_s)
				malignant_cells<-rownames(kmeans_df_s[kmeans_df_s$kmeans_class %in% malignant_clusters,])
			}

			print(length(malignant_cells))
			head(pbmc@meta.data)
			pbmc@meta.data$malignant_type<-'Normal'
			pbmc@meta.data[rownames(pbmc@meta.data) %in% malignant_cells,'malignant_type']<-'Malignant'
			table(pbmc@meta.data$malignant_type)
			
			df_CNV_score<-read.delim(paste0(infercnv_dir,'/CNV_score.tsv'))
			head(df_CNV_score)
			pbmc@meta.data$CNV_score<-0
			for(i in 1:nrow(df_CNV_score))
			{
				cell<-rownames(df_CNV_score)[i]
				pbmc@meta.data[cell,'CNV_score']<-df_CNV_score[i,'CNV_score']
			}
			dim(pbmc@meta.data[pbmc@meta.data$CNV_score>0,])
			#pbmc@meta.data$CNV_score_used<-exp(pbmc@meta.data$CNV_score)
			pbmc@meta.data$CNV_score_used<-10^(pbmc@meta.data$CNV_score*10)
			rescale <- function(x) (x-min(x))/(max(x) - min(x)) * 100
			#ddd <- rescale(aaa)
			pbmc@meta.data$CNV_score_used<-rescale(pbmc@meta.data$CNV_score_used)
			pbmc@meta.data$CNV_score_used1<-pbmc@meta.data$CNV_score_used
			pbmc@meta.data[!pbmc@meta.data$celltype_manual == 'Epithelial_cells','CNV_score_used1']<-0
			
			CNV_score_used<-pbmc@meta.data[pbmc@meta.data$CNV_score>0,'CNV_score_used']
			summary(CNV_score_used)
			
			gene='CNV_score_used'
			plota<-NULL
			plota<-FeaturePlot(object = pbmc, features = gene) & scale_colour_gradientn(colours = c("red",'red',"gray"),values = c(1.0,0.2,0.1,0))
			ggsave(paste0("fig",plotid,"c.FeaturePlot.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,"c.FeaturePlot.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)
			
			gene='CNV_score_used1'
			plota<-NULL
			plota<-FeaturePlot(object = pbmc, features = gene) & scale_colour_gradientn(colours = c("red",'red',"gray"),values = c(1.0,0.2,0.1,0))
			ggsave(paste0("fig",plotid,"d.FeaturePlot.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,"d.FeaturePlot.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)

			plota<-NULL
			plota<-DimPlot(pbmc, reduction = "umap", group.by = 'malignant_type',pt.size = 0.5,label=T)
			ggsave(paste0("fig",plotid,"e.umap.malignant_type.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
			ggsave(paste0("fig",plotid,"e.umap.malignant_type.pdf"), plot = plota, width = 8, height = 8)

			plotc<-NULL
			plotc <- DimPlot(pbmc, reduction = "umap", split.by = "malignant_type")
			ggsave(paste0("fig",plotid,"f.umap.malignant_type.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
			ggsave(paste0("fig",plotid,"f.umap.malignant_type.split.pdf"), plot = plotc, width = 12, height = 6)
			return(pbmc)
}



#infercnv_obj = readRDS(paste0(infercnv_dir,"/run.final.infercnv_obj"))
#pbmc<-seurat_infercnv_step2_CNVscore_used(infercnv_obj,pbmc,CNVscore_cutoff=0.0016,cor_cutoff=0.57,infercnv_dir=infercnv_dir,plotid=plotid)
####  阈值的选择，要么就使用默认的0.001和0.5，要么就根据直方图自己调整一下。
seurat_infercnv_step2_CNVscore_used<-function(infercnv_obj,pbmc,CNVscore_cutoff=0.001,cor_cutoff=0.5,infercnv_dir='infercnv_try2',file_gene="GRCh38_ensembl_v107_gene_coor.txt",plotid=1)
{
	cnvScore <- function(data){
    data <- data %>% as.matrix() %>%
        t() %>% 
        scale() %>% 
        rescale(to=c(-1, 1)) %>% 
        t()
    cnv_score <- as.data.frame(colSums(data * data))
    return(cnv_score)
	}
	estimateCNV <- function(obs, ref, score_threshold, cor_threshold){
    cnv_obs <- colMeans((obs-1)^2)
    cnv_ref <- colMeans((ref-1)^2) 

    cell_top <- names(sort(cnv_obs, decreasing=T))[1:length(cnv_obs)*0.05]
    cnv_top <- rowMeans(obs[, cell_top])

    cor_obs <- apply(obs, 2, function(x)cor(x, cnv_top))
    cor_ref <- apply(ref, 2, function(x)cor(x, cnv_top))

    cnv <- data.frame(score=c(cnv_obs, cnv_ref), cor=c(cor_obs, cor_ref), barcode=c(colnames(obs), colnames(ref)))

    cnv$type <- 'Other'
    cnv$type[cnv$score>score_threshold & cnv$cor>cor_threshold] <- 'Malignant'
    cnv$type[cnv$score<score_threshold & cnv$cor<cor_threshold] <- 'Normal'
    return(cnv)
	}

  expr <- infercnv_obj@expr.data
  dim(expr)  
  
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc_vector <- as.vector(unlist(normal_loc))
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  test_loc_vector <- as.vector(unlist(test_loc))
  
  anno_normal<-vector()
  for(i in 1:length(normal_loc))
  {
    namei<-names(normal_loc)[i]
    print(length(normal_loc[[i]]))
    anno_normal<-c(anno_normal,rep(namei,length(normal_loc[[i]])))
  }
  head(anno_normal)
  
  anno_test<-vector()
  for(i in 1:length(test_loc))
  {
    namei<-names(test_loc)[i]
    print(length(test_loc[[i]]))
    anno_test<-c(anno_test,rep(namei,length(test_loc[[i]])))
  }
  head(anno_test)
  
  gn <- rownames(expr)
  geneFile <- read.table(file_gene,header = F,sep = "\t",stringsAsFactors = F)
  rownames(geneFile)=geneFile$V1
  sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
  expr=expr[intersect(gn,geneFile$V1),]
  head(sub_geneFile,4)
  expr[1:4,1:4]
  dim(expr)
  
  expr_observations<-expr[,test_loc_vector]
  expr_observations[1:4,1:4]
  dim(expr_observations)
  
  cnv_score <- cnvScore(expr_observations)
  
  expr_references<-expr[,normal_loc_vector]
  expr_references[1:4,1:4]
  dim(expr_references)

	cnv_res<-estimateCNV(expr_observations, expr_references, CNVscore_cutoff, cor_cutoff)
	head(cnv_res)
	write.table(cnv_res,file=paste0(infercnv_dir,"/CNV_score_cor.tsv"),quote=F,sep="\t", row.names=T)
	
	summary(cnv_res$score)
	plota<-NULL
	plota<-ggplot(data = cnv_res,aes(x=score))+geom_histogram(bins = 1000)+geom_histogram(bins = 1000)+geom_vline(xintercept = CNVscore_cutoff,linewidth=0.1,col='red')
			ggsave(paste0("fig",plotid,"a.histogram.score",".tiff"), plot = plota, width = 10, height = 5,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,"a.histogram.score",".pdf"), plot = plota, width = 10, height = 5,limitsize=F)
	plota<-NULL
	plota<-ggplot(data = cnv_res,aes(x=cor))+geom_histogram(bins = 1000)+geom_histogram(bins = 1000)+geom_vline(xintercept = cor_cutoff,linewidth=0.1,col='red')
			ggsave(paste0("fig",plotid,"b.histogram.cor",".tiff"), plot = plota, width = 10, height = 5,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,"b.histogram.cor",".pdf"), plot = plota, width = 10, height = 5,limitsize=F)

			head(cnv_res)
			malignant_cells<-rownames(cnv_res[cnv_res$type=='Malignant',])
			print(length(malignant_cells))
			head(pbmc@meta.data)
			pbmc@meta.data$malignant_type<-'Normal'
			pbmc@meta.data[rownames(pbmc@meta.data) %in% malignant_cells,'malignant_type']<-'Malignant'
			table(pbmc@meta.data$malignant_type)
			
			
			
			pbmc@meta.data$CNV_score<-0
			for(i in 1:nrow(cnv_res))
			{
				cell<-rownames(cnv_res)[i]
				pbmc@meta.data[cell,'CNV_score']<-cnv_res[i,'score']
				pbmc@meta.data[cell,'CNV_cor']<-cnv_res[i,'cor']
			}
			dim(pbmc@meta.data[pbmc@meta.data$CNV_score>0,])
			pbmc@meta.data$CNV_score_used<-pbmc@meta.data$CNV_score
			pbmc@meta.data[!pbmc@meta.data$celltype_manual == 'Epithelial_cells','CNV_score_used']<-0
			
			gene='CNV_score'
			plota<-NULL
			plota<-FeaturePlot(object = pbmc, features = gene) & scale_colour_gradientn(colours = c("red",'red',"gray"),values = c(1.0,0.2,0.1,0))
			ggsave(paste0("fig",plotid,"c.FeaturePlot.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,"c.FeaturePlot.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)
			
			gene='CNV_score_used'
			plota<-NULL
			plota<-FeaturePlot(object = pbmc, features = gene) & scale_colour_gradientn(colours = c("red",'red',"gray"),values = c(1.0,0.2,0.1,0))
			ggsave(paste0("fig",plotid,"d.FeaturePlot.",gene,".tiff"), plot = plota, width = 6, height = 5,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,"d.FeaturePlot.",gene,".pdf"), plot = plota, width = 6, height = 5,limitsize=F)

			plota<-NULL
			plota<-DimPlot(pbmc, reduction = "umap", group.by = 'malignant_type',pt.size = 0.5,label=T)
			ggsave(paste0("fig",plotid,"e.umap.malignant_type.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
			ggsave(paste0("fig",plotid,"e.umap.malignant_type.pdf"), plot = plota, width = 8, height = 8)

			plotc<-NULL
			plotc <- DimPlot(pbmc, reduction = "umap", split.by = "malignant_type")
			ggsave(paste0("fig",plotid,"f.umap.malignant_type.split.tiff"), plot = plotc, width = 12, height = 6,compression='lzw')
			ggsave(paste0("fig",plotid,"f.umap.malignant_type.split.pdf"), plot = plotc, width = 12, height = 6)
			return(pbmc)
}



seurat_2celltype_deg <- function(pbmc,plotid='04',group.by='celltype',group_test='CD',group_control='control',plot.only=F,logFCfilter=1,adjPvalFilter=0.05,go_analysis=T){
######################################################################################################################
###  计算差异表达
	DefaultAssay(pbmc) <- "RNA"
	Idents(pbmc) <- group.by
	pbmc.used<-subset(x=pbmc,idents=c(group_test,group_control))
	comparison<-paste0(group_test,'_vs_',group_control)
	comparison<-gsub(' ','',comparison)
	
	    if(!plot.only)
		{
		            deg_all <- FindMarkers(pbmc.used,ident.1 = group_test, ident.2 = group_control,group_by=group.by, verbose = FALSE,only.pos = FALSE,min.pct = 0.1,logfc.threshold = 0)
		            deg_all$comparison<-comparison
					deg_all$gene<-rownames(deg_all)		            
				    write.table(deg_all,file=paste0("fig",plotid,'.',comparison,".deg.all.tsv"),sep="\t",row.names=F,quote=F)
				    deg_sig=deg_all[(abs(as.numeric(as.vector(deg_all$avg_log2FC)))>logFCfilter & as.numeric(as.vector(deg_all$p_val_adj))<adjPvalFilter),]
				    write.table(deg_sig,file=paste0("fig",plotid,'.',comparison,".deg.sig.tsv"),sep="\t",row.names=F,quote=F)
		}else
		{
					deg_all<-read.table(paste0("fig",plotid,'.',comparison,".deg.all.tsv"),header=T,sep='\t')
				    deg_sig<-read.table(paste0("fig",plotid,'.',comparison,".deg.sig.tsv"),header=T,sep='\t')
		}

		ggplot2_volcano(deg_all[,c('p_val_adj','avg_log2FC')],comparison=comparison,output=paste0('fig',plotid,'.volcano.',comparison),
		fc_threshold=2^logFCfilter,pvalue_threshold=adjPvalFilter,ylab='p_val_adj')
		top5.celltype <- deg_sig %>% top_n(n = 5, wt = avg_log2FC)
		top5.head<-top5.celltype$gene
		top5.celltype <- deg_sig %>% top_n(n = -5, wt = avg_log2FC)
		top5.tail<-top5.celltype$gene

		genes_used<-unique(append(top5.head,top5.tail))
		#cell_used<-rownames(pbmc@meta.data[pbmc@meta.data[[group.by]] %in% c(group_test,group_control),])		
		#pbmc.used<-subset(x=pbmc,cells=cell_used)		
		pbmc.used<-subset(x=pbmc,idents=c(group_test,group_control))	
		pdf(file=paste0("fig",plotid,"d.",comparison,".top5.Heatmap.pdf"),width=12,height=9)
			print(DoHeatmap(object = pbmc.used, features = genes_used) + NoLegend())
		dev.off()
		plota<-NULL
		plota<-FeaturePlot(object = pbmc.used, features = genes_used, split.by = group.by)
		ggsave(paste0("fig",plotid,"g.",comparison,".top5.FeaturePlot.tiff"), plot = plota, width = 10, height = 10*5,compression='lzw',limitsize=F)
		ggsave(paste0("fig",plotid,"g.",comparison,".top5.FeaturePlot.pdf"), plot = plota, width = 10, height = 10*5,limitsize=F)
		seurat_markers_plot(pbmc.used,genes_used,group.by=group.by,plotid=plotid)
		if(go_analysis)
		{
		aaa<-deg_sig$gene[deg_sig$avg_log2FC>0]
		gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UP'),paste0('kegg_',comparison,'_UP')),species = species,kegg_pvalueCutoff=0.99)
		
		aaa<-deg_sig$gene[deg_sig$avg_log2FC<0]
		gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_DOWN'),paste0('kegg_',comparison,'_DOWN')),species = species,kegg_pvalueCutoff=0.99)
		
		aaa<-deg_sig$gene
		gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UPDOWN'),paste0('kegg_',comparison,'_UPDOWN')),species = species,kegg_pvalueCutoff=0.99)
		
		aaa<-deg_all[,c('gene','avg_log2FC')]
		gsea_yingbai(aaa,dir_output=paste0('gsea_',comparison),species = species)
		}

		
}


seurat_celltype_deg <- function(pbmc,plotid='04',group.by='celltype',plot.only=F,logFCfilter=1,adjPvalFilter=0.05,topn=5){
######################################################################################################################
###  计算差异表达
	DefaultAssay(pbmc) <- "RNA"
	Idents(pbmc) <- group.by
	levels_backup<-levels(pbmc)
	print(levels(pbmc))
	
	    if(!plot.only)
		{
				cellMarkers=FindAllMarkers(object = pbmc,only.pos = FALSE,min.pct = 0.25,logfc.threshold = 0,return.thresh=0.05)
				write.table(cellMarkers,file=paste0("fig",plotid,".cellMarkers.all.tsv"),sep="\t",row.names=F,quote=F)
				sig.cellMarkers=cellMarkers[(abs(as.numeric(as.vector(cellMarkers$avg_log2FC)))>=logFCfilter & as.numeric(as.vector(cellMarkers$p_val_adj))<=adjPvalFilter),]
				write.table(sig.cellMarkers,file=paste0("fig",plotid,".cellMarkers.sig.tsv"),sep="\t",row.names=F,quote=F)
		}else
		{
				sig.cellMarkers<-read.table(paste0("fig",plotid,".cellMarkers.sig.tsv"),header=T,sep='\t')
		}

	df_count<-as.data.frame(table(pbmc@meta.data[[group.by]]))
	colnames(df_count)<-c(group.by,'Freq')
	df_count<-df_count[order(df_count$Freq,decreasing=T),]
	print(levels(pbmc))
	sig.cellMarkers$cluster<-factor(sig.cellMarkers$cluster,levels=levels_backup)
	sig.cellMarkers<-sig.cellMarkers[order(sig.cellMarkers$cluster),]
	#sig.cellMarkers<-read.table("fig04.cellMarkers.sig.tsv",header=T,sep='\t')
	######################################################################################################################
	###  绘图
	ncelltype<-length(levels(pbmc))

	top5.celltype <- sig.cellMarkers %>% group_by(cluster) %>% top_n(n = topn, wt = avg_log2FC)
	top5.head<-top5.celltype$gene
	top5.head<-unique(top5.head)
			seurat_markers_plot(pbmc,top5.head,group.by=group.by,plotid=plotid,dotplot.only=F,plot.heatmap=T)
}

seurat_cluster_deg <- function(pbmc,plotid='03',plot.only=F,logFCfilter=1,adjPvalFilter=0.05,min.pct = 0.25){

		library(ggplot2)
		library(Seurat)
		library(pheatmap)

		DefaultAssay(pbmc) <- "RNA"
		Idents(pbmc) <- "seurat_clusters" 
		if(!plot.only)
		{
				pbmc.markers <- FindAllMarkers(object = pbmc,only.pos = FALSE,min.pct = min.pct,logfc.threshold = 0,return.thresh=0.05)
				write.table(pbmc.markers,file=paste0("fig",plotid,".clusterMarkers.all.tsv"),sep="\t",row.names=F,quote=F)
				sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>=logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<=adjPvalFilter),]
				write.table(sig.markers,file=paste0("fig",plotid,".clusterMarkers.sig.tsv"),sep="\t",row.names=F,quote=F)
		}else
		{
				sig.markers<-read.table(paste0("fig",plotid,".clusterMarkers.sig.tsv"),header=T,sep='\t')
		}
		#sig.markers<-read.table("top10_markers_for_each_cluster_anno.xls",header=T,sep='\t')

		ncluster<-length(levels(pbmc))

		top5 <- sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
		top5.head<-top5$gene
		top5.head<-unique(top5.head)

		seurat_markers_plot(pbmc,top5.head,group.by="seurat_clusters",ngene_cutoff=20000,plotid=plotid,dotplot.only=F,plot.heatmap=T)
}

##  seurat_cluster_deg_single(pbmc,clusterid='12',ngene=10,plotid='24')	
seurat_cluster_deg_single <- function(pbmc,clusterid='0',ngene=5,plotid='24',plot.only=F,logFCfilter=1,adjPvalFilter=0.05){

		library(ggplot2)
		library(Seurat)
		library(pheatmap)

		DefaultAssay(pbmc) <- "RNA"
		Idents(pbmc) <- "seurat_clusters" 
		if(!plot.only)
		{
				deg.sig <- FindMarkers(object = pbmc,ident.1 = clusterid,only.pos = FALSE,min.pct = 0.25,logfc.threshold = 1,return.thresh=0.05)
				write.table(deg.sig,file=paste0("fig",plotid,'.cluster',clusterid,".clusterMarkers.sig.tsv"),sep="\t",row.names=T,quote=F)
		}else
		{
				deg.sig<-read.table(paste0("fig",plotid,'.cluster',clusterid,".clusterMarkers.sig.tsv"),header=T,sep='\t')
		}
		#sig.markers<-read.table("top10_markers_for_each_cluster_anno.xls",header=T,sep='\t')

		ncluster<-length(levels(pbmc))

		top5 <- deg.sig %>% top_n(n = ngene, wt = avg_log2FC)
		top5$gene<-rownames(top5)
		top5.head<-top5$gene
		top5.head<-unique(top5.head)

		#绘制marker在各个cluster的热图
		if(1)
		{
			pdf(file=paste0("fig",plotid,"b.cluster.marker.top5.Heatmap.pdf"),width=ncluster*0.8+2,height=ncluster*5*16/100,onefile=F)
			#print(DoHeatmap(object = pbmc, features = top5.head) + NoLegend())
			print(DoHeatmap(object = pbmc, features = top5.head,label=F))
			dev.off()
			tiff(file=paste0("fig",plotid,"c.cluster.marker.top5.Heatmap.tiff"),width=ncluster*0.8*300+600,height=ncluster*5*16/100*300,compression='lzw',res=300)
			#print(DoHeatmap(object = pbmc, features = top5.head) + NoLegend())
			print(DoHeatmap(object = pbmc, features = top5.head,label=F))
			dev.off()
		}
		if(0)
		{
			pdf(file=paste0("fig",plotid,"c.cluster.marker.top5.Heatmap.pdf"),width=ncluster*0.8+2,height=ncluster*5*16/100)
			print(DoHeatmap(object = pbmc, features = top5.head,label=F)+scale_fill_gradientn(colors = c("#406aa8", "white", "#d91216")))
			#+theme(text=element_text(size = 20))
			dev.off()
			}
		seurat_markers_plot(pbmc,top5$gene,group.by='seurat_clusters',plotid=plotid)
}

#seurat_deg_bygroup(pbmc,plotid='08',group.by='group',group_test='PE',group_control='control',species = "human",go_analysis=T,go_method='combi',logFCfilter=1,adjPvalFilter=0.05)
seurat_deg_bygroup <- function(pbmc,plotid='08',group.by='group',group_test='PE',group_control='control',species = "human",go_analysis=T,go_method='combi',logFCfilter=1,adjPvalFilter=0.05){
######################################################################################################################
###  计算差异表达
		DefaultAssay(pbmc) <- "RNA"
		Idents(pbmc) <- group.by

		plota<-NULL
		plota<-DimPlot(pbmc, reduction = "umap", label = F, pt.size = 0.5,split.by=group.by)+NoLegend()
		ggsave(paste0("fig",plotid,"a.UMAP.",group.by,".tiff"), plot = plota, width = 10, height = 5,compression='lzw',limitsize=F)
		ggsave(paste0("fig",plotid,"a.UMAP.",group.by,".pdf"), plot = plota, width = 10, height = 5,limitsize=F)
		metadata <- pbmc@meta.data
		plot_cellcount_barplot(metadata,plotid=plotid,x=group.by)

		comparison<-paste0(group_test,'_vs_',group_control)
		pbmc.markers.by.group <- FindMarkers(pbmc,ident.1 = group_test, ident.2 = group_control,group_by=group.by, verbose = FALSE,only.pos = FALSE,min.pct = 0.1,logfc.threshold = 0)
		pbmc.markers.by.group$comparison<-comparison
		pbmc.markers.by.group$gene<-rownames(pbmc.markers.by.group)
		write.table(pbmc.markers.by.group,file=paste0("fig",plotid,"c.deg.",comparison,".all.tsv"),sep="\t",row.names=F,quote=F)
		
		sig.pbmc.markers.by.group=pbmc.markers.by.group[(abs(as.numeric(as.vector(pbmc.markers.by.group$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers.by.group$p_val_adj))<adjPvalFilter),]
		write.table(sig.pbmc.markers.by.group,file=paste0("fig",plotid,"c.deg.",comparison,".sig.tsv"),sep="\t",row.names=F,quote=F)
		
		if(0)
		{
				sig.markers<-read.table(paste0("fig",plotid,"c.deg.",comparison,".sig.tsv"),header=T,sep='\t')
		}
		
		ggplot2_volcano(pbmc.markers.by.group[,c('p_val_adj','avg_log2FC')],comparison=paste0(group_test,'_vs_',group_control),
		output=paste0('fig',plotid,'d.volcano.',comparison),fc_threshold=2^logFCfilter,pvalue_threshold=adjPvalFilter,ylab='p_val_adj')
		seurat_cell_density_plot(pbmc,group.by='group',plotid=plotid)
		
		ndeg<-nrow(sig.pbmc.markers.by.group)
		if(ndeg>0)
		{
			top5.celltype <- sig.pbmc.markers.by.group %>% top_n(n = 5, wt = avg_log2FC)
			top5.head<-top5.celltype$gene
			top5.celltype <- sig.pbmc.markers.by.group %>% top_n(n = -5, wt = avg_log2FC)
			top5.tail<-top5.celltype$gene		
			genes_used<-unique(append(top5.head,top5.tail))

			plota<-NULL
			plota<-FeaturePlot(object = pbmc, features = genes_used, split.by = group.by,cols = c("green", "red"))
			ggsave(paste0("fig",plotid,"g.",comparison,".top5.FeaturePlot.tiff"), plot = plota, width = 10, height = 10*5,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,"g.",comparison,".top5.FeaturePlot.pdf"), plot = plota, width = 10, height = 10*5,limitsize=F)
			seurat_markers_plot(pbmc,genes_used,group.by=group.by,ngene_cutoff=20000,plotid=paste0(plotid,'h'))
		}
		if(ndeg>1)
		{
			seurat_markers_plot_heatmap(pbmc,genes_used,group.by=group.by,output=paste0('fig',plotid,'.DEG.top5.Heatmap'))
		}

}


seurat_deg_bygroup_gokegg <- function(pbmc,plotid='08',group.by='group',group_test='PE',group_control='control',species = "human",go_analysis=T,go_method='combi',logFCfilter=1,adjPvalFilter=0.05){
######################################################################################################################
###  计算差异表达
		DefaultAssay(pbmc) <- "RNA"
		Idents(pbmc) <- group.by
		comparison<-paste0(group_test,'_vs_',group_control)
		pbmc.markers.by.group<-read.table(paste0("fig",plotid,"c.deg.",comparison,".all.tsv"),header=T,sep='\t')
		#sig.pbmc.markers.by.group<-read.table(paste0("fig",plotid,"c.deg.",comparison,".sig.tsv"),header=T,sep='\t')
		if(go_analysis)
		{
		sig.pbmc.markers.by.group=pbmc.markers.by.group[(abs(as.numeric(as.vector(pbmc.markers.by.group$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers.by.group$p_val_adj))<adjPvalFilter),]

		aaa<-sig.pbmc.markers.by.group$gene[sig.pbmc.markers.by.group$avg_log2FC>0]
		ngene<-length(aaa)
		print(paste('up genes: ',ngene ))
		gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UP'),paste0('kegg_',comparison,'_UP')),species = species)
		
		aaa<-sig.pbmc.markers.by.group$gene[sig.pbmc.markers.by.group$avg_log2FC<0]
		ngene<-length(aaa)
		print(paste('down genes: ',ngene ))
		gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_DOWN'),paste0('kegg_',comparison,'_DOWN')),species = species)
		
		aaa<-sig.pbmc.markers.by.group$gene
		ngene<-length(aaa)
		print(paste('up&down genes: ',ngene ))
		gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UPDOWN'),paste0('kegg_',comparison,'_UPDOWN')),species = species)
		
		#aaa<-pbmc.markers.by.group[,c('gene','avg_log2FC')]
		#gsea_yingbai(aaa,dir_output=paste0('gsea_',comparison),species = species)
		}
}

seurat_doubletfinder_new <- function(pbmc,pcSelect=20,plotid='11'){
		library(DoubletFinder)
		installed.packages()[c('Seurat','DoubletFinder'),c('Package','Version','LibPath')]
		sample<-unique(pbmc@meta.data$orig.ident)
		pbmc@meta.data$doubletfinder<-NULL
		
paramSweep<-function (seu, PCs = 1:10, sct = FALSE, num.cores = 1)
{
    require(Seurat)
    require(fields)
    require(parallel)
    pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
    pN <- seq(0.05, 0.3, by = 0.05)
    min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
    pK.test <- round(pK * min.cells)
    pK <- pK[which(pK.test >= 1)]
    orig.commands <- seu@commands
    if (nrow(seu@meta.data) > 10000) {
        real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data),
            10000, replace = FALSE)]
        data <- seu@assays$RNA@counts[, real.cells]
        n.real.cells <- ncol(data)
    }
    if (nrow(seu@meta.data) <= 10000) {
        real.cells <- rownames(seu@meta.data)
        data <- seu@assays$RNA@counts
        n.real.cells <- ncol(data)
    }
    if (num.cores > 1) {
        require(parallel)
        cl <- makeCluster(num.cores)
        output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep,
            n.real.cells, real.cells, pK, pN, data, orig.commands,
            PCs, sct, mc.cores = num.cores)
        stopCluster(cl)
    }
    else {
        output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep,
            n.real.cells, real.cells, pK, pN, data, orig.commands,
            PCs, sct)
    }
    sweep.res.list <- list()
    list.ind <- 0
    for (i in 1:length(output2)) {
        for (j in 1:length(output2[[i]])) {
            list.ind <- list.ind + 1
            sweep.res.list[[list.ind]] <- output2[[i]][[j]]
        }
    }
    name.vec <- NULL
    for (j in 1:length(pN)) {
        name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK,
            sep = "_"))
    }
    names(sweep.res.list) <- name.vec
    return(sweep.res.list)
}

		
		sweep.res.list <- paramSweep(pbmc, PCs = 1:pcSelect, sct = F)
		
		
		
		
		sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
		tiff(file=paste0("fig",plotid,'.',sample,".bcmvn.tiff"),width=1800,heigh=1800,units='px',res=300,compression='lzw')
		bcmvn <- find.pK(sweep.stats) 
		dev.off()
		pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
		# DoubletRate = ncol(pbmc)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
		DoubletRate = 0.075    # 直接查表，10000细胞对应的doublets rate是～7.6%
		DoubletRate = ncol(pbmc)*8*1e-6 #更通用
		#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
		homotypic.prop <- modelHomotypic(pbmc$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
		# 计算双细胞比例
		nExp_poi <- round(DoubletRate*ncol(pbmc)) 
		# 使用同源双细胞比例对计算的双细胞比例进行校正 
		nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
		## 使用确定好的参数鉴定doublets
		pbmc <- doubletFinder(pbmc, PCs = 1:pcSelect, pN = 0.25, pK = pK_bcmvn, 
								  nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
		colnames(pbmc@meta.data)[ncol(pbmc@meta.data)]<-'doubletfinder'
		## 结果展示，分类结果在pbmc@meta.data中
		plota<-DimPlot(pbmc, reduction = "umap",group.by = "doubletfinder")
		ggsave(paste0("fig",plotid,'.',sample,".doubletfinder.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
		ggsave(paste0("fig",plotid,'.',sample,".doubletfinder.pdf"), plot = plota, width = 10, height = 10)
		metadata<-pbmc@meta.data
		metadata<-cbind(rownames(metadata),metadata)
		write.table(metadata,file=paste0("fig",plotid,'.',sample,".metadata.tsv"),quote=F,sep="\t", row.names=F)
		pbmc=subset(x = pbmc, subset = doubletfinder =='Singlet')
		return(pbmc)
}



seurat_doubletfinder <- function(pbmc,pcSelect=20,plotid='11'){
		library(DoubletFinder)
		sample<-unique(pbmc@meta.data$orig.ident)
		pbmc@meta.data$doubletfinder<-NULL
		sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:pcSelect, sct = F)
		sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
		tiff(file=paste0("fig",plotid,'.',sample,".bcmvn.tiff"),width=1800,heigh=1800,units='px',res=300,compression='lzw')
		bcmvn <- find.pK(sweep.stats) 
		dev.off()
		pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
		# DoubletRate = ncol(pbmc)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
		DoubletRate = 0.075    # 直接查表，10000细胞对应的doublets rate是～7.6%
		DoubletRate = ncol(pbmc)*8*1e-6 #更通用
		#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
		homotypic.prop <- modelHomotypic(pbmc$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
		# 计算双细胞比例
		nExp_poi <- round(DoubletRate*ncol(pbmc)) 
		# 使用同源双细胞比例对计算的双细胞比例进行校正 
		nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
		## 使用确定好的参数鉴定doublets
		pbmc <- doubletFinder_v3(pbmc, PCs = 1:pcSelect, pN = 0.25, pK = pK_bcmvn, 
								  nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
		colnames(pbmc@meta.data)[ncol(pbmc@meta.data)]<-'doubletfinder'
		## 结果展示，分类结果在pbmc@meta.data中
		plota<-DimPlot(pbmc, reduction = "umap",group.by = "doubletfinder")
		ggsave(paste0("fig",plotid,'.',sample,".doubletfinder.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
		ggsave(paste0("fig",plotid,'.',sample,".doubletfinder.pdf"), plot = plota, width = 10, height = 10)
		metadata<-pbmc@meta.data
		metadata<-cbind(rownames(metadata),metadata)
		write.table(metadata,file=paste0("fig",plotid,'.',sample,".metadata.tsv"),quote=F,sep="\t", row.names=F)
		pbmc=subset(x = pbmc, subset = doubletfinder =='Singlet')
		return(pbmc)
}

seurat_doubletfinder_version5 <- function(pbmd,pcSelect=20,plotid='11'){
		library(DoubletFinder)
		installed.packages()[c('Seurat','DoubletFinder'),c('Package','Version','LibPath')]
		sample<-unique(pbmd@meta.data$orig.ident)
		pbmd@meta.data$doubletfinder<-NULL
		sweep.res.list <- paramSweep(pbmd, PCs = 1:pcSelect, sct = F)
		sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
		
		tiff(file=paste0("fig",plotid,'.',sample,".bcmvn.tiff"),width=1800,heigh=1800,units='px',res=300,compression='lzw')
		bcmvn <- find.pK(sweep.stats) 
		dev.off()
		
		pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
		# DoubletRate = ncol(pbmc)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
		DoubletRate = 0.075    # 直接查表，10000细胞对应的doublets rate是～7.6%
		DoubletRate = ncol(pbmd)*8*1e-6 #更通用
		#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
		homotypic.prop <- modelHomotypic(pbmd$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
		# 计算双细胞比例
		nExp_poi <- round(DoubletRate*ncol(pbmd)) 
		# 使用同源双细胞比例对计算的双细胞比例进行校正 
		nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
		## 使用确定好的参数鉴定doublets
		pbmd <- doubletFinder(pbmd, PCs = 1:pcSelect, pN = 0.25, pK = pK_bcmvn, 
								  nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
		colnames(pbmd@meta.data)[ncol(pbmd@meta.data)]<-'doubletfinder'
		## 结果展示，分类结果在pbmc@meta.data中
		plota<-DimPlot(pbmd, reduction = "umap",group.by = "doubletfinder")
		ggsave(paste0("fig",plotid,'.',sample,".doubletfinder.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
		ggsave(paste0("fig",plotid,'.',sample,".doubletfinder.pdf"), plot = plota, width = 10, height = 10)
		metadata<-pbmd@meta.data
		metadata<-cbind(rownames(metadata),metadata)
		write.table(metadata,file=paste0("fig",plotid,'.',sample,".metadata.tsv"),quote=F,sep="\t", row.names=F)
		pbmd=subset(x = pbmd, subset = doubletfinder =='Singlet')
		return(pbmd)
}


monocle2_step1_orderCells <- function(pbmc,projectid='monocle2',max_components = 2,groupby='celltype')
{
	
library(monocle)
library(ggsci)
		data <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
		monocle.geneAnn=data.frame(gene_short_name = row.names(data), row.names = row.names(data))
		pd<-new("AnnotatedDataFrame", data = pbmc@meta.data)
		fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
		rm(pbmc)
		gc()
		cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
		dim(cds)
		rm(data)
		gc()

		#细胞轨迹分析流程
		### 2估计size factor和离散度
			print(date())
		cds <- estimateSizeFactors(cds)
			print(date())
		cds <- estimateDispersions(cds)
			print(date())
			gc()
		### 3过滤低质量细胞
		cds<-detectGenes(cds,min_expr = 0.1)
			print(head(fData(cds)))#此时有多少个基因
		#过滤掉在小于10个细胞中表达的基因
		expressed_genes<-row.names(subset(fData(cds),num_cells_expressed>=10))
			dim(cds)
			gc()
		###   各增加了一行，统计基因的细胞数，细胞的基因数。

		### 4细胞分类，暂时使用seurat的分类，跳过

		### 5.1特征基因选择：有多种方法，使用monocle计算高变。
		#######  计算数据离散度，根据离散度选择用于计算的基因
		disp_table <- dispersionTable(cds)
		disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
		write.table(disp.genes, paste0("fig05a.",projectid,".disp.genes.tsv"),sep='\t',row.names = F)
				gc()
		cds <- setOrderingFilter(cds, disp.genes)
				gc()

		tiff(file=paste0("fig05a.",projectid,".featureSelected.tiff"),width=1800,heigh=1800,units='px',res=300,compression='lzw')
			plot_ordering_genes(cds)
		dev.off()



		### 5.2降维
		print(date())
		cds <- reduceDimension(cds, max_components = max_components, reduction_method = 'DDRTree')
		print(date())
		gc()
		### 5.3拟时间轨迹构建，排列细胞
		print(date())
		cds <- orderCells(cds)
		print(date())
		gc()
		write.table(pData(cds), paste0("fig05b.",projectid,".pseudotime.tsv"),sep='\t',row.names = T)
		
		#保存细胞state的细胞轨迹图
		plota<-NULL
		plota<-plot_cell_trajectory(cds,color_by = "State",show_state_number =F)
		ggsave(paste0("fig05a.",projectid,".trajectory.State.tiff"), plot = plota, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05a.",projectid,".trajectory.State.pdf"), plot = plota, width = 12, height =12)

		#保存时间的细胞轨迹图
		plotb<-NULL
		plotb<-plot_cell_trajectory(cds,color_by = "Pseudotime")
		ggsave(paste0("fig05b.",projectid,".trajectory.Pseudotime.tiff"), plot = plotb, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05b.",projectid,".trajectory.Pseudotime.pdf"), plot = plotb, width = 12, height =12)
		#保存细胞名称的细胞轨迹图
		plotc<-NULL
		plotc<-plot_cell_trajectory(cds,color_by = groupby)
		ggsave(paste0("fig05c.",projectid,".trajectory.",groupby,".tiff"), plot = plotc, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05c.",projectid,".trajectory.",groupby,".pdf"), plot = plotc, width = 12, height =12)
		
		ntype<-length(unique(pData(cds)[[groupby]]))
		nr<-ceiling(ntype/4)
		nw=4
		if(ntype<4)nw=ntype
		
		plotf<-NULL
		plotf<-plot_cell_trajectory(cds, color_by = groupby) + facet_wrap(paste0("~",groupby), ncol= 4)
		ggsave(paste0("fig05c1.",projectid,".trajectory.",groupby,".faceted.tiff"), plot = plotf, width = nw*6, height = nr*6,compression='lzw')
		ggsave(paste0("fig05c1.",projectid,".trajectory.",groupby,".faceted.pdf"), plot = plotf, width = nw*6, height =nr*6)
		
		#保存聚类的细胞轨迹图
		plotd<-NULL
		plotd<-plot_cell_trajectory(cds, color_by = "seurat_clusters")
		ggsave(paste0("fig05d.",projectid,".trajectory.seurat_clusters.tiff"), plot = plotd, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05d.",projectid,".trajectory.seurat_clusters.pdf"), plot = plotd, width = 12, height =12)
		
		ntype<-length(unique(pData(cds)$seurat_clusters))
		nr<-ceiling(ntype/4)
		nw=4
		if(ntype<4)nw=ntype
		
		plotf<-NULL
		plotf<-plot_cell_trajectory(cds, color_by = 'seurat_clusters') + facet_wrap(paste0("~",'seurat_clusters'), ncol= 4)
		ggsave(paste0("fig05d1.",projectid,".trajectory.",'seurat_clusters',".faceted.tiff"), plot = plotf, width = nw*6, height = nr*6,compression='lzw')
		ggsave(paste0("fig05d1.",projectid,".trajectory.",'seurat_clusters',".faceted.pdf"), plot = plotf, width = nw*6, height =nr*6)

		###### 拆分轨道面
		ntype<-length(unique(pData(cds)$State))
		nr<-ceiling(ntype/4)
		nw=4
		if(ntype<4)nw=ntype
		
		plote<-NULL
		plote<-plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", ncol = 4)
		ggsave(paste0("fig05e.",projectid,".trajectory.state.faceted.tiff"), plot = plote, width = nw*6, height =nr*6,compression='lzw')
		ggsave(paste0("fig05e.",projectid,".trajectory.state.faceted.pdf"), plot = plote, width = nw*6, height =nr*6)
		
		pearplot<-NULL
		pearplot <-CombinePlots(plots = list(plota, plotb,plotc, plotd),ncol=2)
		ggsave(paste0("fig05e1.",projectid,".trajectory.combination.tiff"), plot = pearplot, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05e1.",projectid,".trajectory.combination.pdf"), plot = pearplot, width = 12, height =12)


		p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "seurat_clusters") + 
		  theme(legend.position='none',panel.border = element_blank()) 
		p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,color_by = "seurat_clusters")+
		  theme(legend.title = element_blank()) 
		plotf<-p1|p2
		ggsave(paste0("fig05f.",projectid,".trajectory.tree.tiff"), plot = plotf, width = 16, height = 12,compression='lzw')
		ggsave(paste0("fig05f.",projectid,".trajectory.tree.pdf"), plot = plotf, width = 16, height =12)
		return(cds)
}

monocle2_step1_orderCells_v5 <- function(pbmc,projectid='monocle2',max_components = 2,groupby='celltype')
{
	
library(monocle)
library(ggsci)

		data <-GetAssayData(object = pbmc, layer = "counts")
		#data <- as(as.matrix(data), 'sparseMatrix')
		monocle.geneAnn=data.frame(gene_short_name = row.names(data), row.names = row.names(data))
		pd<-new("AnnotatedDataFrame", data = pbmc@meta.data)
		fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
		rm(pbmc)
		gc()
		cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
		dim(cds)
		rm(data)
		gc()

		#细胞轨迹分析流程
		### 2估计size factor和离散度
			print(date())
		cds <- estimateSizeFactors(cds)
			print(date())
		cds <- estimateDispersions(cds)
			print(date())
			gc()
		### 3过滤低质量细胞
		cds<-detectGenes(cds,min_expr = 0.1)
			print(head(fData(cds)))#此时有多少个基因
		#过滤掉在小于10个细胞中表达的基因
		expressed_genes<-row.names(subset(fData(cds),num_cells_expressed>=10))
			dim(cds)
			gc()
		###   各增加了一行，统计基因的细胞数，细胞的基因数。

		### 4细胞分类，暂时使用seurat的分类，跳过

		### 5.1特征基因选择：有多种方法，使用monocle计算高变。
		#######  计算数据离散度，根据离散度选择用于计算的基因
		disp_table <- dispersionTable(cds)
		disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
		write.table(disp.genes, paste0("fig05a.",projectid,".disp.genes.tsv"),sep='\t',row.names = F)
				gc()
		cds <- setOrderingFilter(cds, disp.genes)
				gc()

		tiff(file=paste0("fig05a.",projectid,".featureSelected.tiff"),width=1800,heigh=1800,units='px',res=300,compression='lzw')
			plot_ordering_genes(cds)
		dev.off()



		### 5.2降维
		print(date())
		cds <- reduceDimension(cds, max_components = max_components, reduction_method = 'DDRTree')
		print(date())
		gc()
		### 5.3拟时间轨迹构建，排列细胞
		print(date())
		cds <- orderCells(cds)
		print(date())
		gc()
		write.table(pData(cds), paste0("fig05b.",projectid,".pseudotime.tsv"),sep='\t',row.names = T)
		
		#保存细胞state的细胞轨迹图
		plota<-NULL
		plota<-plot_cell_trajectory(cds,color_by = "State",show_state_number =F)
		ggsave(paste0("fig05a.",projectid,".trajectory.State.tiff"), plot = plota, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05a.",projectid,".trajectory.State.pdf"), plot = plota, width = 12, height =12)

		#保存时间的细胞轨迹图
		plotb<-NULL
		plotb<-plot_cell_trajectory(cds,color_by = "Pseudotime")
		ggsave(paste0("fig05b.",projectid,".trajectory.Pseudotime.tiff"), plot = plotb, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05b.",projectid,".trajectory.Pseudotime.pdf"), plot = plotb, width = 12, height =12)
		#保存细胞名称的细胞轨迹图
		plotc<-NULL
		plotc<-plot_cell_trajectory(cds,color_by = groupby)
		ggsave(paste0("fig05c.",projectid,".trajectory.",groupby,".tiff"), plot = plotc, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05c.",projectid,".trajectory.",groupby,".pdf"), plot = plotc, width = 12, height =12)
		
		ntype<-length(unique(pData(cds)[[groupby]]))
		nr<-ceiling(ntype/4)
		nw=4
		if(ntype<4)nw=ntype
		
		plotf<-NULL
		plotf<-plot_cell_trajectory(cds, color_by = groupby) + facet_wrap(paste0("~",groupby), ncol= 4)
		ggsave(paste0("fig05c1.",projectid,".trajectory.",groupby,".faceted.tiff"), plot = plotf, width = nw*6, height = nr*6,compression='lzw')
		ggsave(paste0("fig05c1.",projectid,".trajectory.",groupby,".faceted.pdf"), plot = plotf, width = nw*6, height =nr*6)
		
		#保存聚类的细胞轨迹图
		plotd<-NULL
		plotd<-plot_cell_trajectory(cds, color_by = "seurat_clusters")
		ggsave(paste0("fig05d.",projectid,".trajectory.seurat_clusters.tiff"), plot = plotd, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05d.",projectid,".trajectory.seurat_clusters.pdf"), plot = plotd, width = 12, height =12)
		
		ntype<-length(unique(pData(cds)$seurat_clusters))
		nr<-ceiling(ntype/4)
		nw=4
		if(ntype<4)nw=ntype
		
		plotf<-NULL
		plotf<-plot_cell_trajectory(cds, color_by = 'seurat_clusters') + facet_wrap(paste0("~",'seurat_clusters'), ncol= 4)
		ggsave(paste0("fig05d1.",projectid,".trajectory.",'seurat_clusters',".faceted.tiff"), plot = plotf, width = nw*6, height = nr*6,compression='lzw')
		ggsave(paste0("fig05d1.",projectid,".trajectory.",'seurat_clusters',".faceted.pdf"), plot = plotf, width = nw*6, height =nr*6)

		###### 拆分轨道面
		ntype<-length(unique(pData(cds)$State))
		nr<-ceiling(ntype/4)
		nw=4
		if(ntype<4)nw=ntype
		
		plote<-NULL
		plote<-plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", ncol = 4)
		ggsave(paste0("fig05e.",projectid,".trajectory.state.faceted.tiff"), plot = plote, width = nw*6, height =nr*6,compression='lzw')
		ggsave(paste0("fig05e.",projectid,".trajectory.state.faceted.pdf"), plot = plote, width = nw*6, height =nr*6)
		
		pearplot<-NULL
		pearplot <-CombinePlots(plots = list(plota, plotb,plotc, plotd),ncol=2)
		ggsave(paste0("fig05e1.",projectid,".trajectory.combination.tiff"), plot = pearplot, width = 12, height = 12,compression='lzw')
		ggsave(paste0("fig05e1.",projectid,".trajectory.combination.pdf"), plot = pearplot, width = 12, height =12)


		p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "seurat_clusters") + 
		  theme(legend.position='none',panel.border = element_blank()) 
		p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,color_by = "seurat_clusters")+
		  theme(legend.title = element_blank()) 
		plotf<-p1|p2
		ggsave(paste0("fig05f.",projectid,".trajectory.tree.tiff"), plot = plotf, width = 16, height = 12,compression='lzw')
		ggsave(paste0("fig05f.",projectid,".trajectory.tree.pdf"), plot = plotf, width = 16, height =12)
		return(cds)
}


monocle2_step2_re_orderCells <- function(cds,projectid='monocle2',num_state=num_state,groupby='celltype')
{
library(monocle)
library(ggsci)
		if(!is.na(num_state))
		{
			if(num_state>1)
			{
				cds <- orderCells(cds, root_state = num_state)
				gc()
				###############
				####  基本没有改变，只有Pseudotime.tiff中颜色最深（0时刻）的位置改变了。
				write.table(pData(cds), paste0("fig05b.",projectid,".pseudotime.tsv"),sep='\t',row.names = T)
				
				plotb<-NULL
				plotb<-plot_cell_trajectory(cds,color_by = "Pseudotime")
				ggsave(paste0("fig05b1.",projectid,".trajectory.Pseudotime.tiff"), plot = plotb, width = 12, height = 12,compression='lzw')
				ggsave(paste0("fig05b1.",projectid,".trajectory.Pseudotime.pdf"), plot = plotb, width = 12, height =12)

				plota<-NULL
				plota<-plot_cell_trajectory(cds,color_by = "State",show_state_number =F)
				plotc<-NULL
				plotc<-plot_cell_trajectory(cds,color_by = groupby)
				plotd<-NULL
				plotd<-plot_cell_trajectory(cds, color_by = "seurat_clusters")
				pearplot <-CombinePlots(plots = list(plota, plotb,plotc, plotd),ncol=2)
				ggsave(paste0("fig05e2.",projectid,".trajectory.combination.tiff"), plot = pearplot, width = 12, height = 12,compression='lzw')
				ggsave(paste0("fig05e2.",projectid,".trajectory.combination.pdf"), plot = pearplot, width = 12, height =12)
			}
		}

		### 6 拟时相关基因

		fdata<-fData(cds)
		class(fdata)
		dim(fdata)
		ordering.genes<-fdata$gene_short_name[fdata$use_for_ordering=='TRUE']
		class(ordering.genes)
		length(ordering.genes)

		Time_diff <- differentialGeneTest(cds[ordering.genes,], cores = 1,fullModelFormulaStr = "~sm.ns(Pseudotime)")
		Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
			write.table(Time_diff, paste0("fig05.",projectid,".Time_diff_all.tsv"),sep='\t',row.names = F)
		Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

		sig_gene_names <- row.names(subset(Time_diff, qval < 0.05))
		head(Time_diff[sig_gene_names,])
		ngene<-length(sig_gene_names)
		print(ngene)
		ph=min(16/200*ngene,50)
		print(ph)
		
		#tiff(file=paste0("fig05g.",projectid,".pseudotime_heatmap.tiff"),width=3000,heigh=16/200*ngene*300,units='px',res=300,compression='lzw',limitsize=F)
		#plot_pseudotime_heatmap(cds[sig_gene_names,], num_clusters=4, show_rownames=T, return_heatmap=T)
		#dev.off()
		#pdf(file=paste0("fig05g.",projectid,".pseudotime_heatmap.pdf"),width=10,heigh=16/200*ngene)
			plota<-NULL
			plota<-plot_pseudotime_heatmap(cds[sig_gene_names,], num_clusters=4, show_rownames=T, return_heatmap=T)
			#ggsave(paste0("fig05g.",celltypei,".pseudotime_heatmap.tiff"), plot = plota, width = 10, height = 16/200*ngene,compression='lzw')
			ggsave(paste0("fig05g.",celltypei,".pseudotime_heatmap.pdf"), plot = plota, width = 10, height = ph,limitsize=F)
		#dev.off()


		### 7 指定基因的可视化（选择10个基因来测试）
		top10g_pseudotime <- rownames(Time_diff %>% top_n(n = -10, wt = qval))
		top10g_pseudotime<-top10g_pseudotime[1:10]

		cds_subset <- cds[top10g_pseudotime,]
		pearplot<-NULL
		plot1<-NULL
		plot2<-NULL
		plot3<-NULL
		plot1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
		plot2 <- plot_genes_in_pseudotime(cds_subset, color_by = groupby)
		plot3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
		pearplot <-CombinePlots(plots = list(plot1, plot2,plot3),ncol=3)
		ggsave(paste0("fig05h.",projectid,".Genes_pseudotimeplot.tiff"), plot = pearplot, width = 18, height = 8,compression='lzw')
		ggsave(paste0("fig05h.",projectid,".Genes_pseudotimeplot.pdf"), plot = pearplot, width = 18, height = 8)

		pearplot<-NULL
		plot1<-NULL
		plot2<-NULL
		plot3<-NULL
		plot1 <- plot_genes_jitter(cds_subset, grouping = "State", color_by = "State")
		plot2 <- plot_genes_violin(cds_subset, grouping = "State", color_by = "State")
		plot3 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
		pearplot <-CombinePlots(plots = list(plot1, plot2,plot3),ncol=3)
		ggsave(paste0("fig05i.",projectid,".Genes_Jitterplot.tiff"), plot = pearplot, width = 18, height = 8,compression='lzw')
		ggsave(paste0("fig05i.",projectid,".Genes_Jitterplot.pdf"), plot = pearplot, width = 18, height = 8)


		colnames(pData(cds))
		pearplot<-NULL
		plot1<-NULL
		plot2<-NULL
		plot3<-NULL
		name<-top10g_pseudotime[1]
		name<-gsub('-','_',name)
		
		pData(cds)[[name]] = log2( exprs(cds)[top10g_pseudotime[1],]+1)
		plot1=plot_cell_trajectory(cds, color_by = name)+ scale_color_gsea()+labs(color=top10g_pseudotime[1])
		name<-top10g_pseudotime[2]
		name<-gsub('-','_',name)
		pData(cds)[[name]] = log2(exprs(cds)[top10g_pseudotime[2],]+1)
		plot2=plot_cell_trajectory(cds, color_by = name) + scale_color_gsea()+labs(color=top10g_pseudotime[2])
		pearplot <-CombinePlots(plots = list(plot1, plot2))
		ggsave(paste0("fig05j.",projectid,".trajectory.top2.tiff"), plot = pearplot, width = 16, height = 8,compression='lzw')
		ggsave(paste0("fig05j.",projectid,".trajectory.top2.pdf"), plot = pearplot, width = 16, height = 8)
		return(cds)
}

monocle2_step3_beam <- function(cds,projectid='monocle2',branch_point_sel=1)
{
	library(monocle)
	library(ggsci)
		Time_diff<-read.table(paste0("fig05.",projectid,".Time_diff_all.tsv"),header=T)
		sig_gene_names <- Time_diff$gene_short_name[Time_diff$qval<0.05]

		BEAM_res <- BEAM(cds[sig_gene_names,], branch_point = branch_point_sel, cores = 8,progenitor_method = 'duplicate') 
		BEAM_res <- BEAM_res[order(BEAM_res$qval),]
		BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
		write.table(BEAM_res, paste0("fig05k.",projectid,".BEAM_res.branch",branch_point_sel,".tsv"),sep='\t',row.names = F)
        ### BEAM_res<-read.delim(paste0("fig05k.",celltypei,".BEAM_res.branch",branch_point_sel,".tsv"))
        
		#选前100个基因可视化
		BEAM_genes <-  BEAM_res %>% top_n(n = -100, wt = qval) %>% pull(gene_short_name) %>% as.character()
		length(BEAM_genes)
		##101

		tiff(file=paste0("fig05l.",projectid,".BEAM_res.heatmap.branch",branch_point_sel,".tiff"),width=3000,heigh=3000,units='px',res=300,compression='lzw')
			plot_genes_branched_heatmap(cds[sig_gene_names,],branch_point = branch_point_sel,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
		dev.off()
		pdf(file=paste0("fig05l.",projectid,".BEAM_res.heatmap.branch",branch_point_sel,".pdf"),width=10,heigh=10)
			plot_genes_branched_heatmap(cds[sig_gene_names,],branch_point = branch_point_sel,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
		dev.off()


		tiff(file=paste0("fig05m.",projectid,".BEAM_res.heatmap.branch",branch_point_sel,".top100.tiff"),width=1800,heigh=3000,units='px',res=300,compression='lzw')
			plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = branch_point_sel,num_clusters = 3,cores = 1,use_gene_short_name = T, show_rownames = T)
		dev.off()
		pdf(file=paste0("fig05m.",projectid,".BEAM_res.heatmap.branch",branch_point_sel,".top100.pdf"),width=6,heigh=10)
			plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = branch_point_sel,num_clusters = 3,cores = 1,use_gene_short_name = T, show_rownames = T)
		dev.off()

		top2g_beam<-BEAM_genes[1:2]
		tiff(file=paste0("fig05n.",projectid,".BEAM_res.pseudotime.branch",branch_point_sel,".top2.tiff"),width=1800,heigh=1800,units='px',res=300,compression='lzw')
			plot_genes_branched_pseudotime(cds[top2g_beam,],branch_point = branch_point_sel,color_by = "State", ncol = 1)
		dev.off()
		pdf(file=paste0("fig05n.",projectid,".BEAM_res.pseudotime.branch",branch_point_sel,".top2.pdf"),width=6,heigh=6)
			plot_genes_branched_pseudotime(cds[top2g_beam,],branch_point = branch_point_sel,color_by = "State", ncol = 1)
		dev.off()

		plota<-NULL
		plotb<-NULL
		pearplot<-NULL
		name<-top2g_beam[1]
		name<-gsub('-','_',name)
		pData(cds)[[name]] = log2( exprs(cds)[top2g_beam[1],]+1)
		plota=plot_cell_trajectory(cds, color_by = name)+scale_color_gradientn(top2g_beam[1],colors=c("green","red"))
		name<-top2g_beam[2]
		name<-gsub('-','_',name)
		pData(cds)[[name]] = log2(exprs(cds)[top2g_beam[2],]+1)
		plotb=plot_cell_trajectory(cds, color_by = name)+scale_color_gradientn(top2g_beam[2],colors=c("green","red"))
		pearplot <-CombinePlots(plots = list(plota, plotb))
		ggsave(paste0("fig05o.",projectid,".trajectory.branch",branch_point_sel,".top2.tiff"), plot = pearplot, width = 16, height = 8,compression='lzw')
		ggsave(paste0("fig05o.",projectid,".trajectory.branch",branch_point_sel,".top2.pdf"), plot = pearplot, width = 16, height = 8)
}

#monocle2_step3_beam_plot100(cds,projectid=projectid,branch_point_sel=1)
monocle2_step3_beam_plot100 <- function(cds,projectid='monocle2',branch_point_sel=1)
{
		BEAM_res<-read.table(paste0("fig05k.",projectid,".BEAM_res.branch",branch_point_sel,".tsv"),header=T,sep='\t')
		BEAM_genes <-  BEAM_res %>% top_n(n = -100, wt = qval) %>% pull(gene_short_name) %>% as.character()
		length(BEAM_genes)
		for(i in 1:20)
		{
				#i<-2
				top2g_beam<-BEAM_genes[(5*(i-1)+1):(5*i)]
				plota<-NULL
				plota<-plot_genes_branched_pseudotime(cds[top2g_beam,],branch_point = branch_point_sel,color_by = "sub_celltype", ncol = 1)
				ggsave(paste0("fig05n.",projectid,".BEAM_res.pseudotime.branch",branch_point_sel,'.',i,".tiff"), plot = plota, width = 8, height = 8,compression='lzw')
				ggsave(paste0("fig05n.",projectid,".BEAM_res.pseudotime.branch",branch_point_sel,'.',i,".pdf"), plot = plota, width = 8, height = 8)
		}		
		for(i in 1:100)
		{
				genei<-	BEAM_genes[i]
				plota<-NULL
				pData(cds)[[top2g_beam[1]]] = log2( exprs(cds)[top2g_beam[1],]+1)
				plota=plot_cell_trajectory(cds, color_by = top2g_beam[1])+ scale_color_gsea()
				ggsave(paste0("fig05o.",projectid,".trajectory.branch",branch_point_sel,'.',i,".tiff"), plot = plota, width = 8, height = 8,compression='lzw')
				ggsave(paste0("fig05o.",projectid,".trajectory.branch",branch_point_sel,'.',i,".pdf"), plot = plota, width = 8, height = 8)
		}		
}


#monocle2_markers_pseudotime_plot(cds,m6a_gene,group.by='celltype_sub',plotid='20',time_diff='fig05.Fibroblasts.Time_diff_all.tsv')
monocle2_markers_pseudotime_plot<- function(cds,markers,group.by='celltype',plotid='04',time_diff='fig05.Fibroblasts.Time_diff_all.tsv'){
library(monocle)
  cds_genes<-rownames(cds)
  markers_used<-markers[markers %in% cds_genes]
  ngene<-length(markers_used)
  print(ngene)
  if(ngene)
  {
		cds_subset <- cds[markers_used,]
		plot2<-NULL
		plot2 <- plot_genes_in_pseudotime(cds_subset, color_by = group.by)
		ggsave(paste0("fig",plotid,".Genes_pseudotimeplot.tiff"), plot = plot2, width = 8, height = ngene,compression='lzw',limitsize=F)
		ggsave(paste0("fig",plotid,".Genes_pseudotimeplot.pdf"), plot = plot2, width = 8, height = ngene,limitsize=F)
		if(! is.null(time_diff))
		{
			Time_diff<-read.table(time_diff,header=T)
			sig_gene_names <- Time_diff$gene_short_name[Time_diff$qval<0.05]
			markers_used<-markers_used[markers_used %in% sig_gene_names]
			ngene<-length(markers_used)
			print(ngene)
			if(ngene)
			{
					cds_subset <- cds[markers_used,]
					plot2<-NULL
					plot2 <- plot_genes_in_pseudotime(cds_subset, color_by = group.by)
					ggsave(paste0("fig",plotid,"b.Genes_pseudotimeplot.tiff"), plot = plot2, width = 8, height = ngene,compression='lzw',limitsize=F)
					ggsave(paste0("fig",plotid,"b.Genes_pseudotimeplot.pdf"), plot = plot2, width = 8, height = ngene,limitsize=F)
			}
		}
	}
}


monocle2_markers_pseudotime_plot2<- function(cds,markers,group.by='celltype',plotid='04',time_diff='fig05.Fibroblasts.Time_diff_all.tsv'){
library(monocle)
library(ggsci)
  cds_genes<-rownames(cds)
  markers_used<-markers[markers %in% cds_genes]
  ngene<-length(markers_used)
  print(ngene)
  if(ngene)
  {
		cds_subset <- cds[markers_used,]
		pdf(file=paste0('fig',plotid,'g.',projectid,".pseudotime_heatmap.pdf"),width=10,heigh=4/20*ngene)
			plot_pseudotime_heatmap(cds_subset, num_clusters=4, show_rownames=T, return_heatmap=T)
		dev.off()
		
		if(! is.null(time_diff))
		{
			Time_diff<-read.table(time_diff,header=T)
			sig_gene_names <- Time_diff$gene_short_name[Time_diff$qval<0.05]
			markers_used<-markers_used[markers_used %in% sig_gene_names]
			ngene<-length(markers_used)
			print(ngene)
			if(ngene>0)
			{
					cds_subset <- cds[markers_used,]
					pearplot<-NULL
					plot1<-NULL
					plot2<-NULL
					plot3<-NULL
					plot1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
					plot2 <- plot_genes_in_pseudotime(cds_subset, color_by = group.by)
					plot3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
					pearplot <-CombinePlots(plots = list(plot1, plot2,plot3),ncol=3)
					ggsave(paste0('fig',plotid,'h.',projectid,".Genes_pseudotimeplot.tiff"), plot = pearplot, width = 18, height = 8,compression='lzw')
					ggsave(paste0('fig',plotid,'h.',projectid,".Genes_pseudotimeplot.pdf"), plot = pearplot, width = 18, height = 8)

					pearplot<-NULL
					plot1<-NULL
					plot2<-NULL
					plot3<-NULL
					plot1 <- plot_genes_jitter(cds_subset, grouping = "State", color_by = "State")
					plot2 <- plot_genes_violin(cds_subset, grouping = "State", color_by = "State")
					plot3 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
					pearplot <-CombinePlots(plots = list(plot1, plot2,plot3),ncol=3)
					ggsave(paste0('fig',plotid,'i.',projectid,".Genes_Jitterplot.tiff"), plot = pearplot, width = 18, height = 8,compression='lzw')
					ggsave(paste0('fig',plotid,'i.',projectid,".Genes_Jitterplot.pdf"), plot = pearplot, width = 18, height = 8)
					for(gene_used in markers_used)
					{
						name<-gsub('-','_',gene_used)
						pData(cds)[[name]] = log2( exprs(cds)[gene_used,]+1)
						plot1<-NULL
						plot1<-plot_cell_trajectory(cds, color_by = name)+ scale_color_gsea()+labs(color=gene_used)
						ggsave(paste0('fig',plotid,'j.',projectid,'.trajectory.',name,'.tiff'), plot = plot1, width = 8, height = 8,compression='lzw')
						ggsave(paste0('fig',plotid,'j.',projectid,'.trajectory.',name,'.pdf'), plot = plot1, width = 8, height = 8)
					}
			}
		}
	}
}


#gsva_yingbai(tumor_expr,phenotype,group.by='URG_cluster',plotid='14',species = "human")
gsva_yingbai <- function(expr,phenotype,group.by='group',plotid='14',species = "human"){
      library(GSVA)
      library(msigdbr)
      #BiocManager::install('msigdbr')
      library(pheatmap)
      
      if(species=='human')
      {
        DBkeyset='org.Hs.eg.db'
        org='Homo sapiens'
        category='H'
        #library(org.Hs.eg.db)
        
      }else if(species=='mouse')
      {
        category='H'
        #library(org.Mm.eg.db)
        DBkeyset='org.Mm.eg.db'
        org='Mus musculus'
      }else if(species=='rat')
      {
        #library(org.Rn.eg.db)
        #BiocManager::install('org.Rn.eg.db')
        DBkeyset='org.Rn.eg.db'
        org='Mus musculus'
        category='H'
      }else{
        DBkeyset='org.Hs.eg.db'
        #library(org.Hs.eg.db)
        org='Homo sapiens'
        category='H'
      }
      
      expr<-as.matrix(expr)
      genesets <- msigdbr(species = org, category = category) 
      genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
      genesets <- split(genesets$gene_symbol, genesets$gs_name)
      
      gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
      gsva.df <- data.frame(gsva.res, check.names = F)
      gsva.df[1:3,1:3]
      write.table(gsva.df, paste0("fig",plotid,"a.",group.by,".gsva.hallmark.tsv"),sep='\t',row.names = T,quote = FALSE)
      
      
      ggheatmap_yingbio_2groups_for_gsva(gsva.df,phenotype,geneid='genesets',pheight=6,
                                output=paste0('fig',plotid,'a.gsva.hallmark'),save.data=T)

      
      genesets <- msigdbr(species = org, category = "C2") 
      genesets <- subset(genesets, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
      genesets <- split(genesets$gene_symbol, genesets$gs_name)
      
      gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
      gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
      write.table(gsva.df, paste0("fig",plotid,"b.",group.by,".gsva.kegg.tsv"),sep='\t',row.names = F,quote = FALSE)
      
      if(0)
      {
        gsva.df<-read.delim(paste0("fig",plotid,"b.",group.by,".gsva.kegg.tsv"),row.names = 1)
        gsva.df[1:3,1:3]
        
      }
      ggheatmap_yingbio_2groups_for_gsva(gsva.df,phenotype,geneid='genesets',pheight=16,
                                output=paste0('fig',plotid,'b.gsva.kegg'),save.data=T)
    }


###  df_gsva：行是样本名，列是kegg通路。列名最好先去除“KEGG_”前缀
###  group：样本分组的vector，需要是factor类型
gsva_t_plot <- function(gsva.res=gsva.res,group=df_sample$group,group_test='group1',group_control='group2',ntop=NA,output='output',pw=8,ph=10){
        ## limma差异通路分析
        
        #BiocManager::install('limma')
        library(limma)
        library(stringr)
        library(ggplot2)
        library(ggthemes)
        library(ggprism)
        library(dplyr)
        # 设置或导入分组
        design <- model.matrix(~0+group)
        colnames(design) = levels(factor(group))
        rownames(design) = rownames(gsva.res)
        design
        # Tunor VS Normal
        compareName <- paste0(group_test, "-", group_control)
        compare <- makeContrasts(contrasts = compareName, levels=design)
        fit <- lmFit(t(gsva.res), design)
        fit2 <- contrasts.fit(fit, compare)
        fit3 <- eBayes(fit2)
        Diff <- topTable(fit3, coef=1, number=200)
        head(Diff)
        
        ## 发散条形图绘制
        ## barplot
        dat_plot <- data.frame(id = row.names(Diff),
                               t = Diff$t)
        # 去掉"HALLMARK_"
        
        dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
        # 新增一列 根据t阈值分类
        dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
        # 排序
        dat_plot <- dat_plot %>% arrange(t)
        
        write.table(dat_plot,file=paste0(output,".tsv"),sep="\t",row.names=F,quote=F)
        
        if(is.na(ntop))
        {
          dat_plot_used=dat_plot
        }else
        {
          dat_plot_used1<-head(dat_plot,ntop)
          dat_plot_used2<-tail(dat_plot,ntop)			
          dat_plot_used<-rbind(dat_plot_used1,dat_plot_used2)
          dat_plot_used<-unique(dat_plot_used)
        }
        
        # 变成因子类型
        dat_plot_used$id <- factor(dat_plot_used$id,levels = dat_plot_used$id)
        # 绘制
        ##install.packages("ggthemes")
        #install.packages("ggprism")
        p <- ggplot(data = dat_plot_used,aes(x = id,y = t,fill = threshold)) +
          geom_col()+
          coord_flip() +
          scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
          geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
          xlab('') + 
          ylab(paste("t value of GSVA score", group_test, "vs", group_control)) + #注意坐标轴旋转了
          guides(fill=F)+ # 不显示图例
          theme_prism(border = T) +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          )
        p
        # 添加标签  
        # 小于-1的数量
        low1 <- dat_plot_used %>% filter(t < -1) %>% nrow()
        # 小于0总数量
        low0 <- dat_plot_used %>% filter( t < 0) %>% nrow()
        # 小于1总数量
        high0 <- dat_plot_used %>% filter(t < 1) %>% nrow()
        # 总的柱子数量
        high1 <- nrow(dat_plot_used)
        
        # 依次从下到上添加标签
        p <- p + geom_text(data = dat_plot_used[1:low1,],aes(x = id,y = 0.1,label = id),
                           hjust = 0,color = 'black')  # 小于-1的为黑色标签
        if(low0>low1)
        {
          p <- p + geom_text(data = dat_plot_used[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
                    hjust = 0,color = 'grey')  # 灰色标签
        }
        if(high0>low0)
        {
          p <- p + geom_text(data = dat_plot_used[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
                    hjust = 1,color = 'grey')  # 灰色标签
        }
        p <- p + geom_text(data = dat_plot_used[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
                    hjust = 1,color = 'black') # 大于1的为黑色标签  
        ggsave(filename = paste0(output,'.pdf'),p,width = pw,height  = ph)
        ggsave(filename = paste0(output,'.tiff'), plot = p, width = pw, height = ph,compression='lzw') 
      }
      


gsva_yingbai_custom<- function(expr,phenotype,genesets=gsc,group.by='group',title='Value',plotid='14',species = "human"){
      library(GSVA)
      library(msigdbr)
      #BiocManager::install('msigdbr')
      library(pheatmap)
      library(ggplot2)
      library(reshape2)
      
      if(species=='human')
      {
        DBkeyset='org.Hs.eg.db'
        org='Homo sapiens'
        category='H'
        #library(org.Hs.eg.db)
        
      }else if(species=='mouse')
      {
        category='H'
        #library(org.Mm.eg.db)
        DBkeyset='org.Mm.eg.db'
        org='Mus musculus'
      }else if(species=='rat')
      {
        #library(org.Rn.eg.db)
        #BiocManager::install('org.Rn.eg.db')
        DBkeyset='org.Rn.eg.db'
        org='Mus musculus'
        category='H'
      }else{
        DBkeyset='org.Hs.eg.db'
        #library(org.Hs.eg.db)
        org='Homo sapiens'
        category='H'
      }
      
      expr<-as.matrix(expr)      
      gsvaPar <- ssgseaParam(expr, genesets,minSize=2)
		gsva.res <- gsva(gsvaPar, verbose=FALSE)
      # gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
      gsva.df <- data.frame(gsva.res, check.names = F)
      gsva.df[1:3,1:3]
      write.table(gsva.df, paste0("fig",plotid,"a.",group.by,".gsva.tsv"),sep='\t',row.names = T,quote = FALSE)     
      
      ggheatmap_yingbio_2groups_for_gsva(gsva.df,phenotype,geneid='genesets',pheight=6,
                                output=paste0('fig',plotid,'a.gsva.heatmap'),save.data=T)
	  df_gsva<-as.matrix(gsva.df)
	  df_plot<-melt(df_gsva)
	  head(df_plot)
	  colnames(df_plot)<-c('celltype','sample','score')
	  df_plot_used<-merge(df_plot,phenotype,by.x='sample',by.y=0)
	  head(df_plot_used)
	  ggplot2_boxplot_gsva(df_plot_used,groupby='celltype',score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'b.',projectid,".gsva.boxpot")) 
	  ggplot2_boxplot_paired(df_plot_used,groupby=c('group','celltype'),score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'c.',projectid,".gsva.boxplot_paired"))
}


ggheatmap_yingbio_2groups_for_gsva<-function(mydata,phenotype,geneid='miRNA_ID',pheight=NA,output='heatmap',save.data=F){
  phenotype<-phenotype[order(phenotype[,1]),,drop=F]
  mydata<-mydata[,rownames(phenotype)]
  if(save.data)
  {
    mydata_out<-cbind(rownames(mydata),mydata)
    colnames(mydata_out)[1]<-geneid
    write.table(mydata_out,file = paste0(output,'.tsv'),row.names = FALSE,sep='\t',quote = FALSE)
  }
  ngene<-dim(mydata)[1]
  if(!is.na(pheight))
  {
    ph<-pheight
  }else{
    ph<-max(6,ngene/10)
  }
  print(paste0('pic height',ph))
  ngene<-nrow(mydata)
  cluster_rows=F
  if(ngene>=2)
  {
	cluster_rows=T
	}
  library(pheatmap)
  plot.pheatmap<-NULL
  plot.pheatmap<-pheatmap(mydata, 
                          annotation=phenotype, 
                          color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
                          #color = viridis(8, option = "G")
                          cluster_cols =F,
                          cluster_rows =cluster_rows,
                          show_colnames = F,
                          show_rownames = T,
                          scale="row",
                          fontsize = 10,
                          fontsize_row=5,
                          fontsize_col=10)
  dev.off()
  ####  如果右侧的legend空间不够，增加legend空间。
  ####  plot.pheatmap$gtable$widths[6] <- plot.pheatmap$gtable$widths[6] + unit(30, "bigpts")
  pdf(file=paste0(output,'.pdf'), height=ph, width=10)
  print(plot.pheatmap)
  dev.off()
  tiff(file=paste0(output,'.tiff'), height=ph*300, width=10*300,units="px",res=300,compression='lzw')
  print(plot.pheatmap)
  dev.off()
}


#seurat_gsva(pbmc.x,group.by='seurat_clusters',projectid=celltypei,plotid='07')
seurat_gsva <- function(pbmc,group.by='seurat_clusters',projectid='xxxxxx',plotid='14',species = "human"){
		library(GSVA)
		library(msigdbr)
		#BiocManager::install('msigdbr')
		library(pheatmap)
		
		if(species=='human')
		{
			   DBkeyset='org.Hs.eg.db'
			   org='Homo sapiens'
			   category='H'
			   #library(org.Hs.eg.db)
			   
		}else if(species=='mouse')
		{
				category='H'
				#library(org.Mm.eg.db)
				 DBkeyset='org.Mm.eg.db'
				 org='Mus musculus'
		}else if(species=='rat')
		{
				#library(org.Rn.eg.db)
				#BiocManager::install('org.Rn.eg.db')
				DBkeyset='org.Rn.eg.db'
				org='Mus musculus'
				category='H'
		}else{
				DBkeyset='org.Hs.eg.db'
			    #library(org.Hs.eg.db)
			    org='Homo sapiens'
			    category='H'
		}
		
		Idents(pbmc) <- group.by
		expr <- AverageExpression(pbmc, assays = "RNA", slot = "data")[[1]]
		expr <- expr[rowSums(expr)>0,]  #选取非零基因
		expr <- as.matrix(expr)
		expr_out<-cbind(rownames(expr),expr)
		colnames(expr_out)[1]<-'gene'
		write.table(expr_out, paste0("fig",plotid,"a.",projectid,'.',group.by,".AverageExpression.tsv"),sep='\t',row.names = F,quote = FALSE)
		
		genesets <- msigdbr(species = org, category = category) 
		genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
		genesets <- split(genesets$gene_symbol, genesets$gs_name)

		gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
		gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
		write.table(gsva.df, paste0("fig",plotid,"b.",projectid,'.',group.by,".gsva.hallmark.tsv"),sep='\t',row.names = F,quote = FALSE)

		#gsva.res<-gsva.res[,c('3','1','0','2')]

		pdf(file=paste0("fig",plotid,"b.",projectid,'.',group.by,".gsva.hallmark.Heatmap.pdf"),width=12,height=9)
			print(pheatmap(gsva.res, show_colnames = T, scale = "row",cluster_cols = F))
		dev.off()

		tiff(file=paste0("fig",plotid,"b.",projectid,'.',group.by,".gsva.hallmark.Heatmap.tiff"),width=3600,heigh=2700,units='px',res=300,compression='lzw')
			print(pheatmap(gsva.res, show_colnames = T, scale = "row",cluster_cols = F))
		dev.off()

		tiff(file=paste0("fig",plotid,"b.",projectid,'.',group.by,".gsva.hallmark.Heatmap.nonscale.tiff"),width=3600,heigh=2700,units='px',res=300,compression='lzw')
			print(pheatmap(gsva.res, show_colnames = T))
		dev.off()

		genesets <- msigdbr(species = org, category = "C2") 
		genesets <- subset(genesets, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
		genesets <- split(genesets$gene_symbol, genesets$gs_name)

		gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
		gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
		write.table(gsva.df, paste0("fig",plotid,"c.",projectid,'.',group.by,".gsva.kegg.tsv"),sep='\t',row.names = F,quote = FALSE)


		pdf(file=paste0("fig",plotid,"c.",projectid,'.',group.by,".gsva.kegg.Heatmap.pdf"),width=12,height=24)
			print(pheatmap(gsva.res, show_colnames = T, scale = "row"))
		dev.off()
		tiff(file=paste0("fig",plotid,"c.",projectid,'.',group.by,".gsva.kegg.Heatmap.tiff"),width=3600,heigh=7200,units='px',res=300,compression='lzw')
			print(pheatmap(gsva.res, show_colnames = T, scale = "row"))
		dev.off()
		tiff(file=paste0("fig",plotid,"c.",projectid,'.',group.by,".gsva.kegg.Heatmap.nonscale.tiff"),width=3600,heigh=7200,units='px',res=300,compression='lzw')
			print(pheatmap(gsva.res, show_colnames = T))
		dev.off()
}



seurat_gsva_bp <- function(pbmc,group.by='seurat_clusters',projectid='xxxxxx',plotid='14',species = "human"){
		library(GSVA)
		library(msigdbr)
		#BiocManager::install('msigdbr')
		library(pheatmap)
		library(dplyr)
		
		if(species=='human')
		{
			   DBkeyset='org.Hs.eg.db'
			   org='Homo sapiens'
			   #library(org.Hs.eg.db)
			   
		}else if(species=='mouse')
		{
				#library(org.Mm.eg.db)
				 DBkeyset='org.Mm.eg.db'
				 org='Mus musculus'
		}else if(species=='rat')
		{
				#library(org.Rn.eg.db)
				#BiocManager::install('org.Rn.eg.db')
				DBkeyset='org.Rn.eg.db'
				org='Mus musculus'
		}else{
				DBkeyset='org.Hs.eg.db'
			    #library(org.Hs.eg.db)
			    org='Homo sapiens'
		}
		
		Idents(pbmc) <- group.by
		expr <- AverageExpression(pbmc, assays = "RNA", slot = "data")[[1]]
		expr <- expr[rowSums(expr)>0,]  #选取非零基因
		expr <- as.matrix(expr)
		expr_out<-cbind(rownames(expr),expr)
		colnames(expr_out)[1]<-'gene'
		write.table(expr_out, paste0("fig",plotid,"a.",projectid,'.',group.by,".AverageExpression.tsv"),sep='\t',row.names = F)
		
		
		

		genesets <- msigdbr(species = org, category = "C5",subcategory = "GO:BP") 
		genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
		genesets <- split(genesets$gene_symbol, genesets$gs_name)

		gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
		gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
		write.table(gsva.df, paste0("fig",plotid,"d.",projectid,'.',group.by,".gsva.GO.BP.tsv"),sep='\t',row.names = F)
        
        gsva.df<-gsva.df[,-1]
        nlevels<-dim(gsva.df)[2]
        np<-10
        if(nlevels>10)
        {
			np<-5
		}		
		paths<-as.character()
		for(i in 1:ncol(gsva.df))
		{
			top5.celltype <- gsva.df %>% top_n(n = np, wt = gsva.df[[i]])
			paths<-c(paths,rownames(top5.celltype))
			top5.celltype <- gsva.df %>% top_n(n = -np, wt = gsva.df[[i]])
			paths<-c(paths,rownames(top5.celltype))
		}
        paths<-unique(paths)
        
		gsva.res.pheatmap<-gsva.df[paths,]
		rownames(gsva.res.pheatmap)<-gsub('^GOBP_','',rownames(gsva.res.pheatmap))
		gsva.res.pheatmap<-as.matrix(gsva.res.pheatmap)		
		head(gsva.res.pheatmap)
		rownames(gsva.res.pheatmap)<-substring(rownames(gsva.res.pheatmap),1,30)
		
		
		
		ngene<-nrow(gsva.res.pheatmap)
		ntype<-ncol(gsva.res.pheatmap)
		
		pdf(file=paste0("fig",plotid,"d.",projectid,'.',group.by,".gsva.GO.BP.Heatmap.pdf"),width=ntype*2+3,height=ngene*16/100)
			pheatmap(gsva.res.pheatmap, show_colnames = T, scale = "row",cluster_cols = F,fontsize_row=5)
		dev.off()

		tiff(file=paste0("fig",plotid,"d.",projectid,'.',group.by,".gsva.GO.BP.Heatmap.tiff"),width=ntype*2*300+900,heigh=ngene*16/100*300,units='px',res=300,compression='lzw')
			pheatmap(gsva.res.pheatmap, show_colnames = T, scale = "row",cluster_cols = F,fontsize_row=5)
		dev.off()

		tiff(file=paste0("fig",plotid,"d.",projectid,'.',group.by,".gsva.GO.BP.Heatmap.nonscale.tiff"),width=ntype*2*300+900,heigh=ngene*16/100*300,units='px',res=300,compression='lzw')
			pheatmap(gsva.res.pheatmap, show_colnames = T,fontsize_row=5)
		dev.off()
############################################################
		genesets <- msigdbr(species = org, category = "C5",subcategory = "GO:MF") 
		genesets <- subset(genesets, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
		genesets <- split(genesets$gene_symbol, genesets$gs_name)

		gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
		gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
		write.table(gsva.df, paste0("fig",plotid,"e.",projectid,'.',group.by,".gsva.GO.MF.tsv"),sep='\t',row.names = F)

		gsva.df<-gsva.df[,-1]
        nlevels<-dim(gsva.df)[2]
        np<-10
        if(nlevels>10)
        {
			np<-5
		}		
		paths<-as.character()
		for(i in 1:ncol(gsva.df))
		{
			top5.celltype <- gsva.df %>% top_n(n = np, wt = gsva.df[[i]])
			paths<-c(paths,rownames(top5.celltype))
			top5.celltype <- gsva.df %>% top_n(n = -np, wt = gsva.df[[i]])
			paths<-c(paths,rownames(top5.celltype))
		}
        paths<-unique(paths)
        
		gsva.res.pheatmap<-gsva.df[paths,]
		rownames(gsva.res.pheatmap)<-gsub('^GOMF_','',rownames(gsva.res.pheatmap))
		gsva.res.pheatmap<-as.matrix(gsva.res.pheatmap)		
		head(gsva.res.pheatmap)
		rownames(gsva.res.pheatmap)<-substring(rownames(gsva.res.pheatmap),1,30)
		#rownames(gsva.res.pheatmap)<-NULL
		head(gsva.res.pheatmap)

		ngene<-nrow(gsva.res.pheatmap)
		ntype<-ncol(gsva.res.pheatmap)
		
		pdf(file=paste0("fig",plotid,"e.",projectid,'.',group.by,".gsva.GO.MF.Heatmap.pdf"),width=ntype*2+3,height=ngene*16/100)
			pheatmap(gsva.res.pheatmap, show_colnames = T, scale = "row",fontsize_row=3)
		dev.off()
		tiff(file=paste0("fig",plotid,"e.",projectid,'.',group.by,".gsva.GO.MF.Heatmap.tiff"),width=ntype*2*300+900,heigh=ngene*16/100*300,units='px',res=300,compression='lzw')
			pheatmap(gsva.res.pheatmap, show_colnames = T, scale = "row",fontsize_row=3)
		dev.off()
		tiff(file=paste0("fig",plotid,"e.",projectid,'.',group.by,".gsva.GO.MF.Heatmap.nonscale.tiff"),width=ntype*2*300+900,heigh=ngene*16/100*300,units='px',res=300,compression='lzw')
			pheatmap(gsva.res.pheatmap, show_colnames = T,fontsize_row=3)
		dev.off()
}

seurat_gsva_celltype <- function(pbmc,group.by='seurat_clusters',projectid='xxxxxx',plotid='14',species = "human"){
		library(GSVA)
		library(msigdbr)
		#BiocManager::install('msigdbr')
		library(pheatmap)
		
		if(species=='human')
		{
			   DBkeyset='org.Hs.eg.db'
			   org='Homo sapiens'
			   #library(org.Hs.eg.db)
			   
		}else if(species=='mouse')
		{
				#library(org.Mm.eg.db)
				 DBkeyset='org.Mm.eg.db'
				 org='Mus musculus'
		}else if(species=='rat')
		{
				#library(org.Rn.eg.db)
				#BiocManager::install('org.Rn.eg.db')
				DBkeyset='org.Rn.eg.db'
				org='Mus musculus'
		}else{
				DBkeyset='org.Hs.eg.db'
			    #library(org.Hs.eg.db)
			    org='Homo sapiens'
		}
		
		Idents(pbmc) <- group.by
		expr <- AverageExpression(pbmc, assays = "RNA", slot = "data")[[1]]
		expr <- expr[rowSums(expr)>0,]  #选取非零基因
		expr <- as.matrix(expr)
		expr_out<-cbind(rownames(expr),expr)
		colnames(expr_out)[1]<-'gene'
		write.table(expr_out, paste0("fig",plotid,"a.",projectid,'.',group.by,".AverageExpression.tsv"),sep='\t',row.names = F)
		
		
		

		genesets <- msigdbr(species = org, category = "C7",subcategory = "IMMUNESIGDB") 
		genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
		genesets <- split(genesets$gene_symbol, genesets$gs_name)

		gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
		gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
		write.table(gsva.df, paste0("fig",plotid,"f.",projectid,'.',group.by,".gsva.IMMUNESIGDB.tsv"),sep='\t',row.names = F)
		gsva.df<-gsva.df[,-1]
        nlevels<-dim(gsva.df)[2]
        np<-10
        if(nlevels>10)
        {
			np<-5
		}		
		paths<-as.character()
		for(i in 1:ncol(gsva.df))
		{
			top5.celltype <- gsva.df %>% top_n(n = np, wt = gsva.df[[i]])
			paths<-c(paths,rownames(top5.celltype))
			top5.celltype <- gsva.df %>% top_n(n = -np, wt = gsva.df[[i]])
			paths<-c(paths,rownames(top5.celltype))
		}
        paths<-unique(paths)
        
		gsva.res.pheatmap<-gsva.df[paths,]
		#rownames(gsva.res.pheatmap)<-gsub('^GOBP_','',rownames(gsva.res.pheatmap))
		gsva.res.pheatmap<-as.matrix(gsva.res.pheatmap)		
		head(gsva.res.pheatmap)
		#rownames(gsva.res.pheatmap)<-substring(rownames(gsva.res.pheatmap),1,30)
		
		ngene<-nrow(gsva.res.pheatmap)
		ntype<-ncol(gsva.res.pheatmap)

		pdf(file=paste0("fig",plotid,"f.",projectid,'.',group.by,".gsva.IMMUNESIGDB.Heatmap.pdf"),width=ntype*2+3,height=ngene*16/100)
			pheatmap(gsva.res.pheatmap, show_colnames = T, scale = "row",cluster_cols = F,fontsize_row=5)
		dev.off()

		tiff(file=paste0("fig",plotid,"f.",projectid,'.',group.by,".gsva.IMMUNESIGDB.Heatmap.tiff"),width=ntype*2*300+900,heigh=ngene*16/100*300,units='px',res=300,compression='lzw')
			pheatmap(gsva.res.pheatmap, show_colnames = T, scale = "row",cluster_cols = F,fontsize_row=5)
		dev.off()

		tiff(file=paste0("fig",plotid,"f.",projectid,'.',group.by,".gsva.IMMUNESIGDB.Heatmap.nonscale.tiff"),width=ntype*2*300+900,heigh=ngene*16/100*300,units='px',res=300,compression='lzw')
			pheatmap(gsva.res.pheatmap, show_colnames = T,fontsize_row=5)
		dev.off()
############################################################
		genesets <- msigdbr(species = org, category = "C8") 
		genesets <- subset(genesets, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
		genesets <- split(genesets$gene_symbol, genesets$gs_name)

		gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
		gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
		write.table(gsva.df, paste0("fig",plotid,"g.",projectid,'.',group.by,".gsva.cell.type.tsv"),sep='\t',row.names = F,quote=F)

		gsva.res.pheatmap<-gsva.df
		gsva.res.pheatmap<-gsva.res.pheatmap[,-1]
		#rownames(gsva.res.pheatmap)<-gsub('^GOMF_','',rownames(gsva.res.pheatmap))
		gsva.res.pheatmap<-as.matrix(gsva.res.pheatmap)		
		head(gsva.res.pheatmap)
		#rownames(gsva.res.pheatmap)<-substring(rownames(gsva.res.pheatmap),1,30)
		#rownames(gsva.res.pheatmap)<-NULL
		head(gsva.res.pheatmap)
		ngene<-nrow(gsva.res.pheatmap)
		ntype<-ncol(gsva.res.pheatmap)

		pdf(file=paste0("fig",plotid,"g.",projectid,'.',group.by,".gsva.cell.type.Heatmap.pdf"),width=ntype*2+3,height=36)
			pheatmap(gsva.res.pheatmap, show_colnames = T, scale = "row",fontsize_row=3)
		dev.off()
		tiff(file=paste0("fig",plotid,"g.",projectid,'.',group.by,".gsva.cell.type.Heatmap.tiff"),width=ntype*2*300+900,heigh=10800,units='px',res=300,compression='lzw')
			pheatmap(gsva.res.pheatmap, show_colnames = T, scale = "row",fontsize_row=3)
		dev.off()
		tiff(file=paste0("fig",plotid,"g.",projectid,'.',group.by,".gsva.cell.type.Heatmap.nonscale.tiff"),width=ntype*2*300+900,heigh=10800,units='px',res=300,compression='lzw')
			pheatmap(gsva.res.pheatmap, show_colnames = T,fontsize_row=3)
		dev.off()
}


get_M1M2_marker<-function(species='human')
{
Macro_M1<-c('IL23','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','IDO1','KYNU','CCR7')#PMID： 33545035
Macro_M2<-c('IL4R','CCL4','CCL13','CCL20','CCL17','CCL18','CCL22','CCL24','LYVE1','VEGFA','VEGFB','VEGFC','VEGFD','EGF','CTSA','CTSB',
'CTSC','CTSD','TGFB1','TGFB2','TGFB3','MMP14','MMP19','MMP9','CLEC7A','WNT7B','FASL','TNFSF12','TNFSF8','CD276','VTCN1','MSR1','FN1','IRF4')#PMID： 33545035
proliferating<-c('CDK1','FEN1','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','MKI67','PCNA','PLK1','POLA1','PRIM2','RFC2','RFC3','RFC4','RFC5','RPA1','RPA2',
'RPA3','RRM1','RRM2','SMC1A','SMC3','STAG2')
proinflammatory<-c('CXCL10','CXCL9','IFNG','GZMB','CXCL13','STAT1','IRF1','CCL5','GNLY','TBX21','CD8B','PRF1','IL12A','IL12B','CD19')
Angiogenesis<-c("CCND2", "CCNE1", "CD44", "CXCR4", "E2F3", "EDN1", "EZH2", "FGF18", "FGFR1", "FYN", "HEY1", "ITGAV", "JAG1", "JAG2", "MMP9", "NOTCH1", 
                "PDGFA", "PTK2", "SPP1", "STC1", "TNFAIP6", "TYMP", "VAV2", "VCAN", "VEGFA")#PMID： 33545035
Phagocytosis<-c("MRC1","CD163","MERTK","C1QB")#PMID： 33545035


if(species =='mouse')
{
Macro_M1<-gene_human2mouse_vector(Macro_M1)
Macro_M2<-gene_human2mouse_vector(Macro_M2)
proliferating<-gene_human2mouse_vector(proliferating)
proinflammatory<-gene_human2mouse_vector(proinflammatory)
Angiogenesis<-gene_human2mouse_vector(Angiogenesis)
Phagocytosis<-gene_human2mouse_vector(Phagocytosis)
}

df_geneset<-data.frame(geneset=character(),gene=character())
	head(df_geneset)
	dim(df_geneset)
	colnames(df_geneset)
	
for(i in 1:length(Macro_M1))
{
	df_geneset<-rbind(df_geneset,c('M1',Macro_M1[i]))
}
for(i in 1:length(Macro_M2))
{
	df_geneset<-rbind(df_geneset,c('M2',Macro_M2[i]))
}
for(i in 1:length(proliferating))
{
	df_geneset<-rbind(df_geneset,c('Proliferating',proliferating[i]))
}
for(i in 1:length(proinflammatory))
{
	df_geneset<-rbind(df_geneset,c('Pro-inflammatory',proinflammatory[i]))
}


for(i in 1:length(Angiogenesis))
{
	df_geneset<-rbind(df_geneset,c('Angiogenesis',Angiogenesis[i]))
}
for(i in 1:length(Phagocytosis))
{
	df_geneset<-rbind(df_geneset,c('Phagocytosis',Phagocytosis[i]))
}

colnames(df_geneset)<-c('geneset','gene')
return(df_geneset)
}


#library(GSEABase)
#file="M1M2.gmt"
#gmt <- GSEABase::getGmt(file)
#seurat_gsva_custom(pbmc.x,genesets=gmt,group.by='celltype',projectid=celltypei,plotid='08')
seurat_gsva_custom <- function(pbmc,genesets=gmt,group.by='seurat_clusters',projectid='xxxxxx',plotid='14',species = "human",pw=12,ph=9){
		library(GSVA)
		library(msigdbr)
		#BiocManager::install('msigdbr')
		library(pheatmap)
		
		if(species=='human')
		{
			   DBkeyset='org.Hs.eg.db'
			   org='Homo sapiens'
			   #library(org.Hs.eg.db)
			   
		}else if(species=='mouse')
		{
				#library(org.Mm.eg.db)
				 DBkeyset='org.Mm.eg.db'
				 org='Mus musculus'
		}else if(species=='rat')
		{
				#library(org.Rn.eg.db)
				#BiocManager::install('org.Rn.eg.db')
				DBkeyset='org.Rn.eg.db'
				org='Mus musculus'
		}else{
				DBkeyset='org.Hs.eg.db'
			    #library(org.Hs.eg.db)
			    org='Homo sapiens'
		}
		
		Idents(pbmc) <- group.by
		expr <- AverageExpression(pbmc, assays = "RNA", slot = "data")[[1]]
		expr <- expr[rowSums(expr)>0,]  #选取非零基因
		expr <- as.matrix(expr)
		expr_out<-cbind(rownames(expr),expr)
		colnames(expr_out)[1]<-'gene'
		write.table(expr_out, paste0("fig",plotid,".",projectid,'.',group.by,".AverageExpression.tsv"),sep='\t',row.names = F,quote=F)
		
		gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
		gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
		write.table(gsva.df, paste0("fig",plotid,".",projectid,'.',group.by,".gsva.tsv"),sep='\t',row.names = F,quote=F)
		
		ngene<-nrow(gsva.res)
		ntype<-ncol(gsva.res)

		pdf(file=paste0("fig",plotid,".",projectid,'.',group.by,".gsva.Heatmap.pdf"),width=pw,height=ph)
			pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row",cluster_cols = F)
		dev.off()

		tiff(file=paste0("fig",plotid,".",projectid,'.',group.by,".gsva.Heatmap.tiff"),width=pw*300,heigh=ph*300,units='px',res=300,compression='lzw')
			pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row",cluster_cols = F,angle_col=45)
		dev.off()
		
		pdf(file=paste0("fig",plotid,".",projectid,'.',group.by,".gsva.Heatmap.nonscale.pdf"),width=pw,heigh=ph)
			pheatmap::pheatmap(gsva.res, show_colnames = T)
		dev.off()
		
		tiff(file=paste0("fig",plotid,".",projectid,'.',group.by,".gsva.Heatmap.nonscale.tiff"),width=pw*300,heigh=ph*300,units='px',res=300,compression='lzw')
			pheatmap::pheatmap(gsva.res, show_colnames = T)
		dev.off()

}

#library(GSEABase)	
#	file_hubgene<-'8_Mm_vs_Gs.hubgene_all.txt'
#	df_hubgene <- read.delim(file = file_hubgene,sep='\t',header=F)
#	hubgene<-df_hubgene[,1]
#	list_genesets<-list()
#	list_genesets[[1]]<-hubgene
#	names(list_genesets)[1]<-'ferr_hub_gene'
#	
#	  gsc <- GeneSetCollection(
#   mapply(function(geneIds, keggId) {
#   GeneSet(geneIds, geneIdType=SymbolIdentifier(),
#            collectionType=NullCollection(),
#            setName=keggId)
#    }, list_genesets, names(list_genesets)))
#seurat_gsva_custom_violin(pbmc,genesets=gsc,group.by='celltype',projectid='xxxxxx',plotid='14',species = "human",title="Ferroptosis score")


#m6a_gene<-c('METTL3','METTL14','METTL16','YTHDF1','YTHDF2','YTHDF3','YTHDC1','YTHDC2','RBM15','RBM15B','RBMX','IGF2BP1','IGF2BP2','IGF2BP3','KIAA1429','FMR1','LRPPRC','HNRNPA2B1','HNRNPC','ZC3H13','FTO','ALKBH5','WTAP')
#gs<-genelist_2_genesets(genelist=m6a_gene,geneset_name='m6A_gene',species = "human")
genelist_2_genesets <- function(genelist=m6a_gene,geneset_name='m6A_gene')
{
library(GSEABase)	
list_genesets<-list()
list_genesets[[1]]<-genelist
names(list_genesets)[1]<-geneset_name
	  gsc <- GeneSetCollection(
   mapply(function(geneIds, keggId) {
   GeneSet(geneIds, geneIdType=SymbolIdentifier(),
            collectionType=NullCollection(),
            setName=keggId)
    }, list_genesets, names(list_genesets)))
    return(gsc)
}

###   data.frame含两列，第一列geneset名称，第二列，基因名。
###   gsc<-df_2_genesets(df_geneset)
###   file_genesets<-'genesets.tsv'
###   df_geneset<-read.delim(file_genesets)
###   head(df_geneset)
df_2_genesets <- function(df_geneset,species = "human")
{
  library(GSEABase)
  list_genesets<-list()
  genesets<-unique(df_geneset[,1])
  for(i in 1:length(genesets))
  {
    geneset_name<-genesets[i]		
    list_genesets[[i]]<-df_geneset[df_geneset[,1]==geneset_name,2]
    names(list_genesets)[i]<-geneset_name
  }
  gsc <- GeneSetCollection(
    mapply(function(geneIds, keggId) {
      GeneSet(geneIds, geneIdType=SymbolIdentifier(),
              collectionType=NullCollection(),
              setName=keggId)
    }, list_genesets, names(list_genesets)))
  return(gsc)
}

## seurat_gsva_custom_violin(pbmc,genesets=gs,group.by='celltype_manual',projectid=projectid,plotid=plotid,species = "human",title="m6A score")
seurat_gsva_custom_violin <- function(pbmc,genesets=gmt,group.by='seurat_clusters',projectid='xxxxxx',plotid='14',species = "human",title="Ferroptosis score",plot.only=F){
		library(GSVA)
		library(msigdbr)
		#BiocManager::install('msigdbr')
		
		if(species=='human')
		{
			   DBkeyset='org.Hs.eg.db'
			   org='Homo sapiens'
			   #library(org.Hs.eg.db)
			   
		}else if(species=='mouse')
		{
				#library(org.Mm.eg.db)
				 DBkeyset='org.Mm.eg.db'
				 org='Mus musculus'
		}else if(species=='rat')
		{
				#library(org.Rn.eg.db)
				#BiocManager::install('org.Rn.eg.db')
				DBkeyset='org.Rn.eg.db'
				org='Mus musculus'
		}else{
				DBkeyset='org.Hs.eg.db'
			    #library(org.Hs.eg.db)
			    org='Homo sapiens'
		}
		if(plot.only)
		{
		    df_plot<-read.delim(paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.tsv"),header = T)
		}else{
			DefaultAssay(pbmc) <- "RNA"
		Idents(pbmc) <- group.by
		levels(pbmc)
		expr <- GetAssayData(object = pbmc, layer = "data")
		expr <- expr[rowSums(expr)>0,]  #选取非零基因
		expr <- as.matrix(expr)		
		gsvaPar <- ssgseaParam(expr, genesets,minSize=2)
		gsva.res <- gsva(gsvaPar, verbose=FALSE)
		#gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
		#gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
        metadata<-pbmc@meta.data
        df_plot<-as.data.frame(t(gsva.res))
        dim(df_plot) 
        dim(metadata) 
        df_plot[[group.by]]<-metadata[[group.by]]
        head(df_plot)
        write.table(df_plot, paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.tsv"),sep='\t',row.names = T,quote = FALSE)
	}
  colnames(df_plot)[1]<-'score'
  library(ggplot2)  
  ggplot2_boxplot_gsva(df_plot,groupby=group.by,score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'b.',projectid,".ssgsea")) 
}


## seurat_gsva_custom_violin(pbmc,genesets=gs,group.by='celltype_manual',paired_group='group',projectid=projectid,plotid=plotid,species = "human",title="m6A score")
seurat_gsva_custom_violin_paired <- function(pbmc,genesets=gmt,group.by='seurat_clusters',paired_group='group',projectid='xxxxxx',plotid='14',species = "human",title="Ferroptosis score",plot.only=F){
		#hub_gene得分，按照celltype分组，在每个组中按照group再分组，3个group的待验证
    library(GSVA)
		library(msigdbr)
		#BiocManager::install('msigdbr')
		
		if(plot.only)
		{
		    df_plot<-read.delim(paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.tsv"),header = T)
			metadata<-pbmc@meta.data
		    df_plot$celltype<-factor(df_plot$celltype,levels=levels(metadata[[group.by]]))
			#df_plot<-read.delim('fig17.GSE150825.celltype_nmf.cell.gsva.tsv',header=T)
		    
		    
		    
		}else{
			DefaultAssay(pbmc) <- "RNA"
		Idents(pbmc) <- group.by
		Idents(pbmc) <-'group'
		levels(pbmc)
		expr <- GetAssayData(object = pbmc, layer = "data")
		expr <- expr[rowSums(expr)>0,]  #选取非零基因
		expr <- as.matrix(expr)
		#expr_out<-cbind(rownames(expr),expr)
		#colnames(expr_out)[1]<-'gene'
		#write.table(expr_out, paste0("fig",plotid,".",projectid,'.',group.by,".AverageExpression.tsv"),sep='\t',row.names = F)
		gsvaPar <- ssgseaParam(expr, genesets,minSize=2)
		gsva.res <- gsva(gsvaPar, verbose=FALSE)
		#gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
        metadata<-pbmc@meta.data
        df_plot<-as.data.frame(t(gsva.res))
        dim(df_plot) 
        dim(metadata) 
        df_plot[[group.by]]<-metadata[[group.by]]
        df_plot[[paired_group]]<-metadata[[paired_group]]
        head(df_plot)
        colnames(df_plot)<-c('score','celltype','group')
        write.table(df_plot, paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.tsv"),sep='\t',row.names = T,quote = FALSE)
	}
	library(ggplot2)  
	ggplot2_boxplot_gsva(df_plot,groupby='celltype',score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'b.',projectid,".ssgsea"))  
	ggplot2_boxplot_paired(df_plot,groupby=c('group','celltype'),score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'c.',projectid,".ssgsea"))

}


seurat_CRDscore_violin_paired <- function(pbmc,genesets=geneset_genes,group.by='seurat_clusters',paired_group='group',
                                              projectid='xxxxxx',plotid='14',species = "human",title="Ferroptosis score",plot.only=F){
      # devtools::install_github("yixianfan/CRDscore")
      library(CRDscore)

      if(plot.only)
      {
        df_plot<-read.delim(paste0("fig",plotid,".",projectid,'.',group.by,".CRDscore.tsv"),header = T)
        metadata<-pbmc@meta.data
        df_plot$celltype<-factor(df_plot$celltype,levels=levels(metadata[[group.by]]))
        #df_plot<-read.delim('fig17.GSE150825.celltype_nmf.cell.gsva.tsv',header=T)
        
        
        
      }else{
        DefaultAssay(pbmc) <- "RNA"
        Idents(pbmc) <- group.by
        Idents(pbmc) <-'group'
        levels(pbmc)
        expr <- GetAssayData(object = pbmc, layer = "counts")
        expr <- expr[rowSums(expr)>0,]  #选取非零基因
        expr <- as.data.frame(expr)
        class(expr)
        
        score <- cal_CRDscore(expr = expr, n.bins = 50, circadians = genesets, study.type = "scRNAseq")
        score = as.data.frame(score)
        
        metadata<-pbmc@meta.data
        df_plot<-score
        dim(df_plot) 
        dim(metadata) 
        df_plot[[group.by]]<-metadata[[group.by]]
        df_plot[[paired_group]]<-metadata[[paired_group]]
        head(df_plot)
        colnames(df_plot)<-c('score','celltype','group')
        write.table(df_plot, paste0("fig",plotid,".",projectid,'.',group.by,".CRDscore.tsv"),sep='\t',row.names = T,quote = FALSE)
      }
      library(ggplot2)  
      ggplot2_boxplot_gsva(df_plot,groupby='celltype',score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'b.',projectid,".CRDscore"))  
      ggplot2_boxplot_paired(df_plot,groupby=c('group','celltype'),score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'c.',projectid,".CRDscore"))
      
    }


      
## list_geneSets <- list(URGs=URGs)
## seurat_AUCell_custom_violin(pbmc,genesets=list_geneSets,group.by='celltype_manual',projectid=projectid,plotid=plotid,species = "human",title="m6A score")

seurat_AUCell_custom_violin <- function(pbmc,genesets=list_geneSets,group.by='seurat_clusters',projectid='xxxxxx',plotid='14',species = "human",title="Ferroptosis score",plot.only=F){
		library(AUCell)
		if(plot.only)
		{
		    df_plot<-read.delim(paste0("fig",plotid,".",projectid,'.',group.by,".cell.AUCell.tsv"),header = T)
		}else{
			DefaultAssay(pbmc) <- "RNA"
		Idents(pbmc) <- group.by
		levels(pbmc)
		expr <- GetAssayData(object = pbmc, layer = "data")
		#expr <- expr[rowSums(expr)>0,]  #选取非零基因
		expr <- as.matrix(expr)
		#expr_out<-cbind(rownames(expr),expr)
		#colnames(expr_out)[1]<-'gene'
		#write.table(expr_out, paste0("fig",plotid,".",projectid,'.',group.by,".AverageExpression.tsv"),sep='\t',row.names = F)

		cells_AUC <- AUCell_run(expr, genesets,BPPARAM=BiocParallel::MulticoreParam(4))
		mx_auc<-getAUC(cells_AUC)
        metadata<-pbmc@meta.data
        df_plot<-as.data.frame(t(mx_auc))
        dim(df_plot)
        dim(metadata) 
        df_plot[[group.by]]<-metadata[[group.by]]
        head(df_plot)
        write.table(df_plot, paste0("fig",plotid,".",projectid,'.',group.by,".cell.AUCell.tsv"),sep='\t',row.names = T,quote = FALSE)
	}
	colnames(df_plot)[1]<-'score'
	ymin<-min(df_plot$score)
	ymax<-max(df_plot$score)
	yh<-ymax-ymin
	pw=7
	ph=4

	library(ggplot2)
	plotb<-NULL
	plotb<-ggplot(df_plot,aes(x=.data[[group.by]],y=score,fill=.data[[group.by]]))
	plotb<-plotb+geom_violin(scale = "width",size=0.1,trim = FALSE)
	#plotb<-plotb+scale_fill_manual(values=c('CD'='#ce6700', 'Non_IBD'='#1b8c19'))
	#plotb<-plotb+scale_color_manual(values=c('GO'='red', 'Normal'='blue'))
	#plotb<-plotb+stat_compare_means(aes(group = group),label = "p.signif",label.y = 0.85,hide.ns = TRUE,show.legend = FALSE)
	plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
								panel.background = element_blank(),panel.border=element_blank(),
								plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
								legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
								legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
								axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.line = element_line(linewidth=0.5, colour = "black"))
	plotb<-plotb+labs(y =title)+scale_y_continuous(expand = c(0, 0))#limits=c(0,1)
	#plotb
	ggsave(file=paste0('fig',plotid,'a.',projectid,".AUCell",".violin.tiff"),plot=plotb,width = 7, height = 4,compression='lzw')
	ggsave(file=paste0('fig',plotid,'a.',projectid,".AUCell",".violin.pdf"),plot=plotb,width = 7, height = 4)



	plotb<-NULL
	plotb<-ggplot(df_plot,aes(x=.data[[group.by]],y=score,fill=.data[[group.by]]))
	plotb<-plotb+geom_violin(scale = "width",size=0.1,trim = FALSE)
	plotb<-plotb+geom_boxplot(col='black',linewidth=0.2,linetype = 1,na.rm = T,notch = F,width = 0.3,outlier.shape = NA,position=position_dodge(0.9))
	plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
								panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
								plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
								legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
								legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
								axis.text.x=element_text(size=8,color='black',angle=45,hjust=1))
	plotb<-plotb+labs(y =title,x='')+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
	ggsave(file=paste0('fig',plotid,'c.',projectid,".AUCell",".violin.tiff"),plot=plotb,width = pw, height = ph,compression='lzw')
	ggsave(file=paste0('fig',plotid,'c.',projectid,".AUCell",".violin.pdf"),plot=plotb,width = pw, height = ph)

	library(ggpubr)
	ggviolin(df_plot,x=group.by,y='score',fill=group.by,add='boxplot',add.params = list(fill="white"),ylab=title)+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
								panel.background = element_blank(),panel.border=element_blank(),
								plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
								legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
								legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
								axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.line = element_line(linewidth=0.5, colour = "black"))
	ggsave(file=paste0('fig',plotid,'b.',projectid,".AUCell",".violin.tiff"),width = 7, height = 4,compression='lzw')
	ggsave(file=paste0('fig',plotid,'b.',projectid,".AUCell",".violin.pdf"),width = 7, height = 4)
  
  #https://zhuanlan.zhihu.com/p/482523999
  library(ggraph)
  library(viridis)
  
	df_plot_meta<-data.frame(pbmc@meta.data, pbmc@reductions$umap@cell.embeddings)
	df_plot_meta$score<-df_plot[,1]  
	plotd<-NULL
	plotd<-ggplot(df_plot_meta, aes(UMAP_1, UMAP_2, color=score)) + geom_point( size=1.5)
	plotd<-plotd+scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = title)
	plotd<-plotd+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
	plotd<-plotd+theme(plot.title = element_text(hjust = 0.5))
	ggsave(file=paste0('fig',plotid,'d.',projectid,".AUCell",".featureplot.tiff"),width = 8, height = 6,compression='lzw')
	ggsave(file=paste0('fig',plotid,'d.',projectid,".AUCell",".featureplot.pdf"),width = 8, height = 6)

}

## seurat_gsva_custom_violin(pbmc,genesets=gs,group.by='celltype_manual',paired_group='group',projectid=projectid,plotid=plotid,species = "human",title="m6A score")
seurat_AUCell_custom_violin_paired <- function(pbmc,genesets=gmt,group.by='seurat_clusters',paired_group='group',projectid='xxxxxx',plotid='14',species = "human",title="Ferroptosis score",plot.only=F){
		#hub_gene得分，按照celltype分组，在每个组中按照group再分组，3个group的待验证
    library(GSVA)
		library(msigdbr)
		#BiocManager::install('msigdbr')
		
		if(plot.only)
		{
		    df_plot<-read.delim(paste0("fig",plotid,".",projectid,'.',group.by,".cell.AUCell.tsv"),header = T)
		    
		    #df_plot<-read.delim('fig17.GSE150825.celltype_nmf.cell.gsva.tsv',header=T)
		    
		    
		    
		}else{
			DefaultAssay(pbmc) <- "RNA"
		Idents(pbmc) <- group.by
		Idents(pbmc) <-'group'
		levels(pbmc)
		expr <- GetAssayData(object = pbmc, layer = "data")
		expr <- expr[rowSums(expr)>0,]  #选取非零基因
		expr <- as.matrix(expr)
		#expr_out<-cbind(rownames(expr),expr)
		#colnames(expr_out)[1]<-'gene'
		#write.table(expr_out, paste0("fig",plotid,".",projectid,'.',group.by,".AverageExpression.tsv"),sep='\t',row.names = F)
		cells_AUC <- AUCell_run(expr, genesets,BPPARAM=BiocParallel::MulticoreParam(4))
		
		tiff(file=paste0("fig",plotid,".",projectid,'.',group.by,".cell.AUCell.hist.tiff"),width=3000,heigh=2400,units='px',res=300,compression='lzw')
		cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
		dev.off()
		save(cells_AUC, file=paste0("fig",plotid,".",projectid,'.',group.by,".cell.AUCell.rds"))


		mx_auc<-getAUC(cells_AUC)
        metadata<-pbmc@meta.data
        df_plot<-as.data.frame(t(mx_auc))
        dim(df_plot) 
        dim(metadata) 
        df_plot[[group.by]]<-metadata[[group.by]]
        df_plot[[paired_group]]<-metadata[[paired_group]]
        head(df_plot)
        colnames(df_plot)<-c('score','celltype','group')
        write.table(df_plot, paste0("fig",plotid,".",projectid,'.',group.by,".cell.AUCell.tsv"),sep='\t',row.names = T,quote = FALSE)
	}
	library(ggplot2)
	ymin<-min(df_plot$score)
  ymax<-max(df_plot$score)
  yh<-ymax-ymin
  pw=7
  ph=4 

  
	plotb<-NULL
	plotb<-ggplot(df_plot,aes(x=celltype,y=score,fill=celltype))
	plotb<-plotb+geom_violin(scale = "width",size=0.1,trim=F)
	#plotb<-plotb+scale_fill_manual(values=c('red','blue'))
	plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
								panel.background = element_blank(),panel.border=element_blank(),
								plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
								legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
								legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
								axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.line = element_line(linewidth=0.5, colour = "black"))
	plotb<-plotb+labs(x='',y =title)+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
	#plotb
	ggsave(file=paste0('fig',plotid,'a.',projectid,".AUCell",".violin.tiff"),plot=plotb,width = 7, height = 4,compression='lzw')
	ggsave(file=paste0('fig',plotid,'a.',projectid,".AUCell",".violin.pdf"),plot=plotb,width = 7, height = 4)

	library(ggpubr)
	ggviolin(df_plot,x='celltype',y='score',fill='celltype',add='boxplot',add.params = list(fill="white",width=0.1),xlab='',ylab=title)+
	#scale_fill_manual(values = c("red", "blue"))+
	theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
								panel.background = element_blank(),panel.border=element_blank(),
								plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
								legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
								legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
								axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),
								axis.line = element_line(linewidth=0.5, colour = "black"))+
	scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
	ggsave(file=paste0('fig',plotid,'b.',projectid,".AUCell",".violin.tiff"),width = 7, height = 4,compression='lzw')
	ggsave(file=paste0('fig',plotid,'b.',projectid,".AUCell",".violin.pdf"),width = 7, height = 4)

	plotb<-NULL
	plotb<-ggplot(df_plot,aes(x=celltype,y=score,fill=celltype))
	plotb<-plotb+geom_violin(scale = "width",size=0.1,trim = FALSE)
	plotb<-plotb+geom_boxplot(col='black',linewidth=0.2,linetype = 1,na.rm = T,notch = F,width = 0.2,outlier.shape = NA,position=position_dodge(0.9))
	plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
								panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
								plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
								legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
								legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
								axis.text.x=element_text(size=8,color='black',angle=45,hjust=1))
	plotb<-plotb+labs(y =title,x='')+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
	ggsave(file=paste0('fig',plotid,'c.',projectid,".AUCell",".violin.tiff"),plot=plotb,width = pw, height = ph,compression='lzw')
	ggsave(file=paste0('fig',plotid,'c.',projectid,".AUCell",".violin.pdf"),plot=plotb,width = pw, height = ph)
  
	#https://zhuanlan.zhihu.com/p/482523999
	library(ggraph)
	library(viridis)
	df_plot_meta<-data.frame(pbmc@meta.data, pbmc@reductions$umap@cell.embeddings)
	df_plot_meta$score<-df_plot[,1]  
	plotd<-NULL
	plotd<-ggplot(df_plot_meta, aes(UMAP_1, UMAP_2, color=score)) + geom_point( size=1.5)
	plotd<-plotd+scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = title)
	plotd<-plotd+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
	plotd<-plotd+theme(plot.title = element_text(hjust = 0.5))
	ggsave(file=paste0('fig',plotid,'d.',projectid,".AUCell",".featureplot.tiff"),width = 8, height = 6,compression='lzw')
	ggsave(file=paste0('fig',plotid,'d.',projectid,".AUCell",".featureplot.pdf"),width = 8, height = 6)
  
  ggplot2_boxplot_paired(df_plot,groupby=c('group','celltype'),score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'e.',projectid,".AUCell"))

}


## ggheatmap_yingbio(matrix_cor,ph=8,output=paste0('fig',plotid,'b.',comparison_meta,'.heatmap'))
ggheatmap_yingbio<-function(mydata,orderby=NA,geneid='miRNA_ID',row_name=T,columne_name=T,ph=NA,pw=6,output='heatmap',save.data=F){
  ngene<-dim(mydata)[1]
  nsample<-dim(mydata)[2]
  if(!is.na(ph))
  {
    phu<-ph
  }else{
    phu<-max(6,ngene/10)
  }
  if(!is.na(pw))
  {
    pwu<-pw
  }else{
    pwu<-max(6,nsample/10)
  }
  print(paste0('plot height & width: ',phu,' ',pwu))
  library(ggplot2)
  library(pheatmap)
  library(gtable)
  plot.pheatmap<-NULL
  plot.pheatmap<-pheatmap::pheatmap(mydata,border = F,
                                    color = colorRampPalette(c(rep("dodgerblue",2), "white", rep("brown2",2)))(100),
                                    cluster_cols =F,
                                    show_colnames = columne_name,
                                    show_rownames = row_name,
                                    scale='none',
                                    fontsize = 10,
                                    fontsize_row=5,
                                    fontsize_col=10,angle_col = "45")
  ####  如果右侧的legend空间不够，增加legend空间。
  #plot.pheatmap$gtable$widths[1] <- plot.pheatmap$gtable$widths[1] + unit(200, "bigpts")
  #plot.pheatmap$gtable$widths[6] <- plot.pheatmap$gtable$widths[6] + unit(5, "bigpts")
  dev.off()
  pdf(file=paste0(output,'.pdf'), height=phu, width=pwu)
  print(plot.pheatmap)
  dev.off()
  tiff(file=paste0(output,'.tiff'), height=phu*300, width=pwu*300,units="px",res=300,compression='lzw')
  print(plot.pheatmap)
  dev.off()
}

ggheatmap_yingbio_1col<-function(mydata,orderby=NA,geneid='miRNA_ID',pheight=NA,pw=6,output='heatmap',save.data=F){
  ngene<-dim(mydata)[1]
  if(!is.na(pheight))
  {
    ph<-pheight
  }else{
    ph<-max(6,ngene/10)
  }
  print(paste0('pic height',ph))
  library(ggplot2)
  library(pheatmap)
  library(gtable)
  plot.pheatmap<-NULL
  plot.pheatmap<-pheatmap::pheatmap(mydata,border = T,
                                    color = colorRampPalette(c(rep("dodgerblue",2), "white", rep("brown2",2)))(100),
                                    cluster_rows =F,
                                    cluster_cols =F,
                                    show_colnames = F,
                                    show_rownames = T,
                                    fontsize = 10,
                                    fontsize_row=5,
                                    fontsize_col=10,angle_col = "45")
  ####  如果右侧的legend空间不够，增加legend空间。
  #plot.pheatmap$gtable$widths[1] <- plot.pheatmap$gtable$widths[1] + unit(200, "bigpts")
  #plot.pheatmap$gtable$widths[6] <- plot.pheatmap$gtable$widths[6] + unit(5, "bigpts")
  dev.off()
  pdf(file=paste0(output,'.pdf'), height=ph, width=pw)
  print(plot.pheatmap)
  dev.off()
  tiff(file=paste0(output,'.tiff'), height=ph*300, width=pw*300,units="px",res=300,compression='lzw')
  print(plot.pheatmap)
  dev.off()
}

## ggheatmap_yingbio_2groups(df_tpm_log,phenotype,orderby='group',geneid='tsRNA',output=paste0('fig',plotid,'.heatmap.',comparison),save.data=T)
ggheatmap_yingbio_2groups<-function(mydata,phenotype,orderby=NA,show_colname = F,show_rowname = F,geneid='miRNA_ID',pheight=NA,pw=6,output='heatmap',save.data=F){
  phenotype<-phenotype[colnames(mydata),,drop=F]
  if(is.na(orderby))
  {
    phenotype<-phenotype
  }else{
    phenotype<-phenotype[order(phenotype[[orderby]]),,drop=F]
  }
  mydata<-mydata[,rownames(phenotype)]
  if(save.data)
  {
    mydata_out<-cbind(rownames(mydata),mydata)
    colnames(mydata_out)[1]<-geneid
    write.table(mydata_out,file = paste0(output,'.tsv'),row.names = FALSE,sep='\t',quote = FALSE)
  }
  ngene<-dim(mydata)[1]
  if(!is.na(pheight))
  {
    ph<-pheight
  }else{
    ph<-max(6,ngene/10)
  }
  print(paste0('pic height',ph))
  library(ggplot2)
  library(pheatmap)
  library(gtable)
  ngene<-nrow(mydata)
  if(ngene<100)
  {
	  show_rowname=T
	}
  plot.pheatmap<-NULL
  plot.pheatmap<-pheatmap::pheatmap(mydata,border = F,
                                    annotation=phenotype, 
                                    #color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
                                    color = colorRampPalette(c(rep("dodgerblue",2), "white", rep("brown2",2)))(100),
                                    #color = viridis(8, option = "G")
                                    cluster_cols =F,
                                    show_colnames = show_colname,
                                    show_rownames = show_rowname,
                                    scale="row",
                                    fontsize = 10,
                                    fontsize_row=5,
                                    fontsize_col=10,angle_col = "45")
  ####  如果右侧的legend空间不够，增加legend空间。
  #plot.pheatmap$gtable$widths[1] <- plot.pheatmap$gtable$widths[1] + unit(200, "bigpts")
  plot.pheatmap$gtable$widths[6] <- plot.pheatmap$gtable$widths[6] + unit(5, "bigpts")
  dev.off()
  pdf(file=paste0(output,'.pdf'), height=ph, width=pw)
  print(plot.pheatmap)
  dev.off()
  tiff(file=paste0(output,'.tiff'), height=ph*300, width=pw*300,units="px",res=300,compression='lzw')
  print(plot.pheatmap)
  dev.off()
}



####  df_data_heatmap:数据矩阵，行为基因，列为样本。
####  df_left_annotation:行block的注释矩阵，类似phenotype的单列数据框，列属性最好为factor。
heatmap_with_row_block<-function(df_data_heatmap,df_left_annotation,x=70,lab_width=40,ph=6,pw=10,output='output')
{
        library(ComplexHeatmap)
        group<-colnames(df_left_annotation)[1]
        if(class(df_left_annotation[[group]])=='factor')
        {
          groups<-levels(df_left_annotation[[group]])
        }else{
          groups<-unique(df_left_annotation[[group]])
        }
        ngroup<-length(groups)
        df_left_annotation<-df_left_annotation[order(df_left_annotation[[group]]),,drop=F]
        labels= str_wrap(groups, width=lab_width)
        pathways<-rownames(df_left_annotation)
        df_data_heatmap<-df_data_heatmap[pathways,]
        
        df_data_heatmap<-as.matrix(df_data_heatmap)
        
        left_annotation = rowAnnotation(foo = anno_empty(
          border =F,
          width = max_text_width(unlist(labels)) + unit(4, "mm")
        ),
        show_legend = c("foo" = FALSE)
        )
        tiff(file=paste0(output,'.tiff'), height=ph*300, width=pw*300,res=300,compression='lzw')
        print(Heatmap(
          df_data_heatmap,
          cluster_columns = F,
          cluster_rows = F,
          column_names_rot = 45,
          row_split = df_left_annotation[[group]],
          row_title = NULL,
          row_names_gp = gpar(fontsize = 7),
          heatmap_legend_param = list(title = ''),
          left_annotation = left_annotation
        ))
        for(i in 1:ngroup) {
          #i=1
          decorate_annotation(
            "foo",slice = i, 
            {
              grid.rect(x = unit(x, "mm"), width = unit(4, "mm"), 
                        gp = gpar(fill = rainbow(ngroup)[i], col = NA),just = "left")
              grid.text(labels[i], x = unit(x-2, "mm"), rot=10,
                        gp = gpar(col = 'black',fontsize=10), just = c(1,0))
            })
        }
        dev.off()
        pdf(file=paste0(output,'.pdf'), height=ph, width=pw)
        print(Heatmap(
          df_data_heatmap,
          cluster_columns = F,
          cluster_rows = F,
          column_names_rot = 45,
          row_split = df_left_annotation[[group]],
          row_title = NULL,
          row_names_gp = gpar(fontsize = 7),
          heatmap_legend_param = list(title = ''),
          left_annotation = left_annotation
        ))
        for(i in 1:ngroup) {
          #i=1
          decorate_annotation(
            "foo",slice = i, 
            {
              grid.rect(x = unit(x, "mm"), width = unit(4, "mm"), 
                        gp = gpar(fill = rainbow(ngroup)[i], col = NA),just = "left")
              grid.text(labels[i], x = unit(x-2, "mm"), rot=10,
                        gp = gpar(col = 'black',fontsize=10), just = c(1,0))
            })
        }
        dev.off()
}


### seurat_gsva_2group_bycell(pbmc,group.by='m7G_group',projectid=projectid,plotid=plotid,plot.only=F)
seurat_gsva_2group_bycell <- function(pbmc,group.by='seurat_clusters',projectid='xxxxxx',plotid='14',species = "human",plot.only=F){
		library(GSVA)
		library(msigdbr)
		#BiocManager::install('msigdbr')
		
		if(species=='human')
		{
			   DBkeyset='org.Hs.eg.db'
			   org='Homo sapiens'
			   #library(org.Hs.eg.db)
			   
		}else if(species=='mouse')
		{
				#library(org.Mm.eg.db)
				 DBkeyset='org.Mm.eg.db'
				 org='Mus musculus'
		}else if(species=='rat')
		{
				#library(org.Rn.eg.db)
				#BiocManager::install('org.Rn.eg.db')
				DBkeyset='org.Rn.eg.db'
				org='Mus musculus'
		}else{
				DBkeyset='org.Hs.eg.db'
			    #library(org.Hs.eg.db)
			    org='Homo sapiens'
		}
		if(plot.only)
		{
		    df_plot_hallmarker<-read.delim(paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.hallmarker.tsv"),header = T,row.names=1)
		    df_plot_kegg<-read.delim(paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.kegg.tsv"),header = T,row.names=1)
		}else{
		Idents(pbmc) <- group.by
		levels(pbmc)
		expr <- GetAssayData(object = pbmc, layer = "data")
		expr <- expr[rowSums(expr)>0,]  #选取非零基因
		expr <- as.matrix(expr)
		if(1)  ###  计算hallmarker的gsva
		{
			genesets_hallmarker <- msigdbr(species = org, category = "H") 
			genesets_hallmarker <- subset(genesets_hallmarker, select = c("gs_name","gene_symbol")) %>% as.data.frame()
			genesets_hallmarker <- split(genesets_hallmarker$gene_symbol, genesets_hallmarker$gs_name)
			gsva.res_hallmarker <- gsva(expr, genesets_hallmarker, method="ssgsea",parallel.sz=10,min.sz=1) 
			df_plot_hallmarker<-as.data.frame(t(gsva.res_hallmarker))
			write.table(df_plot_hallmarker, paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.hallmarker.tsv"),sep='\t',row.names = T, quote=F)
		}
		if(1)  ###  计算kegg的gsva
		{
			genesets_kegg <- msigdbr(species = org, category = "C2") 
			genesets_kegg <- subset(genesets_kegg, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
			genesets_kegg <- split(genesets_kegg$gene_symbol, genesets_kegg$gs_name)
			gsva.res_kegg <- gsva(expr, genesets_kegg, method="ssgsea",parallel.sz=10,min.sz=1) 
			df_plot_kegg<-as.data.frame(t(gsva.res_kegg))
			write.table(df_plot_kegg, paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.kegg.tsv"),sep='\t',row.names = T, quote=F)
		}
	}
	df_plot_hallmarker_heatmap<-t(df_plot_hallmarker)
	phenotype<-pbmc@meta.data[,group.by,drop=F]
	ggheatmap_yingbio_2groups_for_gsva(df_plot_hallmarker_heatmap,phenotype,geneid='gene',output=paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.hallmarker"))
	df_plot_kegg_heatmap<-t(df_plot_kegg)
	ggheatmap_yingbio_2groups_for_gsva(df_plot_kegg_heatmap,phenotype,geneid='gene',output=paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.kegg"))
}



label_deg_gene_gf<- function(filename='fig05c.deg.CD_vs_control.sig.tsv',gene_column='gene',species = "human"){
	load('/media/nbfs/panzhong/genome_sequence/scrna/cell_factor.YB20230220.Rdata')
	df_deg<-read.table(filename,header=T,sep='\t')
	df_deg$modification<-'NO'
	df_deg[df_deg[[gene_column]] %in% m6a_genes,'modification']<-'YES'

	df_deg$cell_factor<-'NO'
	df_deg[df_deg[[gene_column]] %in% cell_factor_gene,'cell_factor']<-'YES'

	df_deg$growth_factor<-'NO'
	df_deg[df_deg[[gene_column]] %in% growth_factor_gene,'growth_factor']<-'YES'

	df_deg$checkpoint<-'NO'
	df_deg[df_deg[[gene_column]] %in% immune_checkpoint_gene,'checkpoint']<-'YES'

	df_deg$tf<-'NO'
	df_deg[df_deg[[gene_column]] %in% transcription_factor_gene,'tf']<-'YES'
	write.table(df_deg,file=paste0(filename,".labled.tsv"),sep="\t",row.names=F,quote=F)
}

label_deg_gene<- function(filename='fig05c.deg.CD_vs_control.sig.tsv',genes=c('TP53'),gene_column='gene',label='m7G'){
	df_deg<-read.table(filename,header=T,sep='\t')
	df_deg[,label]<-'NO'
	df_deg[df_deg[[gene_column]] %in% genes,label]<-'YES'
	write.table(df_deg,file=paste0(filename,".labled.tsv"),sep="\t",row.names=F,quote=F)
}

#                                                                      paste0('fig',plotid,'.volcano.',comparison)
#ggplot2_volcano(df_data[,c('Pvalue','log2FC')],comparison=comparison,output=paste0('fig0.volcano.',comparison),fc_threshold=2^logFCfilter,ylab='P value')
ggplot2_volcano<-function(mydata,comparison='comparison',fc_threshold=2,pvalue_threshold=0.05,ylab='Pvalue',output='volcano',ymax=300){
  #mydata<-DEG_full[,c('padj','log2FoldChange')]
  library(ggplot2)
  colnames(mydata)<-c('PValue','logFC')
  mydata$PValue[is.na(mydata$PValue)]<-1
  mydata$regulation<-c('normal')
  mydata$PValue<- -log10(mydata$PValue)		
  ymax_used <- max(mydata$PValue)
  if(ymax_used > ymax)
  {
    ymax_used <- ymax
  }
  mydata$PValue[mydata$PValue > ymax_used]<- ymax_used
  
  mydata[mydata$logFC <= -log2(fc_threshold) & mydata$PValue >= -log10(pvalue_threshold),3]<-'down'
  
  mydata[mydata$logFC >= log2(fc_threshold) & mydata$PValue >= -log10(pvalue_threshold),3]<-'up'
  
  cols <- c('up' = 'red', 'normal' = 'gray', 'down' = 'blue')
  
  xmax<-round(max(abs(mydata$logFC)),0)+1
  
  myplot<-NULL
  myplot<-ggplot(mydata,aes(logFC, PValue,col =regulation))+ geom_point(size=1)+
    labs(title=comparison,x=expression(paste('log'[2],'(fold change)')), y=substitute(paste('-log'["10"],'(',ylab,')'),list(ylab=ylab)))+
    #labs(title=substitute(paste("Histogram of random data with",mu,"=",m,",",sigma^2,"=",s2,",","draws =", numdraws,",",bar(x),"=",xbar,",",s^2,"=",sde),list(m=x_mean,xbar=mean(x),s2=x_sd^2,sde=var(x),numdraws=N))
    scale_color_manual(values =cols,limits = c('up', 'down'))+
    scale_x_continuous(limits=c(-xmax,xmax))+
    #expand = c(0, 0)
    scale_y_continuous(limits=c(0,ymax_used+0.5),expand = c(0, 0))+
    geom_hline(yintercept=-log10(pvalue_threshold),linetype=4,color='black',size=1)+
    geom_vline(xintercept=log(fc_threshold,2),linetype=4,color='black',size=1)+
    geom_vline(xintercept=-log(fc_threshold,2),linetype=4,color='black',size=1)+
    theme(text=element_text(size=20),axis.title.x =element_text(size=30,color='black'),axis.title.y =element_text(size=30,color='black'))+
    theme(panel.border=element_rect(linetype='solid',fill=NA,colour = 'black',size=0.7))+
    #theme(panel.grid.major =element_line(colour = 'black', size = 0.25), panel.grid.minor = element_blank())+
    theme(panel.background = element_blank(),plot.title = element_text(hjust = 0.5,size=30),legend.position='none')
  #element_blank()
  #geom_hline(yintercept=1.3,linetype=2)
  #myplot
  ggsave(paste0(output,".tiff"),plot=myplot, width = 8, height = 8,compression='lzw')
  ggsave(paste0(output,".pdf"),plot=myplot, width = 8, height = 8)
}

ggplot2_volcano20240124<-function(mydata,comparison='comparison',fc_threshold=2,pvalue_threshold=0.05,ylab='Pvalue',output='volcano',ymax=300){
  #mydata<-df_data[,c('Pvalue','log2FC')]
  library(ggplot2)
  colnames(mydata)<-c('PValue','logFC')
  mydata$PValue[is.na(mydata$PValue)]<-1
  mydata$regulation<-c('No significant')
  mydata$PValue<- -log10(mydata$PValue)		
  ymax_used <- max(mydata$PValue)
  mydata$PValue[mydata$PValue > ymax_used]<- ymax_used
  
  mydata[mydata$logFC <= -log2(fc_threshold) & mydata$PValue >= -log10(pvalue_threshold),3]<-'Down'
  
  mydata[mydata$logFC >= log2(fc_threshold) & mydata$PValue >= -log10(pvalue_threshold),3]<-'Up'
  
  cols <- c('Up' = '#ed422e', 'No significant' = '#656565', 'Down' = '#52b6d8')
  cols_fill <- c('Up' = '#ff9b8f', 'No significant' = '#bababa', 'Down' = '#9bdeef')
  
  xmax<-round(max(abs(mydata$logFC)),0)+1
  
  table(mydata$regulation)
  
  comparison_title<-gsub('_',' ',comparison)
  
  myplot<-NULL
  myplot<-ggplot(mydata,aes(logFC, PValue))+ geom_point(aes(fill=regulation,col=regulation), size=2,shape=21)+
    labs(title=comparison_title,x=expression(paste('log'[2],'(Fold change)')), y=substitute(paste('-log'["10"],'(',ylab,')'),list(ylab=ylab)))+
    #labs(title=substitute(paste("Histogram of random data with",mu,"=",m,",",sigma^2,"=",s2,",","draws =", numdraws,",",bar(x),"=",xbar,",",s^2,"=",sde),list(m=x_mean,xbar=mean(x),s2=x_sd^2,sde=var(x),numdraws=N))
    scale_color_manual(values =cols)+
    scale_fill_manual(values =cols_fill)+
    scale_x_continuous(limits=c(-xmax,xmax))+
    #expand = c(0, 0)
    scale_y_continuous(limits=c(0,ymax_used+0.5),expand = c(0, 0))+
    geom_hline(yintercept=-log10(pvalue_threshold),linetype=2,color='#656565',size=1)+
    geom_vline(xintercept=log(fc_threshold,2),linetype=2,color='#656565',size=1)+
    geom_vline(xintercept=-log(fc_threshold,2),linetype=2,color='#656565',size=1)+
    theme(text=element_text(size=20),axis.title.x =element_text(size=20,color='black'),axis.title.y =element_text(size=20,color='black'))+
    theme(panel.border=element_blank(),axis.line = element_line(colour = 'black'))+
    #theme(panel.grid.major =element_line(colour = 'black', size = 0.25), panel.grid.minor = element_blank())+
    theme(panel.background = element_blank(),plot.title = element_text(hjust = 0.5,size=20),
          legend.position='right',legend.key = element_rect(colour = NA, fill = NA),legend.title=element_blank())
  #element_blank()
  #geom_hline(yintercept=1.3,linetype=2)
  myplot
  ggsave(paste0(output,".tiff"),plot=myplot, width = 8, height = 6,compression='lzw')
  ggsave(paste0(output,".pdf"),plot=myplot, width = 8, height = 6)
}


# ggplot2_pca(df_count_disease,phenotype,output=paste0('fig',plotid,'.',projectid,".PCA.CC_group"))
ggplot2_pca<-function(mydata,phenotype,groupby='group',output='fig0.project.PCA'){
	  phenotype<-phenotype[order(phenotype[,1]),,drop=F]
	  mydata<-mydata[,rownames(phenotype)]
	  mydata<-mydata[rowSums(mydata>0)>0,]
	  pca.out<-prcomp(t(mydata),scale=T,rank=4,retx=T)
	  summary(pca.out)
	  bioCol=c("blue","red","green","yellow")
	  ngroup<-length(unique(phenotype[,1]))
	  mycols<-bioCol[1:ngroup]
	  
	  pcaPredict=predict(pca.out)
	  PCA = data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2],group=phenotype[,1],sampleid=rownames(pcaPredict))
	  p=ggplot(data = PCA,aes(PC1, PC2)) + geom_point(aes(shape=group,color = group)) +scale_color_manual(groupby,values=mycols)+
	  #scale_colour_manual(name=var,values =col)+
	  theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	  library(ggrepel)
	  p <- p +  geom_text_repel(aes(label = sampleid), size = 3, show.legend = FALSE, box.padding = unit(0.5, 'lines'))
	  ggsave(paste0(output,".tiff"),plot=p,width = 6, height = 4,compression='lzw')
	  ggsave(paste0(output,".pdf"),plot=p,width = 6, height = 4)
}

ggplot2_tsne<-function(mydata,phenotype,groupby='group',output='fig0.project.TSNE'){
        library(ggpubr)
        library(ggthemes)
        library(Rtsne)
        
        phenotype<-phenotype[order(phenotype[,1]),,drop=F]
        mydata<-mydata[,rownames(phenotype)]
        mydata<-mydata[rowSums(mydata>0)>0,]
        
        
        Atsne <- Rtsne(t(mydata), perplexity = 3)
        
        head(Atsne$Y)
        dim(Atsne$Y)
        df_tsne<-as.data.frame(Atsne$Y)
        head(df_tsne)
        colnames(df_tsne)<-c('TSNE1','TSNE2')
        head(df_tsne)
        head(phenotype)
        df_tsne$sample<-colnames(mydata)
        df_tsne$group<-phenotype[,1]
        
        

        bioCol=c("blue","red","green","yellow")
        ngroup<-length(unique(phenotype[,1]))
        mycols<-bioCol[1:ngroup]

        plota=ggplot(data = df_tsne,aes(TSNE1, TSNE2)) + geom_point(aes(color = group)) +scale_color_manual(groupby,values=mycols)+
          #scale_colour_manual(name=var,values =col)+
          theme_bw()+
          theme(plot.margin=unit(rep(1.5,4),'lines'))+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        library(ggrepel)
        plota <- plota +  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE, box.padding = unit(0.5, 'lines'))
        #print(plota)
        ggsave(paste0(output,".tiff"),plot=plota,width = 6, height = 4,compression='lzw')
        ggsave(paste0(output,".pdf"),plot=plota,width = 6, height = 4)
      }
      
ggplot2_umap<-function(mydata,phenotype,groupby='group',output='fig0.project.UMAP'){
        library(ggpubr)
        library(ggthemes)
        library(umap)
        
        phenotype<-phenotype[order(phenotype[,1]),,drop=F]
        mydata<-mydata[,rownames(phenotype)]
        mydata<-mydata[rowSums(mydata>0)>0,]
        
        umap_results <- umap::umap(t(mydata))
        head(umap_results$layout)

        df_umap<-as.data.frame(umap_results$layout)
        head(df_umap)
        colnames(df_umap)<-c('UMAP1','UMAP2')
        head(df_umap)
        head(phenotype)
        df_umap$sample<-colnames(mydata)
        df_umap$group<-phenotype[,1]
        
        bioCol=c("blue","red","green","yellow")
        ngroup<-length(unique(phenotype[,1]))
        mycols<-bioCol[1:ngroup]
        
        plota=ggplot(data = df_umap,aes(UMAP1, UMAP2)) + geom_point(aes(color = group)) +scale_color_manual(groupby,values=mycols)+
          #scale_colour_manual(name=var,values =col)+
          theme_bw()+
          theme(plot.margin=unit(rep(1.5,4),'lines'))+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        library(ggrepel)
        plota <- plota +  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE, box.padding = unit(0.5, 'lines'))
        #print(plota)
        ggsave(paste0(output,".tiff"),plot=plota,width = 6, height = 4,compression='lzw')
        ggsave(paste0(output,".pdf"),plot=plota,width = 6, height = 4)
      }

###  ggplot2_pca_3d(df_logtpm_used,phenotype,output=paste0('fig',3,'.',projectid,".PCA-3dd"))
ggplot2_pca_3d<-function(mydata,phenotype,output='fig0.project.PCA'){
  	  library(scatterplot3d)
  	  library(dplyr)
  	  phenotype<-phenotype[order(phenotype[,1]),,drop=F]
  	  mydata<-mydata[,rownames(phenotype)]
  	  mydata<-mydata[rowSums(mydata>0)>0,]
  	  pca.out<-prcomp(t(mydata),scale=T,rank=4,retx=T)
  	  summary(pca.out)
  	  bioCol=c("blue","red","green","yellow")
  	  group=levels(phenotype$group)
  	  ngroup<-length(group)
  	  mycolors<-rev(get_colors(ngroup))
  	  colors= mycolors[match(phenotype$group,group)]
  	  pcaPredict=predict(pca.out)
  	  
  	  pdf(file=paste0(output,".pdf"), height=5, width=6)
  	  par(oma=c(0.5,0.5,0.5,0.5))
  	  scatterplot3d(pcaPredict[,1:3], pch = 16, color=colors,lty.hide=2)
  	  legend("bottom",group,fill=mycolors)
  	  dev.off()
  	  
  	  tiff(file=paste0(output,".tiff"), height=5*300, width=6*300,res=300,compression='lzw')
  	  par(oma=c(0.5,0.5,0.5,0.5))
  	  scatterplot3d(pcaPredict[,1:3], pch = 16, color=colors,lty.hide=2)
  	  legend("topleft",group,fill=mycolors,box.col=NA)
  	  dev.off()
}
  	

ggplot2_pca_2factor<-function(mydata,phenotype,output='fig0.project.PCA'){
  
  names<-colnames(phenotype)
  phenotype<-phenotype[order(phenotype[,1]),,drop=F]
  mydata<-mydata[,rownames(phenotype)]
  mydata<-mydata[rowSums(mydata>0)>0,]
  pca.out<-prcomp(t(mydata),scale=T,rank=4,retx=T)
  summary(pca.out)
  pcaPredict=predict(pca.out)
  PCA = data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2],group1=phenotype[,1],group2=phenotype[,2],sampleid=rownames(pcaPredict))
  colnames(PCA)[3:4]<-names
  p=ggplot(data = PCA,aes(PC1, PC2)) + geom_point(aes(shape=.data[[names[1]]],color = .data[[names[2]]])) +
    #scale_colour_manual(name=var,values =col)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  library(ggrepel)
  p <- p +  geom_text_repel(aes(label = sampleid), size = 3, show.legend = FALSE, box.padding = unit(0.5, 'lines'))
  ggsave(paste0(output,".tiff"),plot=p,width = 8, height = 6,compression='lzw')
  ggsave(paste0(output,".pdf"),plot=p,width = 8, height = 6)
}


###  ggplot2_sankey_2column(df_group_plot,output='fig67e.sankey')
ggplot2_sankey_2column<-function(df_plot,output='fig67e.sankey')
{
      library(ggalluvial)
        df_plot$Freq=1
        df_lodes <- to_lodes_form(df_plot,key ="x", value = "stratum", id = "alluvium",axes =1:2)
        head(df_lodes)
        plota<-NULL
        plota<-ggplot(df_lodes,aes(x = x, stratum =stratum, alluvium = alluvium,
                              fill = stratum, label = stratum)) +
        scale_x_discrete(expand = c(0, 0)) +
        geom_flow(width = 0.2, knot.pos = 0.1) +
        geom_stratum(alpha = .9,color="grey20",width = 1/7) +
        geom_text(stat = "stratum", size =8,color="black") +
        #scale_fill_manual(values = mycol3) +
        xlab("") + ylab("") +
        theme_bw() +
        theme(panel.grid =element_blank()) +
        theme(panel.border = element_blank()) +
        theme(axis.line = element_blank(),axis.ticks =element_blank(),
              axis.text.y =element_blank(),axis.text.x =element_text(color='black',size=15))+
        guides(fill = FALSE)
        ggsave(paste0(output,".tiff"),plot=plota, width = 8, height = 6,compression='lzw')
        ggsave(paste0(output,".pdf"),plot=plota, width = 8, height = 6)
}

ggplot2_scatterplot_deg<-function(mydata,comparison=comparison,fc_threshold=2,pvalue_threshold=0.05,ylab='Pvalue',plotid='08'){
  library(ggplot2)
  #mydata<-DEG_full[,c('Normal','Tumor','padj','log2FoldChange')]
  #group_names<-c('Normal','Tumor','padj','log2FoldChange')
  group_names<-colnames(mydata)
  colnames(mydata)<-c('ctrl','test','PValue','logFC')
  mydata$regulation<-c('normal')
  #mydata$PValue[mydata$PValue > ymax_used]<- ymax_used
  mydata[mydata$logFC <= -log2(fc_threshold) & mydata$PValue < pvalue_threshold,'regulation']<-'down'
  mydata[mydata$logFC >= log2(fc_threshold) & mydata$PValue < pvalue_threshold,'regulation']<-'up'
  cols <- c('up' = 'red', 'normal' = 'gray', 'down' = 'blue')
  unique(mydata$regulation)
  head(mydata)
  mydata$ctrl<-log2(mydata$ctrl+1)
  mydata$test<-log2(mydata$test+1)
  xmax<-round(max(max(mydata$ctrl),max(mydata$test)),0)+1
  
  myplot<-NULL
  myplot<-ggplot(mydata,aes(ctrl, test,col =regulation))+ geom_point(size=1)+
    labs(x=group_names[1], y=group_names[2])+
    geom_abline(intercept = 0, slope = 1,linetype=5)+
    #labs(title=substitute(paste("Histogram of random data with",mu,"=",m,",",sigma^2,"=",s2,",","draws =", numdraws,",",bar(x),"=",xbar,",",s^2,"=",sde),list(m=x_mean,xbar=mean(x),s2=x_sd^2,sde=var(x),numdraws=N))
    scale_color_manual(values =cols,limits = c('up', 'down'))+
    scale_x_continuous(expand = c(0, 0),limits=c(0,xmax))+
    #expand = c(0, 0)
    scale_y_continuous(expand = c(0, 0),limits=c(0,xmax))+
    #geom_hline(yintercept=-log10(pvalue_threshold),linetype=4,color='black',size=1)+
    #geom_vline(xintercept=log(fc_threshold,2),linetype=4,color='black',size=1)+
    #geom_vline(xintercept=-log(fc_threshold,2),linetype=4,color='black',size=1)+
    theme(text=element_text(size=20),axis.title.x =element_text(size=30,color='black'),axis.title.y =element_text(size=30,color='black'))+
    theme(panel.border=element_rect(linetype='solid',fill=NA,colour = 'black',size=0.7))+
    #theme(panel.grid.major =element_line(colour = 'black', size = 0.25), panel.grid.minor = element_blank())+
    theme(panel.background = element_blank(),plot.title = element_text(hjust = 0.5,size=30),legend.position=c(0.2,0.9),legend.title=element_blank())
  #element_blank()
  #geom_hline(yintercept=1.3,linetype=2)
  #myplot
  ggsave(paste0("fig",plotid,".",comparison,".scatter.tiff"),plot=myplot, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,".",comparison,".scatter.pdf"),plot=myplot, width = 8, height = 8)
}

ggscatter_2columns<-function(data_plot,output='fig.Scatterplot',method="pearson"){
  
  library(ggpubr)
  library(psych)
  genea<-colnames(data_plot)[1]
  geneb<-colnames(data_plot)[2]
  colnames(data_plot)<-c('a','b')
  head(data_plot)
  plota<-NULL
  plota <- ggscatter(data_plot, x = 'a', y = 'b',
                     add = "reg.line",  # Add regressin line
                     add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                     conf.int = TRUE # Add confidence interval
  )
  plota<-plota+theme(plot.margin=unit(c(1,1,3,1.5), 'lines'),
                     axis.title.x=element_text(size=30,color='black',vjust=0.5,hjust=0.5),
                     axis.title.y=element_text(size=30,color='black'),
                     axis.text.x =element_text(size=10,color='black',hjust = 1), 
                     axis.text.y=element_text(size=10,color='black'),legend.position='none',
                     legend.title=element_blank(),legend.text = element_text(size = 20),
                     legend.key.size = unit(1.2, 'lines'))+labs(x=genea,y=geneb)
  mi<-min(data_plot[,1])
  mx<-max(data_plot[,2])
  matrix_cor_p <- cor.test(data_plot[,1],data_plot[,2],method=method)
  r<-matrix_cor_p$estimate
  p<-matrix_cor_p$p.value
  label<-paste0('R=',round(r,3),', p=',format(p, scientific = TRUE,digits = 3))
  print(label)
  # Add correlation coefficient
  #plota<-plota + stat_cor(method = method,label.x = mi, label.y = mx,color = "blue",size=6)
  plota<-plota + annotate(geom = "text", x = mi, y = mx, label = label,color = "blue",size=6,hjust='left')
  plota<-plota+labs(x=genea,y=geneb) 
  ggsave(paste0(output,".tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = 8, height = 8)
}

ggplot2_scatter_with_density<-function(data_plot,density=T,output='ggplot2_scatter',method="pearson")
{
	library(ggplot2)
    library(ggExtra)
    library(ggpubr)
  corT=cor.test(data_plot[,1],data_plot[,2],method=method)
  cor=corT$estimate
  pValue=corT$p.value
  names<-colnames(data_plot)
  colnames(data_plot)<-c('x','y')
  p1=ggplot(data_plot, aes(x, y)) + 
    xlab(names[1])+ylab(names[2])+  #ggtitle(paste0("\nCancer: ",i))+
    theme(title=element_text(size=10))+
    geom_point()+ geom_smooth(method="lm",formula=y ~ x) + theme_bw()+
    stat_cor(method = method, aes(x =x, y =y))
  p2=ggMarginal(p1, type="density", xparams=list(fill="orange"), yparams=list(fill="blue"))
  if(1)
  {
	  ggsave(paste0(output,".tiff"), plot = p2, width = 6, height = 6,compression='lzw')
	  ggsave(paste0(output,".pdf"), plot = p2, width = 6, height = 6)
	}else
  {
	  ggsave(paste0(output,".tiff"), plot = p1, width = 6, height = 6,compression='lzw')
	  ggsave(paste0(output,".pdf"), plot = p1, width = 6, height = 6)
	}
  return(corT)  
}


### mydata表格中，应该至少包含三列：OS、OS_time、group。
###                            group可以是两组或者更多组。
###                            OS_time最好是以month作为单位。
# ggplot2_KM_plot(lassoSigExp,plotid='36',projectid=projectid,sur='status',sur_time='time',group='group',title='TCGA-LIHC',legene_title='Risk_score',legene_labels=c('high','low'))
ggplot2_KM_plot<-function(mydata,projectid='TCGA',sur='OS_STATUS',sur_time='OS_MONTHS',group='group',title='TCGA',y_title='OS',
legene_title=NA,legene_labels=NULL,output='KM')
{
      library(survival)
      library(survminer)
      #mydata<-df_cox_used
      formula <- as.formula(paste0('Surv(',sur_time,', ',sur,')~', group, sep = '', collapse = ' + '))
      diff=survdiff(formula=formula,data = mydata)
      diff$call$formula <- formula
      pValue=1-pchisq(diff$chisq,df=1)
      print(pValue)
      pvalue_string<-pValue
      if(pValue<0.001){
        pvalue_string="Pvalue<0.001"
      }else{
        pvalue_string=paste0("Pvalue=",sprintf("%.03f",pValue))
      }
      cols<-c("blue","red",'green')
      ngroup<-length(unique(mydata[,group]))
      fit <- survfit(formula=formula, data = mydata)
      fit$call$formula <- formula
      surPlot=ggsurvplot(fit,
                         data=mydata,
                         title=title,
                         pval=pvalue_string,
                         pval.size=6,
                         legend.labs=legene_labels,
                         legend.title=NULL,
                         font.legend=12,
                         xlab="Time(Months)",
                         ylab=y_title,
                         break.time.by = 12,
                         palette=cols[1:ngroup],
                         conf.int=F,
                         fontsize=4,
                         risk.table=F,
                         risk.table.title="",
                         risk.table.height=.25)
      #tiff(file=paste0(sur,".",i,".tiff"),width = 6, height =5)
      tiff(paste0(output,".a.tiff"),width=2400,height=2000,units="px",res=300,compression='lzw')
      print(surPlot)
      dev.off()
      pdf(file=paste0(output,".a.pdf"),width=8,height =6.6)
      print(surPlot)
      dev.off()
      
      surPlot=ggsurvplot(fit,
                         data=mydata,
                         title=title,
                         pval=pvalue_string,
                         pval.size=6,
                         legend.labs=legene_labels,
                         legend.title=NULL,
                         font.legend=12,
                         xlab="Time(Months)",
                         ylab=y_title,
                         break.time.by = 12,
                         palette=cols[1:ngroup],
                         conf.int=F,
                         fontsize=4,
                         risk.table=T,
                         risk.table.title="",
                         risk.table.height=.25)
      #tiff(file=paste0(sur,".",i,".tiff"),width = 6, height =5)
      tiff(paste0(output,".b.tiff"),width=2400,height=2400,units="px",res=300,compression='lzw')
      print(surPlot)
      dev.off()
      pdf(file=paste0(output,".b.pdf"),width=8,height =8)
      print(surPlot)
      dev.off()
      return(pValue)
    }


### venn_genes<-ggplot2_venn2(genes,deg_genes,name1='Transcription factors',name2='Markers of Macro_CXCL14',colors='default',output='fig0.venn2.set1_venn_set2')
ggplot2_venn2<-function(set1,set2,name1='set1',name2='set2',colors='default',output='fig0.venn2.set1_venn_set2',text_size=6,set_name_size = 8)
{
  library(ggvenn)
  library(patchwork)
  
  set1<-unique(set1)
  set2<-unique(set2)
  list_venn <- list(`set1` = set1,`set2` = set2)
  names(list_venn)<-c(name1,name2)
  if(colors == 'default')
  {
    plota<-ggvenn(list_venn,c(name1, name2),show_percentage = FALSE,text_size=text_size,set_name_size = set_name_size)+scale_y_continuous(expand = c(0, 0.2))
  }else if(colors == 'pink')
  {
    plota<-ggvenn(list_venn,c(name1, name2),show_percentage = FALSE,text_size=text_size,set_name_size = set_name_size,
                  fill_color = c("#f8bbcf", "#b1bbf9"),stroke_color='grey') +scale_y_continuous(expand = c(0, 0.2))
  }else{
    print('The colors are not set correctly!!!!!!!!!!!!!!!!!!!!!!!!')
    return()
  }
  ggsave(paste0(output,'.tiff'),plot=plota, width = 6, height = 6,compression='lzw')
  ggsave(paste0(output,'.pdf'),plot=plota, width = 6, height = 6)
  venn_genes<-set1[set1 %in% set2]
  write.table(venn_genes,paste0(output,'.tsv'),quote=F,col.names = F,row.names = F)
  return(venn_genes)
}

ggplot2_venn4<-function(list_venn,colors='default',output='fig0.venn4')
{
  library(ggvenn)
  library(patchwork)

  if(colors == 'default')
  {
    plota<-ggvenn(list_venn,show_percentage = FALSE,text_size=4,set_name_size = 4)+
      scale_x_continuous(expand = c(0.3, 0.3))+
      scale_y_continuous(expand = c(0, 0.2))
  }else if(colors == 'pink')
  {
    plota<-ggvenn(list_venn,show_percentage = FALSE,text_size=6,set_name_size = 8,
                  fill_color = c("#f8bbcf", "#b1bbf9"),stroke_color='grey') +scale_y_continuous(expand = c(0, 0.2))
  }else{
    print('The colors are not set correctly!!!!!!!!!!!!!!!!!!!!!!!!')
    return()
  }
  ggsave(paste0(output,'.tiff'),plot=plota, width = 6, height = 6,compression='lzw')
  ggsave(paste0(output,'.pdf'),plot=plota, width = 6, height = 6)
}


### ggplot2_upset(list_meta,ph=6,pw=10,output='fig1.upset')
ggplot2_upset<-function(list_venn,ph=6,pw=6,output='output')
{
  library(UpSetR)
  library(viridis)
  mycolors<-viridis(56, option = "D")
  nlist<-length(list_venn)
  plota<-upset(fromList(list_venn),nsets = nlist,
        sets.bar.color=brewer.pal(nlist,"Set1"),
        main.bar.color=mycolors)
  
  tiff(file=paste0(output,".tiff"), height=ph*300, width=pw*300,res=300,compression='lzw')
  print(plota)
  dev.off()
  pdf(file=paste0(output,".pdf"), height=ph, width=pw)
  print(plota)
  dev.off()
}



#df_data$Taxon<-factor(df_data$Taxon,levels=rev(df_data$Taxon))
#levels(df_data$Taxon)
#df_plot <- melt(df_data)
#head(df_plot)
#colnames(df_plot)<-c('Taxon','sample','value')
#head(df_plot)
# ggplot2_barplot_paired(df_plot,output=paste0('fig04.table1','.barplot_paired'),factorx='Taxon',groupby='sample',value='value',pw=16,ph=8)


ggplot2_count_table<-function(df_data,byX='malignant_type',byY=NA,output='output')
{
    if(is.na(byY))
    {
      df_count<-as.data.frame(table(df_data[[byX]]))
      colnames(df_count)[1]<-byX
      write.table(df_count,file=paste0(output,'.tsv'),row.names = F,sep='\t',quote = FALSE)
      ggplot2_barplot(df_count,group=byX,value='Freq',title=byX,y_title='',output=paste0(output,'.barplot'))
    }else{
          #  df_data<-metadata  byY='patient'
          df_data_used<-df_data[,c(byX,byY)]
          df_count<-as.data.frame(table(df_data_used))
          write.table(df_count,file=paste0(output,'.tsv'),row.names = F,sep='\t',quote = FALSE)
          
          ng<-length(unique(df_count[[byX]]))
          np<-length(unique(df_count[[byY]]))
          ph=max(np*0.25,10)
          pw=max(ng,10)
          pdf(file=paste0(output,'.pdf'), height=ph, width=pw)
          print(gplots::balloonplot(table(df_data_used[,c(byX,byY)])))
          dev.off()
          tiff(file=paste0(output,'.tiff'), height=ph*300, width=pw*300,units="px",res=300,compression='lzw')
          print(gplots::balloonplot(table(df_data_used[,c(byX,byY)])))
          dev.off()
    }
}

### ggplot2_barplot(df_plot,group=NULL,value=sample,title='celltype',y_title='Percentage',output=paste0("fig",plotid,"a.barplot.Ncells"))
ggplot2_barplot <- function(df_plot,group=NULL,value='value',title='catalog',y_title='value',output=paste0('fig04','.barplot'),pw=8,ph=6)
{
  library(ggplot2)
  if(is.null(group))
  {
    df_plot$group<-rownames(df_plot)
    group<-'group'
  }
  plota<-NULL
  plota<-ggplot(data=df_plot,aes(x=.data[[group]],weight=.data[[value]],fill=.data[[group]])) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold"),legend.position = 'none') +
    labs(title=title,y=y_title)
  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)
}


ggplot2_barplot_paired <- function(df_plot,output=paste0('fig04','.barplot_paired'),factorx='Taxon',groupby='sample',value='value',pw=8,ph=6)
{
  library(ggplot2)
  library(reshape2)
  
  library(ggsci)
  ngroup<-length(unique(df_plot[[groupby]]))
  mycolors= pal_npg("nrc")(ngroup) 
  print(levels(df_plot$type))
  plota<-NULL
  plota<- ggplot(data=df_plot)+geom_bar(aes(x=.data[[factorx]],weight=.data[[value]],fill=.data[[groupby]]),
                                        position = position_dodge2(preserve = 'single',reverse=F)) ### +coord_flip()
  plota<-plota+theme(plot.margin=unit(c(1,1,1,1), 'lines'),
                     panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(),axis.line = element_line(colour = 'black'),
                     axis.title.x=element_blank(),                     axis.title.y=element_blank(),
                     axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), 
                     axis.text.y=element_text(size=20,color='black'),
                     legend.position='right',legend.title=element_blank(),
                     legend.text = element_text(size = 20),legend.key.size = unit(1.2, 'lines'))
  plota<-plota+labs(y = '',x='')
  #+labs(y = '',x='')
  plota<-plota+scale_fill_manual(values=mycolors)+scale_y_continuous(trans='log10',expand = c(0, 0.1))		
  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)
}

ggplot2_barplot_paired_16S <- function(df_plot,output=paste0('fig04','.barplot_paired'),factorx='Taxon',groupby='sample',value='value',pw=8,ph=6)
{
  library(ggplot2)
  library(reshape2)
  ngroup<-length(unique(df_plot[[groupby]]))
  plota<-NULL
  plota<- ggplot(data=df_plot)+geom_bar(aes(x=.data[[factorx]],weight=.data[[value]],fill=.data[[groupby]]),
                                        position = position_dodge2(preserve = 'single',reverse=T))+coord_flip()
  plota<-plota+theme(plot.margin=unit(c(1,1,1,1), 'lines'),
                     panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(),axis.line = element_line(colour = 'black'),
                     axis.title.x=element_blank(),                     axis.title.y=element_blank(),
                     axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), 
                     axis.text.y=element_text(size=20,color='black'),
                     legend.position='bottom',legend.title=element_blank(),
                     legend.text = element_text(size = 20),legend.key.size = unit(1.2, 'lines'))
  plota<-plota+labs(y = '',x='')
  #+labs(y = '',x='')
  plota<-plota+scale_fill_manual(values=rainbow(ngroup+2))+scale_y_continuous(expand = c(0, 0))		
  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)
}


### ggplot2_barplot_paired_with_errorbar(df_plot,value='Freq',x='group',y='celltype_sub',y_title='Cells',output=paste0("fig",plotid,".barplot.Ncells"))
ggplot2_barplot_paired_with_errorbar <- function(df_plot,value='Freq',x='group',y='celltype_sub',y_title='Cells',output='output',pw=8,ph=6)
{
	library(ggplot2)
	library(reshape2)

    Data_summary <- summarySE(df_plot, measurevar=value, groupvars=c(x,y))
    head(Data_summary)
    y_max<-max(df_plot[[value]])
    label_y<-max(Data_summary[[value]]+Data_summary[['se']])
    
    
    ngroup<-length(unique(df_plot[[x]]))    
    mycolors<-get_colors(ngroup)
    
    plota<- ggplot(df_plot, aes(x = .data[[y]], y = .data[[value]], fill = .data[[x]])) +
    geom_bar(color='black',stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7,width = 0.7)
    # plota<-plota+stat_compare_means(method = "wilcox.test",label = "..p.signif..",method.args = list(alternative = "two.sided"),label.y=label_y*1.1,hide.ns = T,show.legend = FALSE)
    plota<-plota+geom_errorbar(data=Data_summary,aes(ymin=.data[[value]], ymax=.data[[value]]+.data[['se']]), width=.5,position=position_dodge(.7))
    plota<-plota+theme(plot.margin=unit(c(1,1,1,1), 'lines'),
                       panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                       panel.background = element_blank(),axis.line = element_line(colour = 'black'),
                       axis.title.x=element_blank(),axis.title.y=element_text(size=20), 
                       axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), 
                       axis.text.y=element_text(size=10,color='black'),
                       legend.position='right',legend.title=element_blank(),
                       legend.text = element_text(size = 10),legend.key.size = unit(1.2, 'lines'))
    plota<-plota+labs(y = y_title)
    plota<-plota+scale_fill_manual(values=mycolors)+scale_y_continuous(expand = c(0, 0),limits=c(0,y_max))
    ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
    ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)
    
    plota<-NULL
    plota<- ggplot(df_plot, aes(x = .data[[y]], y = .data[[value]], fill = .data[[x]])) +
      geom_bar(color='black',stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7,width = 0.7)
    # plota<-plota+stat_compare_means(method = "wilcox.test",label = "..p.signif..",method.args = list(alternative = "two.sided"),label.y=label_y*1.1,hide.ns = F,show.legend = FALSE)
    plota<-plota+geom_errorbar(data=Data_summary,aes(ymin=.data[[value]], ymax=.data[[value]]+.data[['se']]), width=.5,position=position_dodge(.7))
    plota<-plota+theme(plot.margin=unit(c(1,1,1,1), 'lines'),
                       panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                       panel.background = element_blank(),axis.line = element_line(colour = 'black'),
                       axis.title.x=element_blank(),axis.title.y=element_text(size=20), 
                       axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), 
                       axis.text.y=element_text(size=10,color='black'),
                       legend.position='right',legend.title=element_blank(),
                       legend.text = element_text(size = 10),legend.key.size = unit(1.2, 'lines'))
    plota<-plota+labs(y = y_title)
    plota<-plota+scale_fill_manual(values=mycolors)+scale_y_continuous(expand = c(0, 0),limits=c(0,y_max))
    ggsave(paste0(output,".b.tiff"), plot = plota, width = pw, height = ph,compression='lzw')
    ggsave(paste0(output,".b.pdf"), plot = plota, width = pw, height = ph)
}


#使用单细胞metadata画带标准误的柱状图
#x轴是celltype，柱子是group时tf = T，反之为F
#ggplot2_barplot_paired_with_errorbar_jzx(metadata,x='celltype_sub',group='group',tf = T,sample='orig.ident',plotid=plotid,pw=8,ph=6)
ggplot2_barplot_paired_with_errorbar_jzx<-function(metadata,x='celltype_manual',group='group',sample='orig.ident',tf=F,plotid=plotid,pw=8,ph=6)
{
  if(tf)
  {
	df_data<-as.data.frame(table(metadata[,c(sample,x)]))
	df_group<-unique(metadata[,c(sample,group)])
	data<-merge(df_data,df_group,by=sample)
  }else{
	df_data<-as.data.frame(table(metadata[,c(sample,group)]))
	df_group<-unique(metadata[,c(sample,x)])
	data<-merge(df_data,df_group,by=sample)
  }
  ggplot2_barplot_paired_with_errorbar(data,value='Freq',x=group,y=x,y_title='Cell count',output=paste0("fig",plotid,".b.",x,'-',group,"_with_errorbar"),pw=8,ph=6)
  
  # 生成数据摘要
  data_summary <- summarySE(data, measurevar='Freq', groupvars=c(x,group))
  
  write.table(data_summary,file=paste0("fig",plotid,".a.",x,'_',group,"_with_errorbar.tsv"),quote=F,sep="\t", row.names=F)
  
  plota<-ggplot(data_summary, aes(x = .data[[x]], y = Freq, fill = .data[[group]],color = .data[[group]]), color = .data[[group]]) +
	geom_errorbar(aes(ymin = Freq - se, ymax = Freq + se), 
				  position = position_dodge(width = 0.83), 
				  width = 0.4,
				  show.legend = FALSE) + 
	geom_bar(stat = "identity", 
			 position = position_dodge(width = 0.83), 
			 width = 0.8) + 
	labs(y = NULL) +            
	theme_bw(base_size = 10) + 
	theme(
	  axis.title = element_blank(),       # 隐藏坐标轴标题
	  panel.border = element_blank(),     # 去掉包围线
	  axis.line = element_line(),         # 显示 x 和 y 轴线
	  panel.grid = element_blank(),       # 去掉网格线
	  axis.text.x = element_text(colour = 'black',angle = 45,hjust = 1),  # 刻度字体颜色为黑色
	  legend.title=element_blank()
	) +
	scale_y_continuous(expand = c(0, 0),limits = c(0,max(data_summary$Freq) + max(data_summary$sd)))   # 柱子底部贴紧 x 轴
  ggsave(paste0("fig",plotid,".a.",x,'-',group,"_with_errorbar.tiff"), plot = plota, width = pw, height = ph,compression='lzw')
  ggsave(paste0("fig",plotid,".a.",x,'-',group,"_with_errorbar.pdf"), plot = plota, width = pw, height = ph)
}
		
### ggplot2_barplot_pn(DEG_sig4,gene='gene',group='change',ph=6,pw=6,output='fig8.barplot_pn')
#### 绘制正负两个方向的柱状图，（一般表示基因的上下调倍数）
ggplot2_barplot_pn<-function(df_deg,gene='gene',group='change',ph=6,pw=6,output='barplot_pn')
{
	plota<-ggplot(df_deg, aes(x = gene, y = logFC,fill=change)) + geom_bar(stat = "identity", show.legend = FALSE,width = .5)+geom_hline(yintercept = 0,size=1)
	plota<-plota+theme(plot.margin=unit(c(1,1,3,1.5), 'lines'),
					  # axis.title.x=element_text(size=30,color='black',vjust=-5,hjust=0.5),
					   axis.title.y=element_text(size=30,color='black'),
		  axis.text.x =element_blank(), 
		  axis.text.y=element_text(size=20,color='black'),
		  axis.ticks.x = element_blank(), 
		  #axis.text.y = element_blank(),
		  legend.position='None')
	plota<-plota+labs(y = 'logFC',x='')
	plota<-plota+theme(panel.grid.major =element_blank(), 
					   panel.grid.minor = element_blank(), 
					   panel.background = element_blank(),
					   axis.line.y = element_line(colour = 'black'),
					   axis.line.x = element_blank())
	plota<-plota+scale_fill_manual(values=c('green','red'))+scale_y_continuous(expand = c(0, 0))
	plota<-plota+geom_text(data = subset(DEG_sig4, logFC < 0),aes(x=gene, y= 0.07, label= gene), color = 'black',size = 5,vjust = 'inward' ) + 
	geom_text(data = subset(DEG_sig4, logFC > 0),aes(x=gene, y= -0.0575, label=gene), color = 'black',size = 5, vjust = 'outward')
	plota
	ggsave(paste0(output,".tiff"),plot=plota,width = pw, height = ph,compression='lzw')
	ggsave(paste0(output,".pdf"),plot=plota,width = pw, height = ph)
}

## ggplot2_density(df,value='score',output='fig7a.density.score')
ggplot2_density<-function(df_plot,value='score',output='output')
{
	plota<-NULL
	plota<-ggplot(df_plot, aes(x=.data[[value]])) + geom_density(alpha=.2, fill="#FF6666")+theme_classic()
	ggsave(paste0(output,".pdf"), plot = plota, width = 6, height = 6,limitsize=F)
	ggsave(paste0(output,".tiff"), plot = plota, width = 6, height = 6,limitsize=F,compression='lzw')
}
#df_plot<-data.frame(group=c(rep('control',2),rep('test',2)),type=c('common','unique','common','unique'),num=c(100,200,100,500))
#df_plot$type<-factor(df_plot$type,c('unique','common'))
#ggplot2_barplot_strack_with_label(df_plot,plotid='04',x='group',y='num',type='type')

### ggplot2_barplot_strack(df_plot,x='sample',y='percentage',type='celltype',pw=8,ph=6,output=paste0("fig",plotid,".",'celltype_percentage.barplot_stack'))
### ggplot2_barplot_strack_with_label(df_plot,x='sample',y='percentage',type='celltype',pw=8,ph=6,output=paste0("fig",plotid,".",'celltype_percentage.barplot_stack.labled'))

ggplot2_barplot_strack<-function(df_plot,x='group',y='num',type='type',pw=10,ph=6,output='output')
{
  ncelltype<-length(unique(df_plot[[x]]))
  ngroup<-length(unique(df_plot[[type]]))
  
  plota<-NULL
  plota<- ggplot(data=df_plot,aes(x=.data[[x]],y=.data[[y]],fill=.data[[type]]))+geom_bar(position='stack',stat='identity',width=0.6)
  #plota<-plota+geom_text(aes(label = .data[[y]]),position=position_stack(0.5),col="black",size=5)
  plota<-plota+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),
                     axis.title.x=element_text(size=20,color='black',vjust=-5,hjust=0.5),
                     axis.title.y=element_text(size=20,color='black'),
                     axis.text.x =element_text(size=15,color='black',angle = 45,hjust = 1), 
                     axis.text.y=element_text(size=15,color='black'),legend.position='right',
                     legend.title=element_blank(),legend.text = element_text(size = 10),
                     legend.key.size = unit(1.2, 'lines'))
  plota<-plota+labs(y =y,x='')+theme(panel.grid.major =element_blank(), 
                                     panel.grid.minor = element_blank(), panel.background = element_blank(),
                                     axis.line = element_line(colour = 'black'))
  if(ngroup == 2)
  {
    plota<-plota+scale_fill_manual(values=c('#BC3C29','#0072b5'))
  }
  plota<-plota+scale_y_continuous(expand = c(0, 0))
  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)
}


ggplot2_barplot_strack_percentage<-function(df_plot,x='group',y='num',type='type',pw=10,ph=6,output='output',colors=NULL)
{
  ncelltype<-length(unique(df_plot[[x]]))
  ngroup<-length(unique(df_plot[[type]]))
  
  plota<-NULL
  plota<- ggplot(data=df_plot,aes(x=.data[[x]],y=.data[[y]],fill=.data[[type]]))+geom_bar(position='stack',stat='identity',width=0.9)
  #plota<-plota+geom_text(aes(label = .data[[y]]),position=position_stack(0.5),col="black",size=5)
  plota<-plota+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),
                     axis.title.x=element_text(size=20,color='black',vjust=-5,hjust=0.5),
                     axis.title.y=element_text(size=20,color='black'),
                     axis.text.x =element_text(size=15,color='black',angle = 45,hjust = 1), 
                     axis.text.y=element_text(size=15,color='black'),legend.position='right',
                     legend.title=element_blank(),legend.text = element_text(size = 10),
                     legend.key.size = unit(1.2, 'lines'))
  plota<-plota+labs(y ='',x='')+theme(panel.grid.major =element_blank(), 
                                     panel.grid.minor = element_blank(), panel.background = element_blank(),
                                     axis.line = element_line(colour = 'black'))
  if(ngroup == 2)
  {
    plota<-plota+scale_fill_manual(values=c('#BC3C29','#0072b5'))
  }
  if(!is.null(colors))
  {
    plota<-plota+scale_fill_manual(values=colors)
  }
  #plota<-plota+scale_y_continuous(expand = c(0, 0))
  plota<-plota+scale_y_continuous(expand = c(0, 0.05),breaks=c(0,0.25,0.5,0.75,1),labels=c('0%','25%','50%','75%','100%'))
  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)
}

ggplot2_barplot_strack_with_label<-function(df_plot,x='group',y='num',type='type',pw=10,ph=6,output='output')
{
  ncelltype<-length(unique(df_plot[[x]]))
  ngroup<-length(unique(df_plot[[type]]))
  
  plota<-NULL
  plota<- ggplot(data=df_plot,aes(x=.data[[x]],y=.data[[y]],fill=.data[[type]]))+geom_bar(position='stack',stat='identity',width=0.6)
  plota<-plota+geom_text(aes(label = round(.data[[y]],2)),position=position_stack(0.5),col="black",size=5)
  plota<-plota+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),
                     axis.title.x=element_text(size=20,color='black',vjust=-5,hjust=0.5),
                     axis.title.y=element_text(size=20,color='black'),
                     axis.text.x =element_text(size=15,color='black',angle = 45,hjust = 1), 
                     axis.text.y=element_text(size=15,color='black'),legend.position='right',
                     legend.title=element_blank(),legend.text = element_text(size = 10),
                     legend.key.size = unit(1.2, 'lines'))
  plota<-plota+labs(y =y,x='')+theme(panel.grid.major =element_blank(), 
                                     panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = 'black'))
  
  if(ngroup == 2)
  {
    plota<-plota+scale_fill_manual(values=c('#BC3C29','#0072b5'))
  }
  plota<-plota+scale_y_continuous(expand = c(0, 0))
  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)
}

#ggplot2_pie(mydata,group='catalog',count='num',output='ggplot2_pie')
ggplot2_pie<-function(mydata,group='catalog',count='num',output='ggplot2_pie'){
	library(ggplot2)
	#colnames(mydata) <- c('catalog', 'num')
	mydata <- mydata[order(mydata[[count]],decreasing=T),]
	mydata[[group]]<-factor(mydata[[group]],levels=mydata[[group]])
	myLabel<-mydata[[group]]
	myLabel = paste0(myLabel,'\n', round(mydata[[count]] / sum(mydata[[count]]) * 100, 2), '%', sep = '') 
	res<-rev(mydata[[count]])
	labs<-rev(myLabel)
	col = rainbow(nrow(mydata))
	plota<-NULL
	plota <-ggplot(data=mydata,aes(x='',y=.data[[count]],fill=.data[[group]]))
	plota <-plota+geom_bar(stat='identity',width=2,color='black')
	plota <-plota+coord_polar(theta = 'y',direction=-1)
	plota <-plota+labs(x='',y='',title='')
	plota <-plota+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.ticks = element_blank(),
	legend.position='none',axis.text.x = element_blank())
	plota <-plota+geom_text(aes(y=res/2+c(0,cumsum(res)[-length(res)]),x=2.3,label = labs), size = 5,hjust=0.5)
	plota <-plota+scale_fill_manual(values=col)
	ggsave(paste0(output,".tiff"), plot = plota, width = 8, height = 8,compression='lzw')
	ggsave(paste0(output,".pdf"), plot = plota, width = 8, height = 8)
}

ggplot2_pie_with_legend<-function(mydata,group=NULL,count='num',title='catalog',output='ggplot2_pie'){
	library(ggplot2)
	if(is.null(group))
	{
		mydata$group<-rownames(mydata)
		group<-'group'
	}
	#colnames(mydata) <- c('catalog', 'num')
	mydata <- mydata[order(mydata[[count]],decreasing=T),]
	mydata[[group]]<-factor(mydata[[group]],levels=mydata[[group]])
	myLabel<-mydata[[group]]
	myLabel = paste0(myLabel,' ',round(mydata[[count]] / sum(mydata[[count]]) * 100, 2), '%', sep = '') 
	res<-rev(mydata[[count]])
	labs<-rev(myLabel)
	col = rainbow(nrow(mydata))
	plota<-NULL
	plota <-ggplot(data=mydata,aes(x='',y=.data[[count]],fill=.data[[group]]))
	plota <-plota+geom_bar(stat='identity',width=2,color='black')
	plota <-plota+coord_polar(theta = 'y',direction=-1)
	plota <-plota+labs(x='',y='',title='')
	plota <-plota+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.ticks = element_blank(),
	legend.position='right',axis.text.x = element_blank())
	plota <-plota+scale_fill_manual(title,labels=myLabel,values=col)
	ggsave(paste0(output,".tiff"), plot = plota, width = 8, height = 6,compression='lzw')
	ggsave(paste0(output,".pdf"), plot = plota, width = 8, height = 6)
}
 

  
## ggplot2_boxplot(log_normalized_count,group='sample',output=paste0('fig2.',projectid,'.DESeq2.log.boxplot'))
ggplot2_boxplot_matrix<-function(matrix_data,group='sample',samplename=F,ytitle='count',pw=8,output='ggplot2_boxplot'){
  library(ggplot2)
  library(reshape2)
  df_plot<-melt(matrix_data)
  colnames(df_plot)<-c(group,'count')
  plota<-NULL
  plota<-ggplot(df_plot,aes(x=.data[[group]],y=count,fill=.data[[group]]))
  plota<-plota+geom_boxplot(size=0.1,outlier.size=0.3)
  plota<-plota+theme_classic()+labs(y=ytitle)
  plota<-plota+theme(legend.position = 'none')
  if(! samplename)
  {
	plota<-plota+theme(axis.text.x=element_blank())
	}else{
		plota<-plota+theme(axis.text.x=element_text(angle=30,hjust=1))
	}
  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = 8,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = 8)
}

## ggplot2_boxplot_rna_seq(data_ref,phenotype_ref,group='group',ytitle='log(normalized_count)',output=paste0('fig10.test.boxplot'))
ggplot2_boxplot_rna_seq<-function(df_data,phenotype,group='group',samplename=F,ytitle='count',pw=8,output='ggplot2_boxplot'){
  #df_data<-df_data_boxplot
  #phenotype<-phenotype_boxplot
  library(ggplot2)
  library(reshape2)
  df_plot<-melt(df_data)
  colnames(df_plot)<-c('sample','count')
  head(df_plot)
  df_plot_used<-merge(df_plot,phenotype,by.x='sample',by.y=0)
  nsample<-dim(df_data)[2]
  
  plota<-NULL
  plota<-ggplot(df_plot_used,aes(x=sample,y=count,fill=.data[[group]]))
  plota<-plota+geom_boxplot(size=0.1,outlier.size=0.3)
  plota<-plota+theme_classic()+labs(y=ytitle)
  plota<-plota+theme(legend.position = 'none')
  pw<-max(8,nsample/5)
  if(! samplename)
  {
    plota<-plota+theme(axis.text.x=element_blank())
  }else{
    plota<-plota+theme(axis.text.x=element_text(angle=30,hjust=1))
  }
  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = 8,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = 8)
}


ggplot2_boxplot<-function(df_plot,group='sample',output='ggplot2_boxplot'){
  library(ggplot2)
  library(reshape2)
  plota<-NULL
  plota<-ggplot(df_plot,aes(x=.data[[group]],y=count,fill=.data[[group]]))
  plota<-plota+geom_boxplot(size=0.1,outlier.size=0.3)
  plota<-plota+theme_classic()
  plota<-plota+theme(axis.text.x=element_text(angle=30,hjust=1),legend.position = 'NA')
  ggsave(paste0(output,".tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0(output,".pdf"), plot = plota, width = 8, height = 8)
}

#####  数据格式为ggplot2的标准格式，数据在一列，group在一列。
## ggplot2_violin(mydata,score='CNV_score',group='sample',output=paste0('fig2.',projectid,'.DESeq2.log.boxplot'))
ggplot2_violin<-function(mydata,score='CNV_score',group='sample',output='ggplot2_boxplot'){
    library(ggplot2)
    library(reshape2)
    df_plot<-mydata[,c(score,group)]
    df_plot[,group]<-as.factor(df_plot[,group])
    ngroup<-length(unique(df_plot[,group]))
    pw<-max(6,ceiling(ngroup/2))
    if(ngroup>40)
    {
		pw<-20
	}
    plota<-NULL
    plota<-ggplot(df_plot,aes(x=.data[[group]],y=.data[[score]],fill=.data[[group]]))
    plota<-plota+geom_violin(scale = "width",trim = F,linewidth=0.5)
    plota<-plota+geom_boxplot(col='black',linewidth=0.2,linetype = 1,na.rm = T,notch = F,width = 0.3,outlier.shape = NA,position=position_dodge(0.9))
    plota<-plota+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                                  panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                                  plot.margin=unit(c(1,1,0.5,1.5), 'lines'),legend.position = 'none',
                                  axis.text.x=element_text(size=10,color='black',angle=30,hjust=1),
                                  axis.title.x=element_blank(),
                                  axis.text.y=element_text(size=10,color='black'),
                                  axis.title.y=element_text(size=15,color='black'))
    ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = 6,compression='lzw')
    ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = 6)
    return(plota)
  }


#####  数据格式为ggplot2的标准格式，数据在一列，group在一列。
# my_comparisons <- list( c("Normal", "LGSC"), c("Normal", "HGSOC"), c("LGSC", "HGSOC"))
# ggplot2_boxplot_Pvalue(df_plot=df_box,group='group',score='Inflammatory_score',output='fig30.Inflammatory_score',comparisons=my_comparisons,pw=4,ph=6)
ggplot2_boxplot_Pvalue<-function(df_plot,group='group',score='score',output='ggplot2_boxplot',comparisons=my_comparisons,pw=4,ph=6)
{
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  plota<-NULL
  plota<-ggplot(df_plot, aes(x=.data[[group]], y=.data[[score]],color=.data[[group]])) + 
  stat_boxplot(geom = "errorbar", width=0.35,lwd = 0.9) +
  geom_boxplot(width=0.35,outlier.shape = NA,lwd = 0.75) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",hide.ns = F,show.legend = FALSE,)+
  labs(title = colnames(df_plot[group]))+
  theme(panel.background = element_blank(),  # 去除背景
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank(),  # 去除次要网格线
        axis.line = element_line(colour = "black"),  # 显示xy轴
        axis.text.x = element_text(angle = 45, hjust = 0.6, vjust = 0.6,size = 10),  # 旋转x轴标签
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 20),
        plot.title=element_text(size=20,hjust=0.5,vjust=0.5),
        legend.position = "none",
        panel.border = element_blank())  # 去除边框
 print('good good!!!')
 ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
 ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)
}

### ggboxplot_add_pvalue(df_plot,my_comparisons,group='change',value='percentage',output='fig50.boxplot.codon_percentage',y_title='Percentage',title='',pw=4,ph=6,label_y=NA)
ggboxplot_add_pvalue<-function(df_plot,my_comparisons,group='group2',value='value',
							 y_title='celltype',title='celltype',output='ggplot2_boxplot',
							 label_y=NA,pw=8,ph=6){
library(ggplot2)
library(reshape2)
library(ggpubr)

ymin<-min(df_plot[[value]])
ymax<-max(df_plot[[value]])
yh<-ymax-ymin

groups<-unique(df_plot[,group])
ngroup<-length(groups)
mycolors<-get_colors(ngroup)
plota<-NULL
plota<-ggboxplot(df_plot,x=group,y=value,fill=group,add = "jitter",add.params = list(size=0.5))
if(is.na(label_y))
{
  plota<-plota+stat_compare_means(comparisons = my_comparisons,label = "p.format",hide.ns = TRUE,show.legend = FALSE)
}else
{
  plota<-plota+stat_compare_means(comparisons = my_comparisons,label.y = label_y,label = "p.format",hide.ns = TRUE,show.legend = FALSE)
}
plota<-plota+scale_fill_manual(values=mycolors)
plota<-plota+scale_color_manual(values=mycolors)
plota<-plota+theme_classic()
plota<-plota+theme(text=element_text(size=10),plot.title = element_text(hjust = 0.5,size=15),
				   axis.text.x=element_text(size=10,color='black',angle=45,hjust=1),
				   axis.text.y=element_text(size=10,color='black'),
				   axis.title.y=element_text(size=15,color='black'),legend.position = 'none',
				   legend.title =element_blank(),legend.text =element_text(size=15,color='black')  )
plota<-plota+labs(x='',y=y_title,title=title)
ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)

plotb<-NULL
plotb<-ggplot(df_plot,aes(x=.data[[group]],y=.data[[value]],fill=.data[[group]]))
plotb<-plotb+geom_violin(scale = "width",size=0.1,trim = F)
plotb<-plotb+geom_boxplot(col='black',fill='white',linewidth=0.2,linetype = 1,na.rm = T,notch = F,width = 0.3,outlier.shape = NA,position=position_dodge(0.9))
#plotb<-plotb+stat_compare_means(comparisons = my_comparisons,label = "p.format",hide.ns = TRUE,show.legend = FALSE)
if(is.na(label_y))
{
  plotb<-plotb+stat_compare_means(comparisons = my_comparisons,label = "p.format",hide.ns = TRUE,show.legend = FALSE)
}else
{
  plotb<-plotb+stat_compare_means(comparisons = my_comparisons,label.y = label_y,label = "p.format",hide.ns = TRUE,show.legend = FALSE)
}
#plotb<-plotb+annotate(geom="text",label=paste0("Pvalue(W_vs_NC): ", round(stat1$p,4)),x=2,y=ymax+yh*0.1,size=5,color='red')
plotb<-plotb+scale_fill_manual(values=mycolors)
plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
							  panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
							  plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
							  text=element_text(size=10),plot.title = element_text(hjust = 0.5,size=15),
							  legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
							  legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
							  axis.text.x=element_text(size=10,color='black',angle=45,hjust=1),
							  axis.text.y=element_text(size=10,color='black'),
							  axis.title.y=element_text(size=15,color='black'))
plotb<-plotb+labs(x='',y=y_title,title=title)
#plotb<-plotb+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
ggsave(file=paste0(output,".b1.violin.tiff"),plot=plotb,width = pw, height = ph,compression='lzw')
ggsave(file=paste0(output,".b1.violin.pdf"),plot=plotb,width = pw, height = ph)



plotc<-NULL
plotc<-ggplot(df_plot,aes(x=.data[[group]],y=.data[[value]],color=.data[[group]]))
plotc<-plotc+geom_boxplot(width = 0.6,outlier.shape = NA)
plotc<-plotc+stat_boxplot(geom = "errorbar", width=0.2,lwd = 0.6) 
#plotc<-plotc+stat_compare_means(comparisons = my_comparisons,label = "p.format",hide.ns = TRUE,show.legend = FALSE)
if(is.na(label_y))
{
  plotc<-plotc+stat_compare_means(comparisons = my_comparisons,label = "p.format",hide.ns = TRUE,show.legend = FALSE)
}else
{
  plotc<-plotc+stat_compare_means(comparisons = my_comparisons,label.y = label_y,label = "p.format",hide.ns = TRUE,show.legend = FALSE)
}
#plotb<-plotb+annotate(geom="text",label=paste0("Pvalue(W_vs_NC): ", round(stat1$p,4)),x=2,y=ymax+yh*0.1,size=5,color='red')
plotc<-plotc+scale_color_manual(values=mycolors)
plotc<-plotc+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
							  panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
							  plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
							  text=element_text(size=10),plot.title = element_text(hjust = 0.5,size=15),
							  legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
							  legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
							  axis.text.x=element_text(size=10,color='black',angle=45,hjust=1),
							  axis.text.y=element_text(size=10,color='black'),
							  axis.title.y=element_text(size=15,color='black'))
plotc<-plotc+labs(x='',y=y_title,title=title)
#plotb<-plotb+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
ggsave(file=paste0(output,".b2.boxplot.tiff"),plot=plotc,width = pw, height = ph,compression='lzw')
ggsave(file=paste0(output,".b2.boxplot.pdf"),plot=plotc,width = pw, height = ph)

}

 
### 单个基因的boxplot，加绘一张蜜蜂图
### 数据包含两列：gene和group
# ggboxplot_2groups_comparison(df_plot,group='group',value=gene,title=gene,output=paste0('fig',plotid,'.boxplot.comparison.',gene),label_y=ymax+1)
ggboxplot_2groups_comparison_beeswarm<-function(df_plot,group='group',value='count',title='celltype',output='ggplot2_boxplot',label_y=0.8){
      library(ggplot2)
      library(reshape2)
      library(ggpubr)
      library(ggbeeswarm)
      if(class(df_plot[,group]) =='factor')
      {
        print('factor')
        groups<-levels(df_plot[,group])
        
      }else{
        groups<-unique(df_plot[,group])
        groups<-as.character(groups)
      }
      
      my_comparisons <- list( c(groups))
      plota<-NULL
      plota<-ggboxplot(df_plot,x=group,y=value,fill=group)
      plota<-plota+stat_compare_means(comparisons = my_comparisons,label = "p.format",label.y = label_y,hide.ns = T,show.legend = FALSE)
      plota<-plota+scale_fill_manual(values=c('orange','red'))
      plota<-plota+scale_color_manual(values=c('orange','red'))
      plota<-plota+theme_classic()+labs(x='',y=title)
      plota<-plota+theme(axis.text.x=element_text(size=15,color='black',angle=45,hjust=1),
                         axis.text.y=element_text(size=15,color='black'),
                         axis.title.y=element_text(size=20,color='black'),legend.position = 'none',
                         legend.title =element_blank(),legend.text =element_text(size=15,color='black')  )
      ggsave(paste0(output,".tiff"), plot = plota, width = 6, height = 6,compression='lzw')
      ggsave(paste0(output,".pdf"), plot = plota, width = 6, height = 6)
      
      plotb<-NULL
      plotb<-ggplot(df_plot,aes(x=group,y=.data[[value]],col=group))+geom_beeswarm(cex=1.5)
      plotb<-plotb+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
                                geom = 'crossbar', width = 0.5, size = 0.5, color = 'black')
      plotb<-plotb+stat_compare_means(comparisons = my_comparisons,label = "p.format",label.y = label_y,hide.ns = T,show.legend = FALSE)
      plotb<-plotb+theme_bw()
      plotb<-plotb+theme(legend.position ='none',panel.grid = element_blank())
      plotb<-plotb+scale_color_manual(values=c("orange","red"),name="Censored",labels=c("Yes","No"))
      plotb<-plotb+labs(x="",y=value)
      plotb<-plotb+theme(axis.text.x=element_text(size=15,color='black',angle=45,hjust=1),
                         axis.text.y=element_text(size=15,color='black'),
                         axis.title.y=element_text(size=20,color='black'),legend.position = 'none',
                         legend.title =element_blank(),legend.text =element_text(size=15,color='black')  )
      ggsave(paste0(output,".beeswarm.tiff"), plot = plotb, width = 6, height = 6,compression='lzw')
      ggsave(paste0(output,".beeswarm.pdf"), plot = plotb, width = 6, height = 6)
    }

### 单个基因的boxplot，加绘一张蜜蜂图
### 数据包含两列：gene和group
# ggboxplot_2groups_comparison(df_plot,group='group',value=gene,title=gene,output=paste0('fig',plotid,'.boxplot.comparison.',gene),label_y=ymax+1)
ggboxplot_2groups_comparison<-function(df_plot,group='group',value='count',title='celltype',output='ggplot2_boxplot',label_y='max',pw=4){
		library(ggplot2)
		library(reshape2)
		library(ggpubr)
		library(ggbeeswarm)
		if(class(df_plot[,group]) =='factor')
		{
		  print('factor')
		  groups<-levels(df_plot[,group])
		  
		}else{
		  groups<-unique(df_plot[,group])
		  groups<-as.character(groups)
		}
		print(groups)
		ngroup<-length(groups)
		print(ngroup)

		if(label_y=='max')
		{
		  y_max<-max(df_plot[[value]])
		  y_label=y_max
		}else{
		  y_label=label_y
		}
		df_plot<-df_plot[,c(group,value)]
		colnames(df_plot)<-c('group','value')
		head(df_plot)
		if(ngroup==2)
		{
		  formula <- as.formula(paste0('value','~', 'group'))
		  df_res<-as.data.frame(compare_means(formula, data=df_plot))
		  my_comparisons <- list( c(groups))
		  
		  plota<-NULL
		  plota<-ggboxplot(df_plot,x='group',y='value',fill='group',width=0.5)
		  plota<-plota+stat_compare_means(comparisons = my_comparisons,label = "p.format",label.y = y_label,hide.ns = T,show.legend = FALSE)
		  plota<-plota+scale_fill_manual(values=c('orange','red'))
		  #plota<-plota+scale_color_manual(values=c('orange','red'))
		  plota<-plota+theme_classic()+labs(x='',y=title)
		  plota<-plota+theme(axis.text.x=element_text(size=15,color='black',angle=45,hjust=1),
							 axis.text.y=element_text(size=15,color='black'),
							 axis.title.y=element_text(size=20,color='black'),legend.position = 'none',
							 legend.title =element_blank(),legend.text =element_text(size=15,color='black')  )
		  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = 6,compression='lzw')
		  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = 6)
		  
		  plotb<-NULL
		  plotb<-ggplot(df_plot,aes(x=group,y=value,col=group))+geom_beeswarm(cex=1.5,corral.width=0.5)
		  plotb<-plotb+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
									geom = 'crossbar', width = 0.5, size = 0.5, color = 'black')
		  plotb<-plotb+stat_compare_means(comparisons = my_comparisons,label = "p.format",label.y = y_label,hide.ns = T,show.legend = FALSE)
		  plotb<-plotb+theme_bw()
		  plotb<-plotb+theme(legend.position ='none',panel.grid = element_blank())
		  plotb<-plotb+scale_color_manual(values=c("orange","red"))
		  plotb<-plotb+labs(x="",y=title)
		  plotb<-plotb+theme(axis.text.x=element_text(size=15,color='black',angle=45,hjust=1),
							 axis.text.y=element_text(size=15,color='black'),
							 axis.title.y=element_text(size=20,color='black'),legend.position = 'none',
							 legend.title =element_blank(),legend.text =element_text(size=15,color='black')  )
		  ggsave(paste0(output,".beeswarm.tiff"), plot = plotb, width = pw, height = 6,compression='lzw')
		  ggsave(paste0(output,".beeswarm.pdf"), plot = plotb, width = pw, height = 6)
		  return(df_res)
		}else
		{
		  mycolors<-get_colors(ngroup)
		  pw=max(ngroup,6)      
		  plota<-NULL
		  plota<-ggboxplot(df_plot,x='group',y='value',fill='group')
		  plota<-plota+stat_compare_means()
		  #plota<-plota+stat_compare_means(comparisons = my_comparisons,label = "p.format",label.y = y_label,hide.ns = T,show.legend = FALSE)
		  plota<-plota+scale_fill_manual(values=mycolors)
		  #plota<-plota+scale_color_manual(values=c('orange','red'))
		  plota<-plota+theme_classic()+labs(x='',y=title)
		  plota<-plota+theme(axis.text.x=element_text(size=15,color='black',angle=45,hjust=1),
							 axis.text.y=element_text(size=15,color='black'),
							 axis.title.y=element_text(size=20,color='black'),legend.position = 'none',
							 legend.title =element_blank(),legend.text =element_text(size=15,color='black')  )
		  ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = 6,compression='lzw')
		  ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = 6)
		  
		}
}



## 数据矩阵的列是列是样本，行是celltype或者基因。
#ggplot2_boxplot_paired_cibersort(expr, df_sample,groupby='group',plotid='fig11a',ymax=16,label_y=0.6,levels=NULL)
ggplot2_boxplot_paired_cibersort<- function(expr, phenotype,groupby='group',plotid='fig11a',ymax=NA,y_title='Proportion',label_y=0.6,levels=NULL)  ######  绘制hub基因的配对boxplot
{
  library(ggpubr)
  library(reshape2)
  #expr<-data
  phenotype_used<-phenotype[colnames(expr),,drop=F]
  #Type<-df_sample_used[[groupby]]
  data=as.data.frame(t(expr))
  data$group<-phenotype_used[[groupby]]
  df_plot<-melt(data)
  head(df_plot)
  colnames(df_plot)<-c('group','gene','expression')
  
  ymin<-min(df_plot$expression)
  ymax_value<-max(df_plot$expression)
  
  if(is.na(ymax))
  {
    ymax=ymax_value
  }else{
    ymax=ymax
  }
  
  if(!is.null(levels))
  {
    df_plot$gene<-factor(df_plot$gene,levels=levels)
  }
  
  ggplot2_boxplot_paired(df_plot,groupby=c('group','gene'),score='expression',y_title=y_title,pw=7,ph=4,output=paste0(plotid,'.boxplot'))
  
}



summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


##  使用4种方式展示配对boxplot，用于展示单细胞中的计算m7GScore的小提琴图
#ggplot2_boxplot_paired(df_plot,groupby=c('group','celltype'),score='expression',y_title='m1A Score',pw=7,ph=4,output='fig11a')
ggplot2_boxplot_paired<- function(df_plot,groupby=c('group','celltype'),score='expression',x_title='',y_title='m1A Score',pw=7,ph=4,output='fig11a',label_y='max')  ######  绘制hub基因的配对boxplot
{
  #devtools::install_github("JanCoUnchained/ggunchained")
  library(ggpubr)
  library(reshape2)
  library(ggunchained)
  head(df_plot) 
  #colnames(df_plot)[1]<-'score'
  df_plot<-df_plot[,c(score,groupby)]
  colnames(df_plot)[1:3]<-c('score','group','celltype')
  ymin<-min(df_plot$score)
  ymax<-max(df_plot$score)
  yh<-ymax-ymin
  
  		if(label_y=='max')
		{
		  y_max<-max(df_plot[[value]])
		  y_label=y_max
		}else{
		  y_label=label_y
		}
  
  mycolors<-rev(get_colors(2))
  
  
  
  plota<-NULL
  plota<-ggplot(df_plot,aes(x=celltype,y=score,fill=group))
  plota<-plota+geom_boxplot(col='black',size=0.1,outlier.size=0.3)
  plota<-plota+scale_fill_manual(values=mycolors)
  #plota<-plota+scale_color_manual(values=c('red','blue'))
  plota<-plota+stat_compare_means(aes(group = group),label = "p.signif",label.y = y_label,hide.ns = TRUE,show.legend = FALSE,size=1.5)
  plota<-plota+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                                panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position=c(0.9,0.95),legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.title.x=element_text(size=10,color='black',face="bold"))
  plota<-plota+labs(y =y_title,x=x_title)+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  ggsave(file=paste0(output,".a.boxplot.tiff"),plot=plota,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".a.boxplot.pdf"),plot=plota,width = pw, height = ph)
  
  plotb<-NULL
  plotb<-ggplot(df_plot,aes(x=celltype,y=score,fill=group))
  plotb<-plotb+geom_violin(scale = "width",size=0.1,trim = F)
  #plotb<-plotb+stat_summary(position=position_dodge(0.9),fun= mean, geom = "point",shape = 23, size = 2,col='black')
  plotb<-plotb+scale_fill_manual(values=mycolors)
  #plotb<-plotb+scale_color_manual(values=c('GO'='red', 'Normal'='blue'))
  plotb<-plotb+stat_compare_means(aes(group = group),label = "p.signif",label.y = y_label,hide.ns = TRUE,show.legend = FALSE,size=1.5)
  plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                                panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position=c(0.9,0.95),legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.title.x=element_text(size=10,color='black',face="bold"))
  plotb<-plotb+labs(y =y_title,x=x_title)+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  ggsave(file=paste0(output,".b.violin.tiff"),plot=plotb,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".b.violin.pdf"),plot=plotb,width = pw, height = ph)
  
  
  
  plotb<-NULL
  plotb<-ggplot(df_plot,aes(x=celltype,y=score,fill=group))
  plotb<-plotb+geom_violin(scale = "width",linewidth=0.2,trim = F)
  plotb<-plotb+geom_boxplot(col='black',linewidth=0.2,linetype = 1,na.rm = T,notch = F,width = 0.2,outlier.shape = NA,position=position_dodge(0.9))
  #plotb<-plotb+stat_boxplot(aes(ymin = ..lower..,ymax = ..upper..),size = 1,width = 0.3,notch = T,outlier.shape = NA,na.rm = T)
  #plotb<-plotb+stat_summary(position=position_dodge(0.9),fun= mean, geom = "point",shape = 23, size = 2,col='black')
  plotb<-plotb+scale_fill_manual(values=mycolors)
  #plotb<-plotb+scale_color_manual(values=c('GO'='red', 'Normal'='blue'))
  plotb<-plotb+stat_compare_means(aes(group = group),label = "p.signif",label.y = ymax,hide.ns = TRUE,show.legend = FALSE,size=1.5)
  plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                                panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position=c(0.9,0.95),legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.title.x=element_text(size=10,color='black',face="bold"))
  plotb<-plotb+labs(y =y_title,x=x_title)+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  ggsave(file=paste0(output,".b1.violin.tiff"),plot=plotb,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".b1.violin.pdf"),plot=plotb,width = pw, height = ph)
  
  
  plotc<-NULL
  plotc<-ggviolin(df_plot,x='celltype',y='score',fill='group',add='boxplot',xlab=x_title,ylab=y_title) 
  plotc<-plotc+scale_fill_manual(values = mycolors)
  plotc<-plotc+stat_compare_means(aes(group = group),label = "p.signif",label.y = ymax,hide.ns = TRUE,size=1.5)
  plotc<-plotc+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                                plot.title = element_text(hjust = 0.5),
                                panel.background = element_blank(),panel.border=element_blank(),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position=c(0.9,0.95),legend.title=element_blank(),
                                legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.title.x=element_text(size=10,color='black',face="bold"),
                                axis.line = element_line(linewidth=0.5, colour = "black"))
  plotc<-plotc+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))#limits=c(0,1)
  ggsave(file=paste0(output,".c.violin_withboxplot.tiff"),plot=plotc,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".c.violin_withboxplot.pdf"),plot=plotc,width = pw, height = ph)
  
  
  Data_summary <- summarySE(df_plot, measurevar="score", groupvars=c("group","celltype"))
  
  plotd<-NULL
  plotd <- ggplot(data=df_plot, aes(x=celltype, y=score,fill=group)) + 
    geom_split_violin(trim=F,color=NA) + #绘制分半的小提琴图
    geom_point(data = Data_summary,aes(x=celltype, y=score),pch=20,position=position_dodge(0.6),size=1)+ #绘制均值为点图
    geom_errorbar(data = Data_summary,aes(ymin = score-ci, ymax=score+ci), #误差条表示95%的置信区间
                  width=0.1, #误差条末端短横线的宽度
                  position=position_dodge(0.6), 
                  color="black",
                  alpha = 0.7,
                  size=0.5) +
    scale_fill_manual(values = mycolors) #设置填充的颜色
  plotd<-plotd+stat_compare_means(aes(group = group),label = "p.signif",label.y = ymax,hide.ns = TRUE,show.legend = FALSE,size=1.5)
  plotd<-plotd+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                                panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position=c(0.9,0.95),legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.title.x=element_text(size=10,color='black',face="bold"))
  plotd<-plotd+labs(y =y_title,x=x_title)+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  ggsave(file=paste0(output,".d.half_violin.tiff"),plot=plotd,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".d.half_violin.pdf"),plot=plotd,width = pw, height = ph)
}

ggplot2_boxplot_paired_flip<- function(df_plot,groupby=c('group','celltype'),score='expression',y_title='m1A Score',pw=7,ph=4,output='fig11a')  ######  绘制hub基因的配对boxplot
{
	  #devtools::install_github("JanCoUnchained/ggunchained")
	  library(ggpubr)
	  library(reshape2)
	  library(ggunchained)
	  head(df_plot) 
	  #colnames(df_plot)[1]<-'score'
	  df_plot<-df_plot[,c(score,groupby)]
	  colnames(df_plot)[1:3]<-c('score','group','celltype')
	  ymin<-min(df_plot$score)
	  ymax<-max(df_plot$score)
	  yh<-ymax-ymin
	  
	  mycolors<-rev(get_colors(2))
	  
	  plota<-NULL
	  plota<-ggplot(df_plot,aes(x=celltype,y=score,fill=group))
	  plota<-plota+geom_boxplot(col='black',size=0.1,outlier.size=0.3)
	  plota<-plota+scale_fill_manual(values=mycolors)
	  #plota<-plota+scale_color_manual(values=c('red','blue'))
	  plota<-plota+stat_compare_means(aes(group = group),label = "p.signif",label.y = ymax,hide.ns = TRUE,show.legend = FALSE,vjust=0.7)
	  plota<-plota+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
									panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
									plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
									legend.position='top',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
									legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
									axis.text.x=element_text(size=8,color='black',angle=45,hjust=1))
	  plota<-plota+labs(y =y_title,x='')+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))+coord_flip()
	  print(plota)
	  ggsave(file=paste0(output,".a.boxplot.tiff"),plot=plota,width = pw, height = ph,compression='lzw')
	  ggsave(file=paste0(output,".a.boxplot.pdf"),plot=plota,width = pw, height = ph)
}
        
##  使用3种方式展示配对boxplot，用于展示单细胞中的计算m7GScore的小提琴图
#ggplot2_boxplot_gsva(df_plot,groupby='celltype',score='expression',y_title='m1A Score',pw=7,ph=4,output='fig11a')
ggplot2_boxplot_gsva<- function(df_plot,groupby='celltype',score='expression',y_title='m1A Score',pw=7,ph=4,output='fig11a')  ######  绘制hub基因的配对boxplot
{
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
  head(df_plot) 
  df_plot<-df_plot[,c(score,groupby)]
  colnames(df_plot)[1:2]<-c('score','celltype')
  ymin<-min(df_plot$score)
  ymax<-max(df_plot$score)
  yh<-ymax-ymin
  
  plotb<-NULL
  plotb<-ggplot(df_plot,aes(x=celltype,y=score,fill=celltype))
  plotb<-plotb+geom_violin(scale = "width",size=0.1,trim=F)
  plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
                                panel.background = element_blank(),panel.border=element_blank(),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.line = element_line(linewidth=0.5, colour = "black"))
  plotb<-plotb+labs(x='',y =y_title)+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  ggsave(file=paste0(output,".a.violin.tiff"),plot=plotb,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".a.violin.pdf"),plot=plotb,width = pw, height = ph)
  
    
  plotb<-NULL
  plotb<-ggplot(df_plot,aes(x=celltype,y=score,fill=celltype))
  plotb<-plotb+geom_violin(scale = "width",linewidth=0.2,trim = F)
  plotb<-plotb+geom_boxplot(col='black',fill='white',linewidth=0.2,linetype = 1,na.rm = T,notch = F,width = 0.2,outlier.shape = NA,position=position_dodge(0.9))
  #plotb<-plotb+stat_boxplot(aes(ymin = ..lower..,ymax = ..upper..),size = 1,width = 0.3,notch = T,outlier.shape = NA,na.rm = T)
  #plotb<-plotb+stat_summary(position=position_dodge(0.9),fun= mean, geom = "point",shape = 23, size = 2,col='black')
  #plotb<-plotb+scale_fill_manual(values=c('red','blue'))
  #plotb<-plotb+scale_color_manual(values=c('GO'='red', 'Normal'='blue'))
  #plotb<-plotb+stat_compare_means(aes(group = group),label = "p.signif",label.y = ymax,hide.ns = TRUE,show.legend = FALSE)
  plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                                panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1))
  plotb<-plotb+labs(y =y_title,x='')+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  ggsave(file=paste0(output,".a1.violin.tiff"),plot=plotb,width =pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".a1.violin.pdf"),plot=plotb,width = pw, height = ph)
 
  plota<-ggviolin(df_plot,x='celltype',y='score',fill='celltype',add='boxplot',add.params = list(fill="white"),xlab='',ylab=y_title)+
  #scale_fill_manual(values = c("red", "blue"))+
  theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
                                panel.background = element_blank(),panel.border=element_blank(),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),
                                axis.line = element_line(linewidth=0.5, colour = "black"))+
  scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  ggsave(file=paste0(output,".a2.violin.tiff"),plot=plota,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".a2.violin.pdf"),plot=plota,width = pw, height = ph)
}


### ggplot2_boxplot_gsva_add_line(metadata,groupby='orig.ident',score='PANoptosis_score',cutoff=cutoff,y_title='PANoptosis_score',pw=7,ph=4,output=paste0("fig",plotid3,".boxplot"))
ggplot2_boxplot_gsva_add_line<- function(df_plot,groupby='celltype',cutoff=cutoff,score='expression',y_title='m1A Score',pw=7,ph=4,output='fig11a')  ######  绘制hub基因的配对boxplot
{
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
  head(df_plot) 
  df_plot<-df_plot[,c(score,groupby)]
  colnames(df_plot)[1:2]<-c('score','celltype')
  ymin<-min(df_plot$score)
  ymax<-max(df_plot$score)
  yh<-ymax-ymin
  
  plotb<-NULL
  plotb<-ggplot(df_plot,aes(x=celltype,y=score,fill=celltype))
  plotb<-plotb+geom_violin(scale = "width",size=0.1,trim=F)
  plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
                                panel.background = element_blank(),panel.border=element_blank(),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),axis.line = element_line(linewidth=0.5, colour = "black"))
  plotb<-plotb+labs(x='',y =y_title)+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  ggsave(file=paste0(output,".a.violin.tiff"),plot=plotb,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".a.violin.pdf"),plot=plotb,width = pw, height = ph)
  
  
  plotb<-NULL
  plotb<-ggplot(df_plot,aes(x=celltype,y=score,fill=celltype))
  plotb<-plotb+geom_violin(scale = "width",linewidth=0.2,trim = F)
  plotb<-plotb+geom_boxplot(col='black',fill='white',linewidth=0.2,linetype = 1,na.rm = T,notch = F,width = 0.2,
                            outlier.shape = NA,position=position_dodge(0.9))
  plotb<-plotb+geom_hline(yintercept = cutoff,linetype = "dashed", color = 'red') 
  plotb<-plotb+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                                panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1))
  plotb<-plotb+labs(y =y_title,x='')+scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  ggsave(file=paste0(output,".a1.violin.tiff"),plot=plotb,width =pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".a1.violin.pdf"),plot=plotb,width = pw, height = ph)
  
  plota<-ggviolin(df_plot,x='celltype',y='score',fill='celltype',add='boxplot',add.params = list(fill="white"),xlab='',ylab=y_title)+
    #scale_fill_manual(values = c("red", "blue"))+
    theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
                     panel.background = element_blank(),panel.border=element_blank(),
                     plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                     legend.position='none',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                     legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                     axis.text.x=element_text(size=8,color='black',angle=45,hjust=1),
                     axis.line = element_line(linewidth=0.5, colour = "black"))+
    scale_y_continuous(expand = c(0, 0),limits=c(ymin-yh*0.2,ymax+yh*0.2))
  plota<-plota+geom_hline(yintercept = cutoff,linetype = "dashed", color = 'red') 
  ggsave(file=paste0(output,".a2.violin.tiff"),plot=plota,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".a2.violin.pdf"),plot=plota,width = pw, height = ph)
}

### fun_to_corr: 
#### input: heat_in_1, heat_in_2: f_ra/g_ra/biochem/kegg_in...you can select the feature first
#### output: list(), t_cor is correlation index(-1~1), t_p is raw p-value, t_p_sig is formatted with significant levels
#### formatted significant levels: 0~0.001: **; 0.001~0.01: *; 0.01~0.1: +; >0.1 nothing
fun_to_corr <- function(heat_in_1, heat_in_2) {
  t <- corr.test(heat_in_1, heat_in_2, use = "pairwise", method = "spearman", adjust = "fdr")
  t_cor <- data.frame(t$r, check.names = FALSE)
  t_p <- data.frame(t$p, check.names = FALSE)
  cut_sig <- function(p) {
    out <- cut(p, breaks = c(0, 0.001,0.01,0.1,1), include.lowest = T, labels = c("**", "*", "+", ""))
    return(out)
  }
  t_p_sig <- apply(t_p, 2, cut_sig)
  rownames(t_p_sig) <- rownames(t_p)
  return(list(t_cor = t_cor, t_p_sig = t_p_sig, t_p = t_p))
}




#ggcorrplot_2matrix(t(df_tumor_logtpm_used),df_tcga_infiltration_used,ph.factor=2,output='fig2.ggcorrplot',save.data=T)
#ggcorrplot_2matrix(t(expra),t(exprb),output='fig2.ggcorrplot',pw=5,ph=20,save.data=T)
ggcorrplot_2matrix<-function(mx_a,mx_b,cor_method='pearson',output='fig.ggcorrplot',pw=NA,ph=NA,save.data=T,ploty=T){
  
  library(ggcorrplot)  
  library(ggthemes)  
  library(psych)  
  
      sds<-apply(mx_a,1,sd)
      mx_a_used<-mx_a[sds>0,]
      dim(mx_a_used)
      
      sds<-apply(mx_b,1,sd)
      mx_b_used<-mx_b[sds>0,]
      dim(mx_b_used)
      
  
  cor <- corr.test(mx_a_used,mx_b_used,use="complete",method = cor_method,adjust = "BH",ci = F)    #### pearson spearman两种方法计算相关性
  matrix_cor<-cor$r
  matrix_cor_p<-cor$p
  
  if(save.data)
  {
    matrix_cor<-as.data.frame(matrix_cor)
    matrix_cor_p<-as.data.frame(matrix_cor_p)
    matrix_cor_out<-cbind(rownames(matrix_cor),matrix_cor)
    matrix_cor_p_out<-cbind(rownames(matrix_cor_p),matrix_cor_p)
    
    colnames(matrix_cor_out)[1]<-'gene'
    colnames(matrix_cor_p_out)[1]<-'gene'
    write.table(matrix_cor_out,file=paste0(output,'.correlation.r.tsv'),row.names = F,sep='\t',quote = FALSE)
    write.table(matrix_cor_p_out,file=paste0(output,'.correlation.pvalue.tsv'),row.names = F,sep='\t',quote = FALSE)
  }
  
  dim(matrix_cor)
  dim(matrix_cor_p)
  class(matrix_cor)
  class(matrix_cor_p)
  matrix_cor<-as.matrix(matrix_cor)
  matrix_cor_p<-as.matrix(matrix_cor_p)
  print(dim(matrix_cor))
  if(is.na(pw))
  {
    pw<-dim(matrix_cor)[1]/5
    pw<-max(pw,5)
  }
  if(is.na(ph))
  {
    ph<-dim(matrix_cor)[2]/5
    ph<-max(ph,5)
  }
  
  
  if(ploty)
  {
    print(paste0('picture height & width: ',ph,' ',pw))
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,method='circle',ggtheme=theme_bw())
    ggsave(paste0(output,".a.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".a.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw())
    ggsave(paste0(output,".b.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".b.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw(),p.mat = matrix_cor_p,insig="pch",pch.col = "black")
    ggsave(paste0(output,".c.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".c.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
  }
  
}

ggcorrplot_2matrix_plot<-function(matrix_cor,matrix_cor_p,row_name=T,columne_name=T,output='fig.ggcorrplot',pw=NA,ph=NA,save.data=F){
  
  library(ggcorrplot)  
  library(ggthemes)  
  library(psych)  
  if(save.data)
  {
    matrix_cor<-as.data.frame(matrix_cor)
    matrix_cor_p<-as.data.frame(matrix_cor_p)
    matrix_cor_out<-cbind(rownames(matrix_cor),matrix_cor)
    matrix_cor_p_out<-cbind(rownames(matrix_cor_p),matrix_cor_p)
    
    colnames(matrix_cor_out)[1]<-'gene'
    colnames(matrix_cor_p_out)[1]<-'gene'
    write.table(matrix_cor_out,file=paste0(output,'.correlation.r.tsv'),row.names = F,sep='\t',quote = FALSE)
    write.table(matrix_cor_p_out,file=paste0(output,'.correlation.pvalue.tsv'),row.names = F,sep='\t',quote = FALSE)
  }
  
  dim(matrix_cor)
  dim(matrix_cor_p)
  class(matrix_cor)
  class(matrix_cor_p)
  matrix_cor<-as.matrix(matrix_cor)
  matrix_cor_p<-as.matrix(matrix_cor_p)
  print(dim(matrix_cor))
  if(is.na(pw))
  {
    pw<-dim(matrix_cor)[1]/5
    pw<-max(pw,5)
  }
  if(is.na(ph))
  {
    ph<-dim(matrix_cor)[2]/5
    ph<-max(ph,5)
  }
  print(paste0('picture height & width: ',ph,' ',pw))
  if(row_name & columne_name)
  {
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,method='circle',ggtheme=theme_bw())
    ggsave(paste0(output,".a.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".a.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw())
    ggsave(paste0(output,".b.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".b.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw(),p.mat = matrix_cor_p,insig="pch",pch.col = "black")
    ggsave(paste0(output,".c.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".c.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
  }else if(row_name)
  {
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,method='circle',ggtheme=theme_bw())+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),axis.text.x=element_blank())
    ggsave(paste0(output,".a.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".a.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw())+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),axis.text.x=element_blank())
    ggsave(paste0(output,".b.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".b.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw(),p.mat = matrix_cor_p,insig="pch",pch.col = "black")+
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),axis.text.x=element_blank())
    ggsave(paste0(output,".c.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".c.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
  }else if(columne_name)
  {
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,method='circle',ggtheme=theme_bw())+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),axis.text.y=element_blank())
    ggsave(paste0(output,".a.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".a.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw())+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),axis.text.y=element_blank())
    ggsave(paste0(output,".b.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".b.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw(),p.mat = matrix_cor_p,insig="pch",pch.col = "black")+
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),axis.text.y=element_blank())
    ggsave(paste0(output,".c.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".c.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
  }else
  {
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,method='circle',ggtheme=theme_bw())+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),axis.text.x=element_blank(),axis.text.y=element_blank())
    ggsave(paste0(output,".a.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".a.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw())+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),axis.text.x=element_blank(),axis.text.y=element_blank())
    ggsave(paste0(output,".b.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".b.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
    plota<-NULL
    plota<-ggcorrplot(matrix_cor,ggtheme=theme_bw(),p.mat = matrix_cor_p,insig="pch",pch.col = "black")+
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),axis.text.x=element_blank(),axis.text.y=element_blank())
    ggsave(paste0(output,".c.png"), plot = plota, width = pw, height = ph,limitsize=F)
    ggsave(paste0(output,".c.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
  }
  
}


#ggcorrplot_1matrix(df_expr,output='fig2.ggcorrplot',save.data=T)
ggcorrplot_1matrix<-function(mx_a,cor_method='pearson',output='fig.ggcorrplot',hc.order = T,pw=NA,ph=NA,save.data=T,plot_if=T){
             #install.packages('ggcorrplot')
             #install.packages('ggthemes')
             library(ggcorrplot)  
             library(ggthemes)  
             library(psych)  
             
             cor <- corr.test(mx_a,method = cor_method,adjust = "BH",ci = F)  #### pearson spearman两种方法计算相关性
             matrix_cor<-cor$r
             matrix_cor_p<-cor$p
             
             if(save.data)
             {
                 matrix_cor<-as.data.frame(matrix_cor)
                 matrix_cor_p<-as.data.frame(matrix_cor_p)
                 matrix_cor_out<-cbind(rownames(matrix_cor),matrix_cor)
                 matrix_cor_p_out<-cbind(rownames(matrix_cor_p),matrix_cor_p)
                 
                 colnames(matrix_cor_out)[1]<-'gene'
                 colnames(matrix_cor_p_out)[1]<-'gene'
                 write.table(matrix_cor_out,file=paste0(output,'.correlation.r.tsv'),row.names = F,sep='\t',quote = FALSE)
                 write.table(matrix_cor_p_out,file=paste0(output,'.correlation.pvalue.tsv'),row.names = F,sep='\t',quote = FALSE)
             }
             
             dim(matrix_cor)
             dim(matrix_cor_p)
             class(matrix_cor)
             class(matrix_cor_p)
             matrix_cor<-as.matrix(matrix_cor)
             matrix_cor_p<-as.matrix(matrix_cor_p)
             
  if(is.na(pw))
  {
    pw<-dim(matrix_cor)[1]/5
    pw<-max(pw,5)
  }
  if(is.na(ph))
  {
    ph<-dim(matrix_cor)[2]/5
    ph<-max(ph,5)
  }
  
  if(plot_if)
  {
			plota<-NULL
			plota<-ggcorrplot(matrix_cor,hc.order = hc.order,method='circle',ggtheme=theme_bw())
			ggsave(paste0(output,".a.png"), plot = plota, width = pw, height = ph,limitsize=F)
			ggsave(paste0(output,".a.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    
             
			plota<-NULL
			plota<-ggcorrplot(matrix_cor,hc.order = hc.order,ggtheme=theme_bw())          
			ggsave(paste0(output,".b.png"), plot = plota, width = pw, height = ph,limitsize=F)
			ggsave(paste0(output,".b.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
             
             #print(ggcorrplot(matrix_cor,method='circle',ggtheme=theme_bw(),insig = "blank",p.mat = matrix_cor_p))
			plota<-NULL
			plota<-ggcorrplot(matrix_cor,method='circle',hc.order = hc.order,ggtheme=theme_bw(),p.mat = matrix_cor_p,insig="pch",pch.col = "black")
			ggsave(paste0(output,".c.png"), plot = plota, width = pw, height = ph,limitsize=F)
			ggsave(paste0(output,".c.pdf"), plot = plota, width = pw, height = ph,limitsize=F) 
     }
}


ggheatmap_yingbio_for_correlation<-function(mydata,orderby=NA,geneid='miRNA_ID',row_name=T,columne_name=T,ph=NA,pw=6,output='heatmap',save.data=F){
  ngene<-dim(mydata)[1]
  nsample<-dim(mydata)[2]
  if(!is.na(ph))
  {
    phu<-ph
  }else{
    phu<-max(6,ngene/10)
  }
  if(!is.na(pw))
  {
    pwu<-pw
  }else{
    pwu<-max(6,nsample/10)
  }
  print(paste0('plot height & width: ',phu,' ',pwu))
  library(ggplot2)
  library(pheatmap)
  library(gtable)
  plot.pheatmap<-NULL
  plot.pheatmap<-pheatmap::pheatmap(mydata,border = F,
                                    color = colorRampPalette(c(rep("dodgerblue",2), "white", rep("brown2",2)))(100),
                                    cluster_cols =F,
                                    show_colnames = columne_name,
                                    show_rownames = row_name,
                                    scale='none',
                                    fontsize = 10,
                                    fontsize_row=5,
                                    fontsize_col=10,angle_col = "45")
  ####  如果右侧的legend空间不够，增加legend空间。
  #plot.pheatmap$gtable$widths[1] <- plot.pheatmap$gtable$widths[1] + unit(200, "bigpts")
  #plot.pheatmap$gtable$widths[6] <- plot.pheatmap$gtable$widths[6] + unit(5, "bigpts")
  dev.off()
  pdf(file=paste0(output,'.pdf'), height=phu, width=pwu)
  print(plot.pheatmap)
  dev.off()
  tiff(file=paste0(output,'.tiff'), height=phu*300, width=pwu*300,units="px",res=300,compression='lzw')
  print(plot.pheatmap)
  dev.off()
}

ggheatmap_yingbio_1col<-function(mydata,orderby=NA,geneid='miRNA_ID',pheight=NA,pw=6,output='heatmap',save.data=F){
  ngene<-dim(mydata)[1]
  if(!is.na(pheight))
  {
    ph<-pheight
  }else{
    ph<-max(6,ngene/10)
  }
  print(paste0('pic height',ph))
  library(ggplot2)
  library(pheatmap)
  library(gtable)
  plot.pheatmap<-NULL
  plot.pheatmap<-pheatmap::pheatmap(mydata,border = T,
                                    color = colorRampPalette(c(rep("dodgerblue",2), "white", rep("brown2",2)))(100),
                                    cluster_rows =F,
                                    cluster_cols =F,
                                    show_colnames = F,
                                    show_rownames = T,
                                    fontsize = 10,
                                    fontsize_row=5,
                                    fontsize_col=10,angle_col = "45")
  ####  如果右侧的legend空间不够，增加legend空间。
  #plot.pheatmap$gtable$widths[1] <- plot.pheatmap$gtable$widths[1] + unit(200, "bigpts")
  #plot.pheatmap$gtable$widths[6] <- plot.pheatmap$gtable$widths[6] + unit(5, "bigpts")
  dev.off()
  pdf(file=paste0(output,'.pdf'), height=ph, width=pw)
  print(plot.pheatmap)
  dev.off()
  tiff(file=paste0(output,'.tiff'), height=ph*300, width=pw*300,units="px",res=300,compression='lzw')
  print(plot.pheatmap)
  dev.off()
}


filter_cor<-function(matrix_cor,cutoff=0.3)
{
    matrix_cor<-as.data.frame(matrix_cor)
    row_max<-apply(matrix_cor,1,max)
    row_min<-apply(matrix_cor,1,min)
    row_used1<-(row_max > cutoff)
    row_used2<-(row_min < -cutoff)
    row_used<-rownames(matrix_cor)[row_used1 | row_used2]
    col_max<-apply(matrix_cor,2,max)
    col_min<-apply(matrix_cor,2,min)
    col_used1<-(col_max > cutoff)
    col_used2<-(col_min < -cutoff)
    col_used<- colnames(matrix_cor)[col_used1 | col_used2]
    matrix_cor_used<-matrix_cor[row_used,col_used,drop=F]
    return(matrix_cor_used)
}
  
filter_cor_pvalue<-function(matrix_cor_p,cutoff=0.05)
{
    matrix_cor_p<-as.data.frame(matrix_cor_p)
    row_min<-apply(matrix_cor_p,1,min)
    row_used2<-(row_min < cutoff)
    row_used<-rownames(matrix_cor_p)[row_used2]
    

    col_min<-apply(matrix_cor_p,2,min)
    col_used2<-(col_min < cutoff)
    col_used<- colnames(matrix_cor_p)[col_used2]
    matrix_cor_p_used<-matrix_cor_p[row_used,col_used,drop=F]
    return(matrix_cor_p_used)
}



#aaa<-DEG_out[,c('gene','logFC')]
#gsea_yingbai(aaa,dir_output=paste0('gsea_',gene),species = "human")
gsea_yingbai<-function(aaa,dir_output='gsea_results',species = "human",nplot=20,pathways=NULL){
            
            if(species=='human')
            {
              DBkeyset='org.Hs.eg.db'
              organism='hsa'
              library(org.Hs.eg.db)
              
            }else if(species=='mouse')
            {
              library(org.Mm.eg.db)
              DBkeyset='org.Mm.eg.db'
              organism='mmu'
            }else if(species=='rat')
            {
              library(org.Rn.eg.db)
              #BiocManager::install('org.Rn.eg.db')
              DBkeyset='org.Rn.eg.db'
              organism='rno'
            }else{
              DBkeyset='org.Hs.eg.db'
              library(org.Hs.eg.db)
              organism='hsa'
            }
            
            library(clusterProfiler)
            library(enrichplot)
            library(dplyr)
            library(stats)
            library(ggplot2)
            #BiocManager::install("clusterProfiler")
            if(!dir.exists(dir_output))
            {
              dir.create(dir_output)
            }
            colnames(aaa)<-c('gene','logFC')
            genename<-aaa$gene
            gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=DBkeyset)
            gene_map <- dplyr::distinct(gene_map,SYMBOL,.keep_all=TRUE)
            aaa<-merge(aaa,gene_map,by.x='gene',by.y='SYMBOL')
            aaa<-aaa[order(aaa[,2],decreasing=T),]
            class(aaa)
            head(aaa)
            
            bbb<-aaa[,2]
            names(bbb)<-aaa$ENTREZID
            class(bbb)
            head(bbb)
            KEGG_gseresult <- gseKEGG(bbb, organism = organism,minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
            dim(KEGG_gseresult)
            ##335
            head(KEGG_gseresult)
            tail(KEGG_gseresult)
            nkegg<-nrow(KEGG_gseresult)
            
            kksorted<-KEGG_gseresult[order(KEGG_gseresult$NES, decreasing = T),]
            head(kksorted)
            tail(kksorted)
            write.table(kksorted,paste0(dir_output,'/gsea_results.tsv'),sep = "\t",quote = F,row.names = F)
            
            
            for(i in unique(c(1:nplot,(nkegg-nplot+1):nkegg)))
            {
              #i<-1
              pathway<-kksorted[i,1]
              plota<-NULL
              plota<-gseaplot2(KEGG_gseresult, pathway,title = kksorted[i,2])
              pdf(file=paste0(dir_output,'/',pathway,".pdf"), height=6, width=8)
              print(plota)
              dev.off()
              tiff(file=paste0(dir_output,'/',pathway,".tiff"), height=6*300, width=8*300,res=300,compression='lzw')
              print(plota)
              dev.off()
              #ggsave(paste0(dir_output,'/',pathway,".tiff"),width = 8, height = 6,compression='lzw')
              #ggsave(paste0(dir_output,'/',pathway,".pdf"),width = 8, height = 6)
            }
            
            if(!is.null(pathways))
            {
               kkused<-kksorted[kksorted$ID %in% pathways,]
              for(i in 1:nrow(kkused))
              {
                pathway<-kkused[i,1]
                plota<-NULL
                plota<-gseaplot2(KEGG_gseresult, pathway,title = kkused[i,2])
                pdf(file=paste0(dir_output,'/',pathway,".pdf"), height=6, width=8)
                print(plota)
                dev.off()
                tiff(file=paste0(dir_output,'/',pathway,".tiff"), height=6*300, width=8*300,res=300,compression='lzw')
                print(plota)
                dev.off()
              }
            }
}


#aaa<-DEG_out[,c('gene','logFC')]
#gsea_msigdbr_yingbai(aaa,dir_output=paste0('gsea_',gene),species = "human",category="C2",gs_subcat="CGP")
gsea_msigdbr_yingbai<-function(aaa,dir_output='gsea_results',species = "human",pathways=NULL,category="C2",gs_subcat="CGP"){

		if(species=='human')
		{
			   DBkeyset='org.Hs.eg.db'
			   organism='hsa'
			   org='Homo sapiens'
			   library(org.Hs.eg.db)
			   
		}else if(species=='mouse')
		{
				library(org.Mm.eg.db)
				 DBkeyset='org.Mm.eg.db'
				 organism='mmu'
				 org='Mus musculus'
		}else if(species=='rat')
		{
				library(org.Rn.eg.db)
				#BiocManager::install('org.Rn.eg.db')
				DBkeyset='org.Rn.eg.db'
				organism='rno'
				org='Rattus norvegicus'
		}else{
				DBkeyset='org.Hs.eg.db'
			    library(org.Hs.eg.db)
			    organism='hsa'
			    org='Homo sapiens'
		}

		library(clusterProfiler)
		library(enrichplot)
		library(dplyr)
		library(stats)
		library(msigdbr)
		library(ggplot2)
		#BiocManager::install("clusterProfiler")
		if(!dir.exists(dir_output))
		{
				dir.create(dir_output)
		}
		colnames(aaa)<-c('gene','logFC')
		genename<-aaa$gene
		aaa<-aaa[order(aaa[,2],decreasing=T),]
		class(aaa)
		head(aaa)

		bbb<-aaa[,2]
		names(bbb)<-aaa$gene
		class(bbb)
		head(bbb)
		
		genesets <- msigdbr(species = org, category = category)
		genesets <- subset(genesets, gs_subcat==gs_subcat, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
		KEGG_gseresult <- clusterProfiler::GSEA(bbb[bbb!=0], TERM2GENE = genesets, minGSSize = 1,
                              pvalueCutoff = 1, nPermSimple = 10000, verbose = F)
		dim(KEGG_gseresult)
		##335
		head(KEGG_gseresult)
		tail(KEGG_gseresult)
		nkegg<-nrow(KEGG_gseresult)

		kksorted<-KEGG_gseresult[order(KEGG_gseresult$NES, decreasing = T),]
		head(kksorted)
		tail(kksorted)
		write.table(kksorted,paste0(dir_output,'/gsea_results.tsv'),sep = "\t",quote = F,row.names = F)


		for(i in c(1:5,(nkegg-5):nkegg))
		{
		   #i<-1
		   pathway<-kksorted[i,1]
		   plota<-NULL
		    plota<-gseaplot2(KEGG_gseresult, pathway,pvalue_table = F,title = pathway)
		    pdf(file=paste0(dir_output,'/',pathway,".pdf"), height=6, width=8)
		    print(plota)
		    dev.off()
		    tiff(file=paste0(dir_output,'/',pathway,".tiff"), height=6*300, width=8*300,res=300,compression='lzw')
		    print(plota)
		    dev.off()
		   #ggsave(paste0(dir_output,'/',pathway,".tiff"),width = 8, height = 6,compression='lzw')
		   #ggsave(paste0(dir_output,'/',pathway,".pdf"),width = 8, height = 6)
		}
		
		if(!is.null(pathways))
		{
			for(i in 1:length(pathways))
			{
				pathway<-pathways[i]
				  plota<-NULL
				  plota<-gseaplot2(KEGG_gseresult, pathway,pvalue_table = F,title = pathway)
				  pdf(file=paste0(dir_output,'/',pathway,".pdf"), height=6, width=8)
				  print(plota)
				  dev.off()
				  tiff(file=paste0(dir_output,'/',pathway,".tiff"), height=6*300, width=8*300,res=300,compression='lzw')
				  print(plota)
				  dev.off()
			}
		}
}


#aaa<-sig.pbmc.markers.by.group$gene
#aaa<-DE_pcg_down_neg
#gokegg_yingbai(aaa,dir_output=c(paste0('go_',group_test,'_vs_',group_control),paste0('kegg_',group_test,'_vs_',group_control)),species = "rat")
#gokegg_yingbai(DE_pcg_down_neg,dir_output=c(paste0('GO_','DE_pcg_down.coexpression_neg'),paste0('kegg_','DE_pcg_down.coexpression_neg')),species = "human",go_pvalueCutoff=0.05,kegg_pvalueCutoff=0.5)
gokegg_yingbai<-function(aaa,dir_output=c('GO_results','KEGG_results'),species = "human",go_pvalueCutoff=1,kegg_pvalueCutoff=1)
{
  go_yingbai(aaa,dir_output=dir_output[1],species = species,go_pvalueCutoff=go_pvalueCutoff)
  kegg_yingbai_new(aaa,dir_output=dir_output[2],species = species,kegg_pvalueCutoff=kegg_pvalueCutoff)
}

go_yingbai<-function(aaa,dir_output='GO_results',species = "human",go_pvalueCutoff=0.1){
  library(ggplot2)
  #BiocManager::install("topGO",force=T)
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("XVector")
  #BiocManager::install('GO.db')
  #BiocManager::install('DBI')
  #BiocManager::install('AnnotationDbi',force=T)
  #options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
  #BiocManager::install('Rgraphviz')
  library(topGO)
  library(clusterProfiler)
  #library(KEGG.db)
  #library(org.Hs.eg.db)
  library(enrichplot)
  #library(pathview)
  #BiocManager::install("pathview")
  
  if(species=='human')
  {
    #BiocManager::install('org.Hs.eg.db')
    DBkeyset='org.Hs.eg.db'
    org='hsa'
    library(org.Hs.eg.db)
    
  }else if(species=='mouse')
  {
    #BiocManager::install('org.Mm.eg.db')
    library(org.Mm.eg.db)
    DBkeyset='org.Mm.eg.db'
    org='mmu'
  }else if(species=='rat')
  {
    library(org.Rn.eg.db)
    #BiocManager::install('org.Rn.eg.db')
    DBkeyset='org.Rn.eg.db'
    org='rno'
  }else{
    DBkeyset='org.Hs.eg.db'
    library(org.Hs.eg.db)
    org='hsa'
  }
  
  go_dir<-dir_output
  if(dir.exists(go_dir))
    {
      unlink(x = go_dir, recursive = TRUE)
    }

  if(!dir.exists(go_dir))
  {
    dir.create(go_dir)
  }

  bbb<-NULL
  bbb = try(bitr(aaa, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=DBkeyset),silent = T)
  if(class(bbb) =='try-error')
	{
    go_bp<-data.frame()
    go_cc<-data.frame()
    go_mf<-data.frame()
	}else if(nrow(bbb)>0)
  {
    go_bp <- enrichGO(gene = bbb$ENTREZID,OrgDb = DBkeyset,ont = "BP",pAdjustMethod = "BH", pvalueCutoff = go_pvalueCutoff,qvalueCutoff = 1,minGSSize=1,readable = TRUE)
    go_cc <- enrichGO(gene = bbb$ENTREZID,OrgDb = DBkeyset,ont = "CC",pAdjustMethod = "BH", pvalueCutoff = go_pvalueCutoff,qvalueCutoff = 1,minGSSize=1,readable = TRUE)
    go_mf <- enrichGO(gene = bbb$ENTREZID,OrgDb = DBkeyset,ont = "MF",pAdjustMethod = "BH", pvalueCutoff = go_pvalueCutoff,qvalueCutoff = 1,minGSSize=1,readable = TRUE)
  }else{
    go_bp<-data.frame()
    go_cc<-data.frame()
    go_mf<-data.frame()
  }
  go_output_yingbai(go_bp,dir_output=go_dir,term = "BP")
  go_output_yingbai(go_cc,dir_output=go_dir,term = "CC")
  go_output_yingbai(go_mf,dir_output=go_dir,term = "MF")
  
  go_bp_combi<-as.data.frame(go_bp)
  #go_bp_combi<-go_bp_combi[go_bp_combi$ID=='xxx',]
  n_bp<-nrow(go_bp_combi)
  if(n_bp)
  {
    go_bp_combi$Ontology<-'Biological process'
  }
  go_cc_combi<-as.data.frame(go_cc)
  #go_cc_combi<-go_cc_combi[go_cc_combi$ID=='xxx',]
  n_cc<-nrow(go_cc_combi)
  if(n_cc)
  {
    go_cc_combi$Ontology<-'Cellular component'
  }
  go_mf_combi<-as.data.frame(go_mf)
  #go_mf_combi<-go_mf_combi[go_mf_combi$ID=='xxx',]
  n_mf<-nrow(go_mf_combi)
  if(n_mf)
  {
    go_mf_combi$Ontology<-'Molecular function'
  }
  if(n_bp | n_cc |n_mf)
  {
    go_all<-rbind(go_bp_combi,go_cc_combi,go_mf_combi)
    write.table(go_all,paste0(go_dir,'/GO_','all',".tsv"),sep='\t',row.names =FALSE,quote=F)
  }
}

kegg_yingbai<-function(aaa,dir_output='KEGG_results',species = "human",term_number=20,kegg_pvalueCutoff=1,pw=7,ph=7.5,name_width=50){
  library(ggplot2)
  #BiocManager::install("topGO")
  #BiocManager::install("clusterProfiler")
  library(topGO)
  library(clusterProfiler)
  #library(KEGG.db)
  #library(org.Hs.eg.db)
  library(enrichplot)
  #library(pathview)
  #BiocManager::install("pathview")
  
  if(species=='human')
  {
    DBkeyset='org.Hs.eg.db'
    org='hsa'
    organism='Homo sapiens'
    library(org.Hs.eg.db)
    
  }else if(species=='mouse')
  {
    library(org.Mm.eg.db)
    DBkeyset='org.Mm.eg.db'
    org='mmu'
    organism='Mus musculus'
  }else if(species=='rat')
  {
    library(org.Rn.eg.db)
    #BiocManager::install('org.Rn.eg.db')
    DBkeyset='org.Rn.eg.db'
    org='rno'
    organism='Rattus norvegicus'
  }else{
    DBkeyset='org.Hs.eg.db'
    library(org.Hs.eg.db)
    org='hsa'
    organism=''
  }    
  
  kegg_dir<-dir_output
  if(!dir.exists(kegg_dir))
  {
    dir.create(kegg_dir)
  }
  bbb = try(bitr(aaa, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=DBkeyset),silent = T)
  if(class(bbb) =='try-error')
	{
	  print('bad')
	  nkk<-0
	}else if(nrow(bbb)>0)
  {
    kk <- enrichKEGG(gene = bbb$ENTREZID,organism = org,pvalueCutoff = kegg_pvalueCutoff,qvalueCutoff=1,minGSSize=1,use_internal_data =F)
    nkk<-nrow(kk)
  }else{
    nkk<-0
  }
  if(is.null(nkk))nkk<-0
  if(nkk>0)
  {
    df_kk<-as.data.frame(kk)
    df_kk$gene_name<-NA
    for(i in 1:nrow(df_kk))
    {
      df_kk[i,'gene_name']<-kegg_geneid_convert(df_kk[i,'geneID'],species = species)
    }
    write.table(df_kk,paste0(kegg_dir,'/',"KEGG-enrich.tsv"),sep='\t',row.names =FALSE,quote=F)
    
    if(1)
    {
        plota<-NULL
        plota<-cnetplot(kk, showCategory = 5,circular = T,colorEdge = T,color_gene='#0166cc',color_category='#FF7F00')
        ggsave(paste0(kegg_dir,'/',"KEGG",".cnetplot.tiff"), plot = plota, width = 8, height = 6,compression='lzw',limitsize=F)
        ggsave(paste0(kegg_dir,'/',"KEGG",".cnetplot.pdf"), plot = plota, width = 8, height = 6,limitsize=F)
        
        ego_pair <- enrichplot::pairwise_termsim(kk)
        plota<-NULL
        plota<-emapplot(ego_pair,  layout="kk", cex_category=1.5,min_edge = 0.8) 
        ggsave(paste0(kegg_dir,'/',"KEGG",".emapplot.tiff"), plot = plota, width = 16, height = 12,compression='lzw',limitsize=F)
        ggsave(paste0(kegg_dir,'/',"KEGG",".emapplot.pdf"), plot = plota, width = 16, height = 12,limitsize=F)
    }
    
    if(0)  ####  clusterProfiler绘制barplot和dotplot
    {
      plota<-NULL
      plota<-dotplot(kk,title="EnrichmentKEGG",showCategory=20,)+theme(plot.title = element_text(hjust=0.5, face="bold"))
      ggsave(paste0(kegg_dir,'/',"KEGG",".dotplot.tiff"), plot = plota, width = 12, height = 12,compression='lzw',limitsize=F)
      ggsave(paste0(kegg_dir,'/',"KEGG",".dotplot.pdf"), plot = plota, width = 12, height = 12,limitsize=F)
      
      
      plota<-NULL
      plota<-barplot(kk, showCategory=20,title="EnrichmentKEGG")+theme(plot.title = element_text(hjust=0.5, face="bold"))
      ggsave(paste0(kegg_dir,'/',"KEGG",".barplot.tiff"), plot = plota, width = 12, height = 12,compression='lzw',limitsize=F)
      ggsave(paste0(kegg_dir,'/',"KEGG",".barplot.pdf"), plot = plota, width = 12, height = 12,limitsize=F)
    }
    
    if(1)   ####  使用ggplot2绘制barplot和dotplot
    {
      library(stringr)
      #library(tidyverse)
      library(ggplot2)
      df_plot<-head(df_kk,term_number)
      head(df_plot)
      df_plot$matched_all<-0
      df_plot$pathway_ngene<-0
      df_plot$kegg_ngene<-0
      for(i in 1:nrow(df_plot))
      {
        #i=1
        df_plot[i,'matched_all']<-as.numeric(strsplit(df_plot[i,'GeneRatio'],'/')[[1]][2])
        df_plot[i,'pathway_ngene']<-as.numeric(strsplit(df_plot[i,'BgRatio'],'/')[[1]][1])
        df_plot[i,'kegg_ngene']<-as.numeric(strsplit(df_plot[i,'BgRatio'],'/')[[1]][2])
      }
      head(df_plot,20)
      df_plot<-df_plot[order(df_plot$pvalue,-df_plot$Count),]
      colnames(df_plot)
      df_plot$GeneRatio<-as.numeric(df_plot$Count/df_plot$matched_all)
      df_plot$BgRatio<-as.numeric(df_plot$pathway_ngene/df_plot$kegg_ngene)
      df_plot$Enrichment<-as.numeric(df_plot$GeneRatio/df_plot$BgRatio)
      df_plot<-df_plot[order(-df_plot$Enrichment),]
      
      
      rownames(df_plot)<-NULL
      df_plot$index<-as.numeric(rownames(df_plot))
      nlab<-max(nchar(df_plot$Description))
      
      df_plot$Description<- gsub(paste0(' - ',organism,'.*'),'',df_plot$Description)
      
      plota<-NULL
      plota <- ggplot(df_plot) + 
        geom_bar(aes(x=reorder(Description,-index),y=Enrichment,fill=-log10(pvalue)),stat='identity',width=0.8)
      plota<-plota+theme(axis.text.x = element_text(size = 12,face='bold'))
      plota<-plota+theme(axis.text.y = element_text(size = 12,face='bold')) 
      plota<-plota+theme(axis.title.y=element_blank(),legend.position='right',axis.title.x=element_text(size = 15,face='bold'))                
      plota<-plota+labs(y='Enrichment',title='Enrichment of KEGG Pathway')+scale_fill_gradient(low="blue",high="red")
      plota<-plota+theme(axis.line = element_line(linewidth=1, colour = "black"))                  
      plota<-plota+theme(panel.border = element_blank(),plot.title = element_text(size = 15,hjust=0.5, face="bold"))                                     
      plota<-plota+theme(panel.grid =element_blank(),panel.background = element_blank())    
      plota<-plota+coord_flip()
      plota<-plota+scale_y_continuous(expand=c(0,0))
      plota<-plota+scale_x_discrete(labels=function(x) str_wrap(x, width=name_width))

      ggsave(paste0(kegg_dir,'/',"KEGG",".barlot.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
      ggsave(paste0(kegg_dir,'/',"KEGG",".barlot.pdf"), plot = plota, width = pw, height = ph)
      
      head(df_plot,20)
      
      df_plot<-df_plot[order(-df_plot$pvalue,df_plot$Count),]
      df_plot<-df_plot[order(df_plot$GeneRatio),]
      rownames(df_plot)<-NULL
      df_plot$index<-as.numeric(rownames(df_plot))
      
      plota<-NULL
      plota <- ggplot(df_plot) + 
        geom_point(aes(x=GeneRatio,y=index,size=Count,col=-log10(pvalue)),alpha = 0.99, position = position_jitter(w = 0.0, h = 0.0))
      plota<-plota+theme(axis.text.x = element_text(size = 12,face='bold'))
      plota<-plota+theme(axis.text.y = element_text(size = 10,face='bold')) 
      plota<-plota+theme(axis.title.y=element_blank(),axis.title.x=element_text(size = 15,face='bold'),legend.position='right')                
      plota<-plota+labs(x='GeneRatio',title='Enrichment of KEGG Pathway')+scale_colour_distiller(palette = "RdYlBu")
      #+scale_color_gradient(low="blue",high="red")
      plota<-plota+theme(axis.line = element_line(linewidth=1, colour = "black"))                  
      plota<-plota+theme(panel.border = element_blank(),plot.title = element_text(size = 15,hjust=0.5, face="bold"))                                     
      plota<-plota+theme(panel.grid =element_blank(),panel.background = element_blank())    
      plota<-plota+scale_y_continuous(breaks=df_plot$index,labels= str_wrap(df_plot$Description, width=name_width))

      #scale_y_discrete(labels=function(x) str_wrap(x, width=ceiling(nlab/2)))
      ggsave(paste0(kegg_dir,'/',"KEGG",".dotplot.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
      ggsave(paste0(kegg_dir,'/',"KEGG",".dotplot.pdf"), plot = plota, width = pw, height = ph)
    }
    
    #pathways<-kk$ID
    #ccc<-as.numeric(bbb$ENTREZID)
    #for(pathway in pathways)
    #{            
    #pathway='hsa04660'
    #kegg_dir='kegg_results'
    # pv.out <- pathview(gene.data = ccc,pathway.id = pathway,kegg.dir=kegg_dir,species = organism,limit = list(gene=max(abs(ccc)), cpd=1))
    #}
  }else{
    write.table("",paste0(kegg_dir,'/',"NO enriched pathway.tsv"),sep='\t',row.names =FALSE,quote=F)
  }
}


kegg_yingbai_new<-function(aaa,dir_output='KEGG_results',species = "human",term_number=20,kegg_pvalueCutoff=0.05,pw=7,ph=7.5,name_width=50)
{
		  library(ggplot2)
		  #BiocManager::install("topGO")
		  #BiocManager::install("clusterProfiler")
		  library(topGO)
		  library(clusterProfiler)
		  #library(KEGG.db)
		  #library(org.Hs.eg.db)
		  library(enrichplot)
		  #library(pathview)
		  #BiocManager::install("pathview")
		  
		  if(species=='human')
		  {
		    DBkeyset='org.Hs.eg.db'
		    org='hsa'
		    organism='Homo sapiens'
		    library(org.Hs.eg.db)
		    
		  }else if(species=='mouse')
		  {
		    library(org.Mm.eg.db)
		    DBkeyset='org.Mm.eg.db'
		    org='mmu'
		    organism='Mus musculus'
		  }else if(species=='rat')
		  {
		    library(org.Rn.eg.db)
		    #BiocManager::install('org.Rn.eg.db')
		    DBkeyset='org.Rn.eg.db'
		    org='rno'
		    organism='Rattus norvegicus'
		  }else{
		    DBkeyset='org.Hs.eg.db'
		    library(org.Hs.eg.db)
		    org='hsa'
		    organism=''
		  }    
		  
		  kegg_dir<-dir_output
		if(dir.exists(kegg_dir))
			{
			  unlink(x = kegg_dir, recursive = TRUE)
			}
    
		  if(!dir.exists(kegg_dir))
		  {
		    dir.create(kegg_dir)
		  }
		  bbb = try(bitr(aaa, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=DBkeyset),silent = T)
		  if(class(bbb) =='try-error')
		  {
		    print('bad')
		    nres<-0
		  }else if(nrow(bbb)>0)
		  {
		    enrich_result <- enrichKEGG(gene = bbb$ENTREZID,organism = org,pvalueCutoff = kegg_pvalueCutoff,qvalueCutoff=1,minGSSize=1,use_internal_data =F)
		    nres<-nrow(enrich_result)
		  }else{
		    nres<-0
		  }
		  if(is.null(nres))nres<-0
		  if(nres>0)
		  {
		    df_result<-as.data.frame(enrich_result)
		    df_result$gene_name<-NA
		    for(i in 1:nrow(df_result))
		    {
		      df_result[i,'gene_name']<-kegg_geneid_convert(df_result[i,'geneID'],species = species)
		    }
		    write.table(df_result,paste0(kegg_dir,'/',"KEGG-enrich.tsv"),sep='\t',row.names =FALSE,quote=F)
		    
		    if(1)
		    {
		      plota<-NULL
		      plota<-cnetplot(enrich_result, showCategory = 5,circular = T,colorEdge = T,color_gene='#0166cc',color_category='#FF7F00')
		      ggsave(paste0(kegg_dir,'/',"KEGG",".cnetplot.tiff"), plot = plota, width = 8, height = 6,compression='lzw',limitsize=F)
		      ggsave(paste0(kegg_dir,'/',"KEGG",".cnetplot.pdf"), plot = plota, width = 8, height = 6,limitsize=F)
		      
		      ego_pair <- enrichplot::pairwise_termsim(enrich_result)
		      plota<-NULL
		      plota<-emapplot(ego_pair,  layout="kk", cex_category=1.5,min_edge = 0.8) 
		      ggsave(paste0(kegg_dir,'/',"KEGG",".emapplot.tiff"), plot = plota, width = 16, height = 12,compression='lzw',limitsize=F)
		      ggsave(paste0(kegg_dir,'/',"KEGG",".emapplot.pdf"), plot = plota, width = 16, height = 12,limitsize=F)
		    }
		    
		    if(0)  ####  clusterProfiler绘制barplot和dotplot
		    {
		      plota<-NULL
		      plota<-dotplot(enrich_result,title="EnrichmentKEGG",showCategory=20,)+theme(plot.title = element_text(hjust=0.5, face="bold"))
		      ggsave(paste0(kegg_dir,'/',"KEGG",".dotplot.tiff"), plot = plota, width = 12, height = 12,compression='lzw',limitsize=F)
		      ggsave(paste0(kegg_dir,'/',"KEGG",".dotplot.pdf"), plot = plota, width = 12, height = 12,limitsize=F)
		      
		      
		      plota<-NULL
		      plota<-barplot(enrich_result, showCategory=20,title="EnrichmentKEGG")+theme(plot.title = element_text(hjust=0.5, face="bold"))
		      ggsave(paste0(kegg_dir,'/',"KEGG",".barplot.tiff"), plot = plota, width = 12, height = 12,compression='lzw',limitsize=F)
		      ggsave(paste0(kegg_dir,'/',"KEGG",".barplot.pdf"), plot = plota, width = 12, height = 12,limitsize=F)
		    }
		    
		    if(1)   ####  使用ggplot2绘制barplot和dotplot
		    {
		      library(stringr)
		      library(tidyverse)
		      library(ggplot2)
		      
		      df_plot<-head(df_result,term_number)
		      head(df_plot)
		      df_plot$matched_all<-0
		      df_plot$pathway_ngene<-0
		      df_plot$kegg_ngene<-0
		      for(i in 1:nrow(df_plot))
		      {
		        #i=1
		        df_plot[i,'matched_all']<-as.numeric(strsplit(df_plot[i,'GeneRatio'],'/')[[1]][2])
		        df_plot[i,'pathway_ngene']<-as.numeric(strsplit(df_plot[i,'BgRatio'],'/')[[1]][1])
		        df_plot[i,'kegg_ngene']<-as.numeric(strsplit(df_plot[i,'BgRatio'],'/')[[1]][2])
		      }
		      head(df_plot,20)
		      df_plot<-df_plot[order(df_plot$pvalue,-df_plot$Count),]
		      colnames(df_plot)
		      df_plot$GeneRatio<-as.numeric(df_plot$Count/df_plot$matched_all)
		      df_plot$BgRatio<-as.numeric(df_plot$pathway_ngene/df_plot$kegg_ngene)
		      df_plot$Enrichment<-as.numeric(df_plot$GeneRatio/df_plot$BgRatio)
		      df_plot<-df_plot[order(-df_plot$Enrichment),]
		      df_plot$Description<- gsub(paste0(' - ',organism,'.*'),'',df_plot$Description)
		      
		      ggplo2_kegg_barplot(df_plot,name_column='Description',value_column='Enrichment',color_column='pvalue',
		                          ylab='Enrichment',title='Enrichment of KEGG Pathway',pw=7,ph=7.5,name_width=50,output_name=paste0(kegg_dir,'/KEGG.barplot'))
		      
		      ggplo2_kegg_dotplot(df_plot,name_column='Description',value_column='GeneRatio',color_column='pvalue',size_column='Count',
		                          xlab='GeneRatio',title='Enrichment of KEGG Pathway',pw=7,ph=7.5,name_width=50,output_name=paste0(kegg_dir,'/KEGG.dotplot'))
		      
		    }
		    
		    #pathways<-kk$ID
		    #ccc<-as.numeric(bbb$ENTREZID)
		    #for(pathway in pathways)
		    #{            
		    #pathway='hsa04660'
		    #kegg_dir='kegg_results'
		    # pv.out <- pathview(gene.data = ccc,pathway.id = pathway,kegg.dir=kegg_dir,species = organism,limit = list(gene=max(abs(ccc)), cpd=1))
		    #}
		  }else{
		    write.table("",paste0(kegg_dir,'/',"NO enriched pathway.tsv"),sep='\t',row.names =FALSE,quote=F)
		  }
		}

		
ggplo2_kegg_barplot<-function(df_plot,name_column='Description',value_column='Enrichment',color_column='pvalue',
							  ylab='Enrichment',title='Enrichment of KEGG Pathway',pw=7,ph=7.5,name_width=50,output_name='KEGG.barplot')
{
  library(stringr)
  library(tidyverse)
  library(ggplot2)
  
  rownames(df_plot)<-NULL
  df_plot$index<-as.numeric(rownames(df_plot))
  nlab<-max(nchar(df_plot[[name_column]]))
  
  plota<-NULL
  plota <- ggplot(df_plot) + 
	geom_bar(aes(x=reorder(.data[[name_column]],-index),y=.data[[value_column]],fill=-log10(.data[[color_column]])),stat='identity',width=0.8)
  plota<-plota+theme(axis.text.x = element_text(size = 12,face='bold'))
  plota<-plota+theme(axis.text.y = element_text(size = 12,face='bold')) 
  plota<-plota+theme(axis.title.y=element_blank(),legend.position='right',axis.title.x=element_text(size = 15,face='bold'))
  plota<-plota+labs(y=ylab,title=title)+scale_fill_gradient(low="blue",high="red")
  plota<-plota+theme(axis.line = element_line(linewidth=1, colour = "black"))
  plota<-plota+theme(panel.border = element_blank(),plot.title = element_text(size = 15,hjust=0.5, face="bold")) 
  plota<-plota+theme(panel.grid =element_blank(),panel.background = element_blank())
  plota<-plota+coord_flip()
  plota<-plota+scale_y_continuous(expand=c(0,0))
  plota<-plota+scale_x_discrete(labels=function(x) str_wrap(df_plot[[name_column]], width=name_width))
  
  ggsave(paste0(output_name,".tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
  ggsave(paste0(output_name,".pdf"), plot = plota, width = pw, height = ph)
}

ggplo2_kegg_dotplot<-function(df_plot,name_column='Description',value_column='GeneRatio',color_column='pvalue',size_column='Count',
		                              xlab='GeneRatio',title='Enrichment of KEGG Pathway',pw=7,ph=7.5,name_width=50,output_name='KEGG.dotplot')
{
	  df_plot<-df_plot[order(df_plot[[value_column]]),]
	  rownames(df_plot)<-NULL
	  df_plot$index<-as.numeric(rownames(df_plot))
	  
	  plota<-NULL
	  plota <- ggplot(df_plot) + 
		geom_point(aes(x=.data[[value_column]],y=index,size=.data[[size_column]],col=-log10(.data[[color_column]])),alpha = 0.99, position = position_jitter(w = 0.0, h = 0.0))
	  plota<-plota+theme(axis.text.x = element_text(size = 12,face='bold'))
	  plota<-plota+theme(axis.text.y = element_text(size = 10,face='bold')) 
	  plota<-plota+theme(axis.title.y=element_blank(),axis.title.x=element_text(size = 15,face='bold'),legend.position='right')
	  plota<-plota+labs(x=xlab,title=title)+scale_colour_distiller(palette = "RdYlBu")
	  #+scale_color_gradient(low="blue",high="red")
	  plota<-plota+theme(axis.line = element_line(linewidth=1, colour = "black"))
	  plota<-plota+theme(panel.border = element_blank(),plot.title = element_text(size = 15,hjust=0.5, face="bold")) 
	  plota<-plota+theme(panel.grid =element_blank(),panel.background = element_blank())
	  plota<-plota+scale_y_continuous(breaks=df_plot$index,labels= str_wrap(df_plot[[name_column]], width=name_width))
	  ggsave(paste0(output_name,".tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
	  ggsave(paste0(output_name,".pdf"), plot = plota, width = pw, height = ph)
}

reactome_yingbai<-function(aaa,Reactome_dir='Reactome_result',species="human",term_number=20,pvalueCutoff = 0.05,pw=7,ph=7.5,name_width=50)
{
		#Reactome_dir<-'Reactome_GO_vs_Normal_UPDOWN'
		#BiocManager::install("ReactomePA")
		library(ReactomePA)
		library(clusterProfiler)
		library(enrichplot)
		
		if(species=='human')
		  {
			DBkeyset='org.Hs.eg.db'
			org='hsa'
			library(org.Hs.eg.db)
			
		  }else if(species=='mouse')
		  {
			library(org.Mm.eg.db)
			DBkeyset='org.Mm.eg.db'
			org='mmu'
		  }else if(species=='rat')
		  {
			library(org.Rn.eg.db)
			#BiocManager::install('org.Rn.eg.db')
			DBkeyset='org.Rn.eg.db'
			org='rno'
		  }else{
			DBkeyset='org.Hs.eg.db'
			library(org.Hs.eg.db)
			org='hsa'
		  }  
		  
		  if(!dir.exists(Reactome_dir))
		  {
			dir.create(Reactome_dir)
		  }
  
		bbb = bitr(aaa, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb='org.Hs.eg.db')
		enrich_result<- enrichPathway(gene=bbb$ENTREZID, pvalueCutoff = 0.05,qvalueCutoff=1,minGSSize=1, readable=TRUE)
		class(enrich_result)
		nr<-nrow(enrich_result)
		if(!is.null(nr))
		{
		  if(nr>0)
		  {
			df_result<-as.data.frame(enrich_result)
			head(df_result)
			df_result$gene_name<-NULL
			#for(i in 1:nrow(df_result))
			#{
			#  df_result[i,'gene_name']<-kegg_geneid_convert(df_result[i,'geneID'],species = "human")
			#}
			write.table(df_result,paste0(Reactome_dir,'/',"Reactome-enrich.tsv"),sep='\t',row.names =FALSE,quote=F)
			
			if(0)
			{
				plota<-NULL
				plota<-dotplot(enrich_result,title="EnrichmentReactome",showCategory=20)+theme(plot.title = element_text(hjust=0.5, face="bold"))
				ggsave(paste0(Reactome_dir,'/',"Reactome",".dotplot.tiff"), plot = plota, width = 9, height = 6,compression='lzw',limitsize=F)
				ggsave(paste0(Reactome_dir,'/',"Reactome",".dotplot.pdf"), plot = plota, width = 9, height = 6,limitsize=F)
				
				
				plota<-NULL
				plota<-barplot(enrich_result, showCategory=20,title="EnrichmentReactome")+theme(plot.title = element_text(hjust=0.5, face="bold"))
				ggsave(paste0(Reactome_dir,'/',"Reactome",".barplot.tiff"), plot = plota, width = 9, height = 6,compression='lzw',limitsize=F)
				ggsave(paste0(Reactome_dir,'/',"Reactome",".barplot.pdf"), plot = plota, width = 9, height = 6,limitsize=F)
			}
			if(1)
			{
				library(stringr)
				library(tidyverse)
				library(ggplot2)
				
				df_plot<-head(df_result,term_number)
				head(df_plot)
				df_plot$matched_all<-0
				df_plot$pathway_ngene<-0
				df_plot$kegg_ngene<-0
				for(i in 1:nrow(df_plot))
				{
				  #i=1
				  df_plot[i,'matched_all']<-as.numeric(strsplit(df_plot[i,'GeneRatio'],'/')[[1]][2])
				  df_plot[i,'pathway_ngene']<-as.numeric(strsplit(df_plot[i,'BgRatio'],'/')[[1]][1])
				  df_plot[i,'kegg_ngene']<-as.numeric(strsplit(df_plot[i,'BgRatio'],'/')[[1]][2])
				}
				head(df_plot,20)
				df_plot<-df_plot[order(df_plot$pvalue,-df_plot$Count),]
				colnames(df_plot)
				df_plot$GeneRatio<-as.numeric(df_plot$Count/df_plot$matched_all)
				df_plot$BgRatio<-as.numeric(df_plot$pathway_ngene/df_plot$kegg_ngene)
				df_plot$Enrichment<-as.numeric(df_plot$GeneRatio/df_plot$BgRatio)
				df_plot<-df_plot[order(-df_plot$Enrichment),]
				
				ggplo2_kegg_barplot(df_plot,name_column='Description',value_column='Count',color_column='p.adjust',
									ylab='Count',title='Enrichment of Reactome Pathway',pw=7,ph=7.5,name_width=50,output_name=paste0(Reactome_dir,'/Reactome.barplot'))
				ggplo2_kegg_dotplot(df_plot,name_column='Description',value_column='GeneRatio',color_column='p.adjust',size_column='Count',
		                              xlab='GeneRatio',title='Enrichment of Reactome Pathway',pw=7,ph=7.5,name_width=50,output_name=paste0(Reactome_dir,'/Reactome.dotplot'))
			}
	
		  }else{
			write.table("",paste0(Reactome_dir,'/',"NO enriched pathway.tsv"),sep='\t',row.names =FALSE,quote=F)
		  }
		  
		}else{
		  write.table("",paste0(Reactome_dir,'/',"NO enriched pathway.tsv"),sep='\t',row.names =FALSE,quote=F)
		}
}



go_output_yingbai<-function(go_bp,dir_output='gokegg_results',term = "BP",term_number=20,pw=7,ph=7.5,name_width=50){
  if(!dir.exists(dir_output))
  {
    dir.create(dir_output)
  }
  #go_bp<-NA
  ngo<-nrow(go_bp)
  if(!is.null(ngo))
  {
    if(ngo>0)
    {
      write.table(as.data.frame(go_bp),paste0(dir_output,'/GO_',term,".tsv"),sep='\t',row.names =FALSE,quote=F)
      df_kk<-as.data.frame(go_bp)
      
      if(0)  ####  使用enrichplot的函数绘制barplot和dotplot
      {
          plota<-NULL
          plota<-dotplot(go_bp,title=paste0("Enrichment of GO_",term),showCategory=20)+theme(plot.title = element_text(hjust=0.5, face="bold"))
          ggsave(paste0(dir_output,'/',"GO_",term,".dotplot.tiff"), plot = plota, width = 12, height = 12,compression='lzw',limitsize=F)
          ggsave(paste0(dir_output,'/',"GO_",term,".dotplot.pdf"), plot = plota, width = 12, height = 12,limitsize=F)
          
          
          plota<-NULL
          plota<-barplot(go_bp, showCategory=20,title=paste0("EnrichmentGO_",term))+theme(plot.title = element_text(hjust=0.5, face="bold"))
          ggsave(paste0(dir_output,'/',"GO_",term,".barplot.tiff"), plot = plota, width = 12, height = 12,compression='lzw',limitsize=F)
          ggsave(paste0(dir_output,'/',"GO_",term,".barplot.pdf"), plot = plota, width = 12, height = 12,limitsize=F)
        }
      
      if(1)   ####  使用ggplot2绘制barplot和dotplot
      {
        library(stringr)
        library(ggplot2)
        df_plot<-head(df_kk,term_number)
        head(df_plot)
        df_plot$matched_all<-0
        df_plot$pathway_ngene<-0
        df_plot$kegg_ngene<-0
        for(i in 1:nrow(df_plot))
        {
          #i=1
          df_plot[i,'matched_all']<-as.numeric(strsplit(df_plot[i,3],'/')[[1]][2])
          df_plot[i,'pathway_ngene']<-as.numeric(strsplit(df_plot[i,4],'/')[[1]][1])
          df_plot[i,'kegg_ngene']<-as.numeric(strsplit(df_plot[i,4],'/')[[1]][2])
        }
        head(df_plot,20)
        df_plot<-df_plot[order(df_plot$pvalue,-df_plot$Count),]
        colnames(df_plot)
        df_plot$GeneRatio<-as.numeric(df_plot$Count/df_plot$matched_all)
        df_plot$BgRatio<-as.numeric(df_plot$pathway_ngene/df_plot$kegg_ngene)
        df_plot$Enrichment<-as.numeric(df_plot$GeneRatio/df_plot$BgRatio)
        df_plot<-df_plot[order(-df_plot$Enrichment),]
        
        rownames(df_plot)<-NULL
        df_plot$index<-as.numeric(rownames(df_plot))
        nlab<-max(nchar(df_plot$Description))
        
        plota<-NULL
        plota <- ggplot(df_plot) + 
          geom_bar(aes(x=reorder(Description,-index),y=Enrichment,fill=-log10(pvalue)),stat='identity',width=0.8)
        plota<-plota+theme(axis.text.x = element_text(size = 10,face='bold'))
        plota<-plota+theme(axis.text.y = element_text(size = 10,face='bold')) 
        plota<-plota+theme(axis.title.y=element_blank(),axis.title.x=element_text(size = 12,face='bold'),legend.position='right')                 
        plota<-plota+labs(y='Enrichment',title=paste0("Enrichment of GO_",term))+scale_fill_gradient(low="blue",high="red")
        plota<-plota+theme(axis.line = element_line(linewidth=1, colour = "black"))                  
        plota<-plota+theme(panel.border = element_blank(),plot.title = element_text(size = 15,hjust=0.5, face="bold"))                                     
        plota<-plota+theme(panel.grid =element_blank(),panel.background = element_blank())    
        plota<-plota+coord_flip()
        plota<-plota+scale_y_continuous(expand=c(0,0))
        plota<-plota+scale_x_discrete(labels=function(x) str_wrap(x, width=name_width))
        ggsave(paste0(dir_output,'/',"GO_",term,".barplot.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
        ggsave(paste0(dir_output,'/',"GO_",term,".barplot.pdf"), plot = plota, width = pw, height = ph)
        
        head(df_plot,20)
        
        df_plot<-df_plot[order(-df_plot$pvalue,df_plot$Count),]
        df_plot<-df_plot[order(df_plot$GeneRatio),]
        rownames(df_plot)<-NULL
        df_plot$index<-as.numeric(rownames(df_plot))
        
        plota<-NULL
        plota <- ggplot(df_plot) + 
          geom_point(aes(x=GeneRatio,y=index,size=Count,col=-log10(pvalue)),alpha = 0.99, position = position_jitter(w = 0.0, h = 0.0))
        plota<-plota+theme(axis.text.x = element_text(size = 10,face='bold'))
        plota<-plota+theme(axis.text.y = element_text(size = 10,face='bold')) 
        plota<-plota+theme(axis.title.y=element_blank(),axis.title.x=element_text(size = 12,face='bold'),legend.position='right')                
        plota<-plota+labs(x='GeneRatio',title=paste0("Enrichment of GO_",term))+scale_colour_distiller(palette = "RdYlBu")
        #+scale_color_gradient(low="blue",high="red")
        plota<-plota+theme(axis.line = element_line(linewidth=1, colour = "black"))                  
        plota<-plota+theme(panel.border = element_blank(),plot.title = element_text(size = 15,hjust=0.5, face="bold"))                                    
        plota<-plota+theme(panel.grid =element_blank(),panel.background = element_blank())    
        plota<-plota+scale_y_continuous(breaks=df_plot$index,labels= str_wrap(df_plot$Description, width=name_width))

        #scale_y_discrete(labels=function(x) str_wrap(x, width=ceiling(nlab/2)))
        ggsave(paste0(dir_output,'/',"GO_",term,".dotplot.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
        ggsave(paste0(dir_output,'/',"GO_",term,".dotplot.pdf"), plot = plota, width = pw, height = ph)
      }
      
      pdf(file=paste0(dir_output,'/',"GO_",term,".GOgraph.pdf"),width=8,height=8)
      plotGOgraph(go_bp)
      dev.off()
      tiff(file=paste0(dir_output,'/',"GO_",term,".GOgraph.tiff"),width=2400,height=2400,res=300,compression='lzw')
      plotGOgraph(go_bp)
      dev.off()
    }else{
      write.table("",paste0(dir_output,'/',"NO enriched GO_",term,".tsv"),sep='\t',row.names =FALSE,quote=F)
    }
  }else{
    write.table("",paste0(dir_output,'/',"NO enriched GO_",term,".tsv"),sep='\t',row.names =FALSE,quote=F)
  }
}

kegg_geneid_convert<-function(string,species = "human")
{
		if(species=='human')
		{
			    DBkeyset='org.Hs.eg.db'
			    organism='hsa'
			    library(org.Hs.eg.db)
			   
		}else if(species=='mouse')
		{
				library(org.Mm.eg.db)
				 DBkeyset='org.Mm.eg.db'
				 organism='mmu'
		}else if(species=='rat')
		{
				library(org.Rn.eg.db)
				#BiocManager::install('org.Rn.eg.db')
				DBkeyset='org.Rn.eg.db'
				organism='rno'
		}else{
				DBkeyset='org.Hs.eg.db'
			    library(org.Hs.eg.db)
			    organism='hsa'
		}
		aaa<-as.numeric(strsplit(string,'/')[[1]])
		bbb = bitr(aaa, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb=DBkeyset)
		string_new<-paste(bbb$SYMBOL,collapse='/')
		return(string_new)
}



plot_celltype_proportion <- function(metadata,plotid='04',x='group',y='celltype',test_group='',control_group='')
{
		library(ggplot2)
		library(reshape2)
		library(ggpubr)
		#install.packages('ggpubr')
		library(grid)
		df_out<-as.data.frame(table(metadata[c(x,y)]))
		colnames(df_out)<-c('x','y','Freq')
		df2_out <- dcast(df_out, x~y)
		colnames(df2_out)[1]<-x
		write.table(df2_out,paste0("fig",plotid,"d.",x,"-",y,".Ncell.tsv"),sep="\t",row.names=F,quote=F)
		mydata1<-melt(df2_out)
		colnames(mydata1)<-c(x,y,'num')
		ncelltype<-ncol(df2_out)-1
		if( test_group!='' & control_group!='')
		{
			print('good')
			mydata1[,x]<-factor(mydata1[,x],levels=c(control_group,test_group))
		}
		plota<-NULL
		plota<- ggplot(data=mydata1)+geom_bar(aes(x=.data[[y]],weight=num,fill=.data[[x]]),position='fill')+coord_flip()
		plota<-plota+theme(plot.margin=unit(c(1,1,3,1.5), 'lines'),axis.title.x=element_text(size=30,color='black',vjust=-5,hjust=0.5),axis.title.y=element_text(size=30,color='black'),
		axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), axis.text.y=element_text(size=20,color='black'),legend.position='None',
		legend.title=element_blank(),legend.text = element_text(size = 20),legend.key.size = unit(1.2, 'lines'))
		plota<-plota+labs(y = expression(paste('Relative\nProportion')),x='')+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = 'black'))
		plota<-plota+scale_fill_manual(values=c('blue','red'))+scale_y_continuous(expand = c(0, 0))
		plotb<-NULL
		plotb<- ggplot(data=mydata1)+geom_bar(aes(x=.data[[y]],weight=num,fill=.data[[x]]),position='stack')+coord_flip()
		plotb<-plotb+theme(plot.margin=unit(c(1,1,3,1.5), 'lines'),axis.title.x=element_text(size=30,color='black',vjust=-5,hjust=0.5),axis.title.y=element_text(size=30,color='black'),
		axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), axis.text.y=element_blank(),legend.position=c(0.8,0.95),
		legend.title=element_blank(),legend.text = element_text(size = 20),legend.key.size = unit(1.2, 'lines'))
		plotb<-plotb+labs(y = expression(paste('Absolute\nnumber')),x='')+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
		axis.line = element_line(colour = 'black'))
		plotb<-plotb+scale_fill_manual(values=c('blue', 'red'))+scale_y_continuous(expand = c(0, 0))
		pearplot<-NULL
		pearplot <- ggarrange(plota, plotb,ncol = 2)                    
		ggsave(paste0("fig",plotid,"d.",x,"-",y,".cell_proportion.tiff"), plot = pearplot, width = 16, height = ncelltype*1+4,compression='lzw')
		ggsave(paste0("fig",plotid,"d.",x,"-",y,".cell_proportion.pdf"), plot = pearplot, width = 16, height = ncelltype*1+4)
}


plot_celltype_proportion_liujie <- function(metadata,plotid='04',x='group',y='celltype',test_group='',control_group='')
{
		library(ggplot2)
		library(reshape2)
		library(ggpubr)
		#install.packages('ggpubr')
		library(grid)
		df_out<-as.data.frame(table(metadata[c(x,y)]))
		colnames(df_out)<-c('x','y','Freq')
		df2_out <- dcast(df_out, x~y)
		colnames(df2_out)[1]<-x
		write.table(df2_out,paste0("fig",plotid,"d.",x,"-",y,".Ncell.tsv"),sep="\t",row.names=F,quote=F)
		mydata1<-melt(df2_out)
		colnames(mydata1)<-c(x,y,'num')
		ncelltype<-ncol(df2_out)-1
		ngroup<-nrow(df2_out)
		if( test_group!='' & control_group!='')
		{
			print('good')
			mydata1$group<-factor(mydata1$group,levels=c(control_group,test_group))
		}
		plota<-NULL
		plota<- ggplot(data=mydata1)+geom_bar(aes(x=.data[[y]],weight=num,fill=.data[[x]]),position='fill')+coord_flip()
		plota<-plota+theme(plot.margin=unit(c(1,1,0.5,1), 'lines'),axis.title.x=element_text(size=10,color='black',vjust=-5,hjust=0.5),
		axis.title.y=element_text(size=10,color='black'),
		axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), axis.text.y=element_text(size=10,color='black'),legend.position='right',
		legend.title=element_blank(),legend.text = element_text(size = 10),legend.key.size = unit(1.2, 'lines'))
		plota<-plota+labs(y = '',x='')+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
		 panel.background = element_blank(),axis.line = element_line(colour = 'black'))
		plota<-plota+scale_fill_manual(values=rainbow(ngroup+2))+scale_y_continuous(expand = c(0, 0),breaks=c(0,0.25,0.5,0.75,1),labels=c('0%','25%','50%','75%','100%'))
		ggsave(paste0("fig",plotid,"d.",x,"-",y,".cell_proportion.ngroup.tiff"), plot = plota, width = 6, height = 4,compression='lzw')
		ggsave(paste0("fig",plotid,"d.",x,"-",y,".cell_proportion.ngroup.pdf"), plot = plota, width = 6, height = 4)
}


plot_celltype_proportion_ngroup <- function(metadata,plotid='04',x='group',y='celltype',test_group='',control_group='',legend_pos=c(0.8,0.95),pw=NA,ph=NA)
{
		  library(ggplot2)
		  library(reshape2)
		  library(ggpubr)
		  #install.packages('ggpubr')
		  library(grid)
		  df_out<-as.data.frame(table(metadata[c(x,y)]))
		  colnames(df_out)<-c('x','y','Freq')
		  df2_out <- dcast(df_out, x~y)
		  colnames(df2_out)[1]<-x
		  write.table(df2_out,paste0("fig",plotid,"d.",x,"-",y,".Ncell.tsv"),sep="\t",row.names=F,quote=F)
		  mydata1<-melt(df2_out)
		  colnames(mydata1)<-c(x,y,'num')
		  ncelltype<-ncol(df2_out)-1
		  ngroup<-nrow(df2_out)
		  if( test_group!='' & control_group!='')
		  {
		    print('good')
		    mydata1$group<-factor(mydata1$group,levels=c(control_group,test_group))
		  }
		  if(is.na(ph))
		  {
			  ph=max(ncelltype,5)		  
		  }
		  if(is.na(pw))
		  {
			  pw=16		  
		  }
			mycolors<-get_colors(ngroup)
		  plota<-NULL
		  plota<- ggplot(data=mydata1)+geom_bar(aes(x=.data[[y]],weight=num,fill=.data[[x]]),position='fill')+coord_flip()
		  plota<-plota+theme(plot.margin=unit(c(1,1,3,1.5), 'lines'),axis.title.x=element_text(size=30,color='black',vjust=-5,hjust=0.5),axis.title.y=element_text(size=30,color='black'),
		                     axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), axis.text.y=element_text(size=20,color='black'),legend.position='None',
		                     legend.title=element_blank(),legend.text = element_text(size = 20),legend.key.size = unit(1.2, 'lines'))
		  plota<-plota+labs(y = expression(paste('Relative\nProportion')),x='')+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = 'black'))
		  plota<-plota+scale_fill_manual(values=mycolors)+scale_y_continuous(expand = c(0, 0))		
		  
		  plotb<-NULL
		  plotb<- ggplot(data=mydata1)+geom_bar(aes(x=.data[[y]],weight=num,fill=.data[[x]]),position='stack')+coord_flip()
		  plotb<-plotb+theme(plot.margin=unit(c(1,1,3,1.5), 'lines'),axis.title.x=element_text(size=30,color='black',vjust=-5,hjust=0.5),axis.title.y=element_text(size=30,color='black'),
		                     axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), axis.text.y=element_blank(),legend.position=legend_pos,
		                     legend.title=element_blank(),legend.text = element_text(size = 20),legend.key.size = unit(1.2, 'lines'))
		  plotb<-plotb+labs(y = expression(paste('Absolute\nnumber')),x='')+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
		                                                                          axis.line = element_line(colour = 'black'))
		  plotb<-plotb+scale_fill_manual(values=mycolors)+scale_y_continuous(expand = c(0, 0))
		  pearplot<-NULL
		  pearplot <- ggarrange(plota, plotb,ncol = 2)                    
		  ggsave(paste0("fig",plotid,"d.",x,"-",y,".cell_proportion.ngroup.tiff"), plot = pearplot,limitsize = FALSE, width = pw, height = ph,compression='lzw')
		  ggsave(paste0("fig",plotid,"d.",x,"-",y,".cell_proportion.ngroup.pdf"), plot = pearplot,limitsize = FALSE, width = pw, height = ph)
}

plot_barplot_proportion_with_label <- function(metadata,plotid='04',x='group',y='celltype',test_group='',control_group='',label_TF=T)
{
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(grid)
  library(dplyr)
  df_out<-as.data.frame(table(metadata[c(x,y)]))
  colnames(df_out)<-c('x','y','Freq')
  
  myresults_pct <- df_out %>% group_by(y) %>% dplyr::mutate(pct=prop.table(Freq))
  ###  bug记录：使用此算法，有时候group_by部分会失效，导致计算的num的比例是占总的num的比例，而不是按y分组后的比例。
  ###  bug记录：并不是算法错误。大多数情况计算是正确的，偶尔会计算失败，然后重启R后，就计算正确了。
  
  myresults_pct<-as.data.frame(myresults_pct)
  
  df2_out <- dcast(df_out, x~y)
  colnames(df2_out)[1]<-x
  write.table(df2_out,paste0("fig",plotid,".",x,"-",y,".count.tsv"),sep="\t",row.names=F,quote=F)
  write.table(myresults_pct,paste0("fig",plotid,".",x,"-",y,".percentage.tsv"),sep="\t",row.names=F,quote=F)
  mydata1<-melt(df2_out)
  colnames(mydata1)<-c(x,y,'num')
  ncelltype<-ncol(df2_out)-1
  ngroup<-nrow(df2_out)
  if( test_group!='' & control_group!='')
  {
    print('good')
    mydata1[,x]<-factor(mydata1[,x],levels=c(control_group,test_group))
  }
  
  pw<-max(6,ncelltype)
  
# "firebrick"
    plota<-NULL
    plota<- ggplot(data=myresults_pct,aes(x=y,y=pct,fill=x))+geom_bar(position='stack',stat='identity',width=0.9)
    if(label_TF)
    {
    plota<-plota+geom_text(aes(label = scales::percent(round(pct,4))),position=position_stack(0.5),col="black",size=5)
	}
    #+scale_y_continuous(label = scales::percent)
    plota<-plota+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), 'lines'),
                       axis.title.x=element_text(size=15,color='black',vjust=-5,hjust=0.5),
                       axis.title.y=element_text(size=15,color='black'),
                       axis.text.x =element_text(size=15,color='black',angle = 45,hjust = 1), 
                       axis.text.y=element_text(size=15,color='black'),legend.position='right',
                       legend.title=element_blank(),legend.text = element_text(size = 10),
                       legend.key.size = unit(1.2, 'lines'))
    plota<-plota+labs(y = '',x='')+theme(panel.grid.major =element_blank(), 
                                         panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = 'black'))
    if(ngroup == 2)
    {
			plota<-plota+scale_fill_manual(values=c('#0072b5','#BC3C29'))
    }
    plota<-plota+scale_y_continuous(expand = c(0, 0.05),breaks=c(0,0.25,0.5,0.75,1),labels=c('0%','25%','50%','75%','100%'))
  ggsave(paste0("fig",plotid,".",x,"-",y,".percentage.barplot.tiff"), plot = plota, width = pw, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,".",x,"-",y,".percentage.barplot.pdf"), plot = plota, width = pw, height = 6)
}

plot_celltype_proportion_bygroup <- function(metadata,plotid='04',x='group',y='celltype')
{
		library(ggplot2)
		library(reshape2)
		library(ggpubr)
		#install.packages('ggpubr')
		library(grid)
		df_out<-as.data.frame(table(metadata[c(x,y)]))
		colnames(df_out)<-c('x','y','Freq')
		df2_out <- dcast(df_out, x~y)
		colnames(df2_out)[1]<-x
		write.table(df2_out,paste0("fig",plotid,"d.",x,"-",y,".Ncell.tsv"),sep="\t",row.names=F,quote=F)
		mydata1<-melt(df2_out)
		colnames(mydata1)<-c(x,y,'num')
		ncelltype<-ncol(df2_out)-1
		ngroup<-nrow(df2_out)
		plotb<-NULL
		plotb<- ggplot(data=mydata1)+geom_bar(aes(x=.data[[y]],weight=num,fill=.data[[x]]),position='stack')
		plotb<-plotb+theme(plot.margin=unit(c(1,1,3,1.5), 'lines'),axis.title.x=element_text(size=30,color='black',vjust=-5,hjust=0.5),axis.title.y=element_text(size=30,color='black'),
		axis.text.x =element_text(size=10,color='black',angle = 45,hjust = 1), axis.text.y=element_text(size=10,color='black'),legend.position='right',
		legend.title=element_blank(),legend.text = element_text(size = 20),legend.key.size = unit(1.2, 'lines'))
		plotb<-plotb+labs(y = expression(paste('cell number')),x='')+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
		axis.line = element_line(colour = 'black'))
		plotb<-plotb+scale_fill_manual(values=rainbow(ngroup+2))+scale_y_continuous(expand = c(0, 0))
		ggsave(paste0("fig",plotid,"d.",x,"-",y,".cell_proportion.ngroup.tiff"), plot = plotb, width = 8, height = 8,compression='lzw')
		ggsave(paste0("fig",plotid,"d.",x,"-",y,".cell_proportion.ngroup.pdf"), plot = plotb, width = 8, height = 8)
}


plot_cellcount_barplot <- function(metadata,plotid='04',x='celltype')
{
		library(ggplot2)
		library(reshape2)
		library(ggpubr)
		#install.packages('ggpubr')
		library(grid)

		df_celltype_count<-as.data.frame(table(metadata[[x]]))
		colnames(df_celltype_count)<-c(x,'Ncells')
		write.table(df_celltype_count, file=paste0("fig",plotid,"c.",x,".Ncells.tsv"), quote=F, sep="\t", row.names=F)
		plota<-NULL
		plota<-ggplot() + 
			geom_bar(data=df_celltype_count,aes(x=.data[[x]],weight=Ncells, fill=.data[[x]])) +
			theme_classic() +
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
			theme(plot.title = element_text(hjust=0.5, face="bold")) +
			ggtitle("Ncells")
		ggsave(paste0("fig",plotid,"c.",x,".Ncells.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
		ggsave(paste0("fig",plotid,"c.",x,".Ncells.pdf"), plot = plota, width = 8, height = 8)
}





#plot_cellcount_pie(metadata,plotid=plotid,x='seurat_clusters',title='celltype')
plot_cellcount_pie <- function(metadata,plotid='04',x='celltype',title='celltype')
{
		library(ggplot2)
		library(reshape2)
		library(ggpubr)
		#install.packages('ggpubr')
		library(grid)
		df_celltype_count<-as.data.frame(table(metadata[[x]]))
		colnames(df_celltype_count)<-c(x,'Ncells')
		write.table(df_celltype_count, file=paste0("fig",plotid,".",x,".Ncells.tsv"), quote=F, sep="\t", row.names=F)
		ggplot2_pie_with_legend(df_celltype_count,group=x,count='Ncells',title=title,output=paste0("fig",plotid,".pie.Ncells.",x))
}



ggplot2_roc_specificity<-function(df_predict_train,title='Train',plotid='13a'){
  library(pROC)
  library(ggplot2)
  #将真实值和预测值整合到一起
  head(df_predict_train)
  obj_roc_train <- roc(df_predict_train$obs,df_predict_train$predict_res)
  auc_train <- round(auc(df_predict_train$obs, df_predict_train$predict_res),4)
  print(auc_train)
  plota<-NULL
  plota<-ggroc(obj_roc_train, colour = 'red', size = 1) +
    ggtitle(paste0(title, ' (AUC = ', auc_train, ')')) +
    theme_minimal()+theme(panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),plot.title = element_text(hjust = 0.5))
  #scale_y_continuous(expand=c(0,0))
  #scale_x_continuous(expand=c(0,0))
  ggsave(file=paste0('fig',plotid,".roc",".randomForest.tiff"),plot=plota,width = 6, height = 6,compression='lzw')
  ggsave(file=paste0('fig',plotid,".roc",".randomForest.pdf"),plot=plota,width = 6, height = 6)
}
ggplot2_roc_fpr<-function(df_predict_train,title='Train',plotid='13a'){
  library(ROCR)
  library(ggplot2)
  #计算AUC(y.values就是AUC)
  obj_rocr_train <- prediction(df_predict_train$predict_res,df_predict_train$obs)
  auc_rocr_train <- performance(obj_rocr_train, "auc")
  auc_rocr_train_value<-slot(auc_rocr_train,"y.values")[[1]]
  auc_rocr_train_value<-round(auc_rocr_train_value,4)
  text_auc<-paste("AUC=", auc_rocr_train_value,sep="")
  df_rocr <- performance(obj_rocr_train, measure = "tpr", x.measure = "fpr")
  
  tiff(file=paste0('fig',plotid,".roc_fpr",".randomForest.tiff"),width=1800,height=1800,res=300,compression='lzw')
  plot(df_rocr,colorize=F,lwd=2,col='red',
       main=title,
       #xlab="False positive rate", ylab="True positive rate",
       xlab="1 - Specificity", ylab="Sensitivity", 
       box.lty=7, box.lwd=2, box.col="gray")
  text(0.25,0.9,text_auc)
  dev.off()
  pdf(file=paste0('fig',plotid,".roc_fpr",".randomForest.pdf"),width=6,height=6)
  plot(df_rocr,colorize=F,lwd=2,col='red',
       main=title,
       #xlab="False positive rate", ylab="True positive rate",
       xlab="1 - Specificity", ylab="Sensitivity", 
       box.lty=7, box.lwd=2, box.col="gray")
  text(0.25,0.9,text_auc)
  dev.off()
}
ggplot2_roc_fpr_color<-function(df_predict_train,title='Train',plotid='13a'){
  library(ROCR)
  library(ggplot2)
  #计算AUC(y.values就是AUC)
  obj_rocr_train <- prediction(df_predict_train$predict_res,df_predict_train$obs)
  auc_rocr_train <- performance(obj_rocr_train, "auc")
  auc_rocr_train_value<-slot(auc_rocr_train,"y.values")[[1]]
  auc_rocr_train_value<-round(auc_rocr_train_value,4)
  text_auc<-paste("AUC=", auc_rocr_train_value,sep="")
  df_rocr <- performance(obj_rocr_train, measure = "tpr", x.measure = "fpr")
  
  tiff(file=paste0('fig',plotid,".roc_fpr",".randomForest.tiff"),width=1800,height=1800,res=300,compression='lzw')
  plot(df_rocr,colorize=T,lwd=2,
       main=title,
       #xlab="False positive rate", ylab="True positive rate",
       xlab="1 - Specificity", ylab="Sensitivity", 
       box.lty=7, box.lwd=2, box.col="gray")
  text(0.25,0.9,text_auc)
  dev.off()
  pdf(file=paste0('fig',plotid,".roc_fpr",".randomForest.pdf"),width=6,height=6)
  plot(df_rocr,colorize=T,lwd=2,
       main=title,
       #xlab="False positive rate", ylab="True positive rate",
       xlab="1 - Specificity", ylab="Sensitivity", 
       box.lty=7, box.lwd=2, box.col="gray")
  text(0.25,0.9,text_auc)
  dev.off()
}





#seurat_ncell_boxplot(pbmc,plotid='08',groupby.x='celltype',groupby.y='group',groupby.z='orig.ident')
seurat_ncell_boxplot <- function(pbmc,plotid='08',groupby.x='celltype',groupby.y='group',groupby.z='orig.ident')
{
		metadata<-pbmc@meta.data
        head(metadata)
        library(ggplot2)
		library(reshape2)
		library(ggpubr)
		#install.packages('ggpubr')
		library(grid)
		df_out<-as.data.frame(table(metadata[c(groupby.x,groupby.y,groupby.z)]))
		colnames(df_out)<-c(groupby.x,groupby.y,groupby.z,'number')
		df_out1<-as.data.frame(table(metadata[[groupby.z]]))
		colnames(df_out1)<-c(groupby.z,'number')
		df_out$sample_cell_number<-0
		for(i in 1:nrow(df_out1))
		{
			df_out$sample_cell_number[df_out[[groupby.z]]==df_out1[i,1]]<-df_out1[i,2]
		}
		df_out$propotion<-df_out$number/df_out$sample_cell_number
		
		
  plota<-ggplot(df_out,aes(x=.data[[groupby.x]],y=propotion,fill=.data[[groupby.y]]))
  plota<-plota+geom_boxplot(col='black',size=0.1,outlier.size=0.3)
  #plota<-plota+scale_fill_manual(values=c('CD'='red', 'control'='blue'))
  #plota<-plota+scale_color_manual(values=c('CD'='red', 'control'='blue'))
  plota<-plota+stat_compare_means(aes(group = .data[[groupby.y]]),label = "p.signif",label.y = 0.4,hide.ns = TRUE,show.legend = FALSE)
  plota<-plota+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                                panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                                plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                                legend.position=c(0.9,0.95),legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                                legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                                axis.text.x=element_text(size=8,color='black',angle=45,hjust=1))
  plota<-plota+labs(y =paste0('cell propotion'),x='')+scale_y_continuous(expand = c(0, 0),limits=c(0,1))
  ggsave(file=paste0('fig',plotid,'.',projectid,".cell-number",".boxplot.tiff"),plot=plota,width = 7, height = 4,compression='lzw')
  ggsave(file=paste0('fig',plotid,'.',projectid,".cell-number",".boxplot.pdf"),plot=plota,width = 7, height = 4)
}




seurat_to_signature_obsoleted1 <- function(pbmc,group.by='celltype',projectid='xxxxxx',plotid='14'){
Idents(pbmc) <-group.by
df_counts <- as.matrix(GetAssayData(object = pbmc, layer = "counts"))
metadata <- pbmc@meta.data
df_counts<-rbind(metadata[,"celltype"],df_counts)
df_counts <- cbind(rownames(df_counts),df_counts)
write.table(df_counts, paste0("fig",plotid,".",projectid,'.',group.by,".signature.tsv"),sep='\t',row.names = F,col.names=F,quote =F)
}

seurat_to_signature_obsoleted2 <- function(pbmc,group.by='celltype',projectid='xxxxxx',plotid='14'){
Idents(pbmc) <-group.by
df_counts <- as.matrix(GetAssayData(object = pbmc, layer = "counts"))
metadata <- pbmc@meta.data
colnames(df_counts)<-metadata[,"celltype"]
#df_counts <- cbind(rownames(df_counts),df_counts)
#write.table(df_counts, paste0("fig",plotid,".",projectid,'.',group.by,".signature.tsv"),sep='\t',row.names = F,col.names=F,quote =F)
#write.table(df_counts, paste0("fig",plotid,".",projectid,'.',group.by,".signature.tsv"),sep='\t',row.names = T,col.names=T,quote =F)
saveRDS(df_counts,file=paste0("fig",plotid,".",projectid,'.',group.by,".signature.rds"))
}

#seurat_to_signature20221202(pbmc.x,group.by='celltype',projectid=projectid,plotid='14')
seurat_to_signature20221202 <- function(pbmc,group.by='celltype',projectid='xxxxxx',plotid='14'){
	#####   https://www.jianshu.com/p/c38c290006e8
Idents(pbmc) <-group.by
markers <- FindAllMarkers(object = pbmc,only.pos=T,min.pct = 0.25,logfc.threshold = 1,return.thresh=0.05)
markers<-unique(markers$gene)
df_counts <- AverageExpression(pbmc, assays = "RNA", layer = "data")[[1]]
df_counts<-df_counts[markers,]
saveRDS(df_counts,file=paste0("fig",plotid,".",projectid,'.',group.by,".signature.rds"))
}

#seurat_to_signature20240305(pbmc.x,group.by='celltype',projectid=projectid,plotid='14')
seurat_to_signature20240305 <- function(pbmc,group.by='celltype',excludes=NULL,projectid='xxxxxx',plotid='14'){
    #####   https://www.jianshu.com/p/c38c290006e8
    Idents(pbmc) <-group.by
    df_counts <- AverageExpression(pbmc, assays = "RNA", layer = "data")[[1]]
    df_counts<-df_counts[,!colnames(df_counts) %in% excludes]
    saveRDS(df_counts,file=paste0("fig",plotid,".",projectid,'.',group.by,".signature.rds"))
    return(df_counts)
  }


#seurat_to_celltype_geneset(pbmc.x,group.by='celltype',projectid=projectid,plotid='14')
seurat_to_celltype_geneset <- function(pbmc,group.by='celltype',projectid='xxxxxx',plotid='14'){
Idents(pbmc) <-group.by
markers <- FindAllMarkers(object = pbmc,only.pos=T,min.pct = 0.25,logfc.threshold = 1,return.thresh=0.05)
genesets<-markers[,c('cluster','gene')]
#list_genesets <- split(genesets$gene, genesets$cluster)
saveRDS(genesets,file=paste0("fig",plotid,".",projectid,'.',group.by,".genesets.rds"))
}


		
seurat_gene_coorelation <- function(pbmc,genea='a',genes=c('b','c'),pvalue_cut=0.01,projectid='xxxxxx',plotid='14'){
	#####   https://www.jianshu.com/p/9850c4d88a67
		exprSet <- pbmc@assays[["RNA"]]@data  
		exprSet<-as.data.frame(t(exprSet))
		expa <- as.numeric(exprSet[,genea])
		cor_data_df <- data.frame(genes)
		for(i in 1:length(genes))
		{
			geneb<-genes[i]
			expb<-as.numeric(exprSet[,geneb])
			test <- cor.test(expa,expb)
			cor_data_df[i,2] <- test$estimate
			cor_data_df[i,3] <- test$p.value
		}
		names(cor_data_df) <- c("symbol","correlation","pvalue")

		cor_data_sig_pos<-cor_data_df[cor_data_df$pvalue<pvalue_cut,]
		cor_data_sig_pos<-na.omit(cor_data_sig_pos)
		library(ggstatsplot)
		if(nrow(cor_data_sig)>1)
		{
			exprSet_used<-exprSet[,c(genea,cor_data_sig$symbol)]
				for(i in 2:ncol(exprSet_used))
				{
					#i=2
					geneb<-colnames(exprSet_used)[i]
					plota<-NULL
					plota<-ggscatterstats(data = exprSet_used,
								   y = !!geneb,
								   x = !!genea,
								   centrality.para = "mean",
								   margins = "both",
								   xfill = "#CC79A7",
								   yfill = "#009E73",
								   marginal.type = "densigram", # #类型可以换成density,boxplot,violin,densigram
								   title = paste0("Relationship between ",genea," and ",geneb))
							ggsave(paste0("fig",plotid,"a.coorelation.",genea,".",geneb,".tiff"), plot = plota, width = 8, height = 8,compression='lzw')
							ggsave(paste0("fig",plotid,"a.coorelation.",genea,".",geneb,".pdf"), plot = plota, width = 8, height = 8)
					plota<-NULL
					plota<-FeatureScatter(pbmc.x, feature1 = genea, feature2 = geneb)
					ggsave(paste0("fig",plotid,"b.FeatureScatter.",genea,".",geneb,".tiff"), plot = plota, width = 8, height = 8,compression='lzw')
					ggsave(paste0("fig",plotid,"b.FeatureScatter.",genea,".",geneb,".pdf"), plot = plota, width = 8, height = 8)		
				}
		}
		write.table(cor_data_df, file=paste0("fig",plotid,".",genea,".coorelation.tsv"), quote=F, sep="\t", row.names=F)
		return(cor_data_df)
}


df_select_min<-function(df_data,keyname='SYMBOL',value_column='pval')
{
  df_data_min<-aggregate(df_data[[value_column]], by = list(df_data[[keyname]]), FUN = min)
  colnames(df_data_min)<-c(keyname,value_column)
  df_data_used<-merge(df_data,df_data_min,by=c(keyname,value_column))
  return(df_data_used)
}


df_select_first<-function(df_data,keyname='SYMBOL',value_column='annotation',value_first='Promoter (<=1kb)')
{
  df_data_a<-df_data[df_data[[value_column]]==value_first,]
  df_data_b<-df_data[df_data[[value_column]]!=value_first,]
  df_data_b0<-df_data_b[! df_data_b[[keyname]] %in% df_data_a[[keyname]],]
  df_data_used<-rbind(df_data_a,df_data_b0)
  return(df_data_used)
}


###   scenic_rds='GSE150825_SCENIC.regulonAUC.rds'
###   seurat_scenic_heatmap(pbmc.x,scenic_rds=scenic_rds,group.by='celltype_sub',projectid=celltypei,plotid='14')
###   seurat_scenic_heatmap_2celltype(pbmc.x,scenic_rds=scenic_rds,group.by='celltype_sub',projectid=celltypei,plotid='15')
seurat_scenic_heatmap_v5<- function(pbmc,scenic_rds='GSE150825_SCENIC.regulonAUC.rds',group.by='seurat_clusters',projectid='xxxxxx',plotid='14'){
  library(AUCell)
  library(SCENIC)
  library(ComplexHeatmap)
  
	Idents(pbmc) <- group.by
	levels_backup<-levels(pbmc)
	
  #devtools::install_github("aertslab/SCENIC")
     regulonAUC<-readRDS(scenic_rds)
	
	if(0)
	{
	sub_regulonAUC <- regulonAUC[,match(colnames(pbmc),colnames(regulonAUC))]
	identical(colnames(sub_regulonAUC), colnames(pbmc))
	}
	metadata<-pbmc@meta.data
	cells_used<-rownames(metadata)[rownames(metadata) %in% colnames(regulonAUC)]
	sub_regulonAUC <- regulonAUC[,cells_used]
	metadata_sub<-metadata[cells_used,]
	print(identical(colnames(sub_regulonAUC), rownames(metadata_sub)))
	
	
	
	cellTypes <- data.frame(row.names = rownames(metadata_sub),celltype = metadata_sub[[group.by]])
	#selectedResolution <- group.by
	cellTypes[,'celltype']<-factor(cellTypes[,'celltype'],levels=levels_backup)

     cellsPerGroup <- split(rownames(cellTypes),cellTypes[,'celltype']) 
     sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] # 去除extened regulons

	  regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells)rowMeans(sub_regulonAUC[,cells]))
	  regulonActivity_byGroup_output<-cbind(rownames(regulonActivity_byGroup),regulonActivity_byGroup)
	  colnames(regulonActivity_byGroup_output)[1] <-'transcription_factor'
	  write.table(regulonActivity_byGroup_output,file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.txt'),row.names = F,sep='\t',quote = FALSE)
	  
	    # Scale expression:
  regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T)) 
  # 同一个regulon在不同cluster的scale处理
  dim(regulonActivity_byGroup_Scaled)
  regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
  regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
  
  pdf(file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.heatmap.scaled.pdf'), height=14, width=4)
  print(pheatmap(regulonActivity_byGroup_Scaled,fontsize_row=3,heatmap_legend_param = list(title = "TF Activity")))
  dev.off()
  ##########################################################

  rss=regulonActivity_byGroup_Scaled
  head(rss)
  
  df = do.call(rbind,
               lapply(1:ncol(rss), function(i){
                 dat= data.frame(
                   path  = rownames(rss),
                   cluster =   colnames(rss)[i],
                   sd.1 = rss[,i],
                   sd.2 = apply(rss[,-i], 1, median)  
                 )
               }))
  df$fc = df$sd.1 - df$sd.2
  
  top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
  top5<-as.data.frame(top5)
  top5[,'cluster']<-factor(top5[,'cluster'],levels=levels_backup)
  rowcn = data.frame(celltype = top5$cluster) 
  rowcn[,'celltype']<-factor(rowcn[,'celltype'],levels=levels_backup)
  #rowcn<-as.matrix(rowcn)
  rss_top5 = rss[top5$path,]
  rss_top5_output<-cbind(rownames(rss_top5),rss_top5)
  colnames(rss_top5_output)[1] <-'transcription_factor'
  write.table(rss_top5_output,file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.top5.txt'),row.names = F,sep='\t',quote = FALSE)
  
  row_labels<-rownames(rss_top5)
  rownames(rss_top5)<-rownames(rowcn)
  class(rowcn)
  #rownames(rowcn) = rownames(rss_top5)
  ngene<-nrow(rss_top5)
  ph<-max(6,ngene/5)
	  pdf(file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.heatmap.scaled.top5.pdf'), height=ph, width=6)
	  print(pheatmap(rss_top5,annotation_row = rowcn,show_rownames = T,cluster_rows=F,heatmap_legend_param = list(title = "TF Activity"),labels_row = row_labels))
	  dev.off()
}

seurat_scenic_heatmap_2celltype_v5<- function(pbmc,scenic_rds='GSE150825_SCENIC.regulonAUC.rds',group.by='seurat_clusters',projectid='xxxxxx',plotid='14'){
  library(AUCell)
  library(SCENIC)
  library(ComplexHeatmap)
  #devtools::install_github("aertslab/SCENIC")
     regulonAUC<-readRDS(scenic_rds)
	
	metadata<-pbmc@meta.data
	cells_used<-rownames(metadata)[rownames(metadata) %in% colnames(regulonAUC)]
	sub_regulonAUC <- regulonAUC[,cells_used]
	metadata_sub<-metadata[cells_used,]
	print(identical(colnames(sub_regulonAUC), rownames(metadata_sub)))
	
	cellTypes <- data.frame(row.names = rownames(metadata_sub),celltype = metadata_sub[[group.by]])
	selectedResolution <- group.by

     cellsPerGroup <- split(rownames(cellTypes),cellTypes[,'celltype']) 
     sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] # 去除extened regulons

	  regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells)rowMeans(sub_regulonAUC[,cells]))
	  regulonActivity_byGroup_output<-cbind(rownames(regulonActivity_byGroup),regulonActivity_byGroup)
	  colnames(regulonActivity_byGroup_output)[1] <-'transcription_factor'
	  write.table(regulonActivity_byGroup_output,file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.txt'),row.names = F,sep='\t',quote = FALSE)
	  
	    # Scale expression:
  regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T)) 
  # 同一个regulon在不同cluster的scale处理
  dim(regulonActivity_byGroup_Scaled)
  regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
  regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
  
  pdf(file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.heatmap.pdf'), height=14, width=4)
  print(pheatmap(regulonActivity_byGroup,fontsize_row=3,heatmap_legend_param = list(title = "TF Activity")))
  dev.off()
  ##########################################################

  rss=regulonActivity_byGroup
  df_plot<-as.data.frame(rss)
  df_plot$fc1<-df_plot[,1]-df_plot[,2]
  df_plot$fc2<-df_plot[,2]-df_plot[,1]
  top5a <- df_plot %>% top_n(5,wt=fc1)
  top5b <- df_plot %>% top_n(5,wt=fc2)
  top5a$cluster<-colnames(top5a)[1]
  top5b$cluster<-colnames(top5b)[2]
  top5<-rbind(top5a,top5b)
  #top5<-rss_top5[!duplicated(top5),]
  rowcn = data.frame(celltype = top5$cluster) 
  rss_top5<-rss[rownames(top5),]  
  rownames(rowcn) = rownames(rss_top5)
  
  rss_top5_output<-cbind(rownames(rss_top5),rss_top5)
  colnames(rss_top5_output)[1] <-'transcription_factor'
  write.table(rss_top5_output,file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.top5.txt'),row.names = F,sep='\t',quote = FALSE)
  
  row_labels<-rownames(rss_top5)
  rownames(rss_top5)<-rownames(rowcn)
  
	  pdf(file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.heatmap.top5.pdf'), height=6, width=6)
	  print(pheatmap(rss_top5,annotation_row = rowcn,show_rownames = T,cluster_rows=F,heatmap_legend_param = list(title = "TF Activity"),labels_row = row_labels))
	  dev.off()

}

seurat_scenic_heatmap_selected_v5<- function(pbmc,scenic_rds='GSE150825_SCENIC.regulonAUC.rds',hub_tf,group.by='seurat_clusters',projectid='xxxxxx',plotid='14',ph=6,pw=6)
  {
    library(AUCell)
    library(SCENIC)
    library(ComplexHeatmap)
    
    Idents(pbmc) <- group.by
    levels_backup<-levels(pbmc)
    
    #devtools::install_github("aertslab/SCENIC")
    regulonAUC<-readRDS(scenic_rds)
    dim(regulonAUC)
    
    if(0)
    {
      sub_regulonAUC <- regulonAUC[,match(colnames(pbmc),colnames(regulonAUC))]
      identical(colnames(sub_regulonAUC), colnames(pbmc))
    }
    metadata<-pbmc@meta.data
    cells_used<-rownames(metadata)[rownames(metadata) %in% colnames(regulonAUC)]
    sub_regulonAUC <- regulonAUC[,cells_used]
    metadata_sub<-metadata[cells_used,]
    print(identical(colnames(sub_regulonAUC), rownames(metadata_sub)))
    
    
    
    cellTypes <- data.frame(row.names = rownames(metadata_sub),celltype = metadata_sub[[group.by]])
    #selectedResolution <- group.by
    cellTypes[,'celltype']<-factor(cellTypes[,'celltype'],levels=levels_backup)
    
    cellsPerGroup <- split(rownames(cellTypes),cellTypes[,'celltype']) 
    sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] # 去除extened regulons
    
    regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells)rowMeans(sub_regulonAUC[,cells]))
    regulonActivity_byGroup_output<-cbind(rownames(regulonActivity_byGroup),regulonActivity_byGroup)
    colnames(regulonActivity_byGroup_output)[1] <-'transcription_factor'
    write.table(regulonActivity_byGroup_output,file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.txt'),row.names = F,sep='\t',quote = FALSE)
    
    
    head(regulonActivity_byGroup)
    dim(regulonActivity_byGroup)
    
    # Scale expression:
    regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T)) 
    # 同一个regulon在不同cluster的scale处理
    dim(regulonActivity_byGroup_Scaled)
    regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
    regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

    head(regulonActivity_byGroup_Scaled)
    hub_tf_used<-hub_tf[hub_tf %in% rownames(regulonActivity_byGroup_Scaled)]
    regulonActivity_byGroup_Scaled_used<-regulonActivity_byGroup_Scaled[hub_tf_used,]
    pdf(file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.heatmap.scaled_selected.pdf'), height=ph, width=pw)
    print(pheatmap(regulonActivity_byGroup_Scaled_used,fontsize_row=3,heatmap_legend_param = list(title = "TF Activity")))
    dev.off()
  }
  

seurat_scenic_heatmap<- function(pbmc,scenic_rds='GSE150825_SCENIC.regulonAUC.rds',group.by='seurat_clusters',projectid='xxxxxx',plotid='14'){
  library(AUCell)
  library(SCENIC)
  library(ComplexHeatmap)
  
	Idents(pbmc) <- group.by
	levels_backup<-levels(pbmc)
	
  #devtools::install_github("aertslab/SCENIC")
     regulonAUC<-readRDS(scenic_rds)
	
	if(0)
	{
	sub_regulonAUC <- regulonAUC[,match(colnames(pbmc),colnames(regulonAUC))]
	identical(colnames(sub_regulonAUC), colnames(pbmc))
	}
	metadata<-pbmc@meta.data
	cells_used<-rownames(metadata)[rownames(metadata) %in% colnames(regulonAUC)]
	sub_regulonAUC <- regulonAUC[,cells_used]
	metadata_sub<-metadata[cells_used,]
	print(identical(colnames(sub_regulonAUC), rownames(metadata_sub)))
	
	
	
	cellTypes <- data.frame(row.names = rownames(metadata_sub),celltype = metadata_sub[[group.by]])
	#selectedResolution <- group.by
	cellTypes[,'celltype']<-factor(cellTypes[,'celltype'],levels=levels_backup)

     cellsPerGroup <- split(rownames(cellTypes),cellTypes[,'celltype']) 
     sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] # 去除extened regulons

	  regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))
	  regulonActivity_byGroup_output<-cbind(rownames(regulonActivity_byGroup),regulonActivity_byGroup)
	  colnames(regulonActivity_byGroup_output)[1] <-'transcription_factor'
	  write.table(regulonActivity_byGroup_output,file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.txt'),row.names = F,sep='\t',quote = FALSE)
	  
	    # Scale expression:
  regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T)) 
  # 同一个regulon在不同cluster的scale处理
  dim(regulonActivity_byGroup_Scaled)
  regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
  regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
  
  pdf(file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.heatmap.scaled.pdf'), height=14, width=4)
  print(pheatmap(regulonActivity_byGroup_Scaled,fontsize_row=3,heatmap_legend_param = list(title = "TF Activity")))
  dev.off()
  ##########################################################

  rss=regulonActivity_byGroup_Scaled
  head(rss)
  
  df = do.call(rbind,
               lapply(1:ncol(rss), function(i){
                 dat= data.frame(
                   path  = rownames(rss),
                   cluster =   colnames(rss)[i],
                   sd.1 = rss[,i],
                   sd.2 = apply(rss[,-i], 1, median)  
                 )
               }))
  df$fc = df$sd.1 - df$sd.2
  
  top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
  top5<-as.data.frame(top5)
  top5[,'cluster']<-factor(top5[,'cluster'],levels=levels_backup)
  rowcn = data.frame(celltype = top5$cluster) 
  rowcn[,'celltype']<-factor(rowcn[,'celltype'],levels=levels_backup)
  #rowcn<-as.matrix(rowcn)
  rss_top5 = rss[top5$path,]
  rss_top5_output<-cbind(rownames(rss_top5),rss_top5)
  colnames(rss_top5_output)[1] <-'transcription_factor'
  write.table(rss_top5_output,file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.top5.txt'),row.names = F,sep='\t',quote = FALSE)
  
  row_labels<-rownames(rss_top5)
  rownames(rss_top5)<-rownames(rowcn)
  class(rowcn)
  #rownames(rowcn) = rownames(rss_top5)
  ngene<-nrow(rss_top5)
  ph<-max(6,ngene/5)
	  pdf(file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.heatmap.scaled.top5.pdf'), height=ph, width=6)
	  print(pheatmap(rss_top5,annotation_row = rowcn,show_rownames = T,cluster_rows=F,heatmap_legend_param = list(title = "TF Activity"),labels_row = row_labels))
	  dev.off()
}

seurat_scenic_heatmap_2celltype<- function(pbmc,scenic_rds='GSE150825_SCENIC.regulonAUC.rds',group.by='seurat_clusters',projectid='xxxxxx',plotid='14'){
  library(AUCell)
  library(SCENIC)
  library(ComplexHeatmap)
  #devtools::install_github("aertslab/SCENIC")
     regulonAUC<-readRDS(scenic_rds)
	
	metadata<-pbmc@meta.data
	cells_used<-rownames(metadata)[rownames(metadata) %in% colnames(regulonAUC)]
	sub_regulonAUC <- regulonAUC[,cells_used]
	metadata_sub<-metadata[cells_used,]
	print(identical(colnames(sub_regulonAUC), rownames(metadata_sub)))
	
	cellTypes <- data.frame(row.names = rownames(metadata_sub),celltype = metadata_sub[[group.by]])
	selectedResolution <- group.by

     cellsPerGroup <- split(rownames(cellTypes),cellTypes[,'celltype']) 
     sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] # 去除extened regulons

	  regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))
	  regulonActivity_byGroup_output<-cbind(rownames(regulonActivity_byGroup),regulonActivity_byGroup)
	  colnames(regulonActivity_byGroup_output)[1] <-'transcription_factor'
	  write.table(regulonActivity_byGroup_output,file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.txt'),row.names = F,sep='\t',quote = FALSE)
	  
	    # Scale expression:
  regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T)) 
  # 同一个regulon在不同cluster的scale处理
  dim(regulonActivity_byGroup_Scaled)
  regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
  regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
  
  pdf(file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.heatmap.pdf'), height=14, width=4)
  print(pheatmap(regulonActivity_byGroup,fontsize_row=3,heatmap_legend_param = list(title = "TF Activity")))
  dev.off()
  ##########################################################

  rss=regulonActivity_byGroup
  df_plot<-as.data.frame(rss)
  df_plot$fc1<-df_plot[,1]-df_plot[,2]
  df_plot$fc2<-df_plot[,2]-df_plot[,1]
  top5a <- df_plot %>% top_n(5,wt=fc1)
  top5b <- df_plot %>% top_n(5,wt=fc2)
  top5a$cluster<-colnames(top5a)[1]
  top5b$cluster<-colnames(top5b)[2]
  top5<-rbind(top5a,top5b)
  #top5<-rss_top5[!duplicated(top5),]
  rowcn = data.frame(celltype = top5$cluster) 
  rss_top5<-rss[rownames(top5),]  
  rownames(rowcn) = rownames(rss_top5)
  
  rss_top5_output<-cbind(rownames(rss_top5),rss_top5)
  colnames(rss_top5_output)[1] <-'transcription_factor'
  write.table(rss_top5_output,file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.top5.txt'),row.names = F,sep='\t',quote = FALSE)
  
  row_labels<-rownames(rss_top5)
  rownames(rss_top5)<-rownames(rowcn)
  
	  pdf(file=paste0('fig',plotid,'.',projectid,'.regulonActivity_byGroup.heatmap.top5.pdf'), height=6, width=6)
	  print(pheatmap(rss_top5,annotation_row = rowcn,show_rownames = T,cluster_rows=F,heatmap_legend_param = list(title = "TF Activity"),labels_row = row_labels))
	  dev.off()

}





seurat_nmf_step1<- function(pbmc,genes=NULL,ranks=2:9,nrun=5,nmf_p='vp5',projectid='projectid',plotid='03')
{
	library(bigmemory)
	library(synchronicity)
	library(NMF)
	library(Biobase)
  ## 高变基因表达矩阵的分解
  #### vm <- pbmc@assays$RNA@data
  vm<-LayerData(pbmc, assay = "RNA", layer = "data") 
  genes_used<-genes[genes %in% rownames(vm)]
  vm_used<-vm[genes_used,]
  vm_used<-nmf_matrix_prepare(vm_used)  
  res_survey <- nmf(vm_used,ranks,nrun=nrun,.opt=nmf_p,method = "snmf/r",seed='nndsvd')
  if(length(ranks)>1)
  {
	  tiff(file=paste0('fig',plotid,'.',projectid,".NMF_survey.tiff"),width=1800,heigh=1800,units='px',res=300,compression='lzw')
	  print(plot(res_survey))
	  dev.off()
	}
  return(res_survey)
}

seurat_nmf_step2<- function(pbmc,genes=NULL,ranks=4,nrun=50,nmf_p='vp5',projectid='projectid',plotid='03',plot_heatmap=F)
{
  library(NMF)
  library(Biobase)
  ## 高变基因表达矩阵的分解
  #### vm <- pbmc@assays$RNA@data
  vm<-LayerData(pbmc, assay = "RNA", layer = "data") 
  genes_used<-genes[genes %in% rownames(vm)]
  vm_used<-vm[genes_used,]
  vm_used<-nmf_matrix_prepare(vm_used)  
  res <- nmf(vm_used,ranks,nrun=nrun,.opt=nmf_p,method = "snmf/r",seed='nndsvd')
    if(plot_heatmap)
    {
	  ####  绘制一致性聚类：
		pdf(file=paste0('fig',plotid,'.',projectid,".NMF_heatmap.pdf"),width=18,heigh=18,onefile=T)
		consensusmap(res)
		dev.off()
	}
	return(res)
}

seurat_nmf_step3<- function(pbmc,nmf.res=NULL,projectid='projectid',plotid='03')
{
  library(NMF)
  library(Biobase)
  ####  导出NMF的分组信息：
		group <- predict(nmf.res) # 提出亚型
		groupdata<-data.frame(group=paste0('C',as.vector(group)),sample=names(group))
		rownames(groupdata)<-groupdata$sample
		head(groupdata)
		dim(groupdata)
		table(groupdata$group)
		write.table(groupdata,paste0('fig',plotid,'.',projectid,".NMF",'.cluster.tsv'),row.names = F,sep = '\t',quote = F)
		
		metadata<-pbmc@meta.data
		metadata$celltype_nmf<-'unknown'
		for(i in 1:nrow(groupdata))
		{
			metadata[rownames(metadata)==rownames(groupdata)[i],'celltype_nmf']<-groupdata[i,'group']
		}
		ng<-length(unique(groupdata$group))
		metadata[metadata$celltype_nmf=='unknown','celltype_nmf']<-paste0('C',ng+1)
		write.table(metadata,paste0('fig',plotid,'.',projectid,".NMF",'.pbmc.celltype_nmf.tsv'),row.names = F,sep = '\t',quote = F)
		
		pbmc@meta.data$celltype_nmf<-metadata$celltype_nmf
		plotid='03'
		plota<-NULL
		plota<-DimPlot(pbmc, reduction = "umap", group.by = 'celltype_nmf',pt.size = 0.5,label=F)
		ggsave(paste0("fig",plotid,"c.umap.celltype_nmf",".tiff"), plot = plota, width = 8, height = 6,compression='lzw')
		ggsave(paste0("fig",plotid,"c.umap.celltype_nmf",".pdf"), plot = plota, width = 8, height = 6)
		return(pbmc)
}

nmf_matrix_prepare<- function(vm)
{
  vm_used<-as.matrix(vm)
  vm_used[vm_used<0]<-0
  vm_used<-vm_used[rowSums(vm_used)>0,]
  ci0<-which(colSums(vm_used) == 0)
  if(length(ci0)>0)
  {
	vm_used<-vm_used[,-ci0]
	}
  return(vm_used)
}


	# file_gsva<-'fig11.YB20220908006.celltype_sub.cell.gsva.tsv'
	# df_gsva<-read.delim(file_gsva)
	# head(df_gsva)
	
	# df_expr<-as.data.frame(t( pbmc.x@assays$RNA@data[markers,]))	
	# head(df_expr)
	
	# seurat_marker_ggcorrplot(df_gsva,df_expr,group.by='celltype',plotid='19')
seurat_marker_ggcorrplot<- function(df_gsva,df_expr,group.by='celltype',plotid='04',pw=NA,ph=NA){
library(psych) 
library(ggcorrplot)  
library(ggthemes) 
	celltypes<-unique(df_gsva$celltype)	
	ncelltype<-length(celltypes)
	ngene<-ncol(df_expr)
	df_cor<-data.frame()
	df_cor_p<-data.frame()
	for(i in 1:ncelltype)
	{
		#i=1
		celltype<-celltypes[i]
		df_gsva_used<-df_gsva[df_gsva$celltype==celltype,'score',drop=F]
		df_expr_used<-df_expr[rownames(df_gsva_used),]
		head(df_gsva_used)
		head(df_expr_used)
		dim(df_gsva_used)
		dim(df_expr_used)
			cor <- corr.test(df_gsva_used,df_expr_used,method = "spearman",adjust = "BH",ci = F)
			matrix_cor<-cor$r
			matrix_cor_p<-cor$p
			rownames(matrix_cor)<-celltype
			rownames(matrix_cor_p)<-celltype
			
			df_cor<-rbind(df_cor,matrix_cor)
			df_cor_p<-rbind(df_cor_p,matrix_cor_p)
	}
	
	df_cor<-as.matrix(df_cor)
	df_cor_p<-as.matrix(df_cor_p)
	
	   if(is.na(pw))
	  {
		pw<-dim(df_cor)[1]
	  }
	  if(is.na(ph))
	  {
		ph<-dim(df_cor)[2]
	  }
	  
	  	  write.table(df_cor, paste0('fig',plotid,".ggcorrplot.r.tsv"), sep='\t', quote=F,row.names=T)
	  write.table(df_cor_p, paste0('fig',plotid,".ggcorrplot.p.tsv"), sep='\t', quote=F,row.names=T)
	
	 pdf(paste0('fig',plotid,".ggcorrplot.a.pdf"),height=ph,width=pw)              #保存图片的文件名称
  par(oma=c(0.5,1,1,1.2))
  print(ggcorrplot(df_cor,method='circle',ggtheme=theme_bw()))
  dev.off()
  tiff(paste0('fig',plotid,".ggcorrplot.a.tiff"),height=ph*300,width=pw*300,res=300,compression='lzw')              #保存图片的文件名称
  par(oma=c(0.5,1,1,1.2))
  print(ggcorrplot(df_cor,method='circle',ggtheme=theme_bw()))
  dev.off()
  
  pdf(paste0('fig',plotid,".ggcorrplot.b.pdf"),height=ph,width=pw)              #保存图片的文件名称
  par(oma=c(0.5,1,1,1.2))
  print(ggcorrplot(df_cor,ggtheme=theme_bw()))
  dev.off()
  tiff(paste0('fig',plotid,".ggcorrplot.b.tiff"),height=ph*300,width=pw*300,res=300,compression='lzw')              #保存图片的文件名称
  par(oma=c(0.5,1,1,1.2))
  print(ggcorrplot(df_cor,ggtheme=theme_bw()))
  dev.off()
  
  pdf(paste0('fig',plotid,".ggcorrplot.c.pdf"),height=ph,width=pw)              #保存图片的文件名称
  par(oma=c(0.5,1,1,1.2))
  #print(ggcorrplot(matrix_cor,method='circle',ggtheme=theme_bw(),insig = "blank",p.mat = matrix_cor_p))
  print(ggcorrplot(df_cor,method='circle',ggtheme=theme_bw(),p.mat = df_cor_p,insig="pch",pch.col = "black"))
  dev.off()
  
  tiff(paste0('fig',plotid,".ggcorrplot.c.tiff"),height=ph*300,width=pw*300,res=300,compression='lzw')              #保存图片的文件名称
  par(oma=c(0.5,1,1,1.2))
  #print(ggcorrplot(matrix_cor,method='circle',ggtheme=theme_bw(),insig = "blank",p.mat = matrix_cor_p))
  print(ggcorrplot(df_cor,method='circle',ggtheme=theme_bw(),p.mat = df_cor_p,insig="pch",pch.col = "black"))
  dev.off()
} 
	
	
###  pbmc<-seurat_vertion4to5(pbmc)
seurat_vertion4to5<-function(pbmc)
  {
    
    pbmc[["RNA5"]] <- as(object = pbmc[["RNA"]], Class = "Assay5")
    DefaultAssay(pbmc)<-'RNA5'
    pbmc[['RNA']]=NULL  
    pbmc[["RNA"]] <- as(object = pbmc[["RNA5"]], Class = "Assay5")
    DefaultAssay(pbmc)<-'RNA'
    pbmc[['RNA5']]=NULL
    return(pbmc)
  }
  
  
