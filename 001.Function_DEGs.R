rawdata_colname_match <- function(data, column_order) {
  # 获取数据框的列名
  original_columns <- colnames(data)
  
  # 检查向量中的列名是否都在数据框的列名中
  missing_columns <- setdiff(column_order, original_columns)
  if (length(missing_columns) > 0) {
    warning(paste("以下列名在数据框中不存在，将被忽略：", paste(missing_columns, collapse = ", ")))
  }
  
  # 保留向量中的列名，并按照向量的顺序排列
  common_columns <- intersect(column_order, original_columns)
  reordered_data_frame <- data[, common_columns]
  
  # 检查是否有不在向量中的列
  extra_columns <- setdiff(original_columns, column_order)
  if (length(extra_columns) > 0) {
    warning(paste("以下列名不在向量中，已被移除：", paste(extra_columns, collapse = ", ")))
  }
  
  return(reordered_data_frame)
}


#快速ID转换
ID_convert<-function(df,fromType = "ENSEMBL",toType="SYMBOL",species='human') #ENSEMBL  ENTREZID
{
  library(dplyr)
  
  if(fromType == toType){
    return(df)
  }
  
  library(clusterProfiler)
  if(species=='human')
  {
    #BiocManager::install('org.Hs.eg.db')
    OrgDb='org.Hs.eg.db'
    library(org.Hs.eg.db)
  }else if(species=='mouse')
  {
    #BiocManager::install('org.Mm.eg.db')
    library(org.Mm.eg.db)
    OrgDb='org.Mm.eg.db'
  }else if(species=='rat')
  {
    library(org.Rn.eg.db)
    #BiocManager::install('org.Rn.eg.db')
    OrgDb='org.Rn.eg.db'
  }else{
    OrgDb='org.Hs.eg.db'
    library(org.Hs.eg.db)
  }
  
  if(class(df) == "character"){
    geneid_to_symbol<- bitr(df, fromType = fromType, toType=toType,OrgDb = OrgDb)
    geneid_to_symbol<-geneid_to_symbol$SYMBOL
    return(geneid_to_symbol)
  }
  
  df_used<-df
  df_used$gene<-rownames(df_used)
  genes_ensembl_id<-df_used$gene
  geneid_to_symbol<- bitr(genes_ensembl_id, fromType = fromType, toType=toType,OrgDb = OrgDb)
  colnames(geneid_to_symbol)<-c('gene','gene_name')
  head(geneid_to_symbol)
  dim(geneid_to_symbol)
  
  geneid_to_symbol <- dplyr::distinct(geneid_to_symbol)
  
  df_used <- merge(df_used,geneid_to_symbol,by='gene')
  
  duplicated_row_names <- duplicated(df_used$gene_name)
  df_used <- df_used[!duplicated_row_names, ]
  rownames(df_used) <- df_used$gene_name
  
  df_used$gene <- NULL
  df_used$gene_name <- NULL
  
  return(df_used)
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


### input主要有两个矩阵：df_count和df_sample。
### 要求:1,df_tpm是log转化后的tpm，通常为log(TPM+1)，通常行名为基因symbol，列名为sampleid
### 要求:2,df_sample是样本分组矩阵，行名为sampleid。其中一列为group。矩阵中可以包含其它列。
### 要求:3,df_count的列名与df_sample的行名一一对应，且顺序一致（此条经优化后，可以不严格遵循）。
###  wilcox_test的fc列pvalue列fdr列：logFC、P.Value、adj.P.Val
rna_seq_tpm_wilcox_test<-function(df_tpm,df_sample,plotid=2,gcontrol='Normal',gtest='Tumor',
                                  fc_cutoff=2,fc_column='log2FoldChange',pvalue_cutoff=0.05,pvalue_column='padj',
                                  go_analysis=F,species='human',recutoff=F,ID_type)
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
    
    if(go_analysis)
    {
      GO_KEGG_jzx(DEG_sig$gene,species = species,pvalueCutoff=1,term_number=20,pw=30,ph=40,txtsize = 10,ID_type = ID_type)
    }
  }
  
  return(DEG_sig$gene)
}

### input主要有两个矩阵：df_count和df_sample。
### 要求:1,df_tpm是log转化后的tpm，通常为log(TPM+1)，通常行名为基因symbol，列名为sampleid
### 要求:2,df_sample是样本分组矩阵，行名为sampleid。其中一列为group。矩阵中可以包含其它列。
### 要求:3,df_count的列名与df_sample的行名一一对应，且顺序一致（此条经优化后，可以不严格遵循）。
###  wilcox_test的fc列pvalue列fdr列：logFC、P.Value、adj.P.Val
rna_seq_tpm_limma_test<-function(df_tpm,df_sample,plotid=2,gcontrol='Normal',gtest='Tumor',
                                  fc_cutoff=2,fc_column='log2FoldChange',pvalue_cutoff=0.05,pvalue_column='padj',
                                  go_analysis=F,species='human',recutoff=F,ID_type)
{
  # 去除表达量过低的基因
  df_tpm[df_tpm < 0] <- 0
  df_tpm <- df_tpm[rowMeans(df_tpm)>0,]
  
  #排序
  df_tpm <- rawdata_colname_match(data = df_tpm,column_order = rownames(df_sample))
  
  df_sample_used<-df_sample[df_sample$group %in% c(gcontrol,gtest),,drop=F]
  df_sample_used$group<-factor(df_sample_used$group,levels=c(gtest,gcontrol))
  
  comparison<-paste0(gtest,'_vs_',gcontrol)
  group_list<-df_sample_used$group
  df_group<-table(group_list)
  df_group<-as.data.frame(df_group)
  colnames(df_group)<-c('group','n_sample')
  
  group_list <- factor(df_sample$group, levels = c(gtest,gcontrol))
  design <- model.matrix(~0 + group_list)
  colnames(design) = levels(group_list)
  rownames(design) = colnames(df_tpm)
  
  if(range(df_tpm)[2] > 20)
  {
    group_mean <- apply(df_tpm,1,function(x) aggregate(x,by=list(group_list),mean)$x)
    v <- voom(df_tpm,design, normalize="quantile",plot=F)
    
    #voomE
    df_tpm_used<-v$E
    voom_matrix_out<-cbind(rownames(df_tpm_used),df_tpm_used)
    colnames(voom_matrix_out)[1]<-"gene_name"
    write.table(voom_matrix_out,file=paste0('fig1.',projectid,'.limma.voomE.tsv'),row.names = F,sep='\t',quote = FALSE)
  }else{
    group_mean <- apply(2^df_tpm-1,1,function(x) aggregate(x,by=list(group_list),mean)$x)
    df_tpm_used <- df_tpm
    write.table(df_tpm_used,file=paste0('fig1.',projectid,'.limma.tpm_used.tsv'),row.names = F,sep='\t',quote = FALSE)
  }
  
  write.table(df_group,file = paste0('fig',plotid,'a.',projectid,".sample_number.tsv"),row.names = F,sep='\t',quote = FALSE)
  
  if(df_group[1,2]>0 & df_group[2,2]>0)
  {
    print(df_group)
  }else{
    print('sample number good!')
    return()
  }

  if(!recutoff)
  {
    if(1)   #####  弱信号筛选和平均TPM计算
    {
      group_mean<-t(group_mean)
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
    
    if(1)  ## 差异分析
    {
      #比较矩阵，进行两两之间的比较
      fit <- lmFit(df_tpm_used, design)
      contrast_expr <- paste0(gtest, "-", gcontrol)
      contrast <- limma::makeContrasts(contrasts = contrast_expr, levels = design)
      fit2 <- contrasts.fit(fit, contrast)
      fit2=eBayes(fit2)
      
      if(1)  ## 筛选出上调和下调基因
      {
        # topTable
        DEG = topTable(fit2, coef=1,number=Inf)
        # 去掉那些NA值
        DEG = na.omit(DEG)
        head(DEG,6)
        
        #筛选出上调和下调基因
        DEG$change = as.factor(
          ifelse(DEG$P.Value < pvalue_cutoff & abs(DEG$logFC) > log2(fc_cutoff),
                 ifelse(DEG$logFC > log2(fc_cutoff) ,'UP','DOWN'),'NOT'))
        DEG_full<-merge(DEG,group_mean,by='row.names')       
        colnames(DEG_full)[1]<-'gene'
        print(dim(DEG_full))
        write.table(DEG_full,file=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.all.tsv'),row.names = F,sep='\t',quote = FALSE)
        
        DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]
        print(dim(DEG_sig))
        write.table(DEG_sig,file=paste0('fig',plotid,'f.',projectid,'.limma.',comparison,'.deg.sig.tsv'),row.names = F,sep='\t',quote = FALSE)
        
      }
      
    }
  }else{
    DEG_full <- read.table(file = paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.all.tsv'),header=T,sep='\t',check.names=FALSE)  
    DEG_full$change = as.factor(
      ifelse(DEG_full[[pvalue_column]] < pvalue_cutoff & abs(DEG_full[[fc_column]]) > log2(fc_cutoff),
             ifelse(DEG_full[[fc_column]] > log2(fc_cutoff) ,'UP','DOWN'),'NOT'))
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
    phenotype<-df_sample_used[,'group',drop=F]
    df_tpm_used_heatmap<-df_tpm_used[rownames(df_tpm_used) %in% DEG_sig$gene,]
    
    ggheatmap_yingbio_2groups(df_tpm_used_heatmap,phenotype,geneid='gene',pheight=10,pw=10,output=paste0('fig',plotid,'h.',projectid,".Heatmap.",comparison),save.data=T)
  }
  
  if(1) ########## 火山图
  {
    head(DEG_full)
    ggplot2_volcano(DEG_full[,c(pvalue_column,fc_column)],comparison=comparison,output=paste0('fig',plotid,'g.',projectid,".volcano.",comparison),
                    fc_threshold=2^log2(fc_cutoff),pvalue_threshold=pvalue_cutoff,ylab=pvalue_column,ymax=300)
  }
  
  if(nrow(DEG_sig)>1)  ######  进行GO和KEGG pathway分析
  {
    library(stats)
    if(0)
    {
      DEG_full <- read.table(file = paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.all.tsv'),header=T,sep='\t',check.names=FALSE)  
      DEG_sig <- read.table(file = paste0('fig',plotid,'f.',projectid,'.tpm.wilcox.',comparison,'.deg.sig.tsv'),header=T,sep='\t',check.names=FALSE)  
    }
    
    if(go_analysis)
    {
      GO_KEGG_jzx(DEG_sig$gene,species = species,pvalueCutoff=1,term_number=20,pw=10,ph=10,txtsize = 20,ID_type = ID_type)
    }
  }
  
  return(DEG_sig$gene)
}

### input主要有两个矩阵：df_count和df_sample。
### 要求:1,df_count是计数矩阵，通常行名为基因symbol，列名为sampleid
### 要求:2,df_sample是样本分组矩阵，行名为sampleid。其中一列为group。矩阵中可以包含其它列。
### 要求:3,df_count的列名与df_sample的行名一一对应，且顺序一致（此条经优化后，可以不严格遵循）。
###  deseq2的fc列pvalue列fdr列：log2FoldChange、P.Value、padj
rna_seq_count_deseq2<-function(df_count,df_sample,plotid=4,gcontrol='Control',gtest='ALS',
                               fc_cutoff=2,fc_column='log2FoldChange',pvalue_cutoff=0.05,pvalue_column='padj',
                               go_analysis=F,species='human',recutoff=F,countfilter=1,ID_type)
{
  library(stringr)
  library(DESeq2)
  library(BiocParallel)
  
  df_sample_used<-df_sample[df_sample$group %in% c(gcontrol,gtest),,drop=F]
  df_sample_used$group<-factor(df_sample_used$group,levels=c(gtest,gcontrol))
  
  samples_df_count<-colnames(df_count)
  samples_df_sample<-rownames(df_sample_used)
  samples_venn<-samples_df_sample[samples_df_sample %in% samples_df_count]
  df_sample_used<-df_sample_used[samples_venn,,drop=F]
  
  group_list<-df_sample_used$group
  df_group<-table(group_list)
  print(group_list)
  
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

      #DEseq2均一化
      colData <- data.frame(row.names=colnames(df_count_used), group_list)
      head(colData)
      dds <- DESeqDataSetFromMatrix(df_count_used, DataFrame(group_list), design= ~ group_list)
      
      dds_norm <- DESeq(dds,parallel = T) 
      sizeFactors(dds_norm)
      head(dds_norm)
      
      normalized_count<-as.data.frame(counts(dds_norm,normalized=TRUE))

      sum_col<-apply(normalized_count,2,sum)

      write.table(sum_col,file = paste0('fig',plotid,'b.',projectid,".count_colsum.tsv"),row.names = T,sep='\t',quote = FALSE)

      
      rld<-vst(dds_norm,blind=F)
      log_normalized_count<-assay(rld)
      
      normalized_count_out<-cbind(rownames(normalized_count),normalized_count)
      colnames(normalized_count_out)[1]<-'gene'

      write.table(normalized_count_out,file = paste0('fig',plotid,'c.',projectid,".DESeq2.ncount.tsv"),row.names = FALSE,sep='\t',quote = FALSE)
      write.table(log_normalized_count,file = paste0('fig',plotid,'c.',projectid,".DESeq2.ncount_vst.tsv"),row.names = T,sep='\t',quote = FALSE)
      
      group_mean=apply(normalized_count,1,function(x) aggregate(x,by=list(group_list),mean)$x)
      group_mean<-t(group_mean)

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
    
    padj_max<-max(res$padj,na.rm=TRUE)
    res$padj[is.na(res$padj)]<-padj_max
    if(1)   ######  输出差异表达基因表格
    {
      res$change = as.factor(
        ifelse(res[[pvalue_column]] < pvalue_cutoff & abs(res[[fc_column]]) > logFC_cutoff,
               ifelse(res[[fc_column]] > logFC_cutoff ,'UP','DOWN'),'NOT'))
      
      DEG_full<-merge(res,group_mean[,c(gtest,gcontrol)],by='row.names')       
      colnames(DEG_full)[1]<-'gene'
      
      write.table(DEG_full,file=paste0('fig',plotid,'g.',projectid,'.count.DESeq2.',comparison,'.deg.all.tsv'),row.names = F,sep='\t',quote = FALSE)
      DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]
      write.table(DEG_sig,file=paste0('fig',plotid,'g.',projectid,'.count.DESeq2.',comparison,'.deg.sig.tsv'),row.names = F,sep='\t',quote = FALSE)
      
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

    log_normalized_count_heatmap<-log_normalized_count[rownames(log_normalized_count) %in% DEG_sig$gene,]
    
    ggheatmap_yingbio_2groups(log_normalized_count_heatmap,phenotype,geneid='gene',pheight=10,pw=10,output=paste0('fig',plotid,'h.',projectid,".Heatmap.",comparison),save.data=T)
  }
  
  if(1) ########## 火山图
  {
    head(DEG_full)
    ggplot2_volcano(DEG_full[,c(pvalue_column,fc_column)],comparison=comparison,output=paste0('fig',plotid,'v.',projectid,".volcano.",comparison),
                    fc_threshold=2^logFC_cutoff,pvalue_threshold=pvalue_cutoff,ylab=pvalue_column,ymax=300)
  }
  
  if(nrow(DEG_sig)>1)  ######  进行GO和KEGG pathway分析
  {
    library(stats)
    if(0)
    {
      DEG_full <- read.table(file = paste0('fig4.',projectid,'.count.DESeq2.',comparison,'.deg.all.tsv'),header=T,sep='\t',check.names=FALSE)  
      DEG_sig <- read.table(file = paste0('fig4.',projectid,'.count.DESeq2.',comparison,'.deg.sig.tsv'),header=T,sep='\t',check.names=FALSE)  
    }

    if(go_analysis)
    {
      GO_KEGG_jzx(DEG_sig$gene,species = species,pvalueCutoff=1,term_number=20,pw=15,ph=20,txtsize = 20,ID_type = ID_type)
    }
  }
  
  return(DEG_sig$gene)
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

ggheatmap_jzx_2groups<-function(mydata,phenotype,orderby=NA,show_colname = F,show_rowname = F,
                                    geneid='miRNA_ID',pheight=NA,pw=6,fontsize_row = 5,output='heatmap',save.data=F)
{
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
                                    fontsize_row=fontsize_row,
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
      plota<-emapplot(ego_pair,  layout="kk", size_category=1.5,min_edge = 0.8) 
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
		}else if(species=='tabacum')
		{
		  DBkeyset <- loadDb("org.Nicotiana_tabacum.eg.sqlite")
		  org='nta'
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

GO_jzx<-function(aaa,dir_output='GO_results',species = "human",go_pvalueCutoff=0.1,term_number=20,pw=30,ph=20,txtsize = 5,ID_type)
{
  library(ggplot2)
  library(topGO)
  library(clusterProfiler)
  library(enrichplot)
  library(AnnotationDbi)
  
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
  }else if(species=='tabacum')
  {
    DBkeyset <- loadDb("org.Nicotiana_tabacum.eg.sqlite")
    org='nta'
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
  
  if(ID_type == 'ENTREZID'){
    bbb <- data.frame(
      ENTREZID = aaa
      )
  }else{
    bbb = try(bitr(aaa, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=DBkeyset),silent = T)
  }
  
  if(class(bbb) =='try-error')
  {
    go_result<-data.frame()
  }else
  {
    go_result <- enrichGO(gene = bbb$ENTREZID,
                          OrgDb = DBkeyset,
                          ont = "ALL",
                          pAdjustMethod = "BH", 
                          pvalueCutoff = go_pvalueCutoff,
                          qvalueCutoff = 1,
                          minGSSize=1,
                          readable = TRUE)
  }
  
  go_output_all(go_result,term_number=term_number,pw=pw,ph=ph,dir_output=go_dir,term = "ALL",txtsize = txtsize)
  
  go_result_combi<-as.data.frame(go_result)
  #go_bp_combi<-go_bp_combi[go_bp_combi$ID=='xxx',]
  n_bp<-nrow(go_result_combi)
  if(n_bp)
  {
    go_result_combi$Ontology<-'Gene Ontology'
  }
}

kegg_yingbai_jzx<-function(aaa,dir_output='KEGG_results',species = "human",term_number=20,kegg_pvalueCutoff=0.05,pw=7,ph=7.5,txtsize = 20,ID_type)
{
  library(ggplot2)
  library(topGO)
  library(clusterProfiler)
  library(enrichplot)
  library(stringr)
  library(AnnotationDbi)
  
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
  }else if(species=='tabacum')
  {
    DBkeyset <- loadDb("org.Nicotiana_tabacum.eg.sqlite")
    org='nta'
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
  
  if(ID_type == 'ENTREZID'){
    bbb <- data.frame(
      ENTREZID = aaa
    )
  }else{
    bbb = try(bitr(aaa, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=DBkeyset),silent = T)
  }
  
  if(class(bbb) =='try-error')
  {
    print('bad')
    nres<-0
  }else
  {
    enrich_result <- enrichKEGG(gene = bbb$ENTREZID,organism = org,pvalueCutoff = kegg_pvalueCutoff,qvalueCutoff=1,minGSSize=1,use_internal_data =F)
    nres<-nrow(enrich_result)
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
    
    if(1)  ####  clusterProfiler绘制barplot和dotplot
    {
      plota<-NULL
      plota<-dotplot(enrich_result,title="Enrichment KEGG",showCategory=term_number,)+theme(plot.title = element_text(hjust=0.5, face="bold"))+
        theme(
          axis.text.y = element_text(size = txtsize),  # y轴标签字体
        ) +
        scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))
      ggsave(paste0(kegg_dir,'/',"KEGG",".dotplot.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
      ggsave(paste0(kegg_dir,'/',"KEGG",".dotplot.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
      
      plota<-NULL
      plota<-barplot(enrich_result, showCategory=term_number,title="Enrichment KEGG")+theme(plot.title = element_text(hjust=0.5, face="bold"))+
        theme(
          axis.text.y = element_text(size = txtsize),  # y轴标签字体
        ) +
        scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))
      ggsave(paste0(kegg_dir,'/',"KEGG",".barplot.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
      ggsave(paste0(kegg_dir,'/',"KEGG",".barplot.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
    }
    
  }else{
    write.table("",paste0(kegg_dir,'/',"NO enriched pathway.tsv"),sep='\t',row.names =FALSE,quote=F)
  }
}

GO_KEGG_jzx<-function(aaa,species = "human",pvalueCutoff=1,term_number=20,pw=30,ph=40,txtsize = 10 , ID_type)
{
  GO_jzx(aaa,species = species,dir_output = 'GO_results', go_pvalueCutoff = pvalueCutoff, term_number = term_number, pw = pw, ph= ph,txtsize = txtsize,ID_type =ID_type)
  kegg_yingbai_jzx(aaa,dir_output = 'KEGG_results',species = species,term_number = term_number, kegg_pvalueCutoff = pvalueCutoff,pw=pw,ph=ph,txtsize = txtsize,ID_type =ID_type)
}

go_output_all<-function(go_bp,dir_output='gokegg_results',term = "BP",term_number=20,pw=30,ph=20,txtsize = 5)
{
  library(stringr)
  
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
      
      plota<-NULL
      plota<-barplot(go_bp, split="ONTOLOGY")+
        facet_grid(ONTOLOGY~., scale = 'free')+
        theme(
          axis.text.y = element_text(size = txtsize),  # y轴标签字体
        ) +
        scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))#柱状图
      ggsave(paste0(dir_output,'/',"GO_",term,".barplot.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
      ggsave(paste0(dir_output,'/',"GO_",term,".barplot.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
      
      plota<-NULL
      plota<-dotplot(go_bp, split="ONTOLOGY")+
        facet_grid(ONTOLOGY~., scale = 'free')+
        theme(
          axis.text.y = element_text(size = txtsize),  # y轴标签字体
        ) +
        scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))#点状图
      ggsave(paste0(dir_output,'/',"GO_",term,".dotplot.tiff"), plot = plota, width = pw, height = ph,compression='lzw',limitsize=F)
      ggsave(paste0(dir_output,'/',"GO_",term,".dotplot.pdf"), plot = plota, width = pw, height = ph,limitsize=F)
      
    }else{
      write.table("",paste0(dir_output,'/',"NO enriched GO_",term,".tsv"),sep='\t',row.names =FALSE,quote=F)
    }
  }else{
    write.table("",paste0(dir_output,'/',"NO enriched GO_",term,".tsv"),sep='\t',row.names =FALSE,quote=F)
  }
}

