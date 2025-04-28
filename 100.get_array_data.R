args <- commandArgs(trailingOnly = TRUE)

get_array_data<-function(GSE_number = 'GSE121372'){
  library(GEOquery)
  library(limma)
  library(affy)
  
  gset <- getGEO(GSE_number, destdir=".",
                 getGPL = T)       ## 平台文件
  
  exp<-exprs(gset[[1]])
  GPL<-fData(gset[[1]])	## 获取平台信息 
  
  symbol_columns <- grep("symbol", colnames(GPL), ignore.case = TRUE)
  
  # 如果匹配到的列数大于1，则只取第一个匹配的列
  if (length(symbol_columns) > 1) {
    symbol_columns <- symbol_columns[1]
  }else if(length(symbol_columns) == 0){
    stop('自动获取数据失效，请手动使用表格进行ID转换！！！')
  }
  
  gpl_id <- GPL[, symbol_columns, drop = FALSE]
  gpl_id <- subset(x = gpl_id, gpl_id[,1] != "")
  gpl_id[,1]<-data.frame(sapply(gpl_id[,1],function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
  gpl_id <- gpl_id[!duplicated(gpl_id[,1]), ,drop = F]
  
  df_arry <- df_common_merge(gpl_id,exp)
  rownames(df_arry) <- df_arry[,1]
  df_arry <- df_arry[,-1]
  
  write.table(df_arry,file = paste0(GSE_number,'.txt'),sep = '\t')
  print('成功获取数据！！！')
}

#快速行名取交集两个数据框
df_common_merge<-function(df_1,df_2,plotid='1',projectid,rm_low=F,exp_tsv=F)
{
  if(rm_low)
  {
    df_1 <- df_1[rowMeans(df_1)>1,]
    df_2 <- df_2[rowMeans(df_2)>1,]
  }
  
  df_merge<-merge(df_1,df_2,by=0)
  
  if(exp_tsv)
  {
    df_row<-data.frame(df_merge[,1],row.names = df_merge[,1])
    df_1_m<-merge(df_row,df_1,by=0)
    rownames(df_1_m)<-df_1_m[,1]
    df_1_m<-df_1_m[,-1:-2]
    write.table(df_1_m, file= paste0("fig",plotid,"A.train.",projectid,".used.tsv"),sep="\t",row.names=T)
    
    df_2_m<-merge(df_row,df_2,by=0)
    rownames(df_2_m)<-df_2_m[,1]
    df_2_m<-df_2_m[,-1:-2]
    write.table(df_2_m, file= paste0("fig",plotid,"B.test.",projectid,".used.tsv"),sep="\t",row.names=T)
  }
  
  rownames(df_merge)<-df_merge[,1]
  df_merge[,1]<-NULL
  return(df_merge)
}


get_array_data(as.character(args[1]))