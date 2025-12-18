# ===============================================  
#                      csv 转 latex         
# =============================================== 
ns = c(100, 225, 400, 625, 900)
errors = 1
indexs = c(0.6,0.7,0.8,0.9,1,1.1)
for(error in errors){
  for(n in ns){
  	data = round(read.csv(paste0('CP2-',error,'-',n,'.csv')),3)

  	Data = matrix(NA, length(indexs), 4)
    colnames(Data) = c(paste0('CP',1:4))
    rownames(Data) = data[,1]
    Data[, 1] = data[,2]
    Data[, 2] = data[,3]
    Data[, 3] = data[,4]
    Data[, 4] = data[,5]
    

    
    # latex版输出
    library(xtable)
    print(xtable(Data,auto=TRUE,
                 label=paste0('tab-N',n)))
    cat('\n')
  }
}



















                   