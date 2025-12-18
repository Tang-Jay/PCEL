# ===============================================  
#                      csv 转 latex         
# =============================================== 
rm(list = ls()) 
errors = 4
ns = c(100, 225, 400, 625, 900)
indexs = c(0,0.1,0.2,0.3,0.4,0.5)
for(error in errors){
  for(n in ns){
    EL = read.csv(paste0('EL-',error,'-',n,'.csv'))
    nsim = dim(EL)[1]
    ps = round(3*n^indexs)
    q0 = qchisq(0.95, ps+3)
    Q0 = matrix(q0,nsim,length(ps),byrow=TRUE)
    CP0 = colSums(EL<=Q0)/nsim
    
    P <- matrix(ps+3,nsim,length(ps),byrow=TRUE)
    C <- matrix(sqrt(2*(ps+3)),nsim,length(ps),byrow=TRUE)
    MEL <- (EL - P) / C
    q1 = qnorm(0.975)
    Q1 = matrix(q1,nsim,length(ps),byrow=TRUE)
    CP1 = colSums(MEL<=Q1)/nsim
    

    # PCEL1 = read.csv(paste0('PCEL1-',error,'-',n,'.csv'))
    # nsim = dim(PCEL1)[1]
    # ps = round(3*n^indexs)
    # q1 = qchisq(0.95, 1)
    # Q1 = matrix(q1,nsim,length(ps),byrow=TRUE)
    # CP1 = colSums(PCEL1 <= Q1)/nsim
    
    PCEL = round(read.csv(paste0('CP-',error,'-',n,'.csv')),3)
    colnames(PCEL) = c('p','EL',paste0('CP',1:4))
    
    Data = matrix(NA, length(indexs), 9)
    colnames(Data) = c('n','index','p','EL','MEL',paste0('CP',1:4))
    rownames(Data) = c(n,rep('',length(indexs)-1))
    Data[, 1] = c(n,rep('',length(indexs)-1))
    Data[, 2] = indexs
    Data[, 3] = ps
    Data[, 4] = round(CP0, 3)
    Data[, 5] = round(CP1, 3)
    Data[, 6] = PCEL[,3]
    Data[, 7] = PCEL[,4]
    Data[, 8] = PCEL[,5]
    Data[, 9] = PCEL[,6]
    
    
    # latex版输出
    library(xtable)
    print(xtable(Data,auto=T,
                 label=paste0('tab-',error,'-N',n)),
          include.rownames = FALSE,
          # include.colnames = FALSE,  # 移除列名
          )
    cat('\n')
  }
}



















                   