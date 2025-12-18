#==============================================
#                线状分位数图
#==============================================
rm(list = ls()) 
SARAR_QQ <- function(i,EL=NA,PCEL=NA,qtimes=200){
  p = ps[i]
  index = indexs[i]
  nsim1 = dim(EL)[1]
  nsim2 = dim(PCEL)[1]
  as = (1:qtimes-0.5)/qtimes
  
  # # 图一
  # setEPS()
  # postscript(paste0('EL-',error,'-',n,'-',index,'.eps'), width=3.3, height=3.3)
  # el = sort(EL[,i])[ceiling(as*nsim1)]
  # x <- seq(0,max(EL[,i]),5)
  # par(mar = c(4,4,2,1.5))
  # plot(x,x,xaxs = 'i', yaxs = 'i',
  #      xlim =c(0, max(el)),
  #      ylim =c(0, max(el)),
  #      yaxt = 'n', ann = F, type = 'l')
  # axis(2, las = 1)
  # title(main = paste0(' n = ',n,'  p = ',p),
  #       xlab = 'ChiS(p+3) quantile', ylab = 'EL quantile')
  # points(qchisq(as,p+3),el,pch=20,cex=0.5,col='blue1')
  # dev.off()
  
  # 图二
  setEPS()
  postscript(paste0('PCEL1-',error,'-',n,'-',index,'.eps'), width=3.3, height=3.3)
  # PCEL
  qas <- qchisq(as,1)
  pcel = sort(PCEL[,i])[ceiling(as*nsim2)]
  x <- seq(0,max(pcel)+10,5)
  par(mar = c(4,4,2,1.5))
  plot(x,x,xaxs = 'i', yaxs = 'i',
       xlim =c(0, max(pcel)),
       ylim =c(0, max(pcel)),
       yaxt = 'n', ann = F, type = 'l')
  axis(2, las = 1)
  title(main = paste0(' n = ',n,'  p = ',p),
        xlab = 'ChiS(1) quantile', ylab = 'PCEL1 quantile')
  # 点状图
  points(qchisq(as,1),pcel,pch=20,cex=0.5,col='blue1')
  dev.off()
}
#===============================================
#                   赋值画图            
#=============================================== 
error = 1
ns = c(100,225,400,625,900)

# name = 'PCEL1'
# indexs = c(0,0.1,0.2,0.3,0.4,0.5)

name = 'PCEL2'
indexs = c(0.6,0.7,0.8,0.9,1) # 100, 225 有index = 1.1
for(n in ns){
  ps = round(3*n^indexs)
  # EL = read.csv(paste0('EL-',error,'-',n,'.csv'))
  PCEL = read.csv(paste0(name,'-',error,'-',n,'.csv'))
  for(i in 1:length(indexs)){
    # SARAR_QQ(i,EL=EL,PCEL=PCEL)
    # SARAR_QQ(i,EL=EL)
    SARAR_QQ(i,PCEL=PCEL)
  }
}


