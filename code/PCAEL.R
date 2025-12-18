#===============================================
#                 SARAR经验似然统计量
#===============================================
rm(list = ls()) 
source('GlambdaChen.R')
library('MASS')
library('sp')
library('sf')
library('spData')
library('terra')
library('spdep')
SARAR_Data <- function(m,error=1,indexs=c(0,0.1),nsim=1000){
  # ==================修改参数====================
  c = 3
  rou1 = 0.85 
  rou2 = 0.15
  n = m*m
  ps = round(c*n^indexs)
  # ==================函数准备====================
  f<-function(Matrix,Vector){
    irow = nrow(Matrix)
    v = c(0)
    for(i in 2:irow){
      v[i] =Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
    }
    return(v)
  }
  tr<-function(Matrix){
    return(sum(diag(Matrix)))
  }
  L2sqare<-function(Matrix){
    return(sum(diag(Matrix)^2))
  }
  # ==================数据准备====================
  Wnb = cell2nb(m,m,type='queen')
  Ws = nb2listw(Wnb)
  Wn = listw2mat(Ws)
  Mn = Wn
  In = diag(n)
  An = In - rou1*Wn
  Bn = In - rou2*Mn
  Ani = solve(An)
  Bni = solve(Bn)
  # 估计方程准备
  Gn = Bn%*%Wn%*%Ani%*%Bni
  Gnn = 1/2*(Gn + t(Gn))
  Hn = Mn%*%Bni
  Hnn = 1/2*(Hn + t(Hn))
  g = diag(Gnn)
  h = diag(Hnn)
  # ==================开始模拟====================
  EL = matrix(NA,length(ps),6)
  colnames(EL)=c('p','EL','PCA1','PCA2','PCA3','PCA4')
  
  cat('nsim =',nsim,'error =',error,'\n')
  cat('n  ','p ','EL ',' PCAEL1 ',' PCAEL2 ',' PCAEL3 ',' PCAEL4 ','\n')
  j = 0
  for(p in ps){
    j  = j + 1
    beta = matrix(rep(1,p),p,1)
    
    f1 = 0
    f2 = 0
    f3 = 0
    f4 = 0
    f5 = 0
    for(i in 1:nsim){
      Xn = mvrnorm(n,mu = rep(0,p),Sigma = diag(p))
      
      if(error == 1){En = rnorm(n,0,1);sigma2=1;sigma4=sigma2^2;mu3=0;mu4=3*sigma4}
      if(error == 2){En = rnorm(n,0,sqrt(0.75));sigma2 = 0.75;sigma4=sigma2^2;mu3=0;mu4=3*sigma4}
      if(error == 3){En = rt(n,5);sigma2=5/(5-2);sigma4=sigma2^2;sigma4=sigma2^2;mu3=0;mu4=25}
      if(error == 4){En = rchisq(n,4)-4;sigma2=8;sigma4=sigma2^2;mu3=32;mu4=384}
      
      e = En
      b = t(Xn)%*%t(Bn)
      s = Bn%*%Wn%*%Ani%*%Xn%*%beta
      # 模拟Yi(程序运行不需要Yi值)
      # Yn = Ani%*%Xn%*%beta + Ani%*%Bni%*%En
      # 估计方程赋值
      z = matrix(NA,nrow=n,ncol=p+3)
      z[,1:p] = b*e
      z[,p+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
      z[,p+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
      z[,p+3] = e*e - rep(sigma2, n)
      
      # 计算EL值
      lam = lambdaChen(z)
      el = 2*sum(log(1+t(lam)%*%t(z)))
      if(el<=qchisq(0.95,p+3)) f1 = f1+1
      
      # 计算Sigma
      Sigma = matrix(0,p+3,p+3)
      # 赋值上三角形
      Sigma[1:p,p+1] = sigma2*b%*%s+mu3*b%*%diag(Gnn)
      Sigma[1:p,p+2] = mu3*b%*%diag(Hnn)
      Sigma[1:p,p+3] = mu3*b%*%rep(1,n)
      Sigma[p+1,p+2] = 2*sigma4*tr(Gnn%*%Hnn)+(mu4-3*sigma4)*t(diag(Gnn))%*%diag(Hnn)
                       +mu3*t(s)%*%diag(Hnn)
      Sigma[p+1,p+3] = (mu4-sigma4)*tr(Gnn)+mu3*t(s)%*%rep(1,n)
      Sigma[p+2,p+3] = (mu4-sigma4)*tr(Hnn)
      # 赋值下三角形
      Sigma = Sigma + t(Sigma)
      # 赋值对角线
      Sigma[1:p,1:p] = sigma2*b%*%t(b)
      Sigma[p+1,p+1] = 2*sigma4*tr(Gnn%*%Gnn)+sigma2*t(s)%*%s
                       +(mu4-3*sigma4)*L2sqare(Gnn)+2*mu3*t(s)%*%diag(Gnn)
      Sigma[p+2,p+2] = 2*sigma4*tr(Hnn%*%Hnn)*(mu4-3*sigma4)*L2sqare(Hnn)
      Sigma[p+3,p+3] = n*(mu4-sigma4)
      # 计算Sigma特征向量
      E = eigen(Sigma)$vec
      # 特征向量标准化
      E = t(t(E)/sqrt(colSums(E^2)))
      # 检验E模长为1 print(colSums(E^2))
      
      # 计算PCAEL-1值
      s = 1
      Es = E[,1:s]
      pcaz=z%*%Es
      lam=lambdaChen(pcaz)
      pacel=2*sum(log(1+t(lam)%*%t(pcaz)))
      if(pacel<qchisq(0.95,s)) f2 = f2 + 1
      
      # 计算PCAEL-2值
      s = 2
      Es = E[,1:s]
      pcaz=z%*%Es
      lam=lambdaChen(pcaz)
      pacel=2*sum(log(1+t(lam)%*%t(pcaz)))
      if(pacel<qchisq(0.95,s)) f3 = f3 + 1
      
      # 计算PCAEL-3值
      s = 3
      Es = E[,1:s]
      pcaz=z%*%Es
      lam=lambdaChen(pcaz)
      pacel=2*sum(log(1+t(lam)%*%t(pcaz)))
      if(pacel<qchisq(0.95,s)) f4 = f4 + 1

      # 计算PCAEL-4值
      s = 4
      Es = E[,1:s]
      pcaz=z%*%Es
      lam=lambdaChen(pcaz)
      pacel=2*sum(log(1+t(lam)%*%t(pcaz)))
      if(pacel<qchisq(0.95,s)) f5 = f5 + 1
      
    }
    cat(n,p,f1/nsim,f2/nsim,f3/nsim,f4/nsim,f5/nsim,'\n')
    
    # 储存数据
    EL[j,]  =  c(p,f1/nsim,f2/nsim,f3/nsim,f4/nsim,f5/nsim)
  }
  # 写入到本地
  # write.csv( EL, file = paste0('PCA-' ,error,'-',n,'.csv'), row.names = FALSE)
}
#===============================================
#                   赋值运算            
#=============================================== 
errors = 1:4
ms = c(10,15,20,25,30)
for(error in errors){
  for(m in ms){
    indexs = c(0,0.1,0.2,0.3,0.4,0.5)
    SARAR_Data(m,error,indexs,nsim=200)
    cat('\n')
  }
}

