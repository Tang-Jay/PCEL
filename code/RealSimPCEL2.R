#===============================================
#             SARAR model Columbus
#===============================================
rm(list = ls())
library(spdep)
library(emplik)
library(expm)
library(maxLik)
library('sp')
library('terra')
library('spdep')
source('GlambdaChen.R')
data("columbus")
f<-function(Matrix,Vector){
  irow = nrow(Matrix)
  v =c(0)
  for(i in 2:irow){
    v[i] = Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
  }
  return(v)
}
tr<-function(Matrix){
  return(sum(diag(Matrix)))
}
L2sqare<-function(Matrix){
  return(sum(diag(Matrix)^2))
}
Yn_hat<- function(theta,p){
  beta = theta[1:p]
  rou1 = theta[p+1]
  rou2 = theta[p+2]
  An = In - rou1*Wn
  Bn = In - rou2*Mn
  Ani = solve(An)
  Bni = solve(Bn)
  En = Bn%*%(An%*%Yn-Xn%*%beta)
  Yn = Ani%*%Xn%*%beta + Ani%*%Bni%*%En
  return(Yn)
}
#===============================================
#                    修改模拟参数
#===============================================
n = 49
p = 3
Ws = nb2listw(col.gal.nb,style = 'W')
Wn = listw2mat(Ws) 
Mn = Wn
beta = c(49,-0.3,-1)
rou1 = 0.35
rou2 = 0.1
sigma2 = 1
theta = c(beta,rou1,rou2,sigma2)
In = diag(n)
Yn = columbus$CRIME
Xn = Matrix(c(rep(1,n),columbus$HOVAL,columbus$INC),n,p)
#===============================================
#                    求解PCEL
#===============================================
myPCEL <- function(theta){
  # 求解参数
  beta = theta[1:p]
  rou1 = theta[p+1]
  rou2 = theta[p+2]
  # 估计方程准备
  An = In - rou1*Wn
  Bn = In - rou2*Mn
  Ani = solve(An)
  Bni = solve(Bn)
  En = Bn%*%(An%*%Yn-Xn%*%beta)
  Gn = Bn%*%Wn%*%Ani%*%Bni
  Gnn = 1/2*(Gn + t(Gn))
  Hn = Mn%*%Bni
  Hnn = 1/2*(Hn + t(Hn))
  # 简化符号
  b = t(Xn)%*%t(Bn)
  g = diag(Gnn)
  h = diag(Hnn)
  s = Bn%*%Wn%*%Ani%*%Xn%*%beta
  e = En
  # 估计方程赋值
  z1 = b[1,]*e
  z2 = b[2,]*e
  z3 = b[3,]*e
  z4 = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
  z5 = h*(e^2-sigma2) + 2*e*f(Hnn,e)
  z6 = e*e - rep(sigma2, n)
  z = cbind(z1,z2,z3,z4,z5,z6)
  # 估计Sigma
  Sigma = t(z) %*% z
   
  # 计算Sigma特征向量
  E = eigen(Sigma)$vec
  evals <- eigen(Sigma)$values
  evals_prop <- evals / sum(evals)
  # # fig6
  # # setEPS()
  # # s = 1
  # # postscript('fig6-a.eps', width=4.5, height=3.1, paper="special", horizontal=FALSE)
  # # s = 2
  # # postscript('fig6-b.eps', width=4.5, height=3.1, paper="special", horizontal=FALSE)
  # par(mar=c(4.5,4.5,0.5,0.5), oma=c(0,0,0,0))
  # plot(1:length(evals_prop), evals_prop, type = "b",
  #      xlab = "Index of eigenvalue", ylab = "Normalized eigenvalue",
  #      ylim = c(0, 1), las = 1)
  # axis(2, las = 2)
  # text(1+0.05, evals_prop[1], labels = sprintf("%.2f", evals_prop[1]), pos = 3, cex = 0.8)
  # text(2+0.15, evals_prop[2], labels = sprintf("%.2f", evals_prop[2]), pos = 3, cex = 0.8)
  # # dev.off()
  
  # 特征向量标准化
  E = t(t(E)/sqrt(colSums(E^2)))
  # 检验E模长为1 print(colSums(E^2))

  # 计算PCEL-s值 
  s = 2 # 保持与 s 赋值相同
  Es = E[,1:s]
  pczs=z%*%Es
  pclams=lambdaChen(pczs)
  pcels=2*sum(log(1+t(pclams)%*%t(pczs)))
  
  return(pcels)
}
s <- 2
cat('s：',s,'\n')
cat('参数初值：',theta,'\n')
PCELresult = nlminb(theta, myPCEL,
                    lower = c(rep(-Inf,p),-0.99,-0.99,0.01), 
                    upper = c(c(Inf,0,0),0.99,0.99,Inf))
par = PCELresult$par
cat('参数估计值：',par,'\n')
cat('pcels<qchisq(0.95,s)：',PCELresult$objective<qchisq(0.95,s),PCELresult$objective,'\n')
cat('p值：', 1 - pchisq(myPCEL(par), df = s), '\n')
mae = sum(abs(Yn-Yn_hat(par,p)))/49
cat('mae：',mae,'\n')

# 检验参数
# s = 1
# par = c(48.99944, -0.3217947, -1.009112, 0.3667696, 0.2104571, 1.000082) # 1.430998e-14
# par = c(48.9994, -0.3218, -1.0091, 0.3668, 0.2105, 1.0001) # 1.347731e-14
# s = 2
par = c(49.1990, -0.3365, -1.1237, 0.4400, -0.4888, 0.9592) # 1.141252e-14
par = c(49.1990, -0.3365, -1.1237, 0.4400,  0.4888, 0.9592) # 1.284743e-14
par = c(49.1990, -0.3365, -1.1237, 0.3568,  0.1468, 0.9592) # 1.284743e-14
cat('pcels<qchisq(0.95,s)：',myPCEL(par)<qchisq(0.95,s),myPCEL(par),'\n')
# cat('p值：', 1 - pchisq(myPCEL(par), df = s), '\n')