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
rou1 = 0.38
rou2 = 0.18
sigma2 = 1
theta = c(beta,rou1,rou2,sigma2)
mu3 = 0
mu4 = 3
overtheta = c(theta,mu3,mu4)
In = diag(n)
Yn = columbus$CRIME
Xn = Matrix(c(rep(1,n),columbus$HOVAL,columbus$INC),n,p)
#===============================================
#                    求解AEL
#===============================================
myPCEL<- function(overtheta){
  # 求解参数
  beta = overtheta[1:p]
  rou1 = overtheta[p+1]
  rou2 = overtheta[p+2]
  sigma2 = overtheta[p+3]
  mu3 = overtheta[p+4]
  mu4 = overtheta[p+5]
  sigma4 = sigma2^2
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
  
  # 计算Sigma
  Sigma = matrix(1:36,p+3,p+3)
  # 赋值上三角形
  Sigma[1:p,p+1] = as.vector(sigma2*b%*%s+mu3*b%*%diag(Gnn))
  Sigma[1:p,p+2] = as.vector(mu3*b%*%diag(Hnn))
  Sigma[1:p,p+3] = as.vector(mu3*b%*%rep(1,n))
  Sigma[p+1,p+2] = as.numeric(2*sigma4*tr(Gnn%*%Hnn)+(mu4-3*sigma4)*t(diag(Gnn))%*%diag(Hnn)
                   +mu3*t(s)%*%diag(Hnn))
  Sigma[p+1,p+3] = as.numeric((mu4-sigma4)*tr(Gnn)+mu3*t(s)%*%rep(1,n))
  Sigma[p+2,p+3] = as.numeric((mu4-sigma4)*tr(Hnn))
  # 赋值下三角形
  Sigma = Sigma + t(Sigma)
  # 赋值对角线
  Sigma[1:p,1:p] = as.matrix(sigma2*b%*%t(b))
  Sigma[p+1,p+1] = as.numeric(2*sigma4*tr(Gnn%*%Gnn)+sigma2*t(s)%*%s
                   +(mu4-3*sigma4)*L2sqare(Gnn)+2*mu3*t(s)%*%diag(Gnn))
  Sigma[p+2,p+2] = as.numeric(2*sigma4*tr(Hnn%*%Hnn)+(mu4-3*sigma4)*L2sqare(Hnn))
  Sigma[p+3,p+3] = n*(mu4-sigma4)
  # 计算Sigma特征向量
  E = eigen(Sigma)$vec
  # 特征向量标准化
  E = t(t(E)/sqrt(colSums(E^2)))
  # 检验E模长为1 print(colSums(E^2))

  # 计算PCEL-s值
  s = 1
  Es = E[,1:s]
  pczs=z%*%Es
  pclams=lambdaChen(pczs)
  pcels=2*sum(log(1+t(pclams)%*%t(pczs)))
  
  return(pcels)
}
s = 1
cat('s：',s,'\n')
cat('参数初值：',overtheta,'\n')
PCELresult = nlminb(overtheta, myPCEL, lower = c(rep(-Inf,p),-0.99,-0.99,0.01,rep(-Inf,2)), upper = c(c(Inf,0,0),0.99,0.99,Inf,rep(Inf,2)))
par = PCELresult$par
cat('参数估计值：',par,'\n')
# cat('pcels<qchisq(0.95,s)：',PCELresult$objective<qchisq(0.95,s),'\n')
cat('p值：', 1 - pchisq(myPCEL(par), df = s), '\n')
mae = sum(abs(Yn-Yn_hat(par[1:6],p)))/49
cat('mae：',mae,'\n')

# 检验参数
# s = 1
# par = c(48.99963, -0.3129058, -1.004385, 0.3678953, 0.200798, 0.9999725, -4.105493e-05, 3.000015) # 1.634124e-14
# par = c(48.9996, -0.3129, -1.0044, 0.3679, 0.2008, 1.0000, -0.0000, 3.0000) # 1.477333e-14
# s = 2
# par = c(48.99966, -0.3117704, -1.004343, 0.3666539, 0.179707, 1.000014, 4.310436e-05, 3.000038) # 1.479939e-14 
# par = c(49.0000, -0.3118, -1.0043, 0.3667, 0.1797, 1.0000, 0.0000, 3.0000) # 1.273414e-14
# cat('pcels<qchisq(0.95,s):',myPCEL(par)<qchisq(0.95,s),'myPCEL(par):',myPCEL(par),'\n')
# cat('p值：', 1 - pchisq(myPCEL(par), df = s), '\n')

# rmse = sqrt(mean((Yn-Yn_hat(par,p))^2))
# cat('rmse：',rmse,'\n')




