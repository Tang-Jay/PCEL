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
Yn_hat<- function(theta,k){
  beta = theta[1:k]
  rou1 = theta[k+1]
  rou2 = theta[k+2]
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
k = 3
Ws = nb2listw(col.gal.nb,style = 'W')
Wn = listw2mat(Ws) 
Mn = Wn
beta = c(0,-0.5,-0.5)
rou1 = 0.5
rou2 = 0.5
sigma2 = 1
theta = c(beta,rou1,rou2,sigma2)
In = diag(n)
Yn = columbus$CRIME
Xn = Matrix(c(rep(1,n),columbus$HOVAL,columbus$INC),n,k)
#===============================================
#                    求解EL参数
#===============================================
myEL<- function(theta){
  # 求解参数
  beta = theta[1:k]
  rou1 = theta[k+1]
  rou2 = theta[k+2]
  sigma2 = theta[k+3]

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
  # 计算EL值
  lam = lambdaChen(z)
  el = 2*sum(log(1+t(lam)%*%t(z)))

  return(el)
}
cat('参数初值：',theta,'\n')
result = nlminb(theta,myEL,lower = c(rep(-Inf,k),-0.99,-0.99,0.01),upper = c(c(Inf,0,0),0.99,0.99,Inf))
par = result$par
cat('参数估计值：',par,'\n')
cat('el<qchisq(0.95,s)：',ELresult$objective<qchisq(0.95,k+3),'\n')
mae = sum(abs(Yn-Yn_hat(par,k)))/49
cat('mae：',mae,'\n')

# 检验参数
par = c(-104.24, -0.42, -0.30, -0.55, 0.99, 93.89) # mae: 8.90801e-14
par = c(49.05148, -0.283115, -1.068778, 0.3532607, 0.1319959, 99.42319) # 1.335043e-14 
par = c(49.0519,-0.2832,-1.0688,0.3533,0.1320,99.4232) # 1.413892e-14
par = c(  49.05, -0.28, -1.07,  0.35, 0.13, 99.42) # mae: 1.324054e-14
# par = c(0.048, -2.26, -0.38, 0.99, 0.99, 1.33) # PCEL s=1 参数，显示拒绝
cat('el<qchisq(0.95,k+3):',myEL(par)<qchisq(0.95,k+3),'myEL(par):',myEL(par),'\n')
cat('p值：', 1 - pchisq(myEL(par), df = k+3), '\n')
# rmse = sqrt(mean((Yn-Yn_hat(par,k))^2))
# cat('rmse：',rmse,'\n')




