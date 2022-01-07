library(spatstat)
library(RandomFields)
#library(Rcplex)
#library(MASS)
#library(slam)
#library(Matrix)
library(matlab)
#library(LICORS)
#library(VGAM)
library(hypergeo)

Plos = function(htx,hrx,r,alpha,beta,gamma){
  n = floor((r/1000)*sqrt(alpha*beta))
  a = 0:max(0,n-1)
  P = prod(1-exp(-((htx-(a+1/2)*(htx-hrx)/(n))^2)/(2*gamma^2)))
  if(n==0 && is.na(P))
  {P = 1}
  return(P)
}

BSvoidProb = function(hu,hm,r1,angle,alpha,beta,gamma,density){
  if(hu<hm){
  if(r1>=(hm-hu)/cos(angle/2)){r1=(hm-hu)/cos(angle/2)}
  h = seq(from=hu,to=cos(angle/2)*(r1)+hu,by=(cos(angle/2)*r1)/100)
  cov1 = vector(length=length(h))
  for (j in 1:length(h)){
    max = floor((tan(angle/2)*(h[j]-hu)/1000)*sqrt(beta*alpha))
    sum = 0
    for (k in 0:max){
      lo = max(0,(k)*1000/sqrt(beta*alpha))
      up = min((tan(angle/2)*(h[j]-hu)),(k+1)*1000/sqrt(beta*alpha))
      P=Plos(htx=h[j],hrx=hu,r=(up+lo)/2,alpha=alpha,beta=beta,gamma=gamma)
      sum = sum + P*(up^2-lo^2)/2
    }
    cov1[j] = sum*(h[2]-h[1])
  }
  }
  else{cov1=0}
  
  if((cos(angle/2)*(r1)+hu)<hm){
  h = seq(from=cos(angle/2)*(r1)+hu,to=min(r1+hu,hm),by=(min(r1+hu,hm)-(cos(angle/2)*(r1)+hu))/100)
  cov2 = vector(length=length(h))
  for (j in 1:length(h)){
    max = floor((sqrt(max(0,(r1)^2-(h[j]-hu)^2))/1000)*sqrt(beta*alpha))
    sum = 0
    for (k in 0:max){
      lo = max(0,(k)*1000/sqrt(beta*alpha))
      up = min((sqrt(max(0,(r1)^2-(h[j]-hu)^2))),(k+1)*1000/sqrt(beta*alpha))
      P=Plos(htx=h[j],hrx=hu,r=(up+lo)/2,alpha=alpha,beta=beta,gamma=gamma)
      sum = sum + P*(up^2-lo^2)/2
    }
    cov2[j] = sum*(h[2]-h[1])
  }
  }
  else{cov2 = 0}
  
  return(exp(-2*pi*density*(sum(cov1[!is.na(cov1)])+sum(cov2[!is.na(cov2)]))))
}

BSdistPDF = function(hu,hm,r1,angle,alpha,beta,gamma,density){
  if((cos(angle/2)*(r1)+hu)<hm){
  h = seq(from=cos(angle/2)*(r1)+hu,to=min(r1+hu,hm),by=(min(r1+hu,hm)-(cos(angle/2)*(r1)+hu))/500)
  P = 0
  for (j in 1:length(h)){
    P = P + Plos(htx=h[j],hrx=hu,r=(sqrt(max(0,(r1)^2-(h[j]-hu)^2))),alpha=alpha,beta=beta,gamma=gamma)*(h[2]-h[1])/(min(r1+hu,hm)-(cos(angle/2)*(r1)+hu))#  }
  }
  
  # pdf[i,l] = 2*pi*density*(r1[l]^2)*(1-cos(angle[i]/2))*exp(-2*pi*density*(r1[l]^3)*(1-cos(angle[i]/2))/3)*(r1[2]-r1[1])
  return(2*pi*density*P*(r1)*(min(r1+hu,hm)-(cos(angle/2)*(r1)+hu))*BSvoidProb(hu,hm,r1,angle,alpha,beta,gamma,density))
  }
  else{return(0)}
}

BSLaplaceI = function(hu,hm,r1,s,angle,uBW,alpha,beta,gamma,density,al){
  if((cos(angle/2)*(r1)+hu)<hm){
  h = seq(from=cos(angle/2)*r1+hu,to=min((r1+hu),hm),by=(min((r1+hu),hm)-(cos(angle/2)*r1+hu))/100)
 # h = seq(from=sqrt(r1^2/(1+tan(angle/2)^2))+hu,to=min((r1+hu),hm),by=(min((r1+hu),hm)-(sqrt(r1^2/(1+tan(angle/2)^2))+hu))/100)
  cov1 = vector(length=length(h))
  for (j in 1:length(h)){
      max = floor((tan(angle/2)*(h[j]-hu)/1000)*sqrt(beta*alpha))
      min = floor((sqrt(max(0,(r1)^2-(h[j]-hu)^2))/1000)*sqrt(beta*alpha))
      sum = 0
      for (k in min:max){
        lo = max(sqrt(max(0,(r1)^2-(h[j]-hu)^2)),(k)*1000/sqrt(beta*alpha))
        up = min((tan(angle/2)*(h[j]-hu)),(k+1)*1000/sqrt(beta*alpha))
       # P=Plos(htx=h[j],hrx=hu,r=(up+lo)/2,alpha=alpha,beta=beta,gamma=gamma)
         range = 0:max(k-1,0)
         P = prod(1-exp(-((h[j]-(range+1/2)*(h[j]-hu)/(k))^2)/(2*gamma^2)))
   
        if(al==2){
          sum = sum + P*(s*log(1+(up^2+(h[j]-hu)^2)/s)-s*log(1+(lo^2+(h[j]-hu)^2)/s))
        }
        else{
        sum = sum + P*((up^2+(h[j]-hu)^2)*hypergeo(1,2/al,1+2/al,-((up^2+(h[j]-hu)^2)^(al/2))/s)) 
        sum = sum - P*(lo^2+(h[j]-hu)^2)*hypergeo(1,2/al,1+2/al,-((lo^2+(h[j]-hu)^2)^(al/2))/s)
        }
      }
      cov1[j] = sum*(h[2]-h[1])
  }
  }
  else{cov1=0}
  
    if((r1+hu)<hm){
    h = seq(from=(r1+hu),to=hm,by=(hm-(r1+hu))/100)
    cov2 = vector(length=length(h))
    for (j in 1:length(h)){
      max = floor((tan(angle/2)*(h[j]-hu)/1000)*sqrt(beta*alpha))
      sum = 0
      for (k in 0:max){
        lo = max(0,(k)*1000/sqrt(beta*alpha))
        up = min((tan(angle/2)*(h[j]-hu)),(k+1)*1000/sqrt(beta*alpha))
      #  P=Plos(htx=h[j],hrx=hu,r=(up+lo)/2,alpha=alpha,beta=beta,gamma=gamma)
         range = 0:max(k-1,0)
         P = prod(1-exp(-((h[j]-(range+1/2)*(h[j]-hu)/(k))^2)/(2*gamma^2)))
      
        if(al==2){
          sum = sum + P*(s*log(1+(up^2+(h[j]-hu)^2)/s)-s*log(1+(lo^2+(h[j]-hu)^2)/s))
        }
        else{
          sum = sum + P*((up^2+(h[j]-hu)^2)*hypergeo(1,2/al,1+2/al,-((up^2+(h[j]-hu)^2)^(al/2))/s)) 
          sum = sum - P*(lo^2+(h[j]-hu)^2)*hypergeo(1,2/al,1+2/al,-((lo^2+(h[j]-hu)^2)^(al/2))/s)
        }
      }
      cov2[j] = sum*(h[2]-h[1])
    }
    }
  else{cov2 = 0}
  
  
  return(exp(-(uBW/2)*density*(sum(cov1[!is.na(cov1)])+sum(cov2[!is.na(cov2)]))))
}
