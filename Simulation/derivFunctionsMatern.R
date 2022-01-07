##new version of derivFunction that allows for both transmitter and receiver non-zero heights

library(spatstat)
library(RandomFields)
#library(Rcplex)
library(MASS)
#library(slam)
library(Matrix)
library(matlab)
#library(LICORS)
#library(VGAM)
library(hypergeo)

source("derivFunctionsv5.R")

##like the strongestLaplaceLOSSINR, but with a cluster interferer
strongestLaplaceClusterInt = function(max,min,type,N,al,an,alpha,beta,gamma,rxheight,txheight,mal,man,l,const,s,clusterType,clusterRad,clusterTypeProb){
  ints = vector(mode="integer",length=4)
  answer = 0
  loop = TRUE
  while(loop==TRUE){
    
    if(sum(ints)==l){
      
      foo = (factorial(l)/(prod(factorial(ints))))
      prod=1
      
      rn = sqrt(max(0,(min^2+(rxheight-txheight)^2)^(al/an) - (rxheight-txheight)^2))
      rl = min(max,sqrt((min^2+(rxheight-txheight)^2)^(an/al) - (rxheight-txheight)^2))
      
      if(type==al){
      prod = prod*deriveLaplaceLOS(max=max,min=min,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
      prod = prod*deriveLaplaceNLOS(max=max,min=rn,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)
      if(clusterType==al){
      if(min>=min(max,clusterRad)){print("Error")}
      prod = prod*ClusterLaplaceLOS(max=min(max,clusterRad),min=min,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[3],b=2/al,c=1+2/al,s=s,clusterRad,clusterTypeProb)
      }
      else{
      if(rn>=min(max,clusterRad)){print("Error")}  
      prod = prod*ClusterLaplaceNLOS(max=min(max,clusterRad),min=rn,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[3],b=2/an,c=1+2/an,s=s,clusterRad,clusterTypeProb)  
      }
      }
      else{
      prod = prod*deriveLaplaceLOS(max=max,min=rl,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
      prod = prod*deriveLaplaceNLOS(max=max,min=min,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)
      if(clusterType==al){
        if(rl>=min(max,clusterRad)){print("Error")}
        prod = prod*ClusterLaplaceLOS(max=min(max,clusterRad),min=rl,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[3],b=2/al,c=1+2/al,s=s,clusterRad,clusterTypeProb)
      }
      else{
        if(min>=min(max,clusterRad)){print("Error")}
        prod = prod*ClusterLaplaceNLOS(max=min(max,clusterRad),min=min,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[3],b=2/an,c=1+2/an,s=s,clusterRad,clusterTypeProb)  
      }
      }
      prod = prod*((-1)^ints[4])*((N)^(ints[4]))*exp(-N*s)
      answer = answer+foo*prod
    }
    
    
    ints[1]=ints[1]+1
    for(i in 1:3){
      if(ints[i]>l){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    if(ints[4]>l){
      loop = FALSE
    }  
  }
  return(answer)
}


strongestLaplaceClusterIntv2 = function(max,min,type,N,al,an,alpha,beta,gamma,rxheight,txheight,mal,man,l,const,s,clusterType,clusterRad,clusterTypeProb){
  ints = vector(mode="integer",length=5)
  answer = 0
  loop = TRUE
  while(loop==TRUE){
    
    if(sum(ints)==l){
      
      foo = (factorial(l)/(prod(factorial(ints))))
      prod=1
      
      rn = sqrt(max(0,(min^2+(rxheight-txheight)^2)^(al/an) - (rxheight-txheight)^2))
      rl = min(max,sqrt((min^2+(rxheight-txheight)^2)^(an/al) - (rxheight-txheight)^2))
      
      if(type==al){
        prod = prod*deriveLaplaceLOS(max=max,min=min,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
        prod = prod*deriveLaplaceNLOS(max=max,min=rn,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)
          
        if(min<min(max,clusterRad)){
          prod = prod*ClusterLaplaceLOS(max=min(max,clusterRad),min=min,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[3],b=2/al,c=1+2/al,s=s,clusterRad,clusterTypeProb)
        }
        if(rn<min(max,clusterRad)){
     #     prod = prod*ClusterLaplaceNLOS(max=min(max,clusterRad),min=rn,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[4],b=2/an,c=1+2/an,s=s,clusterRad,clusterTypeProb)  
        }
      }
      else{
        prod = prod*deriveLaplaceLOS(max=max,min=rl,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
        prod = prod*deriveLaplaceNLOS(max=max,min=min,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)
        
          if(rl<min(max,clusterRad)){
          prod = prod*ClusterLaplaceLOS(max=min(max,clusterRad),min=rl,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[3],b=2/al,c=1+2/al,s=s,clusterRad,clusterTypeProb)
          }
         if(min<min(max,clusterRad)){
        #  prod = prod*ClusterLaplaceNLOS(max=min(max,clusterRad),min=min,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[4],b=2/an,c=1+2/an,s=s,clusterRad,clusterTypeProb)
         }
      }
      prod = prod*((-1)^ints[5])*((N)^(ints[5]))*exp(-N*s)
      answer = answer+foo*prod
    }
    
    
    ints[1]=ints[1]+1
    for(i in 1:4){
      if(ints[i]>l){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    if(ints[5]>l){
      loop = FALSE
    }  
  }
  return(answer)
}

strongestLaplaceClusterIntv3 = function(max,min,type,N,al,an,alpha,beta,gamma,rxheight,txheight,mal,man,l,const,s,clusterType,clusterRad,rangeProb,clusterLOSProb,clusterNLOSProb){
  ints = vector(mode="integer",length=4)
  answer = 0
  loop = TRUE
  while(loop==TRUE){
    
    if(sum(ints)==l){
      
      foo = (factorial(l)/(prod(factorial(ints))))
      prod=1
      
      rn = sqrt(max(0,(min^2+(rxheight-txheight)^2)^(al/an) - (rxheight-txheight)^2))
      rl = min(max,sqrt((min^2+(rxheight-txheight)^2)^(an/al) - (rxheight-txheight)^2))
      
      if(type==al){
        prod = prod*deriveLaplaceLOS(max=max,min=min,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
        prod = prod*deriveLaplaceNLOS(max=max,min=rn,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)
        prod = prod*ClusterLaplace(max=min(max,clusterRad),lmin=min,nmin=rn,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,mal=mal,man=man,l=ints[3],s=s,clusterRad=clusterRad,rangeProb=rangeProb,clusterLOSProb=clusterLOSProb,clusterNLOSProb=clusterNLOSProb)
      }
      else{
        prod = prod*deriveLaplaceLOS(max=max,min=rl,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
        prod = prod*deriveLaplaceNLOS(max=max,min=min,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)
        prod = prod*ClusterLaplace(max=min(max,clusterRad),lmin=rl,nmin=min,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,mal=mal,man=man,l=ints[3],s=s,clusterRad=clusterRad,rangeProb=rangeProb,clusterLOSProb=clusterLOSProb,clusterNLOSProb=clusterNLOSProb)
      }
      prod = prod*((-1)^ints[4])*((N)^(ints[4]))*exp(-N*s)
      answer = answer+foo*prod
    }
    
    
    ints[1]=ints[1]+1
    for(i in 1:3){
      if(ints[i]>l){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    if(ints[4]>l){
      loop = FALSE
    }  
  }
  return(answer)
}

strongestLaplaceClusterIntv4 = function(max,min,type,N,al,an,alpha,beta,gamma,rxheight,txheight,mal,man,l,const,s,clusterType,clusterRad,rangeProb,clusterLOSProb,clusterNLOSProb,isClusterProb){
  ints = vector(mode="integer",length=4)
  answer = 0
  loop = TRUE
  while(loop==TRUE){
    
    if(sum(ints)==l){
      
      foo = (factorial(l)/(prod(factorial(ints))))
      prod=1
      
      rn = sqrt(max(0,(min^2+(rxheight-txheight)^2)^(al/an) - (rxheight-txheight)^2))
      rl = min(max,sqrt((min^2+(rxheight-txheight)^2)^(an/al) - (rxheight-txheight)^2))
      
      if(type==al){
        prod = prod*deriveLaplaceLOS(max=max,min=min,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
        prod = prod*deriveLaplaceNLOS(max=max,min=rn,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)
        if(ints[3]==0){di=1}
        else{di=0}
        prod = prod*(isClusterProb*di + (1-isClusterProb)*ClusterLaplace(max=max,lmin=min,nmin=rn,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,mal=mal,man=man,l=ints[3],s=s,clusterRad=clusterRad,rangeProb=rangeProb,clusterLOSProb=clusterLOSProb,clusterNLOSProb=clusterNLOSProb))
      }
      else{
        prod = prod*deriveLaplaceLOS(max=max,min=rl,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
        prod = prod*deriveLaplaceNLOS(max=max,min=min,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)
        if(ints[3]==0){di=1}
        else{di=0}
        prod = prod*(isClusterProb*di+(1-isClusterProb)*ClusterLaplace(max=max,lmin=rl,nmin=min,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,mal=mal,man=man,l=ints[3],s=s,clusterRad=clusterRad,rangeProb=rangeProb,clusterLOSProb=clusterLOSProb,clusterNLOSProb=clusterNLOSProb))
      }
      prod = prod*((-1)^ints[4])*((N)^(ints[4]))*exp(-N*s)
      answer = answer+foo*prod
    }
    
    
    ints[1]=ints[1]+1
    for(i in 1:3){
      if(ints[i]>l){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    if(ints[4]>l){
      loop = FALSE
    }  
  }
  return(answer)
}

ClusterLaplace = function(max,lmin,nmin,al,an,alpha,beta,gamma,rxheight,txheight,mal,man,l,s,clusterRad,rangeProb,clusterLOSProb,clusterNLOSProb){
  # ints = vector(mode="integer",length=l)
  # answer = 0
  # loop = TRUE
  
  lfoo = 0
  nfoo = 0
  nofoo=0
  
  Prange = rangeProb
  
  if((clusterRad-lmin)>0.001){
  if(min(max,clusterRad)>lmin){
  for (k in 0:(mal)){
    lfoo = lfoo+Cluster2F1LOS(max=min(max,clusterRad),min=lmin,m=mal,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=l,const=choose(mal,k),b=2/al,c=1+2/al,s=s)
  }
  }
    if(l==0 && (clusterRad>max)){
      
      foo = seq(from=max(max,lmin),to=clusterRad,by=(clusterRad-max(max,lmin))/200)
      Plos = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Plos[i] = prod(1-exp(-((max(rxheight,txheight)-(a+1/2)*(abs(rxheight-txheight))/(n))^2)/(2*gamma^2)))
      }
      
      lfoo = lfoo+ 2*sum((Plos)*foo*(foo[2]-foo[1]))
      
    }
 #   lfoo = lfoo/((max^2-lmin^2)*min(1,(clusterLOSProb+clusterNLOSProb)))
    lfoo = lfoo/((clusterRad^2)*min(1,(clusterLOSProb+clusterNLOSProb)))
  }
  if((clusterRad-nmin)>0.001){
  if(min(max,clusterRad)>nmin){
  for (k in 0:(man)){
    nfoo = nfoo+Cluster2F1NLOS(max=min(max,clusterRad),min=nmin,m=man,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=l,const=choose(man,k),b=2/an,c=1+2/an,s=s)
  }
  }
    
    if(l==0 && (clusterRad>max)){
      
      foo = seq(from=max(max,nmin),to=clusterRad,by=(clusterRad-max(max,nmin))/200)
      Plos = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Plos[i] = prod(1-exp(-((max(rxheight,txheight)-(a+1/2)*(abs(rxheight-txheight))/(n))^2)/(2*gamma^2)))
      }
      
      nfoo = nfoo+ 2*sum((1-Plos)*foo*(foo[2]-foo[1]))
      
    }
    
    
  #  nfoo = nfoo/((max^2-nmin^2)*min(1,(clusterLOSProb+clusterNLOSProb)))
 #   nfoo = nfoo/((clusterRad^2)*min(1,(clusterLOSProb+clusterNLOSProb)))
    nfoo = nfoo/((clusterRad^2)*Prange)
  }
  
 # if((lmin-nmin)>0.001){
#    for (k in 0:(man)){
#      nofoo = nofoo+Cluster2F1NLOS(max=lmin,min=nmin,m=man,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=l,const=choose(man,k),b=2/an,c=1+2/an,s=s)
#    }
#    nofoo = nofoo/(max^2-nmin^2)
#  }
#  answer = Prange*(lfoo+nfoo+nofoo)
  answer = (lfoo+nfoo+nofoo)
  #  return(answer*(1/((clusterRad^2))))
  #   return(answer*(1/((clusterRad^2-min^2))))
  #  return(answer*(1/((max^2-min^2))))
  if(l==0){
#  return(answer+(1-Prange))
    return(answer)
  }
  else{
    return(answer)
  }
  #  return(answer*(1/((clusterRad^2)*clusterTypeProb)))
  #  return(answer*(1/((clusterRad^2-min^2)*clusterTypeProb)))
}


ClusterLaplaceLOS = function(max,min,al,alpha,beta,gamma,rxheight,txheight,m,l,b,c,s,clusterRad,clusterTypeProb){
 # ints = vector(mode="integer",length=l)
 # answer = 0
 # loop = TRUE
 
    foo = 0
    for (k in 0:(m)){
      foo = foo+Cluster2F1LOS(max=max,min=min,m=m,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=l,const=choose(m,k),b=b,c=c,s=s)
    }
    answer = foo
#  return(answer*(1/((clusterRad^2))))
#   return(answer*(1/((clusterRad^2-min^2))))
#  return(answer*(1/((max^2-min^2))))
  return(answer*(1/((max^2-min^2)*clusterTypeProb)))
#  return(answer*(1/((clusterRad^2)*clusterTypeProb)))
#  return(answer*(1/((clusterRad^2-min^2)*clusterTypeProb)))
}

Cluster2F1LOS = function(max,min,m,al,alpha,beta,gamma,rxheight,txheight,k,l,const,b,c,s){
  n1 = floor((min/1000)*sqrt(beta*alpha))
  nmax = floor((max/1000)*sqrt(beta*alpha))
  
  sum = 0
  for (i in n1:nmax){
    lo = max(min,(i)*1000/sqrt(beta*alpha))
    up = min(max,(i+1)*1000/sqrt(beta*alpha))  
    
    range = 0:max(i-1,0)
    Plos = prod(1-exp(-((max(rxheight,txheight)-(range+1/2)*abs(rxheight-txheight)/(i))^2)/(2*gamma^2)))
    if(is.na(Plos)&&i==0){Plos=1}
    
    lo = lo^2+(rxheight-txheight)^2
    up = up^2+(rxheight-txheight)^2
    
    if(l==0){
      if(k==0){
        sum = sum+Plos*(up - lo)
      }
      else{
      sum = sum + Plos*const*up*(-1)^(k)*hypergeo(k,b,c,-m*up^(al/2)/s)
      sum = sum - Plos*const*lo*(-1)^(k)*hypergeo(k,b,c,-m*lo^(al/2)/s)
      }
    }
    else{
      if(k>0){
      sum = sum + derive2F1s(const=Plos*const*up*(-1)^(k),l=l,a=k,b=b,c=c,m=m*up^(al/2),s=s)
      sum = sum - derive2F1s(const=Plos*const*lo*(-1)^(k),l=l,a=k,b=b,c=c,m=m*lo^(al/2),s=s)
      }
      else{
        sum=sum+0
      }
    }
  }
  
  return(sum)
}

ClusterLaplaceNLOS = function(max,min,an,alpha,beta,gamma,rxheight,txheight,m,l,b,c,s,clusterRad,clusterTypeProb){
  ints = vector(mode="integer",length=l)
  answer = 0
  loop = TRUE

    foo = 0
    for (k in 0:(m)){
      foo = foo+Cluster2F1NLOS(max=max,min=min,m=m,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=l,const=choose(m,k),b=b,c=c,s=s)
    }
    answer = foo
  
#  return(answer*(1/((clusterRad^2))))
#  return(answer*(1/((clusterRad^2-min^2))))
 # return(answer*(1/((max^2-min^2))))
  return(answer*(1/(((max^2-min^2))*clusterTypeProb)))
#  return(answer*(1/((clusterRad^2)*clusterTypeProb)))
#  return(answer*(1/((clusterRad^2-min^2)*clusterTypeProb)))
}

Cluster2F1NLOS = function(max,min,m,an,alpha,beta,gamma,rxheight,txheight,k,l,const,b,c,s){
  n1 = floor((min/1000)*sqrt(beta*alpha))
  nmax = floor((max/1000)*sqrt(beta*alpha))
  
  sum = 0
  for (i in n1:nmax){
    lo = max(min,(i)*1000/sqrt(beta*alpha))
    up = min(max,(i+1)*1000/sqrt(beta*alpha))  
    
    range = 0:max(i-1,0)
    PNlos = 1-prod(1-exp(-((max(rxheight,txheight)-(range+1/2)*abs(rxheight-txheight)/(i))^2)/(2*gamma^2)))
    if(is.na(PNlos)&&i==0){PNlos=0}
    
    lo = lo^2+(rxheight-txheight)^2
    up = up^2+(rxheight-txheight)^2
    
    if(l==0){
      if(k==0){
        sum = sum+PNlos*(up - lo)
      }
      else{
      sum = sum + PNlos*const*up*(-1)^(k)*hypergeo(k,b,c,-m*up^(an/2)/s)
      sum = sum - PNlos*const*lo*(-1)^(k)*hypergeo(k,b,c,-m*lo^(an/2)/s)
      }
    }
    else{
      if(k>0){
      sum = sum + derive2F1s(const=PNlos*const*up*(-1)^(k),l=l,a=k,b=b,c=c,m=m*up^(an/2),s=s)
      sum = sum - derive2F1s(const=PNlos*const*lo*(-1)^(k),l=l,a=k,b=b,c=c,m=m*lo^(an/2),s=s)
      }
      else{sum=sum+0}
    }
  }
  
  return(sum)
}

Cluster2F1NLOSonly = function(max,min,m,an,alpha,beta,gamma,rxheight,txheight,k,l,const,b,c,s){
  n1 = floor((min/1000)*sqrt(beta*alpha))
  nmax = floor((max/1000)*sqrt(beta*alpha))
  
  sum = 0
  for (i in n1:nmax){
    lo = max(min,(i)*1000/sqrt(beta*alpha))
    up = min(max,(i+1)*1000/sqrt(beta*alpha))  
    
    range = 0:max(i-1,0)
    PNlos = 1#-prod(1-exp(-((max(rxheight,txheight)-(range+1/2)*abs(rxheight-txheight)/(i))^2)/(2*gamma^2)))
    if(is.na(PNlos)&&i==0){PNlos=0}
    
    lo = lo^2+(rxheight-txheight)^2
    up = up^2+(rxheight-txheight)^2
    
    if(l==0){
      if(k==0){
        sum = sum+PNlos*(up - lo)
      }
      else{
        sum = sum + PNlos*const*up*(-1)^(k)*hypergeo(k,b,c,-m*up^(an/2)/s)
        sum = sum - PNlos*const*lo*(-1)^(k)*hypergeo(k,b,c,-m*lo^(an/2)/s)
      }
    }
    else{
      if(k>0){
        sum = sum + derive2F1s(const=PNlos*const*up*(-1)^(k),l=l,a=k,b=b,c=c,m=m*up^(an/2),s=s)
        sum = sum - derive2F1s(const=PNlos*const*lo*(-1)^(k),l=l,a=k,b=b,c=c,m=m*lo^(an/2),s=s)
      }
      else{sum=sum+0}
    }
  }
  
  return(sum)
}
