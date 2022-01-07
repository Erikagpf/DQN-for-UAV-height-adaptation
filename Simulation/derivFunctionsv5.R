##new version of derivFunction that allows for both transmitter and receiver non-zero heights

library(spatstat)
library(RandomFields)
#library(Rcplex)
#library(MASS)
#library(slam)
library(Matrix)
library(matlab)
#library(LICORS)
#library(VGAM)
library(hypergeo)

risingFac = function(x,n){
  foo = choose(x+n-1,n)*factorial(n)
    return(foo)
}

#get the l-th derivative of (-b/s)^n wrt s
deriveZ = function(l,n,b,s){
  result = 0
  for (A in 0:l){
    for (i in 0:A){
      result = result + ((-1)^(i))*((-b)^A)*((-b/s)^(n-A))*(s^(-l-A))*risingFac(1+n-A,A)*risingFac(1+i-l-A,l)/(factorial(i)*factorial(A-i))
    }
  }
  return(result)
}

#get the l-th derivative of 2F1(a,b,c,z) wrt z
derive2F1z = function(l,a,b,c,z){
  foo = (risingFac(a,l)*risingFac(b,l)/(risingFac(c,l)))*hypergeo(a+l,b+l,c+l,z)
  return(foo)
}

#get the l-th derivative of const*2F1(a,b,c,-m/s) wrt s, using 0.430.1 in Ryzhik's tables for the composite function
derive2F1s = function(const,l,a,b,c,m,s){
  foo = 0 
  for(q in 1:l){
    U = 0
    for(t in 0:(q-1)){
      U = U+((-1)^t)*choose(q,t)*((-m/s)^t)*deriveZ(l=l,n=q-t,b=m,s=s)
    }
    foo = foo+(U/factorial(q))*const*derive2F1z(l=q,a=a,b=b,c=c,z=-m/s)
  }
  return(foo)
}

#use Faa Di Bruno's formula to get l-th derivative of exp(const*2F1(a,b,c,-m/s)) wrt s
deriveEy = function(const,l,a,b,c,m,s){
  ints = vector(mode="integer",length=l)
  answer = 0
  loop = TRUE
  while(loop==TRUE){
    test = 0
    for(i in 1:l){
      test = test+i*ints[i]
    }
    if(test==l){
    
    foo = (factorial(l)/(prod(factorial(ints))))*exp(const*hypergeo(a,b,c,-m/s))
    prod=1
    for(i in 1:l){
      prod=prod*(derive2F1s(const,i,a,b,c,m,s)/(factorial(i)))^(ints[i])
    }
    answer = answer+foo*prod
      
    }  
    ints[1]=ints[1]+1
    if(l>1){
    for(i in 1:(l-1)){
      if(ints[i]>floor(l/i)){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    }
    if(ints[l]>1){
      loop = FALSE
    } 
    
  }
  return(answer)
}

#get l-th order derivative of Laplace transform in the form of 
#prod_(k=1)^(2m)exp(-pi*density*(m choose ceil(k/2))*(-1)^(ceil(k/2)+1)*(-1)^(k+1)*(max or min)*2F1(a,b,c,-(max or min)^(al/2)/s))
###Very computationally intensive, made obsolete by better-performing function below
deriveLaplaceObs1 = function(max,min,alpha,m,l,const,b,c,s){
  ints = vector(mode="integer",length=(2*m))
  answer = 0
  loop = TRUE
  while(loop==TRUE){
  
    if(sum(ints)==l){
      
      foo = (factorial(l)/(prod(factorial(ints))))
      prod=1
      for(k in 1:(2*m)){
        if(k%%2==1){
          B=max
        }
        else{
          B=min
        }
        if(ints[k]>0){
        prod=prod*deriveEy(const=const*choose(m,ceil(k/2))*B*(-1)^(ceil(k/2)+1)*(-1)^(k+1),l=ints[k],a=ceil(k/2),b=b,c=c,m=B^(alpha/2),s=s)
        }
        else{
        prod=prod*exp(const*choose(m,ceil(k/2))*B*(-1)^(ceil(k/2)+1)*(-1)^(k+1)*hypergeo(ceil(k/2),b,c,-B^(alpha/2)/s)) 
        }
        }
      answer = answer+foo*prod
    }
      
    
    ints[1]=ints[1]+1
    for(i in 1:(2*m-1)){
      if(ints[i]>l){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    if(ints[2*m]>l){
      loop = FALSE
    }  
  }
  return(answer)
}


#get l-th order derivative of Laplace transform in the form of 
#exp(sum_(k=1)^(2m)-pi*density*(m choose ceil(k/2))*(-1)^(ceil(k/2)+1)*(-1)^(k+1)*(max or min)*2F1(a,b,c,-(max or min)^(al/2)/s))
##getting strange errors in some specific cases, like when SINR threshold = 0
deriveLaplace = function(max,min,alpha,m,l,const,b,c,s){
  ints = vector(mode="integer",length=l)
  answer = 0
  loop = TRUE
  if(l>0){
  while(loop==TRUE){
    test = 0
    for(i in 1:l){
      test = test+i*ints[i]
    }
    if(test==l){
      
      foo = (factorial(l)/(prod(factorial(ints))))
      for (k in 1:(2*m)){
        if(k%%2==1){
          B=max
        }
        else{
          B=min
        }
        foo = foo*exp(const*choose(m,ceil(k/2))*B*(-1)^(ceil(k/2)+1)*(-1)^(k+1)*hypergeo(ceil(k/2),b,c,-m*B^(alpha/2)/s)) 
      }
      prod=1
      for(i in 1:l){
        der = 0
        for (k in 1:(2*m)){
          if(k%%2==1){
            B=max
          }
          else{
            B=min
          }
          der = der + derive2F1s(const=const*choose(m,ceil(k/2))*B*(-1)^(ceil(k/2)+1)*(-1)^(k+1),l=i,a=ceil(k/2),b=b,c=c,m=m*B^(alpha/2),s=s)
        }
        prod=prod*(der/(factorial(i)))^(ints[i])
      }
      answer = answer+foo*prod
      
    }  
    ints[1]=ints[1]+1
    if(l>1){
      for(i in 1:(l-1)){
        if(ints[i]>floor(l/i)){
          ints[i]=0
          ints[i+1]=ints[i+1]+1
        }
      }
    }
    if(ints[l]>1){
      loop = FALSE
    } 
    
  }
  }
  else{
    foo = (factorial(l)/(prod(factorial(ints))))
    for (k in 1:(2*m)){
      if(k%%2==1){
        B=max
      }
      else{
        B=min
      }
      foo = foo*exp(const*choose(m,ceil(k/2))*B*(-1)^(ceil(k/2)+1)*(-1)^(k+1)*hypergeo(ceil(k/2),b,c,-m*B^(alpha/2)/s)) 
    }
    answer = foo
  }
  return(answer)
}

deriveLaplaceLOSSINR = function(max,min,N,al,an,alpha,beta,gamma,rxheight,txheight,mal,man,l,const,s){
  ints = vector(mode="integer",length=3)
  answer = 0
  loop = TRUE
  while(loop==TRUE){
    
    if(sum(ints)==l){
      
      foo = (factorial(l)/(prod(factorial(ints))))
      prod=1
      
      prod = prod*deriveLaplaceLOS(max=max,min=min,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
      prod = prod*deriveLaplaceNLOS(max=max,min=min,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)
      prod = prod*((-1)^ints[3])*((N)^(ints[3]))*exp(-N*s)
      answer = answer+foo*prod
    }
    
    
    ints[1]=ints[1]+1
    for(i in 1:2){
      if(ints[i]>l){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    if(ints[3]>l){
      loop = FALSE
    }  
  }
  return(answer)
}

##variant of the above for Noise only
deriveLaplaceLOSSNR = function(N,s,l){
  ints = vector(mode="integer",length=3)
  answer = 0
  loop = TRUE
  while(loop==TRUE){
    
    if(sum(ints)==l){
      
      foo = (factorial(l)/(prod(factorial(ints))))
      prod=1
      prod = prod*((-1)^ints[3])*((N)^(ints[3]))*exp(-N*s)
      answer = answer+foo*prod
    }
    
    
    ints[1]=ints[1]+1
    for(i in 1:2){
      if(ints[i]>l){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    if(ints[3]>l){
      loop = FALSE
    }  
  }
  return(answer)
}

strongestLaplaceLOSSINR = function(max,min,type,N,al,an,alpha,beta,gamma,rxheight,txheight,mal,man,l,const,s){
  ints = vector(mode="integer",length=3)
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
      }
      else{
      prod = prod*deriveLaplaceLOS(max=max,min=rl,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=mal,l=ints[1],const=const,b=2/al,c=1+2/al,s=s)
      prod = prod*deriveLaplaceNLOS(max=max,min=min,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,m=man,l=ints[2],const=const,b=2/an,c=1+2/an,s=s)  
      }
      prod = prod*((-1)^ints[3])*((N)^(ints[3]))*exp(-N*s)
      answer = answer+foo*prod
    }
    
    
    ints[1]=ints[1]+1
    for(i in 1:2){
      if(ints[i]>l){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    if(ints[3]>l){
      loop = FALSE
    }  
  }
  return(answer)
}

##variant of the above function for the noise-limited case
strongestLaplaceLOSSNR = function(max,min,type,N,al,an,alpha,beta,gamma,rxheight,txheight,mal,man,l,const,s){
  ints = vector(mode="integer",length=3)
  answer = 0
  loop = TRUE
  while(loop==TRUE){
    
    if(sum(ints)==l){
      
      foo = (factorial(l)/(prod(factorial(ints))))
      prod=1
      prod = prod*((-1)^ints[3])*((N)^(ints[3]))*exp(-N*s)
      answer = answer+foo*prod
    }
    
    
    ints[1]=ints[1]+1
    for(i in 1:2){
      if(ints[i]>l){
        ints[i]=0
        ints[i+1]=ints[i+1]+1
      }
    }
    if(ints[3]>l){
      loop = FALSE
    }  
  }
  return(answer)
}

deriveLaplaceLOS = function(max,min,al,alpha,beta,gamma,rxheight,txheight,m,l,const,b,c,s){
  ints = vector(mode="integer",length=l)
  answer = 0
  loop = TRUE
  if(l>0){
    while(loop==TRUE){
      test = 0
      for(i in 1:l){
        test = test+i*ints[i]
      }
      if(test==l){
        
        foo = (factorial(l)/(prod(factorial(ints))))
        for (k in 1:(m)){
         foo = foo*exp(Derive2F1LOS(max=max,min=min,m=m,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=0,const=const*choose(m,k),b=b,c=c,s=s))
        }
        prod=1
        for(i in 1:l){
          der = 0
          for (k in 1:(m)){
            der = der + Derive2F1LOS(max=max,min=min,m=m,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=i,const=const*choose(m,k),b=b,c=c,s=s)
          }
          prod=prod*(der/(factorial(i)))^(ints[i])
        }
        answer = answer+foo*prod
        
      }  
      ints[1]=ints[1]+1
      if(l>1){
        for(i in 1:(l-1)){
          if(ints[i]>floor(l/i)){
            ints[i]=0
            ints[i+1]=ints[i+1]+1
          }
        }
      }
      if(ints[l]>1){
        loop = FALSE
      } 
      
    }
  }
  else{
    foo = (factorial(l)/(prod(factorial(ints))))
    for (k in 1:(m)){
      foo = foo*exp(Derive2F1LOS(max=max,min=min,m=m,al=al,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=0,const=const*choose(m,k),b=b,c=c,s=s))
    }
    answer = foo
  }
  return(answer)
}

Derive2F1LOS = function(max,min,m,al,alpha,beta,gamma,rxheight,txheight,k,l,const,b,c,s){
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
      sum = sum + Plos*const*up*(-1)^(k+1)*hypergeo(k,b,c,-m*up^(al/2)/s)
      sum = sum - Plos*const*lo*(-1)^(k+1)*hypergeo(k,b,c,-m*lo^(al/2)/s)
    }
    else{
      sum = sum + derive2F1s(const=Plos*const*up*(-1)^(k+1),l=l,a=k,b=b,c=c,m=m*up^(al/2),s=s)
      sum = sum - derive2F1s(const=Plos*const*lo*(-1)^(k+1),l=l,a=k,b=b,c=c,m=m*lo^(al/2),s=s)
    }
  }
  
  return(sum)
}

deriveLaplaceNLOS = function(max,min,an,alpha,beta,gamma,rxheight,txheight,m,l,const,b,c,s){
  ints = vector(mode="integer",length=l)
  answer = 0
  loop = TRUE
  if(l>0){
    while(loop==TRUE){
      test = 0
      for(i in 1:l){
        test = test+i*ints[i]
      }
      if(test==l){
        
        foo = (factorial(l)/(prod(factorial(ints))))
        for (k in 1:(m)){
          foo = foo*exp(Derive2F1NLOS(max=max,min=min,m=m,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=0,const=const*choose(m,k),b=b,c=c,s=s))
        }
        prod=1
        for(i in 1:l){
          der = 0
          for (k in 1:(m)){
            der = der + Derive2F1NLOS(max=max,min=min,m=m,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=i,const=const*choose(m,k),b=b,c=c,s=s)
          }
          prod=prod*(der/(factorial(i)))^(ints[i])
        }
        answer = answer+foo*prod
        
      }  
      ints[1]=ints[1]+1
      if(l>1){
        for(i in 1:(l-1)){
          if(ints[i]>floor(l/i)){
            ints[i]=0
            ints[i+1]=ints[i+1]+1
          }
        }
      }
      if(ints[l]>1){
        loop = FALSE
      } 
      
    }
  }
  else{
    foo = (factorial(l)/(prod(factorial(ints))))
    for (k in 1:(m)){
      foo = foo*exp(Derive2F1NLOS(max=max,min=min,m=m,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxheight,txheight=txheight,k=k,l=0,const=const*choose(m,k),b=b,c=c,s=s))
    }
    answer = foo
  }
  return(answer)
}

Derive2F1NLOS = function(max,min,m,an,alpha,beta,gamma,rxheight,txheight,k,l,const,b,c,s){
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
      sum = sum + PNlos*const*up*(-1)^(k+1)*hypergeo(k,b,c,-m*up^(an/2)/s)
      sum = sum - PNlos*const*lo*(-1)^(k+1)*hypergeo(k,b,c,-m*lo^(an/2)/s)
    }
    else{
      sum = sum + derive2F1s(const=PNlos*const*up*(-1)^(k+1),l=l,a=k,b=b,c=c,m=m*up^(an/2),s=s)
      sum = sum - derive2F1s(const=PNlos*const*lo*(-1)^(k+1),l=l,a=k,b=b,c=c,m=m*lo^(an/2),s=s)
    }
  }
  
  return(sum)
}



#get l-th order derivative of Laplace transform in the form of 
#exp(sum_(k=1)^(2m)-pi*density*(m choose ceil(k/2))*(-1)^(ceil(k/2)+1)*(-1)^(k+1)*(max or min)*2F1(a,b,c,-(max or min)^(al/2)/s))
##getting strange errors in some specific cases, like when SINR threshold = 0
kmuLaplace = function(max,min,alpha,kappa,mu,m,l,const,b,c,s){
  ints = vector(mode="integer",length=l)
  answer = 0
  loop = TRUE
  if(l>0){
    while(loop==TRUE){
      test = 0
      for(i in 1:l){
        test = test+i*ints[i]
      }
      if(test==l){
        foo = (factorial(l)/(prod(factorial(ints))))
        foo = foo*exp(const*(max-min))
        for (j in 0:(m-mu)){
          B = choose(m-mu,j)*((m/(mu*kappa+m))^j)*(mu*kappa/(mu*kappa+m))^(m-mu-j)
          wb = (mu*kappa+m)*(m-j)/(m*mu*(1+kappa))
          for (k in 0:((m-j))){
            
            foo = foo*exp(-const*B*choose(m-j,k)*max*(-1)^(k)*hypergeo(k,b,c,-(m-j)*max^(alpha/2)/(wb*s)))
            foo = foo*exp(+const*B*choose(m-j,k)*min*(-1)^(k)*hypergeo(k,b,c,-(m-j)*min^(alpha/2)/(wb*s))) 
          }
        }
  
        prod=1
        for(i in 1:l){
          der = 0
          for (j in 0:(m-mu)){
          B = choose(m-mu,j)*((m/(mu*kappa+m))^j)*(mu*kappa/(mu*kappa+m))^(m-mu-j)
          wb = (mu*kappa+m)*(m-j)/(m*mu*(1+kappa))
          for (k in 0:(m-j)){
            der = der - derive2F1s(const=const*B*choose(m-j,k)*max*(-1)^(k),l=i,a=k,b=b,c=c,m=(m-j)*max^(alpha/2)/wb,s=s)
            der = der + derive2F1s(const=const*B*choose(m-j,k)*min*(-1)^(k),l=i,a=k,b=b,c=c,m=(m-j)*min^(alpha/2)/wb,s=s)
          }
          }
          prod=prod*(der/(factorial(i)))^(ints[i])
        }
        answer = answer+foo*prod
        
      }  
      ints[1]=ints[1]+1
      if(l>1){
        for(i in 1:(l-1)){
          if(ints[i]>floor(l/i)){
            ints[i]=0
            ints[i+1]=ints[i+1]+1
          }
        }
      }
      if(ints[l]>1){
        loop = FALSE
      } 
      
    }
  }
  else{
    foo = (factorial(l)/(prod(factorial(ints))))
    foo = foo*exp(const*(max-min))
    for (j in 0:(m-mu)){
      B = choose(m-mu,j)*((m/(mu*kappa+m))^j)*(mu*kappa/(mu*kappa+m))^(m-mu-j)
      wb = (mu*kappa+m)*(m-j)/(m*mu*(1+kappa))
      for (k in 0:((m-j))){
       
        foo = foo*exp(-const*B*choose(m-j,k)*max*(-1)^(k)*hypergeo(k,b,c,-(m-j)*max^(alpha/2)/(wb*s)))
        foo = foo*exp(+const*B*choose(m-j,k)*min*(-1)^(k)*hypergeo(k,b,c,-(m-j)*min^(alpha/2)/(wb*s))) 
      }
    }
    answer = foo
  }
  return(answer)
}

#function to calculate the Laplace of the interference using the Sigmoid function. Note that this has to be done numerically and therefore only works for m=1 (Rayleigh Fading
sigmoidLaplace = function(max,min,A,B,height,density,alpha,s,steps,LOS){
  foo = 0
  r = seq(from=min,to=max,by=(max-min)/steps)
  for (i in 1:length(r)){
    angle = atan(height/r[i])*(180/pi)
    if(LOS==TRUE){
      P = 1/(1+A*exp(-B*(angle-A)))
    }
    else{
      P = 1-(1/(1+A*exp(-B*(angle-A))))
    }
    foo = foo + P*(r[i]/(1+((r[i]^2+height^2)^(alpha/2))/s))*(r[2]-r[1])
  }
  return(exp(-2*density*pi*foo))
}