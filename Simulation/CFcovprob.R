##define the analytical coverage prob calculation as a function


source("derivFunctionsMatern.R")

####UAV Coverage Probability calculation####
CFcov = function(rxh,txh,angle,density,T,gain,al,an,mal,man,alpha,beta,gamma,dx){
    Tlin = 10^(T/10)
    rmax = abs(rxh-txh)*tan(angle/2)
    M = (rmax^2+(rxh-txh)^2)
    
    r = seq(from=rmax/dx,to=rmax,by=rmax/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
    
    Cov = zeros(nrow=length(r),ncol=4)
    for (q in 1:length(r)){
    ds = (r[q]^2+(rxh-txh)^2)
    
    rn = sqrt(max(0,(r[q]^2+(rxh-txh)^2)^(al/an) - (rxh-txh)^2))
    rl = min(rmax,sqrt((r[q]^2+(rxh-txh)^2)^(an/al) - (rxh-txh)^2))
    
    foo = seq(from=0,to=r[q],by=r[q]/500)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
    n = floor((foo[i]/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
  
    ilP = 2*pi*density*sum(Plos*foo*(foo[2]-foo[1]))
    inP = 2*pi*density*sum((1-Plos)*foo*(foo[2]-foo[1]))
  #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      if(p==1){
       s = mal*Tlin*ds^(al/2)
       ms=mal
       type=al
      }
      else{
        s = man*Tlin*ds^(an/2) 
        ms=man
        type=an
      }
    Cov[q,p] = 0
    for (i in 0:(ms-1)){
      L = strongestLaplaceLOSSINR(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s)
        Cov[q,p]= Cov[q,p] + (s^i)/(factorial(i))*((-1)^i)*Re(L)
    }
    }
    
    if(rn==0){
      Prob0n = 0
    }
    else{
      foo = seq(from=0,to=rn,by=rn/500)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0n = 2*pi*density*sum((1-Ptemp)*foo*(foo[2]-foo[1]))
    }
    
    #  rl = min(rmax[j],(r[k]^2+h^2)^(an/al) - h^2)
    if(rl==0){
      Prob0l = 0
    }
    else{
      foo = seq(from=0,to=rl,by=rl/500)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0l = 2*pi*density*sum(Ptemp*foo*(foo[2]-foo[1]))
    } 
  
     n = floor((r[q]/1000)*sqrt(alpha*beta))
     a = 0:max(0,n-1)
     P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))

    Cov[q,1] = Cov[q,1]*2*P*density*pi*r[q]*exp(-ilP)*(r[2]-r[1])*exp(-Prob0n)
    Cov[q,2] = Cov[q,2]*2*(1-P)*density*pi*r[q]*exp(-inP)*(r[2]-r[1])*exp(-Prob0l)
   # Cov[q,] = P*Cov[q,]
   Cov[q,1] = Cov[q,1]+Cov[q,2]
  #  Cov[q,1] = Cov[q,1]*y[q]
    }
    
    return(sum(Cov[,1]))
}


CFBHcov = function(rxh,txh,angle,rmax,density,alignprob,T,horgain,gain,al,an,mal,man,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
 # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=2)
  for (q in 1:length(r)){
    if(angle<pi/2){
    ang = atan2(abs(rxh-txh),r[q])
    if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
      threshold = abs(rxh-txh)/tan(ang-(angle/2))  
    }
    else if(ang>(pi/2-(angle/2))){
       threshold = abs(rxh-txh)/tan(pi/2 - angle)  
    }
    else{
      threshold=rmax
    }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
   
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      if(p==1){
        s = mal*Tlin*ds^(al/2)
        ms=mal
        type=al
      }
      else{
        s = man*Tlin*ds^(an/2) 
        ms=man
        type=an
      }
      Cov[q,p] = 0
      for (i in 0:(ms-1)){
        L = deriveLaplaceLOSSINR(max=threshold,min=r[q],N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-(angle/2)*density*alignprob,s=s)
        Cov[q,p]= Cov[q,p] + (s^i)/(factorial(i))*((-1)^i)*Re(L)
      }
    }
    
   
    
  #  n = floor((r[q]/1000)*sqrt(alpha*beta))
  #  a = 0:max(0,n-1)
  #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
  
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    #  Cov[q,1] = Cov[q,1]*y[q]
  }
  
  return(sum(Cov[,1]))
}

MCCFBHcov = function(rxh,txh,angle,rmax,density,alignprob,T,horgain,gain,al,an,mal,man,alpha,beta,gamma,dx,mc){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
 # Cov = zeros(nrow=length(r),ncol=2)
  iteration =function(q){
    Cov = zeros(nrow=1,ncol=2)
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      if(p==1){
        s = mal*Tlin*ds^(al/2)
        ms=mal
        type=al
      }
      else{
        s = man*Tlin*ds^(an/2) 
        ms=man
        type=an
      }
      Cov[1,p] = 0
      for (i in 0:(ms-1)){
        L = deriveLaplaceLOSSINR(max=threshold,min=r[q],N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-(angle/2)*density*alignprob,s=s)
        Cov[1,p]= Cov[1,p] + (s^i)/(factorial(i))*((-1)^i)*Re(L)
      }
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[1,1] = Cov[1,1]*P
    Cov[1,2] = Cov[1,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[1,1] = (Cov[1,1]+Cov[1,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    #  Cov[q,1] = Cov[q,1]*y[q]
    return(Cov[1,1])
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=mc)
  
  covP = 0
  for(k in 1:length(r)){
    covP = covP+c[[k]]
  }
  
  
  return(covP)
}

##Closed form expression for backhaul calculation for LTE BSs, involving numerically integrating the Laplace transform of interference
CFLTEBHcov = function(rxh,txh,angle,rmax,density,uptilt,T,gain,horgain,al,an,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=2)
  for (q in 1:length(r)){
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
      P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
      ang = atan2((rxh-txh),dr[i])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
      Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
      g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
      Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    #  Cov[q,1] = Cov[q,1]*y[q]
  }
  
  return(sum(Cov[,1]))
}

##As above, but where the Nth closest BS is the serving BS, and the closest (N-1) BSs are interfering
CFNthLTEBHcov = function(rxh,txh,angle,rlim,rmax,skipN,Nint,density,uptilt,T,gain,horgain,al,an,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = rlim
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=2)
  for (q in 1:length(r)){
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
      
      if(ang<(pi/2-angle/2)){
        lthreshold = abs(rxh-txh)/tan(ang+(angle/2))  
      }else{lthreshold=0}
      
    }else{threshold=rmax
    lthreshold = 0
    }
    
    
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/500) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      if(Nint==TRUE && skipN>1){
        probIn = (((r[q]^2-lthreshold^2)*angle/2)/(pi*r[q]^2))
        dr = seq(from=lthreshold,to=r[q],by=(r[q]-lthreshold)/500)
        LIN = vector(length=length(dr))
        for(i in 1:length(dr)){
          P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
          ang = atan2((rxh-txh),dr[i])*(180/pi)
          bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
          gn = max(10^(horgain/10)*bv,10^(-2.5))
          g1 = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
          g2 = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
          LIN[i] = ((1/(g1+1))*P+(1/(g2+1))*(1-P))*2*(dr[i]/(r[q]^2-lthreshold^2))*(dr[2]-dr[1])
        }
      }
      else{LIN=1
      probIn = 0
      }
      
      Cov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)*((probIn*sum(LIN)+(1-probIn))^(skipN-1))
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*((pi*density)^(skipN))*(r[q]^(2*skipN-1))*exp(-pi*density*r[q]^2)*(r[2]-r[1])/factorial(skipN-1)
    #  Cov[q,1] = Cov[q,1]*y[q]
  }
  
  return(sum(Cov[,1]))
}


##As above, but where the Nth closest BS is the serving BS, and the closest (N-1) BSs are interfering (2 antenna variant)
CFNthLTEBHcov2ant = function(rxh,txh,angle,rlim,rmax,skipN,Nint,selfInt,density,tilt1,tilt2,T,gain,horgain,al,an,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = rlim
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=2)
  for (q in 1:length(r)){
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
      
      if(ang<(pi/2-angle/2)){
        lthreshold = abs(rxh-txh)/tan(ang+(angle/2))  
      }else{lthreshold=0}
      
    }else{threshold=rmax
    lthreshold = 0
    }
    
    
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-tilt1)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/500) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr))
      Ls =vector(length=length(dr))
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv1 = 10^(-min(12*((ang-tilt1)/10)^2,20)/10)
        gn1 = max(10^(horgain/10)*bv1,10^(-2.5))
        g1 = s*gn1*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        bv2 = 10^(-min(12*((ang-tilt2)/10)^2,20)/10)
        gn2 = max(10^(horgain/10)*bv2,10^(-2.5))
        g2 = s*gn2*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-(1/(g1+1))*(1/(g2+1)))*P*dr[i]*(dr[2]-dr[1])
        g1 = s*gn1*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        g2 = s*gn2*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-(1/(g1+1))*(1/(g2+1)))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      if(Nint==TRUE && skipN>1){
        probIn = (((r[q]^2-lthreshold^2)*angle/2)/(pi*r[q]^2))
        dr = seq(from=lthreshold,to=r[q],by=(r[q]-lthreshold)/500)
        LIN = vector(length=length(dr))
        for(i in 1:length(dr)){
          P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
          ang = atan2((rxh-txh),dr[i])*(180/pi)
          bv1 = 10^(-min(12*((ang-tilt1)/10)^2,20)/10)
          g1 = max(10^(horgain/10)*bv1,10^(-2.5))
          bv2 = 10^(-min(12*((ang-tilt2)/10)^2,20)/10)
          g2 = max(10^(horgain/10)*bv2,10^(-2.5))
          gl1 = s*g1*(dr[i]^2+(rxh-txh)^2)^(-al/2)
          gn1 = s*g1*(dr[i]^2+(rxh-txh)^2)^(-an/2)
          gl2 = s*g2*(dr[i]^2+(rxh-txh)^2)^(-al/2)
          gn2 = s*g2*(dr[i]^2+(rxh-txh)^2)^(-an/2)
          LIN[i] = ((1/(gl1+1))*(1/(gl2+1))*P+(1/(gn1+1))*(1/(gn2+1))*(1-P))*2*(dr[i]/(r[q]^2-lthreshold^2))*(dr[2]-dr[1])
        }
      }
      else{LIN=1
      probIn = 0
      }
      
      if(selfInt==TRUE){
      bv2 = 10^(-min(12*((ang-tilt2)/10)^2,20)/10)
      gn2 = max(10^(horgain/10)*bv2,10^(-2.5))
      g2 = s*gn2*(r[q]^2+(rxh-txh)^2)^(-al/2)
      Ls=(1/(g2+1))
      }
      else{Ls=1}
      
      Cov[q,p] = Ls*exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)*((probIn*sum(LIN)+(1-probIn))^(skipN-1))
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*((pi*density)^(skipN))*(r[q]^(2*skipN-1))*exp(-pi*density*r[q]^2)*(r[2]-r[1])/factorial(skipN-1)
    #  Cov[q,1] = Cov[q,1]*y[q]
  }
  
  return(sum(Cov[,1]))
}


##Closed form expression for backhaul calculation for LTE BSs, given another network of LTE BSs is causing interference
CFsameFreqLTEBHcov = function(rxh,txh,angle,rmax,density,idensity,uptilt,iuptilt,T,gain,horgain,al,an,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=2)
  for (q in 1:length(r)){
    ang = atan2(abs(rxh-txh),r[q])
    if(angle<pi/2){
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    if(ang<(pi/2-angle/2)){
      lthreshold = abs(rxh-txh)/tan(ang+(angle/2))  
    }else{lthreshold=0}
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      dr = seq(from=lthreshold,to=threshold,by=(threshold-lthreshold)/200) 
      iLl = vector(length=length(dr))
      iLn = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-iuptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        iLl[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        iLn[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-idensity*angle*sum(iLl))*exp(-idensity*angle*sum(iLn))*exp(-s*N/gain)
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    #  Cov[q,1] = Cov[q,1]*y[q]
  }
  
  return(sum(Cov[,1]))
}

##As above, but UAV is within a certain radius of its serving BS 
CFClusterLTEBHcov = function(rxh,txh,angle,rmax,rrad,density,uptilt,T,gain,horgain,al,an,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = rrad
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=2)
  for (q in 1:length(r)){
    ang = atan2(abs(rxh-txh),r[q])
    if(angle<pi/2){
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    
    if(ang<(pi/2-angle/2)){
      lthreshold = abs(rxh-txh)/tan(ang+(angle/2))  
    }else{lthreshold=0}
    
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=lthreshold,to=threshold,by=(threshold-lthreshold)/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*r[q]*(r[2]-r[1])/(rrad^2)
    #  Cov[q,1] = Cov[q,1]*y[q]
  }
  
  return(sum(Cov[,1]))
}



##Closed form expression for backhaul calculation for LTE BSs, involving numerically integrating the Laplace transform of interference
CFLTEBHcovOmni = function(rxh,txh,rmax,density,uptilt,T,gain,horgain,al,an,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = min(rmax,sqrt(-log(0.001)/(pi*density)))
  angle = 2*pi
  r = seq(from=rm/dx,to=(rm-0.1),by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=2)
  for (q in 1:length(r)){
    threshold = min(rmax,10000)
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    #  Cov[q,1] = Cov[q,1]*y[q]
  }
  
  return(sum(Cov[,1]))
}


####Closed form expression for backhaul calculation for Mmwave BSs, assuming no interference

CFMMBHcovNoInt = function(rxh,txh,angle,rmax,density,alignprob,T,N,horgain,gain,al,an,mal,man,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=2)
  for (q in 1:length(r)){
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      if(p==1){
        s = mal*Tlin*ds^(al/2)
        ms=mal
        type=al
      }
      else{
        s = man*Tlin*ds^(an/2) 
        ms=man
        type=an
      }
      Cov[q,p] = 0
      for (i in 0:(ms-1)){
        L = ((-1)^i)*((N/gain)^i)*exp(-(N/gain)*s)
        Cov[q,p]= Cov[q,p] + (s^i)/(factorial(i))*((-1)^i)*Re(L)
      }
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    #  Cov[q,1] = Cov[q,1]*y[q]
  }
  
  return(sum(Cov[,1]))
}


#######Matern UAV Distribution calculation#############


###get PDF of serving UAV distance under matern distribution
MaternDistPDF = function(r,Rmax,rRad,rxh,txh,al,an,mal,man,alpha,beta,gamma){
  
    
    rn = sqrt(max(0,(r^2+(rxh-txh)^2)^(al/an) - (rxh-txh)^2))
    rl = min(rRad,sqrt((r^2+(rxh-txh)^2)^(an/al) - (rxh-txh)^2))
    
    foo = seq(from=0,to=r,by=r/500)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    ilP = 2*pi*density*sum(Plos*foo*(foo[2]-foo[1]))
    inP = 2*pi*density*sum((1-Plos)*foo*(foo[2]-foo[1]))
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    if(rn==0){
      Prob0n = 0
      Prob0nc = 1
    }
    else{
      foo = seq(from=0,to=rn,by=rn/500)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0n = 2*pi*density*sum((1-Ptemp)*foo*(foo[2]-foo[1]))
      
      foo = seq(from=0,to=min(rn,Rmax),by=min(rn,Rmax)/500)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      
      Prob0nc = 1-2*sum((1-Ptemp)*foo*(foo[2]-foo[1]))/(Rmax^2)
    }
    
    #  rl = min(rmax[j],(r[k]^2+h^2)^(an/al) - h^2)
    if(rl==0){
      Prob0l = 0
      Prob0lc = 1
    }
    else{
      foo = seq(from=0,to=rl,by=rl/500)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0l = 2*pi*density*sum(Ptemp*foo*(foo[2]-foo[1]))
      
      foo = seq(from=0,to=min(rl,Rmax),by=min(rl,Rmax)/500)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0lc = 1-2*sum((Ptemp)*foo*(foo[2]-foo[1]))/(Rmax^2)
    } 
    
    n = floor((r/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    
    PPPLOS = 2*P*density*pi*r*exp(-ilP)*exp(-Prob0n)*Prob0nc*Prob0lc
    PPPNLOS = 2*(1-P)*density*pi*r*exp(-inP)*exp(-Prob0l)*Prob0nc*Prob0lc
    if(r<Rmax){
    CLOS = 2*P*r*exp(-ilP)*exp(-Prob0n)/(Rmax^2)
    CNLOS = 2*(1-P)*r*exp(-inP)*exp(-Prob0l)/(Rmax^2)
    }else{
      CLOS = 0
      CNLOS = 0
    }
    # Cov[q,] = P*Cov[q,]
    #  Cov[q,1] = Cov[q,1]*y[q]
  
  return(PPPLOS+PPPNLOS+CLOS+CNLOS)
}

###get CDF of serving UAV distance under matern distribution
MaternDistCDF = function(rlim,Rmax,rRad,rxh,txh,al,an,mal,man,alpha,beta,gamma){
  
  r = seq(from=0.1,to=rlim,by=(rlim-0.1)/200)
  pdf=vector(length=length(r))
  for (q in 1:length(r)){
  rn = sqrt(max(0,(r[q]^2+(rxh-txh)^2)^(al/an) - (rxh-txh)^2))
  rl = min(rRad,sqrt((r[q]^2+(rxh-txh)^2)^(an/al) - (rxh-txh)^2))
  
  foo = seq(from=0,to=r[q],by=r[q]/200)
  Plos = vector(mode="numeric",length=length(foo))
  for(i in 1:length(foo)){
    n = floor((foo[i]/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
  }
  
  ilP = 2*pi*density*sum(Plos*foo*(foo[2]-foo[1]))
  inP = 2*pi*density*sum((1-Plos)*foo*(foo[2]-foo[1]))
  
  foo = seq(from=0,to=min(r[q],Rmax),by=min(r[q],Rmax)/200)
  Plos = vector(mode="numeric",length=length(foo))
  for(i in 1:length(foo)){
    n = floor((foo[i]/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
  }
  
  cilP = 2*sum(Plos*foo*(foo[2]-foo[1]))/(Rmax^2)
  cinP = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(Rmax^2)
  
  #  P = c(Plos,(1-Plos))
  #ProbCov[l,j]= 0
  
  if(rn==0){
    Prob0n = 0
    Prob0nc = 0
  }
  else{
    foo = seq(from=0,to=rn,by=rn/200)
    Ptemp = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    }
    Prob0n = 2*pi*density*sum((1-Ptemp)*foo*(foo[2]-foo[1]))
    
    foo = seq(from=0,to=min(rn,Rmax),by=min(rn,Rmax)/200)
    Ptemp = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    }
    
    Prob0nc = (2*sum((1-Ptemp)*foo*(foo[2]-foo[1]))/(Rmax^2))
  }
  
  #  rl = min(rmax[j],(r[q][k]^2+h^2)^(an/al) - h^2)
  if(rl==0){
    Prob0l = 0
    Prob0lc = 0
  }
  else{
    foo = seq(from=0,to=rl,by=rl/200)
    Ptemp = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    }
    Prob0l = 2*pi*density*sum(Ptemp*foo*(foo[2]-foo[1]))
    
    foo = seq(from=0,to=min(rl,Rmax),by=min(rl,Rmax)/200)
    Ptemp = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    }
    Prob0lc = 2*sum((Ptemp)*foo*(foo[2]-foo[1]))/(Rmax^2)
  } 
  
  n = floor((r[q]/1000)*sqrt(alpha*beta))
  a = 0:max(0,n-1)
  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
  
  PPPLOS = 2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*(1-(Prob0nc+cilP))
  PPPNLOS = 2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*(1-(cinP+Prob0lc))
  if(r[q]<Rmax){
    CLOS = 2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(Rmax^2)
    CNLOS = 2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(Rmax^2)
  }else{
    CLOS = 0
    CNLOS = 0
  }
  pdf[q] = PPPLOS+PPPNLOS+CLOS+CNLOS
  pdf[q] = pdf[q]*(r[2]-r[1])
  }
  
  return(sum(pdf))
}

MCMaternDistCDF = function(rlim,Rmax,rRad,rxh,txh,al,an,mal,man,alpha,beta,gamma,dx,cores){
  
  r = seq(from=0.1,to=rlim,by=(rlim-0.1)/dx)
  iteration=function(q){
    rn = sqrt(max(0,(r[q]^2+(rxh-txh)^2)^(al/an) - (rxh-txh)^2))
    rl = min(rRad,sqrt((r[q]^2+(rxh-txh)^2)^(an/al) - (rxh-txh)^2))
    
    foo = seq(from=0,to=r[q],by=r[q]/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    ilP = 2*pi*density*sum(Plos*foo*(foo[2]-foo[1]))
    inP = 2*pi*density*sum((1-Plos)*foo*(foo[2]-foo[1]))
    
    foo = seq(from=0,to=min(r[q],Rmax),by=min(r[q],Rmax)/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    cilP = 2*sum(Plos*foo*(foo[2]-foo[1]))/(Rmax^2)
    cinP = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(Rmax^2)
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    if(rn==0){
      Prob0n = 0
      Prob0nc = 0
    }
    else{
      foo = seq(from=0,to=rn,by=rn/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0n = 2*pi*density*sum((1-Ptemp)*foo*(foo[2]-foo[1]))
      
      foo = seq(from=0,to=min(rn,Rmax),by=min(rn,Rmax)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      
      Prob0nc = (2*sum((1-Ptemp)*foo*(foo[2]-foo[1]))/(Rmax^2))
    }
    
    #  rl = min(rmax[j],(r[q][k]^2+h^2)^(an/al) - h^2)
    if(rl==0){
      Prob0l = 0
      Prob0lc = 0
    }
    else{
      foo = seq(from=0,to=rl,by=rl/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0l = 2*pi*density*sum(Ptemp*foo*(foo[2]-foo[1]))
      
      foo = seq(from=0,to=min(rl,Rmax),by=min(rl,Rmax)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0lc = 2*sum((Ptemp)*foo*(foo[2]-foo[1]))/(Rmax^2)
    } 
    
    n = floor((r[q]/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    
    PPPLOS = 2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*(1-(Prob0nc+cilP))
    PPPNLOS = 2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*(1-(cinP+Prob0lc))
    if(r[q]<Rmax){
      CLOS = 2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(Rmax^2)
      CNLOS = 2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(Rmax^2)
    }else{
      CLOS = 0
      CNLOS = 0
    }
    pdf = PPPLOS+PPPNLOS+CLOS+CNLOS
    pdf = pdf*(r[2]-r[1])
    return(pdf)
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  PDF = 0
  for(k in 1:length(r)){
    PDF = PDF+c[[k]]
  }
  
  return(PDF)
}

#analytical calculation of the coverage probability
MaternCFCov = function(angle,density,T,gain,rRad,rxh,txh,al,an,mal,man,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  r = seq(from=min(rmax,rRad*1.5)/dx,to=min(rmax,rRad),by=min(rmax,rRad*1.5)/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=4)
  Coverage = vector(length=length(r))
  for (q in 1:length(r)){
      if(q%%10==0){
        print(q)
      }
    ds = (r[q]^2+(rxh-txh)^2)
    
    rn = sqrt(max(0,(r[q]^2+(rxh-txh)^2)^(al/an) - (rxh-txh)^2))
    rl = min(rmax,sqrt((r[q]^2+(rxh-txh)^2)^(an/al) - (rxh-txh)^2))
    
    foo = seq(from=0,to=r[q],by=r[q]/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    ilP = 2*pi*density*sum(Plos*foo*(foo[2]-foo[1]))
    inP = 2*pi*density*sum((1-Plos)*foo*(foo[2]-foo[1]))
    
    foo = seq(from=0,to=min(r[q],rRad),by=min(r[q],rRad)/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    cilP = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
    cinP = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
    
    if(rn==0){
      Prob0n = 0
      Prob0nc = 0
    }
    else{
      foo = seq(from=0,to=rn,by=rn/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0n = 2*pi*density*sum((1-Ptemp)*foo*(foo[2]-foo[1]))
      
      foo = seq(from=0,to=min(rn,rRad),by=min(rn,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      
      Prob0nc = (2*sum((1-Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2))
    }
    
    #  rl = min(rmax[j],(r[q][k]^2+h^2)^(an/al) - h^2)
    if(rl==0){
      Prob0l = 0
      Prob0lc = 0
    }
    else{
      foo = seq(from=0,to=rl,by=rl/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0l = 2*pi*density*sum(Ptemp*foo*(foo[2]-foo[1]))
      
      foo = seq(from=0,to=min(rl,rRad),by=min(rl,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0lc = 2*sum((Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2)
    } 
    
    
    ##1- PPP LOS
    ##2- PPP NLOS
    ##3- cluster center LOS
    ##4- cluster center NLOS
    for (p in 1:4){
      if(p==1 | p==3){
        s = mal*Tlin*ds^(al/2)
        ms=mal
        type=al
      }
      else{
        s = man*Tlin*ds^(an/2) 
        ms=man
        type=an
      }
      Cov[q,p] = 0
      if(p == 3 | p==4){
      if(r[q]<rRad){
      for (i in 0:(ms-1)){
        L = strongestLaplaceLOSSINR(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s)
        Cov[q,p]= Cov[q,p] + (s^i)/(factorial(i))*((-1)^i)*Re(L)
      }
      }
      else{
        Cov[q,p]=0
      }
      }
      ##user serviced by PPP UAV
      else{
        
        ##get probability of interfering cluster UAV being LOS
        if(type==al){
          lower=r[q]
        }else{lower=rl}
        
        if((min(rmax,rRad)-lower)<10^(-4)){
          intisLOS=0 
        }
        else{
        foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
        Plos = vector(mode="numeric",length=length(foo))
        for(i in 1:length(foo)){
          n = floor((foo[i]/1000)*sqrt(alpha*beta))
          a = 0:max(0,n-1)
          Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
        }
#        intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
        intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
#        intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/((rRad^2)*max(0,(1-(Prob0nc+cilP))))
  #      intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/((rRad^2-lower^2)*max(0,(1-(Prob0nc+cilP))))
        intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
  #      intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
            if(is.na(intisLOS)){intisLOS=0}
        }
      
        ##get probability of interfering cluster UAV being NLOS
        if(type==an){
          lower=r[q]
        }else{lower=rn}
        
        if((min(rmax,rRad)-lower)<10^-4){
        intisNLOS=0 
        }
        else{
        foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
        Plos = vector(mode="numeric",length=length(foo))
        for(i in 1:length(foo)){
          n = floor((foo[i]/1000)*sqrt(alpha*beta))
          a = 0:max(0,n-1)
          Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
        }
#        intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
        intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
#      intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/((rRad^2)*max(0,(1-(cinP+Prob0lc))))
#        intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/((rRad^2-lower^2)*max(0,(1-(cinP+Prob0lc))))
        intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
 #       intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
        if(is.na(intisNLOS)){intisNLOS=0}
        }
        
   #     if(rRad<rmax){
  #        intis0 = 0
  #      }else{
  #      intis0 = (rRad^2-rmax^2)/(rRad^2)
      #  intis0 = (rRad^2-rmax^2)/((rRad^2)*(1-(Prob0nc+cilP)))
  #      }
        
        intisLOS = max(0,min(1,intisLOS))
        intisNLOS = max(0,min(1,intisNLOS))
        
        intis0 = max(0,(1-(intisLOS+intisNLOS)))
        
      #  if(p==1){
#        print(intisLOS+intisNLOS)
       # }
        for (i in 0:(ms-1)){
          L= 0
          if(intisLOS>0){
          L = strongestLaplaceClusterInt(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s,clusterType=al,clusterRad=rRad,clusterTypeProb=intisLOS2)
          Cov[q,p]= Cov[q,p] + intisLOS*((s^i)/(factorial(i))*((-1)^i)*Re(L))
          }
          if(intisNLOS>0){
            L = strongestLaplaceClusterInt(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s,clusterType=an,clusterRad=rRad,clusterTypeProb=intisNLOS2)
            Cov[q,p]= Cov[q,p] + intisNLOS*((s^i)/(factorial(i))*((-1)^i)*Re(L))
          }
          if(intis0>0){
          L = strongestLaplaceLOSSINR(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s)
          Cov[q,p]= Cov[q,p] + intis0*((s^i)/(factorial(i))*((-1)^i)*Re(L))
          }
        }
      }
    }
    
    
    n = floor((r[q]/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    
   # PPPLOS = 2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*(1-(Prob0nc+cilP))
   # PPPNLOS = 2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*(1-(cinP+Prob0lc))
    Cov[q,1] = Cov[q,1]*2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP)))
 #   cat(2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP)))/(r[2]-r[1])," ")
    Cov[q,2] = Cov[q,2]*2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc)))
 #   cat(2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc)))/(r[2]-r[1])," ")
    if(r[q]<rRad){
      Cov[q,3]=Cov[q,3]*2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
 #     cat(2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)/(r[2]-r[1])," ")
      Cov[q,4]=Cov[q,4]*2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
 #     cat(2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)/(r[2]-r[1]), "\n")
  #    CLOS = 2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
 #     CNLOS = 2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
    }else{
      Cov[q,3] = 0
 #     cat(0, " ")
      Cov[q,4] = 0
  #    cat(0, "\n")
#      CLOS = 0
#      CNLOS = 0
    }
 #   cat(Cov[q,1]," ",Cov[q,2]," ",Cov[q,3]," ",Cov[q,4],"\n")
    Coverage[q] = sum(Cov[q,])*(r[2]-r[1])
  }
  
  return(sum(Coverage))
}


MaternCFCov2 = function(angle,density,T,gain,rRad,rxh,txh,al,an,mal,man,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  r = seq(from=min(rmax,rRad)/dx,to=min(rmax,rRad),by=min(rmax,rRad)/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=4)
  Coverage = vector(length=length(r))
  for (q in 1:length(r)){
    if(q%%10==0){
      print(q)
    }
    ds = (r[q]^2+(rxh-txh)^2)
    
    rn = sqrt(max(0,(r[q]^2+(rxh-txh)^2)^(al/an) - (rxh-txh)^2))
    rl = min(rmax,sqrt((r[q]^2+(rxh-txh)^2)^(an/al) - (rxh-txh)^2))
    
    foo = seq(from=0,to=r[q],by=r[q]/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    ilP = 2*pi*density*sum(Plos*foo*(foo[2]-foo[1]))
    inP = 2*pi*density*sum((1-Plos)*foo*(foo[2]-foo[1]))
    
    foo = seq(from=0,to=min(r[q],rRad),by=min(r[q],rRad)/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    cilP = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
    cinP = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
    
    cilP = max(0,min(1,cilP))
    cinP = max(0,min(1,cinP))
    
    if(rn==0){
      Prob0n = 0
      Prob0nc = 0
    }
    else{
      foo = seq(from=0,to=rn,by=rn/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0n = 2*pi*density*sum((1-Ptemp)*foo*(foo[2]-foo[1]))
      Prob0n = max(0,min(1,Prob0n))
      
      foo = seq(from=0,to=min(rn,rRad),by=min(rn,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      
      Prob0nc = (2*sum((1-Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2))
      Prob0nc = max(0,min(1,Prob0nc))
    }
    
    #  rl = min(rmax[j],(r[q][k]^2+h^2)^(an/al) - h^2)
    if(rl==0){
      Prob0l = 0
      Prob0lc = 0
    }
    else{
      foo = seq(from=0,to=rl,by=rl/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0l = 2*pi*density*sum(Ptemp*foo*(foo[2]-foo[1]))
      Prob0l = max(0,min(1,Prob0l))
      
      foo = seq(from=0,to=min(rl,rRad),by=min(rl,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0lc = 2*sum((Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2)
      Prob0lc = max(0,min(1,Prob0lc))
    } 
    
    
    ##1- PPP LOS
    ##2- PPP NLOS
    ##3- cluster center LOS
    ##4- cluster center NLOS
    for (p in 1:4){
      if(p==1 | p==3){
        s = mal*Tlin*ds^(al/2)
        ms=mal
        type=al
      }
      else{
        s = man*Tlin*ds^(an/2) 
        ms=man
        type=an
      }
      Cov[q,p] = 0
      if(p == 3 | p==4){
        if(r[q]<rRad){
          for (i in 0:(ms-1)){
            L = strongestLaplaceLOSSINR(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s)
            Cov[q,p]= Cov[q,p] + (s^i)/(factorial(i))*((-1)^i)*Re(L)
          }
        }
        else{
          Cov[q,p]=0
        }
      }
      ##user serviced by PPP UAV
      else{
        
        ##get probability of interfering cluster UAV being LOS
        if(type==al){
          lower=r[q]
        }else{lower=rl}
        
        if((min(rmax,rRad)-lower)<10^(-4)){
          intisLOS=0
          intisLOS2=0
        }
        else{
          foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
          Plos = vector(mode="numeric",length=length(foo))
          for(i in 1:length(foo)){
            n = floor((foo[i]/1000)*sqrt(alpha*beta))
            a = 0:max(0,n-1)
            Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
          }
      #            intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
          intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
      #    intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
       #   intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2)
          #      intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
          if(is.na(intisLOS)){intisLOS=0}
        }
        
        ##get probability of interfering cluster UAV being NLOS
        if(type==an){
          lower=r[q]
        }else{lower=rn}
        
        if((min(rmax,rRad)-lower)<10^-4){
          intisNLOS=0 
          intisNLOS2 = 0
        }
        else{
          foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
          Plos = vector(mode="numeric",length=length(foo))
          for(i in 1:length(foo)){
            n = floor((foo[i]/1000)*sqrt(alpha*beta))
            a = 0:max(0,n-1)
            Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
          }

      #            intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
          intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
      #    intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
      #    intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2)
          #       intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
          if(is.na(intisNLOS)){intisNLOS=0}
        }
        
        
        
        
        
        
        for (i in 0:(ms-1)){
          L= 0
            L = strongestLaplaceClusterIntv3(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s,clusterType=al,clusterRad=rRad,rangeProb=min(1,intisLOS+intisNLOS),clusterLOSProb=intisLOS,clusterNLOSProb=intisNLOS)
            Cov[q,p]= Cov[q,p] + ((s^i)/(factorial(i))*((-1)^i)*Re(L))
        }
      }
      if(is.na(Cov[q,p])){Cov[q,p]=0} 
      if(Cov[q,p]==Inf){Cov[q,p]=0}
    }
    
 #   cat(Cov[q,1]," ",Cov[q,2]," ",Cov[q,3]," ",Cov[q,4],"\n")
    
    n = floor((r[q]/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    
    # PPPLOS = 2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*(1-(Prob0nc+cilP))
    # PPPNLOS = 2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*(1-(cinP+Prob0lc))
    Cov[q,1] = Cov[q,1]*2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP)))
      cat(Cov[q,1]*2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP)))," ")
    Cov[q,2] = Cov[q,2]*2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc)))
       cat(Cov[q,2]*2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc)))," ")
    if(r[q]<rRad){
      Cov[q,3]=Cov[q,3]*2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
           cat(Cov[q,3]*2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)," ")
      Cov[q,4]=Cov[q,4]*2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
           cat(Cov[q,4]*2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2), "\n")
      #    CLOS = 2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
      #     CNLOS = 2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
    }else{
      Cov[q,3] = 0
           cat(0, " ")
      Cov[q,4] = 0
          cat(0, "\n")
      #      CLOS = 0
      #      CNLOS = 0
    }
    #   cat(Cov[q,1]," ",Cov[q,2]," ",Cov[q,3]," ",Cov[q,4],"\n")
    Coverage[q] = sum(Cov[q,])*(r[2]-r[1])
  }
  
  return(sum(Coverage))
}

MaternCFCov3 = function(angle,rlim,density,T,gain,rRad,rxh,txh,al,an,mal,man,alpha,beta,gamma,dx){
  Tlin = 10^(T/10)
  rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
 # r = seq(from=min(rmax,rRad)/dx,to=min(rmax,rRad),by=min(rmax,rRad)/dx)
#  r = seq(from=1.5*min(rmax,rRad)/dx,to=1.5*min(rmax,rRad),by=1.5*min(rmax,rRad)/dx)
  r = seq(from=rlim/dx,to=rlim-0.1,by=(rlim-0.1)/dx)
 # r = seq(from=rmax/dx,to=rmax,by=rmax/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=2)
  Coverage = vector(length=length(r))
  for (q in 1:length(r)){
    if(q%%10==0){
      print(q)
    }
    ds = (r[q]^2+(rxh-txh)^2)
    
    rn = sqrt(max(0,(r[q]^2+(rxh-txh)^2)^(al/an) - (rxh-txh)^2))
    rl = min(rmax,sqrt((r[q]^2+(rxh-txh)^2)^(an/al) - (rxh-txh)^2))
    
    foo = seq(from=0,to=r[q],by=r[q]/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    ilP = 2*pi*density*sum(Plos*foo*(foo[2]-foo[1]))
    inP = 2*pi*density*sum((1-Plos)*foo*(foo[2]-foo[1]))
    
    foo = seq(from=0,to=min(r[q],rRad),by=min(r[q],rRad)/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    cilP = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
    cinP = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
    
    cilP = max(0,min(1,cilP))
    cinP = max(0,min(1,cinP))
    
    if(rn==0){
      Prob0n = 0
      Prob0nc = 0
    }
    else{
      foo = seq(from=0,to=rn,by=rn/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0n = 2*pi*density*sum((1-Ptemp)*foo*(foo[2]-foo[1]))
      Prob0n = max(0,min(1,Prob0n))
      
      foo = seq(from=0,to=min(rn,rRad),by=min(rn,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      
      Prob0nc = (2*sum((1-Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2))
      Prob0nc = max(0,min(1,Prob0nc))
    }
    
    #  rl = min(rmax[j],(r[q][k]^2+h^2)^(an/al) - h^2)
    if(rl==0){
      Prob0l = 0
      Prob0lc = 0
    }
    else{
      foo = seq(from=0,to=rl,by=rl/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0l = 2*pi*density*sum(Ptemp*foo*(foo[2]-foo[1]))
      Prob0l = max(0,min(1,Prob0l))
      
      foo = seq(from=0,to=min(rl,rRad),by=min(rl,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0lc = 2*sum((Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2)
      Prob0lc = max(0,min(1,Prob0lc))
    } 
    
    
    ##1-LOS
    ##2-NLOS
   
    for (p in 1:2){
      if(p==1){
        s = mal*Tlin*ds^(al/2)
        ms=mal
        type=al
        
        if(r[q]<rRad){
        isClusterProb = 2*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
        isClusterProb = isClusterProb/(2*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)+2*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP))))
        }
        else{isClusterProb=0}
        probOutRange = max(0,(1-(Prob0nc+cilP)))
      }
      else{
        s = man*Tlin*ds^(an/2) 
        ms=man
        type=an
      
        if(r[q]<rRad){
          isClusterProb = 2*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
          isClusterProb = isClusterProb/(2*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)+2*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc))))
        }
        else{isClusterProb=0}
        probOutRange = max(0,(1-(cinP+Prob0lc)))
      }

        
        ##get probability of interfering cluster UAV being LOS
        if(type==al){
          lower=r[q]
        }else{lower=rl}
        
        if((min(rmax,rRad)-lower)<10^(-4)){
          intisLOS=0
          intisLOS2=0
        }
        else{
   #       foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
          foo = seq(from=lower,to=rRad,by=(rRad-lower)/200)
          Plos = vector(mode="numeric",length=length(foo))
          for(i in 1:length(foo)){
            n = floor((foo[i]/1000)*sqrt(alpha*beta))
            a = 0:max(0,n-1)
            Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
          }
            intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
    #      intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
     #         intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
          #   intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2)
          #      intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
          if(is.na(intisLOS)){intisLOS=0}
        }
        
        ##get probability of interfering cluster UAV being NLOS
        if(type==an){
          lower=r[q]
        }else{lower=rn}
        
        if((min(rmax,rRad)-lower)<10^-4){
          intisNLOS=0 
          intisNLOS2 = 0
        }
        else{
       #   foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
          foo = seq(from=lower,to=rRad,by=(rRad-lower)/200)
          Plos = vector(mode="numeric",length=length(foo))
          for(i in 1:length(foo)){
            n = floor((foo[i]/1000)*sqrt(alpha*beta))
            a = 0:max(0,n-1)
            Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
          }
          
           intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
     #     intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
              intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
          #    intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2)
          #       intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
          if(is.na(intisNLOS)){intisNLOS=0}
        }
        for (i in 0:(ms-1)){
          L= 0
          L = strongestLaplaceClusterIntv4(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s,clusterType=al,clusterRad=rRad,rangeProb=probOutRange,clusterLOSProb=intisLOS,clusterNLOSProb=intisNLOS,isClusterProb=isClusterProb)
          Cov[q,p]= Cov[q,p] + ((s^i)/(factorial(i))*((-1)^i)*Re(L))
        }
      if(is.na(Cov[q,p])){Cov[q,p]=0} 
      if(Cov[q,p]==Inf){Cov[q,p]=0}
    }
    
    #   cat(Cov[q,1]," ",Cov[q,2]," ",Cov[q,3]," ",Cov[q,4],"\n")
    
    n = floor((r[q]/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    
    # PPPLOS = 2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*(1-(Prob0nc+cilP))
    # PPPNLOS = 2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*(1-(cinP+Prob0lc))
    PPPlos = 2*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP)))
    PPPnlos = 2*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc)))
    if(r[q]<rRad){
   clusterLOS = 2*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
    clusterNLOS = 2*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
      #    CLOS = 2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
      #     CNLOS = 2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
    }else{
      clusterLOS = 0
      clusterNLOS = 0
      #      CLOS = 0
      #      CNLOS = 0
    }
    #   cat(Cov[q,1]," ",Cov[q,2]," ",Cov[q,3]," ",Cov[q,4],"\n")
    Coverage[q] = P*Cov[q,1]*(PPPlos+clusterLOS)+(1-P)*Cov[q,2]*(PPPnlos+clusterNLOS)
  }
  
  return(sum(Coverage*(r[2]-r[1])))
}

MCMaternCFCov = function(angle,rlim,density,T,gain,rRad,rxh,txh,al,an,mal,man,alpha,beta,gamma,dx,cores){
  Tlin = 10^(T/10)
  rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  r = seq(from=rlim/dx,to=rlim,by=rlim/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  Cov = zeros(nrow=length(r),ncol=4)
 # Coverage = vector(length=length(r))
  iteration=function(q){
    if(q%%10==0){
      cat(q,"\n")
    }
    ds = (r[q]^2+(rxh-txh)^2)
    
    rn = sqrt(max(0,(r[q]^2+(rxh-txh)^2)^(al/an) - (rxh-txh)^2))
    rl = min(rmax,sqrt((r[q]^2+(rxh-txh)^2)^(an/al) - (rxh-txh)^2))
    
    foo = seq(from=0,to=r[q],by=r[q]/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    ilP = 2*pi*density*sum(Plos*foo*(foo[2]-foo[1]))
    inP = 2*pi*density*sum((1-Plos)*foo*(foo[2]-foo[1]))
    
    foo = seq(from=0,to=min(r[q],rRad),by=min(r[q],rRad)/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    cilP = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
    cinP = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
    
    if(rn==0){
      Prob0n = 0
      Prob0nc = 0
    }
    else{
      foo = seq(from=0,to=rn,by=rn/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0n = 2*pi*density*sum((1-Ptemp)*foo*(foo[2]-foo[1]))
      
      foo = seq(from=0,to=min(rn,rRad),by=min(rn,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      
      Prob0nc = (2*sum((1-Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2))
    }
    
    #  rl = min(rmax[j],(r[q][k]^2+h^2)^(an/al) - h^2)
    if(rl==0){
      Prob0l = 0
      Prob0lc = 0
    }
    else{
      foo = seq(from=0,to=rl,by=rl/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0l = 2*pi*density*sum(Ptemp*foo*(foo[2]-foo[1]))
      
      foo = seq(from=0,to=min(rl,rRad),by=min(rl,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0lc = 2*sum((Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2)
    } 
    
    
    ##1- PPP LOS
    ##2- PPP NLOS
    ##3- cluster center LOS
    ##4- cluster center NLOS
    for (p in 1:4){
      if(p==1 | p==3){
        s = mal*Tlin*ds^(al/2)
        ms=mal
        type=al
      }
      else{
        s = man*Tlin*ds^(an/2) 
        ms=man
        type=an
      }
      Cov[1,p] = 0
      if(p == 3 | p==4){
        if(r[q]<rRad){
          for (i in 0:(ms-1)){
            L = strongestLaplaceLOSSINR(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s)
            Cov[1,p]= Cov[1,p] + (s^i)/(factorial(i))*((-1)^i)*Re(L)
          }
        }
        else{
          Cov[1,p]=0
        }
      }
      ##user serviced by PPP UAV
      else{
        
        ##get probability of interfering cluster UAV being LOS
        if(type==al){
          lower=r[q]
        }else{lower=rl}
        
        if((min(rmax,rRad)-lower)<10^(-4)){
          intisLOS=0 
        }
        else{
          foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
          Plos = vector(mode="numeric",length=length(foo))
          for(i in 1:length(foo)){
            n = floor((foo[i]/1000)*sqrt(alpha*beta))
            a = 0:max(0,n-1)
            Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
          }
  #        intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
                  intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
          #        intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/((rRad^2)*max(0,(1-(Prob0nc+cilP))))
          #      intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/((rRad^2-lower^2)*max(0,(1-(Prob0nc+cilP))))
                 intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
       #   intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
          if(is.na(intisLOS)){intisLOS=0}
        }
        
        ##get probability of interfering cluster UAV being NLOS
        if(type==an){
          lower=r[q]
        }else{lower=rn}
        
        if((min(rmax,rRad)-lower)<10^-4){
          intisNLOS=0 
        }
        else{
          foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
          Plos = vector(mode="numeric",length=length(foo))
          for(i in 1:length(foo)){
            n = floor((foo[i]/1000)*sqrt(alpha*beta))
            a = 0:max(0,n-1)
            Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
          }
          #        intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
          intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
          #      intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/((rRad^2)*max(0,(1-(cinP+Prob0lc))))
          #        intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/((rRad^2-lower^2)*max(0,(1-(cinP+Prob0lc))))
                  intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
       #   intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
          if(is.na(intisNLOS)){intisNLOS=0}
        }
        
        #     if(rRad<rmax){
        #        intis0 = 0
        #      }else{
        #      intis0 = (rRad^2-rmax^2)/(rRad^2)
        #  intis0 = (rRad^2-rmax^2)/((rRad^2)*(1-(Prob0nc+cilP)))
        #      }
        
        intisLOS = max(0,min(1,intisLOS))
        intisNLOS = max(0,min(1,intisNLOS))
        
        intis0 = max(0,(1-(intisLOS+intisNLOS)))
        
        #  if(p==1){
        #        print(intisLOS+intisNLOS)
        # }
        for (i in 0:(ms-1)){
          L= 0
          if(intisLOS>0){
            L = strongestLaplaceClusterInt(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s,clusterType=al,clusterRad=rRad,clusterTypeProb=intisLOS2)
            Cov[1,p]= Cov[1,p] + intisLOS*((s^i)/(factorial(i))*((-1)^i)*Re(L))
          }
          if(intisNLOS>0){
            L = strongestLaplaceClusterInt(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s,clusterType=an,clusterRad=rRad,clusterTypeProb=intisNLOS2)
            Cov[1,p]= Cov[1,p] + intisNLOS*((s^i)/(factorial(i))*((-1)^i)*Re(L))
          }
          if(intis0>0){
            L = strongestLaplaceLOSSINR(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s)
            Cov[1,p]= Cov[1,p] + intis0*((s^i)/(factorial(i))*((-1)^i)*Re(L))
          }
        }
      }
    }
    
    
    n = floor((r[q]/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    
    # PPPLOS = 2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*(1-(Prob0nc+cilP))
    # PPPNLOS = 2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*(1-(cinP+Prob0lc))
    Cov[1,1] = Cov[1,1]*2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP)))
    #   cat(2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP)))/(r[2]-r[1])," ")
    Cov[1,2] = Cov[1,2]*2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc)))
    #   cat(2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc)))/(r[2]-r[1])," ")
    if(r[q]<rRad){
      Cov[1,3]=Cov[1,3]*2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
      #     cat(2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)/(r[2]-r[1])," ")
      Cov[1,4]=Cov[1,4]*2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
      #     cat(2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)/(r[2]-r[1]), "\n")
      #    CLOS = 2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
      #     CNLOS = 2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
    }else{
      Cov[1,3] = 0
      #     cat(0, " ")
      Cov[1,4] = 0
      #    cat(0, "\n")
      #      CLOS = 0
      #      CNLOS = 0
    }
    #   cat(Cov[q,1]," ",Cov[q,2]," ",Cov[q,3]," ",Cov[q,4],"\n")
    Coverage = sum(Cov[1,])*(r[2]-r[1])
    return(Coverage)
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  covP = 0
  for(k in 1:length(r)){
    covP = covP+c[[k]]
  }
  
  
  return(covP)
 
}

MCMaternCFCov3 = function(angle,rlim,density,T,gain,rRad,rxh,txh,al,an,mal,man,alpha,beta,gamma,dx,cores){
  Tlin = 10^(T/10)
  rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  # r = seq(from=min(rmax,rRad)/dx,to=min(rmax,rRad),by=min(rmax,rRad)/dx)
  #  r = seq(from=1.5*min(rmax,rRad)/dx,to=1.5*min(rmax,rRad),by=1.5*min(rmax,rRad)/dx)
#  r = seq(from=min(rRad*3,rmax)/dx,to=min(rRad*3,rmax),by=min(rRad*3,rmax)/dx)
  # r = seq(from=rmax/dx,to=rmax,by=rmax/dx)
  r = seq(from=rlim/dx,to=rlim-0.1,by=(rlim-0.1)/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
 # Cov = zeros(nrow=length(r),ncol=2)
 # Coverage = vector(length=length(r))
  iteration=function(q){
    Cov = zeros(nrow=1,ncol=2)
    if(q%%10==0){
      print(q)
    }
    ds = (r[q]^2+(rxh-txh)^2)
    
    rn = sqrt(max(0,(r[q]^2+(rxh-txh)^2)^(al/an) - (rxh-txh)^2))
    rl = min(rmax,sqrt((r[q]^2+(rxh-txh)^2)^(an/al) - (rxh-txh)^2))
    
    foo = seq(from=0,to=r[q],by=r[q]/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    ilP = 2*pi*density*sum(Plos*foo*(foo[2]-foo[1]))
    inP = 2*pi*density*sum((1-Plos)*foo*(foo[2]-foo[1]))
    
    foo = seq(from=0,to=min(r[q],rRad),by=min(r[q],rRad)/200)
    Plos = vector(mode="numeric",length=length(foo))
    for(i in 1:length(foo)){
      n = floor((foo[i]/1000)*sqrt(alpha*beta))
      a = 0:max(0,n-1)
      Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
    }
    
    cilP = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
    cinP = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
    
    cilP = max(0,min(1,cilP))
    cinP = max(0,min(1,cinP))
    
    if(rn==0){
      Prob0n = 0
      Prob0nc = 0
    }
    else{
      foo = seq(from=0,to=rn,by=rn/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0n = 2*pi*density*sum((1-Ptemp)*foo*(foo[2]-foo[1]))
      Prob0n = max(0,min(1,Prob0n))
      
      foo = seq(from=0,to=min(rn,rRad),by=min(rn,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      
      Prob0nc = (2*sum((1-Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2))
      Prob0nc = max(0,min(1,Prob0nc))
    }
    
    #  rl = min(rmax[j],(r[q][k]^2+h^2)^(an/al) - h^2)
    if(rl==0){
      Prob0l = 0
      Prob0lc = 0
    }
    else{
      foo = seq(from=0,to=rl,by=rl/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0l = 2*pi*density*sum(Ptemp*foo*(foo[2]-foo[1]))
      Prob0l = max(0,min(1,Prob0l))
      
      foo = seq(from=0,to=min(rl,rRad),by=min(rl,rRad)/200)
      Ptemp = vector(mode="numeric",length=length(foo))
      for(i in 1:length(foo)){
        n = floor((foo[i]/1000)*sqrt(alpha*beta))
        a = 0:max(0,n-1)
        Ptemp[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
      }
      Prob0lc = 2*sum((Ptemp)*foo*(foo[2]-foo[1]))/(rRad^2)
      Prob0lc = max(0,min(1,Prob0lc))
    } 
    
    
    ##1-LOS
    ##2-NLOS
    
    for (p in 1:2){
      if(p==1){
        s = mal*Tlin*ds^(al/2)
        ms=mal
        type=al
        
        if(r[q]<rRad){
          isClusterProb = 2*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
          isClusterProb = isClusterProb/(2*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)+2*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP))))
        }
        else{isClusterProb=0}
      probOutRange = max(0,(1-(Prob0nc+cilP)))
      }
      else{
        s = man*Tlin*ds^(an/2) 
        ms=man
        type=an
        
        if(r[q]<rRad){
          isClusterProb = 2*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
          isClusterProb = isClusterProb/(2*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)+2*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc))))
        }
        else{isClusterProb=0}
        probOutRange = max(0,(1-(cinP+Prob0lc)))
      }
      
      
      ##get probability of interfering cluster UAV being LOS
      if(type==al){
        lower=r[q]
      }else{lower=rl}
      
      if((min(rmax,rRad)-lower)<10^(-4)){
        intisLOS=0
        intisLOS2=0
      }
      else{
        #       foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
        foo = seq(from=lower,to=rRad,by=(rRad-lower)/200)
        Plos = vector(mode="numeric",length=length(foo))
        for(i in 1:length(foo)){
          n = floor((foo[i]/1000)*sqrt(alpha*beta))
          a = 0:max(0,n-1)
          Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
        }
        intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2)
        #      intisLOS = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
        #         intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
        #   intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2)
        #      intisLOS2 = 2*sum(Plos*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
        if(is.na(intisLOS)){intisLOS=0}
      }
      
      ##get probability of interfering cluster UAV being NLOS
      if(type==an){
        lower=r[q]
      }else{lower=rn}
      
      if((min(rmax,rRad)-lower)<10^-4){
        intisNLOS=0 
        intisNLOS2 = 0
      }
      else{
        #   foo = seq(from=lower,to=min(rmax,rRad),by=(min(rmax,rRad)-lower)/200)
        foo = seq(from=lower,to=rRad,by=(rRad-lower)/200)
        Plos = vector(mode="numeric",length=length(foo))
        for(i in 1:length(foo)){
          n = floor((foo[i]/1000)*sqrt(alpha*beta))
          a = 0:max(0,n-1)
          Plos[i] = prod(1-exp(-((max(rxh,txh)-(a+1/2)*(abs(rxh-txh))/(n))^2)/(2*gamma^2)))
        }
        
        intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2)
        #     intisNLOS = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
        intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2-lower^2)
        #    intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(min(rmax,rRad)^2)
        #       intisNLOS2 = 2*sum((1-Plos)*foo*(foo[2]-foo[1]))/(rRad^2-lower^2)
        if(is.na(intisNLOS)){intisNLOS=0}
      }

      for (i in 0:(ms-1)){
        L= 0
        L = strongestLaplaceClusterIntv4(max=rmax,min=r[q],type=type,N=N/gain,al=al,an=an,alpha=alpha,beta=beta,gamma=gamma,rxheight=rxh,txheight=txh,mal=mal,man=man,l=i,const=-pi*density,s=s,clusterType=al,clusterRad=rRad,rangeProb=probOutRange,clusterLOSProb=intisLOS,clusterNLOSProb=intisNLOS,isClusterProb=isClusterProb)
        Cov[1,p]= Cov[1,p] + ((s^i)/(factorial(i))*((-1)^i)*Re(L))
      }
      if(is.na(Cov[1,p])){Cov[1,p]=0} 
      if(Cov[1,p]==Inf){Cov[1,p]=0}
    }
    
    #   cat(Cov[q,1]," ",Cov[q,2]," ",Cov[q,3]," ",Cov[q,4],"\n")
    
    n = floor((r[q]/1000)*sqrt(alpha*beta))
    a = 0:max(0,n-1)
    P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    
    # PPPLOS = 2*P*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*(1-(Prob0nc+cilP))
    # PPPNLOS = 2*(1-P)*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*(1-(cinP+Prob0lc))
    PPPlos = 2*density*pi*r[q]*exp(-ilP)*exp(-Prob0n)*max(0,(1-(Prob0nc+cilP)))
    PPPnlos = 2*density*pi*r[q]*exp(-inP)*exp(-Prob0l)*max(0,(1-(cinP+Prob0lc)))
    if(r[q]<rRad){
      clusterLOS = 2*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
      clusterNLOS = 2*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
      #    CLOS = 2*P*r[q]*exp(-ilP)*exp(-Prob0n)/(rRad^2)
      #     CNLOS = 2*(1-P)*r[q]*exp(-inP)*exp(-Prob0l)/(rRad^2)
    }else{
      clusterLOS = 0
      clusterNLOS = 0
      #      CLOS = 0
      #      CNLOS = 0
    }
    #   cat(Cov[q,1]," ",Cov[q,2]," ",Cov[q,3]," ",Cov[q,4],"\n")
    Coverage = P*Cov[1,1]*(PPPlos+clusterLOS)+(1-P)*Cov[1,2]*(PPPnlos+clusterNLOS)
    return(Coverage*(r[2]-r[1]))
  }
  
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  covP = 0
  for(k in 1:length(r)){
    covP = covP+c[[k]]
  }
  
  
  
  return(covP)
}


#######Mobile UAV Probability calculation#############

handoverProbabilityCond = function(r0,angle,velocity, density){
  
      r1 = sqrt(r0^2+velocity^2-2*r0*velocity*cos(angle))
      angle1 = acos((velocity^2+r1^2-r0^2)/(2*velocity*r1))
      
      if(angle<pi/2 & angle1<pi/2){
        lens0 = ((r0^2)/2)*(angle-sin(angle))
        lens1 = ((r1^2)/2)*(angle1-sin(angle1))
        lensarea = lens0+lens1
        foo = 1-exp(-density*(pi*r1^2-lensarea))
      }else if(angle>=pi/2 & angle1<pi/2){
        lens0 = ((r0^2)/2)*((pi-angle)-sin((pi-angle)))
        lens1 = ((r1^2)/2)*(angle1-sin(angle1))
        lensarea = (pi*r0^2-abs(lens0-lens1))
        foo = 1-exp(-density*(pi*r1^2-lensarea))
      }
      else{
        lens0 = ((r0^2)/2)*(angle-sin(angle))
        lens1 = ((r1^2)/2)*((pi-angle1)-sin((pi-angle1)))
        foo = 1-exp(-density*abs(lens0-lens1))
      }
      
      #  a = min(r0[i],r1)
      #  b= max(r0[i],r1)
      #  c= velocity
      #  delta = (1/4)*sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c))
      #lunearea = 2*delta+(a^2)*acos((b^2-a^2-c^2)/(2*a*c))-(b^2)*acos((b^2+c^2-a^2)/(2*b*c))
      #foo[j] = 1-exp(-density*lunearea)
      #  delta = (1/2)*sqrt((-c+a+b)*(c+a-b)*(c-a+b)*(c+a+b))
      #  lensarea = (a^2)*acos((c^2+a^2-b^2)/(2*c*a))+(b^2)*acos((c^2+b^2-a^2)/(2*c*b))+delta
      # foo[j] = 1-exp(-density*(pi*r1^2-lensarea))
      if(is.na(foo)){foo=0}

  return(foo)
}

#as above, but in a multi-network scenario. Cond probability of handing over to network 1 and not network 2
handoverProbabilityCondMulti = function(r0,angle,velocity, density1, density2){
  
  r1 = sqrt(r0^2+velocity^2-2*r0*velocity*cos(angle))
  angle1 = acos((velocity^2+r1^2-r0^2)/(2*velocity*r1))
  
  if(angle<pi/2 & angle1<pi/2){
    lens0 = ((r0^2)/2)*(angle-sin(angle))
    lens1 = ((r1^2)/2)*(angle1-sin(angle1))
    lensarea = lens0+lens1
    foo = (1-exp(-(density1+density2)*(pi*r1^2-lensarea)))*(density1/(density1+density2))#*(exp(-density2*(pi*r1^2-lensarea)))
  }else if(angle>=pi/2 & angle1<pi/2){
    lens0 = ((r0^2)/2)*((pi-angle)-sin((pi-angle)))
    lens1 = ((r1^2)/2)*(angle1-sin(angle1))
    lensarea = (pi*r0^2-abs(lens0-lens1))
    foo = (1-exp(-(density1+density2)*(pi*r1^2-lensarea)))*(density1/(density1+density2))#*exp(-density2*(pi*r1^2-lensarea))
  }
  else{
    lens0 = ((r0^2)/2)*(angle-sin(angle))
    lens1 = ((r1^2)/2)*((pi-angle1)-sin((pi-angle1)))
    foo = (1-exp(-(density1+density2)*abs(lens0-lens1)))*(density1/(density1+density2))#*exp(-density2*abs(lens0-lens1))
  }
  
  #  a = min(r0[i],r1)
  #  b= max(r0[i],r1)
  #  c= velocity
  #  delta = (1/4)*sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c))
  #lunearea = 2*delta+(a^2)*acos((b^2-a^2-c^2)/(2*a*c))-(b^2)*acos((b^2+c^2-a^2)/(2*b*c))
  #foo[j] = 1-exp(-density*lunearea)
  #  delta = (1/2)*sqrt((-c+a+b)*(c+a-b)*(c-a+b)*(c+a+b))
  #  lensarea = (a^2)*acos((c^2+a^2-b^2)/(2*c*a))+(b^2)*acos((c^2+b^2-a^2)/(2*c*b))+delta
  # foo[j] = 1-exp(-density*(pi*r1^2-lensarea))
  if(is.na(foo[j])){foo[j]=0}
  
  return(foo)
}

handoverProbabilityMulti = function(velocity, density1,density2,dx){
  
  rm = sqrt(-log(0.001)/(pi*density1))
  
  r0 = seq(from=rm/dx,to=rm,by=rm/dx)
  #  r0 = 100
  angle = seq(from=(pi)/dx,to=(pi*(dx-1)/dx),by=(pi)/dx)
  
  handoverProb = vector(length=length(r0))
  for(i in 1:length(r0)){
    foo = vector(length=length(angle))
    for(j in 1:length(angle)){
    foo[j] = handoverProbabilityCondMulti(r0=r0[i],angle=angle[j],velocity=velocity, density1=density1,density2=density2)
    }
    handoverProb[i]= sum(foo*(angle[2]-angle[1])/(pi))*2*pi*density1*r0[i]*exp(-pi*density1*r0[i]^2)*(r0[2]-r0[1]) 
  }
  return(sum(handoverProb))
}


handoverProbability = function(velocity,density,dx){
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r0 = seq(from=rm/dx,to=rm,by=rm/dx)
  #  r0 = 100
  angle = seq(from=(pi)/dx,to=(pi*(dx-1)/dx),by=(pi)/dx)
  
  handoverProb = vector(length=length(r0))
  for(i in 1:length(r0)){
    foo = vector(length=length(angle))
    for(j in 1:length(angle)){
      foo[j] = handoverProbabilityCond(r0=r0[i],angle=angle[j],velocity=velocity, density=density)
    }
    handoverProb[i]= sum(foo*(angle[2]-angle[1])/(pi))*2*pi*density*r0[i]*exp(-pi*density*r0[i]^2)*(r0[2]-r0[1]) 
  }
  return(sum(handoverProb))
}

#probability that a handover will happen to a BS which is in the current BSs neighbour list (ie, is one of the B closest BSs to it)
NeighbourHandoverProbabilityCond  = function(r0,angle, velocity, B, density){
  
  
  
  r1 = sqrt(r0^2+velocity^2-2*r0*velocity*cos(angle))
  angle1 = acos((velocity^2+r1^2-r0^2)/(2*velocity*r1))
  
  if(r0>r1/2){
  R = max(r0,r1)
  r = min(r0,r1)
  d=r0
  delta = (1/4)*sqrt((-d+r+R)*(d-r+R)*(d+r-R)*(d+r+R))
  
  area = (r^2)*acos((d^2+r^2-R^2)/(2*d*r))+(R^2)*acos((d^2+R^2-r^2)/(2*d*R))-2*delta
  
  area = pi*r1^2 - area
  }
  else{
    area = pi*r1^2 - pi*r0^2
  }
  
#  P1 = 1-exp(-density*area)
  n = 1:B
  P2 = sum(exp(-density*area)*((density*area)^n)/(factorial(n)))
  
  if(angle<pi/2 & angle1<pi/2){
    lens0 = ((r0^2)/2)*(angle-sin(angle))
    lens1 = ((r1^2)/2)*(angle1-sin(angle1))
    lensarea = lens0+lens1
    foo = (1-exp(-(density)*(pi*r1^2-lensarea)))*(P2)#*(exp(-density2*(pi*r1^2-lensarea)))
  }else if(angle>=pi/2 & angle1<pi/2){
    lens0 = ((r0^2)/2)*((pi-angle)-sin((pi-angle)))
    lens1 = ((r1^2)/2)*(angle1-sin(angle1))
    lensarea = (pi*r0^2-abs(lens0-lens1))
    foo = (1-exp(-(density)*(pi*r1^2-lensarea)))*(P2)#*exp(-density2*(pi*r1^2-lensarea))
  }
  else{
    lens0 = ((r0^2)/2)*(angle-sin(angle))
    lens1 = ((r1^2)/2)*((pi-angle1)-sin((pi-angle1)))
    foo = (1-exp(-(density)*abs(lens0-lens1)))*(P2)#*exp(-density2*abs(lens0-lens1))
  }
  
  if(is.na(foo)){foo=0}
  
  return(foo)
}

NeighbourHandoverProbability = function(velocity,B,density,dx){
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r0 = seq(from=rm/dx,to=rm,by=rm/dx)
  #  r0 = 100
  angle = seq(from=(pi)/dx,to=(pi*(dx-1)/dx),by=(pi)/dx)
  
  neighbourhandoverProb = vector(length=length(r0))
  for(i in 1:length(r0)){
    foo = vector(length=length(angle))
    for(j in 1:length(angle)){
  #    r1 = sqrt(r0[i]^2+velocity^2-2*r0[i]*velocity*cos(angle[j]))
      foo[j] = NeighbourHandoverProbabilityCond(r0=r0[i],angle=angle[j],velocity=velocity, B=B,density=density)
    }
    neighbourhandoverProb[i]= sum(foo*(angle[2]-angle[1])/(pi))*2*pi*density*r0[i]*exp(-pi*density*r0[i]^2)*(r0[2]-r0[1]) 
  }
  return(sum(neighbourhandoverProb))
}


DualConnectHandoverProbabilityCond = function(r0,r1,angle0,angle1,velocity, density){
  
  r2 = sqrt(r0^2+velocity^2-2*r0*velocity*cos(angle0))
  r3 = sqrt(r1^2+velocity^2-2*r1*velocity*cos(angle1))
  
  if(r2>r3){
  r1=r2
  angle= angle0
  }else{
  r0=r1
  r1=r3
  angle=angle1
  }
  
  angle1 = acos((velocity^2+r1^2-r0^2)/(2*velocity*r1))
  
  if(angle<pi/2 & angle1<pi/2){
    lens0 = ((r0^2)/2)*(angle-sin(angle))
    lens1 = ((r1^2)/2)*(angle1-sin(angle1))
    lensarea = lens0+lens1
    foo = 1-exp(-density*(pi*r1^2-lensarea))
  }else if(angle>=pi/2 & angle1<pi/2){
    lens0 = ((r0^2)/2)*((pi-angle)-sin((pi-angle)))
    lens1 = ((r1^2)/2)*(angle1-sin(angle1))
    lensarea = (pi*r0^2-abs(lens0-lens1))
    foo = 1-exp(-density*(pi*r1^2-lensarea))
  }
  else{
    lens0 = ((r0^2)/2)*(angle-sin(angle))
    lens1 = ((r1^2)/2)*((pi-angle1)-sin((pi-angle1)))
    foo = 1-exp(-density*abs(lens0-lens1))
  }
  
  #  a = min(r0[i],r1)
  #  b= max(r0[i],r1)
  #  c= velocity
  #  delta = (1/4)*sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c))
  #lunearea = 2*delta+(a^2)*acos((b^2-a^2-c^2)/(2*a*c))-(b^2)*acos((b^2+c^2-a^2)/(2*b*c))
  #foo[j] = 1-exp(-density*lunearea)
  #  delta = (1/2)*sqrt((-c+a+b)*(c+a-b)*(c-a+b)*(c+a+b))
  #  lensarea = (a^2)*acos((c^2+a^2-b^2)/(2*c*a))+(b^2)*acos((c^2+b^2-a^2)/(2*c*b))+delta
  # foo[j] = 1-exp(-density*(pi*r1^2-lensarea))
  if(is.na(foo[j])){foo[j]=0}
  
  return(foo)
}

DualConnectHandoverProbability = function(velocity,density,dx){
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r0 = seq(from=rm/dx,to=rm,by=rm/dx)
  #  r0 = 100
  angle0 = seq(from=(pi)/dx,to=(pi*(dx-1)/dx),by=(pi)/dx)
  angle1= angle0
  
  handoverProb = vector(length=length(r0))
  for(i in 1:length(r0)){
    r1 = seq(from=r0[i]+1,to=rm*2,by=(rm*2-r0[i]-1)/dx)
    foo1 = vector(length=length(r1))
    for(j in 1:length(r1)){
      fooA = matrix(nrow=length(angle0),ncol=length(angle1))
      for(k in 1:length(angle0)){
        for(l in 1:length(angle1)){
      fooA[k,l] = DualConnectHandoverProbabilityCondMulti(r0=r0[i],r1=r1[j],angle0=angle0[k],angle1=angle1[l],velocity=velocity,density=density)
      }
      }
    foo1 = sum(fooA)*((2*(pi*density)^2)/factorial(1))*(r1[j]^3)*exp(-pi*density*r1[j]^2)*(r1[2]-r1[1]) 
    }
    handoverProb[i]= sum(foo*((angle0[2]-angle0[1])/(pi))^2)*2*pi*density*r0[i]*exp(-pi*density*r0[i]^2)*(r0[2]-r0[1]) 
  }
  return(sum(handoverProb))
}

MCDualConnectHandoverProbability = function(velocity,density,dx,cores){
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r0 = seq(from=rm/dx,to=rm,by=rm/dx)
  #  r0 = 100
  angle0 = seq(from=(pi)/dx,to=(pi*(dx-1)/dx),by=(pi)/dx)
  angle1= angle0
  
  handoverProb = vector(length=length(r0))
  iteration = function(i){
  #for(i in 1:length(r0)){
    r1 = seq(from=r0[i]+1,to=rm*2,by=(rm*2-r0[i]-1)/dx)
    foo1 = vector(length=length(r1))
    for(j in 1:length(r1)){
      fooA = matrix(nrow=length(angle0),ncol=length(angle1))
      for(k in 1:length(angle0)){
        for(l in 1:length(angle1)){
          fooA[k,l] = DualConnectHandoverProbabilityCondMulti(r0=r0[i],r1=r1[j],angle0=angle0[k],angle1=angle1[l],velocity=velocity,density=density)
        }
      }
      foo1 = sum(fooA)*((2*(pi*density)^2)/factorial(1))*(r1[j]^3)*exp(-pi*density*r1[j]^2)*(r1[2]-r1[1]) 
    }
     return(sum(foo*((angle0[2]-angle0[1])/(pi))^2)*2*pi*density*r0[i]*exp(-pi*density*r0[i]^2)*(r0[2]-r0[1])) 
  }
  
  X=1:length(r0)
  opt = mclapply(X=X,FUN=iteration,cores=cores)
  
  
  for(k in 1:length(r0)){
    handoverProb[k]=opt[[k]][1]
  }
  return(sum(handoverProb))
}



##Coverage probability, when taking into account mobility and handover
MCMobileCovProb = function(rxh,txh,angle,rmax,density,velocity,handPenalty,uptilt,T,gain,horgain,al,an,alpha,beta,gamma,dx,cores){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  iteration=function(q){
    Cov = zeros(nrow=length(r),ncol=2)
    HCov = zeros(nrow=length(r),ncol=2)
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      HCov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      HCov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    HCov[q,1] = HCov[q,1]*P
    HCov[q,2] = HCov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    horangle = seq(from=(pi)/dx,to=(pi*(dx-1)/dx),by=(pi)/dx)
    foo= vector(length=length(horangle))
    for(i in 1:length(horangle)){
    foo[i]=(HCov[q,1]+HCov[q,2])*(1-handoverProbabilityCond(r0=r[q],angle=horangle[i],velocity=velocity, density=density))*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])*(horangle[2]-horangle[1])  
    }
    HCov[q,1] = sum(foo)/pi#(Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    
   return(c(Cov[q,1],HCov[q,1]))
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  Cov = 0
  HCov = 0
  for(k in 1:length(r)){
    Cov = Cov+c[[k]][1]
    HCov = HCov + c[[k]][2]
  }
  
  return((1-handPenalty)*sum(Cov)+(handPenalty)*sum(HCov))
}


####Cov Prob functions using Ramy's antenna model####
##as above, but using Ramy's BS antenna model for the radiation pattern
MCMobileCovProbRamy = function(rxh,txh,angle,rmax,density,velocity,handPenalty,uptilt,T,gain,Nt,al,an,alpha,beta,gamma,dx,cores){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  iteration=function(q){
    Cov = zeros(nrow=length(r),ncol=2)
    HCov = zeros(nrow=length(r),ncol=2)
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])
      #gn=(1/Nt)*((sin(Nt*pi*(sin(ang))/2)^2)/(sin(pi*(sin(ang))/2)^2))
      gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2))
      
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])
       # gn=(1/Nt)*((sin(Nt*pi*(sin(ang))/2)^2)/(sin(pi*(sin(ang))/2)^2))
        gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2))
        
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])
   #   gn=(1/Nt)*((sin(Nt*pi*(sin(ang))/2)^2)/(sin(pi*(sin(ang))/2)^2))
      gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2))
      
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      HCov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])
     #   gn=(1/Nt)*((sin(Nt*pi*(sin(ang))/2)^2)/(sin(pi*(sin(ang))/2)^2))
        gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      HCov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    HCov[q,1] = HCov[q,1]*P
    HCov[q,2] = HCov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    horangle = seq(from=(pi)/dx,to=(pi*(dx-1)/dx),by=(pi)/dx)
    foo= vector(length=length(horangle))
    for(i in 1:length(horangle)){
      foo[i]=(HCov[q,1]+HCov[q,2])*(1-handoverProbabilityCond(r0=r[q],angle=horangle[i],velocity=velocity, density=density))*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])*(horangle[2]-horangle[1])  
    }
    HCov[q,1] = sum(foo)/pi#(Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    
    return(c(Cov[q,1],HCov[q,1]))
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  Cov = 0
  HCov = 0
  for(k in 1:length(r)){
    Cov = Cov+c[[k]][1]
    HCov = HCov + c[[k]][2]
  }
  
  return((1-handPenalty)*sum(Cov)+(handPenalty)*sum(HCov))
}


##Coverage probability using Ramy's BS antenna model, and the Nakagami CDF approximation
MCCovProbRamyNakApprox = function(rxh,txh,angle,rmax,density,uptilt,T,gain,Nt,al,an,mal,man,alpha,beta,gamma,dx,cores){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  iteration=function(q){
    Cov = zeros(nrow=1,ncol=2)
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
  
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])
#      gn=(1/Nt)*((sin(Nt*pi*(sin(ang))/2)^2)/(sin(pi*(sin(ang))/2)^2))
      gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2))
      
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
        m=mal
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2)
        m=man
      }
      Cov[1,p] = 0
      
      bx= 0
      if(m==1){
        bx=1  
      }else if(m==2){
        bx=1.487
      }else if(m==3){
        bx = 1.81
      }else if(m==10){
        bx = 2.872
      }
      
      for(j in 1:m){
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])
      #  gn=(1/Nt)*((sin(Nt*pi*(sin(ang))/2)^2)/(sin(pi*(sin(ang))/2)^2))
        gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2))
        g = bx*j*s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = bx*j*s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[1,p] = Cov[1,p]+choose(m,j)*((-1)^(j+1))*exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-bx*j*s*N/gain)
      }
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[1,1] = Cov[1,1]*P
    Cov[1,2] = Cov[1,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[1,1] = (Cov[1,1]+Cov[1,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    
    return(Cov[1,1])
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  Cov = 0
  for(k in 1:length(r)){
    Cov = Cov+c[[k]]
  }
  
  return(Cov)
}


##as above, but using Ramy's BS antenna model for the radiation pattern
MCCovProbRamyCond = function(d,phi,rxh,txh,angle,rmax,density,velocity,handPenalty,uptilt,T,gain,Nt,al,an,alpha,beta,gamma,dx,cores){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  rd = sqrt(r^2+d^2-2*r*d*cos(phi))
  #  y = 2*density*pi*r*exp(-density*pi*r^2)*(r[2]-r[1])
  
  iteration=function(q){
    Cov = zeros(nrow=length(r),ncol=2)
    HCov = zeros(nrow=length(r),ncol=2)
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),rd[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (rd[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),rd[q])*(180/pi)
      gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*rd[q])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*rd[q])))/2)^2))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=rd[q],to=threshold,by=(threshold-rd[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[q,p] = exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    
    
    #  n = floor((rd[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=rd[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    
    
    return(Cov[q,1])
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  Cov = 0
  for(k in 1:length(r)){
    Cov = Cov+c[[k]][1]
  }
  
  return(sum(Cov))
}


#as above, but in a two-network scenario.
MCMobileCovProbMulti = function(rxh,txh,angle,rmax,density1,density2,velocity,handPenalty1,handPenalty2,uptilt,T,gain,horgain,al,an,alpha,beta,gamma,dx,cores){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  densityTotal=(density1+density2)
  rm = sqrt(-log(0.001)/(pi*densityTotal))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*densityTotal*pi*r*exp(-densityTotal*pi*r^2)*(r[2]-r[1])
  
  iteration=function(q){
    Cov = zeros(nrow=length(r),ncol=2)
    noHCov = zeros(nrow=length(r),ncol=2)
    HsCov = zeros(nrow=length(r),ncol=2)
    HhCov = zeros(nrow=length(r),ncol=2)
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[q,p] = exp(-density1*angle*sum(Ll))*exp(-density1*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*densityTotal*r[q]*exp(-pi*densityTotal*r[q]^2)*(r[2]-r[1])
    
    horangle = seq(from=(pi)/dx,to=(pi*(dx-1)/dx),by=(pi)/dx)
    ProbNoH = vector(length=length(horangle))
    ProbHs = vector(length=length(horangle))
    ProbHh = vector(length=length(horangle))
    for(i in 1:length(horangle)){
    ProbNoH[i] = (1-handoverProbabilityCond(r0=r[q],angle=horangle[i],velocity=velocity, density=densityTotal))
    ProbHs[i] = handoverProbabilityCondMulti(r0=r[q],angle=horangle[i],velocity=velocity, density1=density1, density2=density2)
    ProbHh[i] = handoverProbabilityCondMulti(r0=r[q],angle=horangle[i],velocity=velocity, density1=density2, density2=density1)
    }
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      noHCov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      noHCov[q,p] = exp(-density1*angle*sum(Ll))*exp(-density1*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    noHCov[q,1] = noHCov[q,1]*P
    noHCov[q,2] = noHCov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    foo= vector(length=length(horangle))
    for(i in 1:length(horangle)){
      foo[i]=(noHCov[q,1]+noHCov[q,2])*ProbNoH[i]*2*pi*densityTotal*r[q]*exp(-pi*densityTotal*r[q]^2)*(r[2]-r[1])*(horangle[2]-horangle[1])  
    }
    noHCov[q,1] = sum(foo)/pi#(Cov[q,1]+Cov[q,2])*2*pi*density1*r[q]*exp(-pi*density1*r[q]^2)*(r[2]-r[1])
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      HsCov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      HsCov[q,p] = exp(-density1*angle*sum(Ll))*exp(-density1*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    HsCov[q,1] = HsCov[q,1]*P
    HsCov[q,2] = HsCov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    foo= vector(length=length(horangle))
    for(i in 1:length(horangle)){
      foo[i]=(HsCov[q,1]+HsCov[q,2])*ProbHs[i]*2*pi*densityTotal*r[q]*exp(-pi*densityTotal*r[q]^2)*(r[2]-r[1])*(horangle[2]-horangle[1])  
    }
    HsCov[q,1] = sum(foo)/pi#(Cov[q,1]+Cov[q,2])*2*pi*density1*r[q]*exp(-pi*density1*r[q]^2)*(r[2]-r[1])
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
      gn = max(10^(horgain/10)*bv,10^(-2.5))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      HhCov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        bv = 10^(-min(12*((ang-uptilt)/10)^2,20)/10)
        gn = max(10^(horgain/10)*bv,10^(-2.5))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      HhCov[q,p] = exp(-density1*angle*sum(Ll))*exp(-density1*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    HhCov[q,1] = HhCov[q,1]*P
    HhCov[q,2] = HhCov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]

    foo= vector(length=length(horangle))
    for(i in 1:length(horangle)){
      foo[i]=(HhCov[q,1]+HhCov[q,2])*ProbHh[i]*2*pi*densityTotal*r[q]*exp(-pi*densityTotal*r[q]^2)*(r[2]-r[1])*(horangle[2]-horangle[1])  
    }
    HhCov[q,1] = sum(foo)/pi#(Cov[q,1]+Cov[q,2])*2*pi*density1*r[q]*exp(-pi*density1*r[q]^2)*(r[2]-r[1])
    
    
    return(c(Cov[q,1],noHCov[q,1],HsCov[q,1],HhCov[q,1]))
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  Cov = 0
  noHCov = 0
  HsCov = 0
  HhCov = 0
  for(k in 1:length(r)){
    Cov = Cov+c[[k]][1]
    noHCov = noHCov + c[[k]][2]
    HsCov = HsCov + c[[k]][3]
    HhCov = HhCov + c[[k]][4]
  }
  
#  return((2-handPenalty1-handPenalty2)*sum(Cov)+(handPenalty1+handPenalty2-1)*sum(noHCov)+(handPenalty2-1)*sum(HsCov)+(handPenalty1-1)*sum(HhCov))
  return((noHCov+(1-handPenalty1)*HsCov+(1-handPenalty2)*HhCov))
}

#as above, but in a two-network scenario.
MCMobileCovProbMultiRamy = function(rxh,txh,angle,rmax,density1,density2,velocity,handPenalty1,handPenalty2,uptilt,T,gain,Nt,al,an,alpha,beta,gamma,dx,cores){
  Tlin = 10^(T/10)
  # rmax = abs(rxh-txh)*tan(angle/2)
  M = (rmax^2+(rxh-txh)^2)
  densityTotal=(density1+density2)
  rm = sqrt(-log(0.001)/(pi*densityTotal))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  #  y = 2*densityTotal*pi*r*exp(-densityTotal*pi*r^2)*(r[2]-r[1])
  
  iteration=function(q){
    Cov = zeros(nrow=length(r),ncol=2)
    noHCov = zeros(nrow=length(r),ncol=2)
    HsCov = zeros(nrow=length(r),ncol=2)
    HhCov = zeros(nrow=length(r),ncol=2)
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    #  P = c(Plos,(1-Plos))
    #ProbCov[l,j]= 0
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      Cov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[q,p] = exp(-density1*angle*sum(Ll))*exp(-density1*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    
    
    #  n = floor((r[q]/1000)*sqrt(alpha*beta))
    #  a = 0:max(0,n-1)
    #  P = prod(1-exp(-((max(rxh,txh)-(a+1/2)*abs(rxh-txh)/(n))^2)/(2*gamma^2)))
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[q,1] = Cov[q,1]*P
    Cov[q,2] = Cov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    Cov[q,1] = (Cov[q,1]+Cov[q,2])*2*pi*densityTotal*r[q]*exp(-pi*densityTotal*r[q]^2)*(r[2]-r[1])
    
    horangle = seq(from=(pi)/dx,to=(pi*(dx-1)/dx),by=(pi)/dx)
    ProbNoH = vector(length=length(horangle))
    ProbHs = vector(length=length(horangle))
    ProbHh = vector(length=length(horangle))
    for(i in 1:length(horangle)){
      ProbNoH[i] = (1-handoverProbabilityCond(r0=r[q],angle=horangle[i],velocity=velocity, density=densityTotal))
      ProbHs[i] = handoverProbabilityCondMulti(r0=r[q],angle=horangle[i],velocity=velocity, density1=density1, density2=density2)
      ProbHh[i] = handoverProbabilityCondMulti(r0=r[q],angle=horangle[i],velocity=velocity, density1=density2, density2=density1)
    }
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      noHCov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      noHCov[q,p] = exp(-density1*angle*sum(Ll))*exp(-density1*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    noHCov[q,1] = noHCov[q,1]*P
    noHCov[q,2] = noHCov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    foo= vector(length=length(horangle))
    for(i in 1:length(horangle)){
      foo[i]=(noHCov[q,1]+noHCov[q,2])*ProbNoH[i]*2*pi*densityTotal*r[q]*exp(-pi*densityTotal*r[q]^2)*(r[2]-r[1])*(horangle[2]-horangle[1])  
    }
    noHCov[q,1] = sum(foo)/pi#(Cov[q,1]+Cov[q,2])*2*pi*density1*r[q]*exp(-pi*density1*r[q]^2)*(r[2]-r[1])
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      HsCov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      HsCov[q,p] = exp(-density1*angle*sum(Ll))*exp(-density1*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    HsCov[q,1] = HsCov[q,1]*P
    HsCov[q,2] = HsCov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    foo= vector(length=length(horangle))
    for(i in 1:length(horangle)){
      foo[i]=(HsCov[q,1]+HsCov[q,2])*ProbHs[i]*2*pi*densityTotal*r[q]*exp(-pi*densityTotal*r[q]^2)*(r[2]-r[1])*(horangle[2]-horangle[1])  
    }
    HsCov[q,1] = sum(foo)/pi#(Cov[q,1]+Cov[q,2])*2*pi*density1*r[q]*exp(-pi*density1*r[q]^2)*(r[2]-r[1])
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])*(180/pi)
      gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*r[q])))/2)^2))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2) 
      }
      HhCov[q,p] = 0
      
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])*(180/pi)
        gn=(1/Nt)*((sin(Nt*pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2)/(sin(pi*(sin(pi*(rxh-txh)/(4*dr[i])))/2)^2))
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      HhCov[q,p] = exp(-density1*angle*sum(Ll))*exp(-density1*angle*sum(Ln))*exp(-s*N/gain)
      
    }
    
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    HhCov[q,1] = HhCov[q,1]*P
    HhCov[q,2] = HhCov[q,2]*(1-P)
    # Cov[q,] = P*Cov[q,]
    
    foo= vector(length=length(horangle))
    for(i in 1:length(horangle)){
      foo[i]=(HhCov[q,1]+HhCov[q,2])*ProbHh[i]*2*pi*densityTotal*r[q]*exp(-pi*densityTotal*r[q]^2)*(r[2]-r[1])*(horangle[2]-horangle[1])  
    }
    HhCov[q,1] = sum(foo)/pi#(Cov[q,1]+Cov[q,2])*2*pi*density1*r[q]*exp(-pi*density1*r[q]^2)*(r[2]-r[1])
    
    
    return(c(Cov[q,1],noHCov[q,1],HsCov[q,1],HhCov[q,1]))
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  Cov = 0
  noHCov = 0
  HsCov = 0
  HhCov = 0
  for(k in 1:length(r)){
    Cov = Cov+c[[k]][1]
    noHCov = noHCov + c[[k]][2]
    HsCov = HsCov + c[[k]][3]
    HhCov = HhCov + c[[k]][4]
  }
  
  #  return((2-handPenalty1-handPenalty2)*sum(Cov)+(handPenalty1+handPenalty2-1)*sum(noHCov)+(handPenalty2-1)*sum(HsCov)+(handPenalty1-1)*sum(HhCov))
  return((noHCov+(1-handPenalty1)*HsCov+(1-handPenalty2)*HhCov))
}


