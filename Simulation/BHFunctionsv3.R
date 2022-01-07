##Functions for the UAV backhaul
#library(LICORS)
#library(distr)
#library(VGAM)
#library(expint)

#Calculate the Sallid approach
#sallidMinMax = function()




##Get the achieveable spectral efficiency of UAV-BS link, but instead of treating the UAV gain as being rectangular we consider a smooth gain lobe
getDualRateAntennaGain = function(x,BHBS,BSh,withFading,LOS,BHtilt,iBSBH,iBSBHHeight,al,an,mal,man,antgain,htilt,vtilt,beamwidth,PWRgain,N,alpha,beta,gamma){
  dist = sqrt((BHBS[1]-x[1])^2+(BHBS[2]-x[2])^2+(BSh-x[3])^2)
  
  rdist = sqrt((BHBS[1]-x[1])^2+(BHBS[2]-x[2])^2)
  angle = (atan2(x[3]-BSh,rdist)*180/pi)
  hangle = atan2((BHBS[2]-x[2]),(BHBS[1]-x[1]))
  Gv = -min(12*((angle-BHtilt)/10)^2,20)
  #Gv = 10^(Gv/10)
  g = -min(-(Gv+antgain),25)
  g = 10^(g/10)
  g=g*PWRgain
  
  hanglediff = abs(hangle-htilt)
  vanglediff = abs((angle*pi/180)-vtilt)
  g = g*dnorm(x=hanglediff,mean=0,sd=beamwidth)*dnorm(x=vanglediff,mean=0,sd=beamwidth)
  
  if(LOS==TRUE){
    PL = dist^(-al)
    if(withFading==TRUE){
      PL = PL*rgamma(n=1,rate=mal,shape=mal)
    }
  }
  else{
    PL = dist^(-an)
    if(withFading==TRUE){
      PL = PL*rgamma(n=1,rate=man,shape=man) 
    }
  }
  
  I=0
  
  if(length(iBSBH[,1])>0){
    for(i in 1:length(iBSBH[,1])){
      idist = sqrt((iBSBH[i,1]-x[1])^2+(iBSBH[i,2]-x[2])^2+(BSh-x[3])^2)
      
      rdist =sqrt((iBSBH[i,1]-x[1])^2+(iBSBH[i,2]-x[2])^2)
      angle = (atan2(x[3]-BSh,rdist)*180/pi)
      hangle = atan2((iBSBH[i,2]-x[2]),(iBSBH[i,1]-x[1]))
      Gv = -min(12*((angle-BHtilt)/10)^2,20)
      ig = -min(-(Gv+antgain),25)
      ig = 10^(ig/10)
      ig=ig*PWRgain
      
      hanglediff = abs(hangle-htilt)
      vanglediff = abs((angle*pi/180)-vtilt)
      ig = ig*dnorm(x=hanglediff,mean=0,sd=beamwidth)*dnorm(x=vanglediff,mean=0,sd=beamwidth)
      
      test = (iBSBHHeight[i]<=x[3])
      if(test==TRUE){
        iLOS = TRUE
        iPL = idist^(-al) 
        if(withFading==TRUE){
          iPL =iPL*rgamma(n=1,rate=mal,shape=mal)
        }
      }
      else{
        iLOS = FALSE
        iPL = idist^(-an) 
        if(withFading==TRUE){
          iPL =iPL*rgamma(n=1,rate=man,shape=man)
        }
      }
      I = I+ig*iPL
    } 
  }
  
  
  SINR = (g*PL/(N+I))
  return(log2(1+SINR))
  
} 

##Get the achieveable rate, when the base stations have antennas that use Ramy's directional model
getDualRateRamyAntenna = function(x,BHBS,BSh,withFading,LOS,BHtilt,iBSBH,iBSBHHeight,al,an,mal,man,Nt,PWRgain,N,alpha,beta,gamma){
  dist = sqrt((BHBS[1]-x[1])^2+(BHBS[2]-x[2])^2+(BSh-x[3])^2)
  
  rdist = sqrt((BHBS[1]-x[1])^2+(BHBS[2]-x[2])^2)
  angle = (atan2(x[3]-BSh,rdist))
  hangle = atan2((BHBS[2]-x[2]),(BHBS[1]-x[1]))
  g=(1/Nt)*((sin(Nt*pi*(sin(angle))/2)^2)/(sin(pi*(sin(angle))/2)^2))
#  g=(1/Nt)*((sin(Nt*pi*(sin(pi*(x[3]-BSh)/(4*rdist)))/2)^2)/(sin(pi*(sin(pi*(x[3]-BSh)/(4*rdist)))/2)^2))
  
  g=g*PWRgain
  
  
  if(LOS==TRUE){
    PL = dist^(-al)
    if(withFading==TRUE){
      PL = PL*rgamma(n=1,rate=mal,shape=mal)
    }
  }
  else{
    PL = dist^(-an)
    if(withFading==TRUE){
      PL = PL*rgamma(n=1,rate=man,shape=man) 
    }
  }
  
  I=0
  
  if(length(iBSBH[,1])>0){
    for(i in 1:length(iBSBH[,1])){
      idist = sqrt((iBSBH[i,1]-x[1])^2+(iBSBH[i,2]-x[2])^2+(BSh-x[3])^2)
      
      rdist =sqrt((iBSBH[i,1]-x[1])^2+(iBSBH[i,2]-x[2])^2)
      angle = (atan2(x[3]-BSh,rdist))
      hangle = atan2((iBSBH[i,2]-x[2]),(iBSBH[i,1]-x[1]))
      ig=(1/Nt)*((sin(Nt*pi*(sin(angle))/2)^2)/(sin(pi*(sin(angle))/2)^2))
   #   ig=(1/Nt)*((sin(Nt*pi*(sin(pi*(x[3]-BSh)/(4*rdist)))/2)^2)/(sin(pi*(sin(pi*(x[3]-BSh)/(4*rdist)))/2)^2))
      ig=ig*PWRgain
      
      
      test = (iBSBHHeight[i]<=x[3])
      if(test==TRUE){
        iLOS = TRUE
        iPL = idist^(-al) 
        if(withFading==TRUE){
          iPL =iPL*rgamma(n=1,rate=mal,shape=mal)
        }
      }
      else{
        iLOS = FALSE
        iPL = idist^(-an) 
        if(withFading==TRUE){
          iPL =iPL*rgamma(n=1,rate=man,shape=man)
        }
      }
      I = I+ig*iPL
    } 
  }
  
  
  SINR = (g*PL/(N+I))
  return(log2(1+SINR))
  
} 


getInt = function(x,int,BHBS,grid,UAVBHBW){
  intBSx = grid$x[int]-x[1]
  intBSy = grid$y[int]-x[2]
  BSx = BHBS[1]-x[1]
  BSy = BHBS[2]-x[2]
  
  di = acos((BSx*intBSx+BSy*intBSy)/(sqrt(BSx^2+BSy^2)*sqrt(intBSx^2+intBSy^2)))
  if(length(find(is.na(di)))>0){
    di[is.na(di)] = Inf
  }
  int = int[abs(di)<UAVBHBW/2]
  return(cbind(grid$x[int],grid$y[int]))
} 

getInt2 = function(x,int,BHBS,grid,UAVBHBW){
  intBSx = grid$x[int]-x[1]
  intBSy = grid$y[int]-x[2]
  BSx = BHBS[1]-x[1]
  BSy = BHBS[2]-x[2]
  
  di = acos((BSx*intBSx+BSy*intBSy)/(sqrt(BSx^2+BSy^2)*sqrt(intBSx^2+intBSy^2)))
  if(length(find(is.na(di)))>0){
    di[is.na(di)] = Inf
  }
  int = int[abs(di)<UAVBHBW/2]
  return(int)
} 


