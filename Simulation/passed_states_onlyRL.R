#This code generates topologies with different number of BS and Building densities and for each topology we run 4 different RL aproaches and 3 baselines solutions.
#Height_model_1.R saves 1 file for each RL solution and in the file Noint, it saves the values for the baselines.
# Archiaveblerate saves the last 100 rewards (spectrum efficiency) from the last 100 steps (last episode). Height saves the last 100 heights the UAV had for 100 steps (1 episode).
#Heights are saved each 100 epochs. Archieaveblerates each 10.

##Given a certain BS to connect to and a certain starting height, have the UAV pick the best height as quickly as possible
rm(list=ls())
library(spatstat)
library(VGAM)
library(hypergeo)
library(keras)
#install_keras(version = "2.1.3")
library(tensorflow)
#install_tensorflow(gpu=TRUE)
library(parallel)
library(RhpcBLASctl)
source("HetFunctionsv1.R")
source("BHFunctionsv3.R")
source("CFcovprob.R")

try(k_constant(1), silent=TRUE)
try(k_constant(1), silent=TRUE)
options(bitmapType='cairo')

#number of MC trials, episodes and steps
MCtrials = 500
cores=10
blas_set_num_threads(cores)
steps <<- 101
BHdensity = 2.5/(1000^2)
buildDens = 300/(1000^2)
BScand = 5
BScandi = 5
velocity = 10
d = 10
#which BS to be connected to (at the moment the UAV only picks height, not BS association)
whichBS = 1
alpha = 0.5
beta = 300
gamma = 20
buildWidth = 40
buildR = buildWidth/(2*sin(pi/4))
heightParam = 20
#number of antenna elements
Nt= 8
#max and min UAV heights
minH = 40
maxH = 240
BShparam = 25
BSh=30
windowWidth = 3000#800 for salid paper
Freq = 3.5*10^9
BStilt=-30
BHtilt=-30
UAVBHBW = pi*1/4
Tb = -6
Tu = 0
al = 2.1
an = 4
N=10^(-9)
mal = 1
man = 1
BStx = 40 # in watts
N = -174+10*log10(20*10^6)+10
N = 10^(N/10)/1000
#Building and BS parameters
#To study BS density, I have to vary it the same number of MCtrials
tot= c()
#tot <- c(1,5,10,100,500)*1/(1000^2)#
for (i in 1:MCtrials){
  #for (j in 4:length(seq(from=4,to=10,by=1)))
  tot <- c(tot,c(1,2.5,5,100,300)*1/(1000^2))
}

save_append <- function(..., list = character(), file,episodes) {
  #  add objects to existing Rdata file. Original code written by "flodel"
  # on StackOverflow (http://www.linkedin.com/in/florentdelmotte)  . 
  #file="result/Noint1_1e-06_0.00025_.Rdata"
  #file = "BS1_1e-06_0.00025_.Rdata"
  previous  <- load(file, envir = (ev = new.env()) )
   # if (episodes == 1){
    #history1 <-append(ev$history1,history1)
    AchieveableRate <-append(ev$AchieveableRate,AchieveableRate)
    Height <- append(ev$Height, Height)
    maxAchieveableRate <- append(ev$maxAchieveableRate,maxAchieveableRate)
    randAchieveableRate <-append(ev$randAchieveableRate,randAchieveableRate)
    constAchieveableRate <- append(ev$constAchieveableRate,constAchieveableRate)
    const1AchieveableRate <- append(ev$const1AchieveableRate,const1AchieveableRate)
    const2AchieveableRate <- append(ev$const2AchieveableRate,const2AchieveableRate)
    # const3AchieveableRate <- append(ev$const3AchieveableRate,const3AchieveableRate)
    # const4AchieveableRate <- append(ev$const4AchieveableRate,const4AchieveableRate)
    # const5AchieveableRate <- append(ev$const5AchieveableRate,const5AchieveableRate)
    optmAchieveableRate <-append(ev$optmAchieveableRate,optmAchieveableRate)
    hminAchieveableRate <-append(ev$hminAchieveableRate,hminAchieveableRate)
    hmaxAchieveableRate <-append(ev$hmaxAchieveableRate,hmaxAchieveableRate)
    hmedAchieveableRate <-append(ev$hmedAchieveableRate,hmedAchieveableRate)
    hminrealAchieveableRate <-append(ev$hminrealAchieveableRate,hminrealAchieveableRate)
    hmaxrealAchieveableRate <-append(ev$hmaxrealAchieveableRate,hmaxrealAchieveableRate)
    hmedrealAchieveableRate <-append(ev$hmedrealAchieveableRate,hmedrealAchieveableRate)
    maxheight <-append(ev$maxheight,maxheight)
    randheight <-append(ev$randheight,randheight)
    constheight <-append(ev$constheight,constheight)
    const1height <-append(ev$const1height,const1height)
    const2height <-append(ev$const2height,const2height)
    # const3height <-append(ev$const3height,const3height)
    # const4height <-append(ev$const4height,const4height)
    # const5height <-append(ev$const5height,const5height)
    optmheight <-append(ev$optmheight,optmheight)
    hminheight <-append(ev$hminheight,hminheight)
    hmaxheight <-append(ev$hmaxheight,hmaxheight)
    hmedheight <-append(ev$hmedheight,hmedheight)
    hminrealheight <-append(ev$hminrealheight,hminrealheight)
    hmaxrealheight <-append(ev$hmaxrealheight,hmaxrealheight)
    hmedrealheight <-append(ev$hmedrealheight,hmedrealheight)
    save(history1,maxAchieveableRate,randAchieveableRate,constAchieveableRate,const1AchieveableRate,const2AchieveableRate,optmAchieveableRate,hminAchieveableRate,hmaxAchieveableRate,hmedAchieveableRate,hmaxrealAchieveableRate,hminrealAchieveableRate,hmedrealAchieveableRate,AchieveableRate,maxheight,randheight,constheight,const1height,const2height,optmheight,hminheight,hmaxheight,hmedheight,hmaxrealheight,hminrealheight,hmedrealheight,Height, file=file)
    # }else{
    #   history1 <-append(ev$history1,history1)
    #   achieveableRate <-append(ev$achieveableRate,achieveableRate)
    #   height <- append(ev$height, height)
    #   maxAchieveableRate <- ev$maxAchieveableRate
    #   randAchieveableRate <-append(ev$randAchieveableRate,randAchieveableRate)
    #   constAchieveableRate <- ev$constAchieveableRate
    #   const1AchieveableRate <- ev$const1AchieveableRate
    #   const2AchieveableRate <- ev$const2AchieveableRate
    #   # const3AchieveableRate <- ev$const3AchieveableRate
    #   # const4AchieveableRate <- ev$const4AchieveableRate
    #   # const5AchieveableRate <- ev$const5AchieveableRate
    #   optmAchieveableRate <-ev$optmAchieveableRate
    #   hminAchieveableRate <-ev$hminAchieveableRate
    #   hmaxAchieveableRate <-ev$hmaxAchieveableRate
    #   hmedAchieveableRate <-ev$hmedAchieveableRate
    #   hminrealAchieveableRate <-ev$hminrealAchieveableRate
    #   hmaxrealAchieveableRate <-ev$hmaxrealAchieveableRate
    #   hmedrealAchieveableRate <-ev$hmedrealAchieveableRate
    #   maxheight <-ev$maxheight
    #   randheight <-append(ev$randheight, randheight)
    #   constheight <-ev$constheight
    #   const1height <-ev$const1height
    #   const2height <-ev$const2height
    #   # const3height <-ev$const3height
    #   # const4height <-ev$const4height
    #   # const5height <-ev$const5height
    #   optmheight <-ev$optmheight
    #   hminheight <-ev$hminheight
    #   hmaxheight <-ev$hmaxheight
    #   hmedheight <-ev$hmedheight
    #   hminrealheight <-ev$hminrealheight
    #   hmaxrealheight <-ev$hmaxrealheight
    #   hmedrealheight <<-ev$hmedrealheight
    #   save(history1,maxAchieveableRate,randAchieveableRate,constAchieveableRate,const1AchieveableRate,const2AchieveableRate,optmAchieveableRate,hminAchieveableRate,hmaxAchieveableRate,hmedAchieveableRate,hmaxrealAchieveableRate,hminrealAchieveableRate,hmedrealAchieveableRate,achieveableRate,maxheight,randheight,constheight,const1height,const2height,optmheight,hminheight,hmaxheight,hmedheight,hmaxrealheight,hminrealheight,hmedrealheight,height, file=file)
    #   
    # }
  
}

#Raytracing function for checking LOS
isLOS = function(buildings,buildR,buildH,x0,y0,x,y,h,BSh){
  angle = atan2((y-y0),(x-x0))
  dist = sqrt((x-x0)^2+(y-y0)^2)
  
  build = buildings
  build = shift(build,c(-x0,-y0))
  build = rotate(build,angle=-angle)
  
  buildX = build$x
  buildY = build$y
  foo = which(buildX<dist)
  buildX = buildX[foo]
  buildY = buildY[foo]
  buildH = buildH[foo]
  
  foo = which(buildX>0)
  buildX = buildX[foo]
  buildY = buildY[foo]
  buildH = buildH[foo]
  
  foo = which(abs(buildY)<=buildR)
  buildX = buildX[foo]
  buildY = buildY[foo]
  buildH = buildH[foo]
  
  foo = buildH>((abs(h-BSh)*(buildX/dist)+min(BSh,h)))
  if(length(which(foo==TRUE))>0){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}


antennaGain = function(r,BSh,h,Nt){
  angle = atan2(BSh-h,r)
  g=(1/Nt)*((sin(Nt*pi*(sin(angle))/2)^2)/(sin(pi*(sin(angle))/2)^2))
  g=g*BStx*K
  return(g)
}



#normalise the neural network input data

normaliseData = function(SINR,h,maxSINR){#normaliseData(SINR = mcSINR[foo,j,],h=h[foo],maxSINR = mcSINR)
  
  #normalise observed power
  h =(h-minH)/(maxH-minH)
  
  mnfoo=min(maxSINR)
  mxfoo=max(maxSINR)
  if(mxfoo>0){
    SINR=(SINR)/(mxfoo)#(SINR - mean(SINR)) / sd(SINR)#(P-mnfoo)/(mxfoo-mnfoo)
  }#returning absolut values
  
  return(c(SINR,h))
}
#Inicialize variables; it is need before run each of the models
EPSILON_DECAY <<- 0.9#9#.9
DISCOUNT <<- 0.99#}else{DISCOUNT <<- 0.9}
#DISCOUNT <<- 0.99
REPLAY_MEMORY_SIZE <<- 100  # How many last steps to keep for model training
MINIBATCH_SIZE <<- 16  #64 How many steps (samples) to use for training
epoch <<-200
# Exploration settings
epsilon <<- 1
MIN_EPSILON <<- 0.05
i=1
while (epsilon > MIN_EPSILON) {
  epsilon = epsilon * EPSILON_DECAY
  i=i+1
}
print(i)
inicialize = function(m){
  Tu <<- 1
  epsilon <<- 1
  # Initializing vectors that will be need to each mc trial
  
  #our performance metric. We measure the rolling average of the cumulative rate of the UAV link
  achieveableRate <<- c()
  AchieveableRate <<- c(0)
  history1 <<-c()
  #Vector to save the heights trhough the path
  height <<-c()
  Height <<-c()
  #the replay memory for the training
  action_replay <<-vector(length=REPLAY_MEMORY_SIZE)
  reward_replay <<-vector(length=REPLAY_MEMORY_SIZE)
  replay_index <<-1
  ##UAV antenna gain
  UAVBHgain <<- 4*pi/((UAVBHBW/2)^2)
  K <<- (((3*10^8)/Freq)/(4*pi))^2
  episodes <<- 10#30 # 
  
  maxAchieveableRate <<-  c()
  randAchieveableRate <<- c()
  
  constAchieveableRate <<- c()
  const2AchieveableRate <<- c()
  const1AchieveableRate <<- c()
  # const3AchieveableRate <<- c()
  # const4AchieveableRate <<- c()
  # const5AchieveableRate <<- c()
  optmAchieveableRate <<-c()
  hminAchieveableRate <<-c()
  hmaxAchieveableRate <<-c()
  hmedAchieveableRate <<-c()
  hminrealAchieveableRate <<-c()
  hmaxrealAchieveableRate <<-c()
  hmedrealAchieveableRate <<-c()
  maxheight <<-c()
  randheight <<-c()
  
  constheight <<-c()
  const1height <<-c()
  const2height <<-c()
  # const3height <<-c()
  # const4height <<-c()
  # const5height <<-c()
  optmheight <<-c()
  hminheight <<-c()
  hmaxheight <<-c()
  hmedheight <<-c()
  hminrealheight <<-c()
  hmaxrealheight <<-c()
  hmedrealheight <<-c()
  neurons <<-38
  model <<- keras_model_sequential()
  model %>%
    layer_dense(units=neurons+1,input_shape = c(neurons)) %>%
    layer_dense(units = 200, activation = 'relu') %>%
    layer_dense(units = 200, activation = 'relu') %>%
    layer_dense(units = 200, activation = 'relu') %>%
    layer_dense(units = 3, activation = 'softplus')
  #model <- load_model_tf("model_Height_")
  if (m>1){model <- load_model_tf("model_Height_")}
  
  model %>% compile(
    optimizer = 'adam', 
    loss = 'sparse_categorical_crossentropy',
    metrics = c('accuracy'),
    # learning_rate=0.01
  )
  
  current_state_replay<<-zeros(nrow=REPLAY_MEMORY_SIZE,ncol=neurons)
  next_state_replay<<-zeros(nrow=REPLAY_MEMORY_SIZE,ncol=neurons)
 
  #copy the model
  targetModel <<-model
  targetModel  %>% set_weights(model %>% get_weights())
}

##UAV antenna gain
UAVBHgain = 4*pi/((UAVBHBW/2)^2)
K = (((3*10^8)/Freq)/(4*pi))^2
#interactions considers all the needed interactions to run all the densities. 
#Creating the files to save - so append can use it
#for(m in 1:5){
  
   # BHdensity = tot[m]
   # buildDens = 200/(1000^2)
   # if(m==4|| m==5){
   #   print("building change")
   #   buildDens= tot[m]
   #   BHdensity = 5/(1000^2)
   # }

   history1=c()
   AchieveableRate=c()
   Height= c()
   maxAchieveableRate= c()
   randAchieveableRate= c()
   constAchieveableRate= c()
   const1AchieveableRate= c()
   const2AchieveableRate= c()
   # const3AchieveableRate= c()
   # const4AchieveableRate= c()
   # const5AchieveableRate= c()
   optmAchieveableRate=c()
   hminAchieveableRate=c()
   hmaxAchieveableRate=c()
   hmedAchieveableRate=c()
   hminrealAchieveableRate=c()
   hmaxrealAchieveableRate=c()
   hmedrealAchieveableRate=c()
   maxheight=c()
   randheight=c()
   constheight=c()
   const1height=c()
   const2height=c()
   # const3height=c()
   # const4height=c()
   # const5height=c()
   optmheight=c()
   hminheight=c()
   hmaxheight=c()
   hmedheight=c()
   hminrealheight=c()
   hmaxrealheight=c()
   hmedrealheight=c()
   base_name = "result/Final_trying_1"
   local = paste0(base_name, BHdensity,"_",buildDens,"_",".Rdata")
   save(history1,maxAchieveableRate,randAchieveableRate,constAchieveableRate,const1AchieveableRate,const2AchieveableRate,optmAchieveableRate,hminAchieveableRate,hmaxAchieveableRate,hmedAchieveableRate,hmaxrealAchieveableRate,hminrealAchieveableRate,hmedrealAchieveableRate,AchieveableRate,maxheight,randheight,constheight,const1height,const2height,optmheight,hminheight,hmaxheight,hmedheight,hmaxrealheight,hminrealheight,hmedrealheight,Height, file=local)
 #}


#multiple MC trials. In each trial we generate the environment and then run the algorithm over multiple episodes
for(m in 1:MCtrials){
 # m=1
  # BHdensity = tot[order(tot)[m]]#tot[m]
  # buildDens = 200/(1000^2)
  # 
  # if(tot[order(tot)[m]] == 3e-4 || tot[order(tot)[m]]==1e-4){
  #   print("building change")
  #   buildDens= tot[order(tot)[m]]
  #   BHdensity = 5/(1000^2)
  # } 
   # BHdensity = 5/(1000^2)
   # buildDens = 200/(1000^2)
    print(m)
    h = seq(from=minH,to=maxH,by=1)#the possible heights, sequency from min UAV height to maximum
    local = paste0("enviroment/BSh_fixed_", BHdensity,"_",buildDens,"_",m,".Rdata")
    load(local)
    print("after interacton")
    mP = 1:(BScand)
    Int = (mP[BScand]+1):(mP[BScand]+BScand*BScandi)
    Distances = (Int[BScand*BScandi]+1):(Int[BScand*BScandi]+BScand)
    C = (Distances[BScand]+1):(Distances[BScand]+BScand)
    sSINR = (C[BScand]+1):(C[BScand]+BScand)
    # nSINR = C[BScand]+2
    #the S,Int,Dist,Cov and SNIR for all possible h
    mcS =array(dim=c(length(h),steps,BScand))
    mcInt =array(dim=c(length(h),steps,BScand*BScandi))
    mcDist =array(dim=c(length(h),steps,BScand))
    mcCov =array(dim=c(length(h),steps,BScand))
    mcSINR = array(dim=c(length(h),steps,BScand))
    hmin = array(dim=c(length(h),steps))
    hmax = array(dim=c(length(h),steps))
    uavy_change = array(dim=c(length(h),steps))
    yBScand = array(dim=c(length(h),steps))
      
    for(k in 0:(length(h)*steps-1)){
      j = floor(k/steps)+1
      w = k%%steps+1
      mcS[j,w,]=opt[[k+1]][mP]
      mcInt[j,w,]=opt[[k+1]][Int]
      mcDist[j,w,]=opt[[k+1]][Distances]
      mcCov[j,w,] =opt[[k+1]][C]
      mcSINR[j,w,] = opt[[k+1]][sSINR]
      hmin[j,w] = opt[[k+1]][46]
      hmax[j,w] = opt[[k+1]][47]
      uavy_change[j,w] = opt[[k+1]][48]
      yBScand[j,w] = opt[[k+1]][49]
    }

    #which height is optimum
    optH=vector(length=steps)
    for(i in 1:steps){
      #this is the optimal mean height for the city, the UAV will be at this height
      optH[i] = h[which.max(mcSINR[,i,whichBS])]  
    }
    
   
    #Inicializing the simplest solution, here I will also calculate the benchmarks, in the other proposed solution we dont need to calculate the benchmarks anymore
    inicialize(m)
    print("after inicialize noint")
    model %>% compile(
      optimizer = 'adam', 
      loss = 'sparse_categorical_crossentropy',
      metrics = c('accuracy'),
      #learning_rate=0.01
    )
    
    for(k in 1:episodes){
      #starting height
      #k=1
      #print(k)
      
      #time1=Sys.time()
      #All the solutions start at the same height
      UAVh = minH#floor(runif(n=1,min=minH,max=maxH))
     
      episode_reward = 0
      previous_reward = 0
    
      rand_episode_reward=0
      
     # if (k==1){#k=1
        randUAVh = minH
        genieUAVh = minH
        constUAVh = minH
        const1UAVh = minH
        const2UAVh = minH
        # const3UAVh = minH
        # const4UAVh = minH
        # const5UAVh = minH
        optmUAVh = minH
        hminUAVh = minH
        hmaxUAVh = minH
        hmedUAVh = minH
        hminrealUAVh = minH
        hmaxrealUAVh = minH
        hmedrealUAVh = minH
        rand_episode_reward=0
        genie_episode_reward=0
        const_episode_reward=0
        const1_episode_reward=0
        const2_episode_reward=0
        #const3_episode_reward=0
        # const4_episode_reward=0
        # const5_episode_reward=0
        optm_episode_reward=0
        hmin_episode_reward = 0
        hmax_episode_reward = 0
        hmed_episode_reward = 0
        hminreal_episode_reward = 0
        hmaxreal_episode_reward = 0
        hmedreal_episode_reward = 0
     # }else{epsilon=0.3}
      #print("new steps")
      state_4 = 0
      state_3 = 0
      state_2 = 0
      state_1 = 0
      action_4 = 0
      action_3 = 0
      action_2 = 0
      action_1 = 0
      reward_4 = 0
      reward_3 = 0
      reward_2 = 0
      reward_1 = 0
      for(j in 1:(steps-1)){
        #j=1
        if (j>1){
          state_4 = state_3
          state_3 = state_2
          state_2 = state_1
          state_1 = current_state
          action_4 = action_3
          action_3 = action_2
          action_2 = action_1
          action_1 = action
          reward_4 = reward_3
          reward_3 = reward_2
          reward_2 = reward_1
          reward_1 = reward
        }
        foo = which(h==UAVh)
        # if ((UAVh-d)>minH){foo1 = which(h==(UAVh-d))}else{foo1=foo}
        # if ((UAVh+d)<maxH){foo2 = which(h==(UAVh+d))}else{foo2=foo}
        current_state = normaliseData(SINR = mcSINR[foo,j,],h=h[foo],maxSINR = mcSINR)
        
        
        #take an action
        if (j>4){
          state = c(current_state, state_1,action_1,reward_1,state_2,action_2,reward_2,state_3,action_3,reward_3,state_4,action_4,reward_4)
          foo = rbind(state,vector(length=neurons))
        #1- do nothing, 2- move down a meter, 3- move up a meter
          if(runif(n=1,min=0,max=1)>epsilon){
            action = NULL
            #print("action = 0")
            while(is.null(action) | length(action)==0){        
              action = model %>% predict(foo)#state)
              action = which.max(action[1,])
              #print("inside while action")
            }
          }else{
            action = sample(1:3, 1)#floor(runif(n=1,min=1,max=4))
            #print("else action")
          }
        #print(action)
        }else{ action = sample(1:3, 1)}
        #update UAV height
        if (!is.null(action)){
          if(action==1){
            if(UAVh<=(d+minH)){UAVh=minH}
            else {UAVh=UAVh-d}
          }else if(action==3){
            if((UAVh+d)>=maxH){UAVh = maxH
            }else{UAVh=UAVh+d}
          }
        }
        
        if(j+1<=length(optH)){subs = j+1
        }else {subs = j}      
        #get the new state and reward
        foo = which(h==UAVh)
        reward = log2(1+mcSINR[foo,subs,whichBS]) - previous_reward 
        # if ((UAVh-d)>minH){foo1 = which(h==(UAVh-d))}else{foo1=foo}
        # if ((UAVh+d)<maxH){foo2 = which(h==(UAVh+d))}else{foo2=foo}
        
        new_state = normaliseData(SINR = mcSINR[foo,j,],h=h[foo],maxSINR = mcSINR)#normaliseDataComplete(BSdensity = BHdensity ,Builddensity =buildDens ,P=mcS[foo,j,],Dist=mcDist[foo,j,],h=h[foo],SINR = mcSINR[foo,j,])#,SINRmin = mcSINR[foo1,j+1,],SINRmax = mcSINR[foo2,j+1,] )
        new_state = c(new_state,current_state,action,reward, state_1,action_1,reward_1,state_2,action_2,reward_2,state_3,action_3,reward_3)
      #print("replay")
        
        ##add to the replay memory
        if (j>MINIBATCH_SIZE){
          if(replay_index<REPLAY_MEMORY_SIZE){
            current_state_replay[replay_index,]=state
            action_replay[replay_index]=action
            reward_replay[replay_index]= reward
            next_state_replay[replay_index,]=new_state
            replay_index = replay_index+1
          }else{
            current_state_replay=rbind(current_state_replay[2:REPLAY_MEMORY_SIZE,],state)
            action_replay = c(action_replay[2:REPLAY_MEMORY_SIZE],action)
            reward_replay = c(reward_replay[2:REPLAY_MEMORY_SIZE],reward)
            next_state_replay=rbind(next_state_replay[2:REPLAY_MEMORY_SIZE,],new_state)
          }

          ##training
          if(replay_index>=MINIBATCH_SIZE){

            #randomly sample from the replay memory for a minibatch of data
            minibatch = sample(1:replay_index,size=MINIBATCH_SIZE,replace=FALSE)#(j-MINIBATCH_SIZE):j#
            curr = current_state_replay[minibatch,]
            act = action_replay[minibatch]
            rev = reward_replay[minibatch]
            future = next_state_replay[minibatch,]

            #get the Q values for the current states from the neural network
            curr_q_list =  model %>% predict(curr)

            #get the Q values for the future states from the target neural network
            future_q_list = targetModel  %>% predict(future)

            X = zeros(nrow=MINIBATCH_SIZE,ncol=neurons)
            Y = zeros(nrow=MINIBATCH_SIZE,ncol=1)

            #for each entry in the minibatch do training
            for(i in 1:MINIBATCH_SIZE){

              #max Q value from the future states  i=1
              max_future_Q = max(future_q_list[i,])
              #new Q value given reward and next max Q value
              new_Q = rev[i]+DISCOUNT*max_future_Q

              #Q values for the current state
              Qs = curr_q_list[i,]
              #update the Q value for the action taken
              Qs[act[i]]=new_Q
              X[i,]=curr[i,]
              #Y[i,]=Qs
              Y[i]  = which.max(Qs)-1
            }

            #fit the model
            history1<- model  %>% fit(X, Y, epochs = epoch, verbose=0)#
          }
        }
        
        #Add other solutions
       # if (k==1){
          
          if (const1UAVh < 120){const1UAVh = const1UAVh+d}
          if (const2UAVh < 240){const2UAVh = const2UAVh+d}
          # if (const3UAVh < 70){const3UAVh = const3UAVh+d}
          # if (const4UAVh < 95){const4UAVh = const4UAVh+d}
          # if (const5UAVh < 120){const5UAVh = const5UAVh+d}
          randAction = sample(1:3, 1)#floor(runif(n=1,min=1,max=4))
          #random walk
          if(randAction==2){
            randUAVh=max(randUAVh-d,minH)
          }else if(randAction==3){
            randUAVh=min(randUAVh+d,maxH)
          }
          
          #optimum move
          if(j+1<=length(optH)){subs = j+1
          }else {subs = j}
          
          if(genieUAVh>optH[subs]){
            if(genieUAVh<=(d+minH)){genieUAVh=minH
            }else {genieUAVh=genieUAVh-d}
          }else if(genieUAVh<optH[subs]){
            if((genieUAVh+d)>=maxH){genieUAVh=maxH
            }else {genieUAVh=genieUAVh+d}
          }
          optmUAVh = optH[subs] 
          
          #hmin  move# j=1
          #print(hmin[hminUAVh - minH,j])
          #this is to be sure the values calculated in Sallid solution will be inside the min max height
          if (hminUAVh - minH+1>0){if (hmin[hminUAVh - minH +1,j] != 0){
            hminUAVh =hmin[hminUAVh - minH +1,j]
          }
          }
          #this epositioning is because the vector has max-min indexes and we have to realocate the sallid solution to deslocated position
          
          if (hminUAVh < minH ){hminUAVh = minH}
          if (hminUAVh > maxH){hminUAVh = maxH}
          #this is because the possible values of UAV height is from 20 to 200, so we have to move the values to the right position
          #hmax
          if (hmaxUAVh - minH +1>0){
            hmaxUAVh = hmax[hmaxUAVh - minH+1,j]
          }
          if (hmaxUAVh > maxH){hmaxUAVh = maxH}
          if (hmaxUAVh < minH){hmaxUAVh = minH}
          #print(hmax[hmaxUAVh - minH,j])
          
          hmedUAVh = as.integer((hmaxUAVh+hminUAVh)/2)
          
          if (hmaxUAVh > hmaxrealUAVh) {hmaxrealUAVh = min(hmaxrealUAVh+d,maxH)
          }else if (hmaxUAVh < hmaxrealUAVh)  {hmaxrealUAVh = max(hmaxrealUAVh-d,minH)}
          
          if (hminUAVh > hminrealUAVh) {hminrealUAVh = min(hminrealUAVh+d,maxH)
          }else if (hminUAVh < hminrealUAVh)  {hminrealUAVh = max(hminrealUAVh-d,minH)}
          
          if (hmedUAVh > hmedrealUAVh) {hmedrealUAVh = min(hmedrealUAVh+d,maxH)
          }else if (hmedUAVh < hmedrealUAVh)  {hmedrealUAVh = max(hmedrealUAVh-d,minH)}
          
          
          foo = which(h==randUAVh)
          #if(randUAVh<maxH & randUAVh>mi
          
          randReward = log2(1+mcSINR[foo,subs,whichBS])
          #}else{randReward=0}
          
          foo = which(h==genieUAVh)
          genieReward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==as.integer(constUAVh))
          constReward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==as.integer(const1UAVh))
          const1Reward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==as.integer(const2UAVh))
          const2Reward = log2(1+mcSINR[foo,subs,whichBS])
          
          # foo = which(h==as.integer(const3UAVh))
          # const3Reward = log2(1+mcSINR[foo,subs,whichBS])
          # 
          # foo = which(h==as.integer(const4UAVh))
          # const4Reward = log2(1+mcSINR[foo,subs,whichBS])
          # 
          # foo = which(h==as.integer(const5UAVh))
          # const5Reward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==as.integer(optmUAVh))
          optmReward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==hminUAVh)
          hminReward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==hmaxUAVh)
          hmaxReward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==hmedUAVh)
          hmedReward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==hmaxrealUAVh)
          hmaxrealReward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==hminrealUAVh)
          hminrealReward = log2(1+mcSINR[foo,subs,whichBS])
          
          foo = which(h==hmedrealUAVh)
          hmedrealReward = log2(1+mcSINR[foo,subs,whichBS])
          rand_episode_reward=randReward+rand_episode_reward
          genie_episode_reward = genieReward+ genie_episode_reward 
          const_episode_reward = constReward+ const_episode_reward
          const1_episode_reward = const1Reward+ const1_episode_reward
          const2_episode_reward = const2Reward+ const2_episode_reward
          # const3_episode_reward = const3Reward+ const3_episode_reward
          # const4_episode_reward = const4Reward+ const4_episode_reward
          # const5_episode_reward = const5Reward+ const5_episode_reward
          optm_episode_reward = optmReward+ optm_episode_reward
          hmin_episode_reward=hminReward+hmin_episode_reward
          hmax_episode_reward=hmaxReward+hmax_episode_reward
          hmed_episode_reward=hmedReward+hmed_episode_reward
          hminreal_episode_reward=hminReward+hminreal_episode_reward
          hmaxreal_episode_reward=hmaxReward+hmaxreal_episode_reward
          hmedreal_episode_reward=hmedReward+hmedreal_episode_reward
          randAchieveableRate[j] = randReward
          maxAchieveableRate[j] = genieReward
          constAchieveableRate[j] = constReward
          const1AchieveableRate[j] = const1Reward
          const2AchieveableRate[j] = const2Reward
          # const3AchieveableRate[j] = const3Reward
          # const4AchieveableRate[j] = const4Reward
          # const5AchieveableRate[j] = const5Reward
          optmAchieveableRate[j] = optmReward
          hminAchieveableRate[j] = hminReward
          hmaxAchieveableRate[j] = hmaxReward
          hmedAchieveableRate[j] = hmedReward
          hminrealAchieveableRate[j] = hminrealReward
          hmaxrealAchieveableRate[j] = hmaxrealReward
          hmedrealAchieveableRate[j] = hmedrealReward
          randheight[j] = randUAVh
          maxheight[j] = genieUAVh
          constheight[j] = constUAVh
          const1height[j] = const1UAVh
          const2height[j] = const2UAVh
          # const3height[j] = const3UAVh
          # const4height[j] = const4UAVh
          # const5height[j] = const5UAVh
          optmheight[j] =optmUAVh
          hminheight[j] =hminUAVh
          hmaxheight[j] =hmaxUAVh
          hmedheight[j] =hmedUAVh
          hminrealheight[j] =hminrealUAVh
          hmaxrealheight[j] =hmaxrealUAVh
          hmedrealheight[j] =hmedrealUAVh
             
        #}#other solutions
        
      
         #store the rolling averages
        foo = which(h==UAVh)
        achieveableRate[j] = log2(1+mcSINR[foo,subs,whichBS])
        episode_reward= achieveableRate[j]+episode_reward
        height[j] = UAVh
      
        #decay the exploration factor
        if(epsilon > MIN_EPSILON){
          epsilon = epsilon*EPSILON_DECAY
          epsilon = max(MIN_EPSILON, epsilon)
        }
        
        if(j%%5==0){
          #once every 5 episodes update the targetModel to have the same weights as the model
          targetModel  %>% set_weights(model %>% get_weights())
        }
        previous_reward = reward
      }#steps
     
      print(mean(achieveableRate))
      print(mean(constAchieveableRate))
      print(mean(const1AchieveableRate))
      print(mean(maxAchieveableRate))
     
     
        model %>% save_model_tf("model_Height_")
        print("dentro")
        cat("\n salvando rl part.", length(mean(AchieveableRate)))
        local = paste0(base_name, BHdensity,"_",buildDens,"_",".Rdata")#as.numeric(Sys.time())
        #targetModel %>% save_model_tf("model_Height_")
        save_append(history1,maxAchieveableRate,randAchieveableRate,constAchieveableRate,const1AchieveableRate,const2AchieveableRate,optmAchieveableRate,hminAchieveableRate,hmaxAchieveableRate,hmedAchieveableRate,hmaxrealAchieveableRate,hminrealAchieveableRate,hmedrealAchieveableRate,AchieveableRate,maxheight,randheight,constheight,const1height,const2height,optmheight,hminheight,hmaxheight,hmedheight,hmaxrealheight,hminrealheight,hmedrealheight,Height, file=local,episodes=k)
      
      
    }#episodes
}#tot

