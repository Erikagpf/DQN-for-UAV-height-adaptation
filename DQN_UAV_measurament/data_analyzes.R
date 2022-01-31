#This file apply the simple solution that adapts the height to the real world data 

rm(list=ls())
library("readxl")
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
#source("CFcovprob.R")
library(geosphere)

try(k_constant(1), silent=TRUE)
try(k_constant(1), silent=TRUE)
options(bitmapType='cairo')

#Reading the data in order
#Steps on the real world measurament files
steps = 171
cores=6
blas_set_num_threads(cores)
MCtrials = 10
#Reading rate
file_list =list.files(pattern = "n.xls" )
mcRate = array(dim=c(length(file_list),steps))
for (i in 1:length(file_list)){
  my_data =c(data.frame(read_excel(file_list[i])) )
  if (i <4){ mcRate[i+8,] = my_data$Rate#to order because it reads the files starting with 1 first
  #print(file_list[i])
       # print(i+8)
  }else{mcRate[i-3,] = my_data$Rate
  #print(file_list[i])
  #print(i-3)
  }
}
#Reading SINR
#steps = 161
file_list =list.files(pattern = "msinr.xls")
mcSINR = array(dim=c(length(file_list),steps))
for (i in 1:length(file_list)){
  my_data =c(data.frame(read_excel(file_list[i])) )
  if (i <4){ mcSINR[i+8,] = my_data$SINR#to order because it reads the files starting with 1 first
  #print(file_list[i])
 # print(i+8)
  }else{mcSINR[i-3,] = my_data$SINR
  #print(file_list[i])
  #print(i-3)
  }
}
#for the Gran cannal experiment, the ue was connected all the timein the Trinity Enterprise
latitude_cell = 53.341667
longitude_cell = -6.2225

save_append <- function(..., list = character(), file) {
  #  add objects to existing Rdata file. Original code written by "flodel"
  # on StackOverflow (http://www.linkedin.com/in/florentdelmotte)  . 
  #file="result/Noint1_1e-06_0.00025_.Rdata"
  #file = "BS1_1e-06_0.00025_.Rdata"
  previous  <- load(file, envir = (ev = new.env()) )
    history1 <-append(ev$history1,history1)
    achieveableRate <-append(ev$achieveableRate,achieveableRate)
    height <- append(ev$height, height)
    maxAchieveableRate <- ev$maxAchieveableRate
    randAchieveableRate <-append(ev$randAchieveableRate,randAchieveableRate)
    constAchieveableRate <- ev$constAchieveableRate
    const1AchieveableRate <- ev$const1AchieveableRate
    const2AchieveableRate <- ev$const2AchieveableRate
    optmAchieveableRate <-ev$optmAchieveableRate
    maxheight <-ev$maxheight
    randheight <-append(ev$randheight, randheight)
    constheight <-ev$constheight
    const1height <-ev$const1height
    const2height <-ev$const2height
    save(history1,maxAchieveableRate,randAchieveableRate,constAchieveableRate,const1AchieveableRate,const2AchieveableRate,optmAchieveableRate,achieveableRate,maxheight,randheight,constheight,const1height,const2height,optmheight,height, file=file)
  
}

normaliseData = function(SINR,h,maxSINR,minSINR){#normaliseData(SINR = mcSINR[foo,j,],h=h[foo],maxSINR = mcSINR)
  
  #normalise observed power
  h =(h-minH)/(maxH-minH)
  
  mnfoo=minSINR#0#min(SINR)
  mxfoo=maxSINR#30#max(SINR)
  if(mxfoo>0){
    #    assP=(assP-mnfoo)/(mxfoo-mnfoo)
    SINR=(SINR-mnfoo)/(mxfoo-mnfoo)#(SINR - mean(SINR)) / sd(SINR)#(P-mnfoo)/(mxfoo-mnfoo)
  }#returning absolut values
  
  return(c(SINR,h))
}
inicialize = function(tipo){
  DISCOUNT <<- 0.99#}else{DISCOUNT <<- 0.9}
  #DISCOUNT <<- 0.99
  REPLAY_MEMORY_SIZE <<- 100  # How many last steps to keep for model training
  MINIBATCH_SIZE <<- 16  #64 How many steps (samples) to use for training
  # Exploration settings
  epsilon <<- 1
  EPSILON_DECAY <<- 0.9#99#.9
  MIN_EPSILON <<- 0.05
  # Initializing vectors that will be need to each mc trial
  epoch <<- 200#500
  #our performance metric. We measure the rolling average of the cumulative rate of the UAV link
  achieveableRate <<- c()
  history1 <<-c()
  #Vector to save the heights trhough the path
  height <<-c()
  
  #the replay memory for the training
  action_replay <<-vector(length=REPLAY_MEMORY_SIZE)
  reward_replay <<-vector(length=REPLAY_MEMORY_SIZE)
  replay_index <<-1
  ##UAV antenna gain
  episodes <<- 1 # 
  
  maxAchieveableRate <<-  c()
  randAchieveableRate <<- c()
  constAchieveableRate <<- c()
  const2AchieveableRate <<- c()
  const1AchieveableRate <<- c()
  optmAchieveableRate <<-c()
  maxheight <<-c()
  randheight <<-c()
  constheight <<-c()
  const1height <<-c()
  const2height <<-c()
  optmheight <<-c()
  
  neurons <<-18
  model <<- keras_model_sequential()
  model %>%
    layer_dense(units=neurons+1,input_shape = c(neurons)) %>%
    layer_dense(units = 200, activation = 'relu') %>%
    layer_dense(units = 200, activation = 'relu') %>%
    layer_dense(units = 200, activation = 'relu') %>%
    layer_dense(units = 3, activation = 'softplus')
  
  #model <- load_model_tf("model_Height_supervised")
  #if (m>1){model <- load_model_tf("model_Height_supervised")}
  
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



minH=1
maxH=11
h=seq(from=minH,to=maxH,by=1)

#Inicializing the simplest solution, here I will also calculate the benchmarks, in the other proposed solution we dont need to calculate the benchmarks anymore
for(k in 1:MCtrials){
  print(k)
inicialize("noint")
print("after inicialize noint")
model %>% compile(
  optimizer = 'adam',
  loss = 'sparse_categorical_crossentropy',
  metrics = c('accuracy'),
  #learning_rate=0.01
)
#calculation the optimal move
for (j in 1:steps){
  #optimum move
  optmheight[j] = which(mcSINR[,j]== max(mcSINR[,j]))
  optmAchieveableRate[j] = max(mcRate[,j])#log2(1+max(mcSINR[,j]))#mcSINR[foo,j]
}
optm_episode_reward = sum(optmAchieveableRate)
maxreward = max(optmAchieveableRate)
d=1

  #starting height
  #k=1
  
  #time1=Sys.time()
  #All the solutions start at the same height
  UAVh = minH#floor(runif(n=1,min=minH,max=maxH))

  episode_reward = 0


  

  if (k==1){#k=1
    randUAVh = minH
    genieUAVh = minH
    constUAVh = minH
    const1UAVh = minH
    const2UAVh = minH
    hminUAVh = minH
    optmUAVh = minH
    rand_episode_reward=0
    genie_episode_reward=0
    const_episode_reward=0
    const1_episode_reward=0
    const2_episode_reward=0
    hmin_episode_reward=0
    optm_episode_reward=0
  }
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
    #print(j)
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
    current_state = normaliseData(SINR = mcSINR[foo,j],h=h[foo],maxSINR = max(mcSINR),minSINR=min(mcSINR))

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
        if(UAVh<=(1+minH)){UAVh=minH}
        else {UAVh=UAVh-1}
      }else if(action==3){
        if((UAVh+1)>=maxH){UAVh = maxH
        }else{UAVh=UAVh+1}
      }
    }

    #get the new state and reward
    foo = which(h==UAVh)
    # if (mcSINR[foo,j]>0){
    #   mnfoo=0
    #   mxfoo=maxreward
      reward = mcRate[foo,j]/maxreward#(mcRate[foo,j]-mnfoo)/mxfoo#(log2(1+mcSINR[foo,j])-mnfoo)/mxfoo
    # }else{reward = 0}
    #
    new_state = normaliseData(SINR = mcSINR[foo,j],h=h[foo],maxSINR = max(mcSINR),minSINR=min(mcSINR))#normaliseDataComplete(BSdensity = BHdensity ,Builddensity =buildDens ,P=mcS[foo,j,],Dist=mcDist[foo,j,],h=h[foo],SINR = mcSINR[foo,j,])#,SINRmin = mcSINR[foo1,j+1,],SINRmax = mcSINR[foo2,j+1,] )
    new_state = c(new_state,current_state,action,reward, state_1,action_1,reward_1,state_2,action_2,reward_2,state_3,action_3,reward_3)

    ##add to the replay memory
    if (j>4){
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
        minibatch = sample(1:replay_index,size=MINIBATCH_SIZE,replace=FALSE)
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
    episode_reward=reward+episode_reward

    #store the rolling averages
    foo = which(h==UAVh)
    achieveableRate[j] = reward*maxreward
    height[j] = UAVh

    #decay the exploration factor
    if(epsilon > MIN_EPSILON){
      epsilon = epsilon*EPSILON_DECAY
      epsilon = max(MIN_EPSILON, epsilon)
    }
    if (k==1){

      if (const1UAVh < 5){const1UAVh = const1UAVh+1}
      if (const2UAVh < 10){const2UAVh = const2UAVh+1}
      randAction = sample(1:3, 1)#floor(runif(n=1,min=1,max=4))
      #random walk
      if(randAction==2){
        randUAVh=max(randUAVh-1,minH)
      }else if(randAction==3){
        randUAVh=min(randUAVh+1,maxH)
      }



      if(genieUAVh>optmheight[j+1]){
        if(genieUAVh<=(1+minH)){genieUAVh=minH
        }else {genieUAVh=genieUAVh-1}
      }else if(genieUAVh<optmheight[j+1]){
        if((genieUAVh+1)>=maxH){genieUAVh=maxH
        }else {genieUAVh=genieUAVh+1}
      }
      
      
      # #######Implementing the Wallid min solution#############################
      # #Assumption: 12 rbs
      # distance = distm(c(longitude_cell,latitude_cell), c( my_data$Longitude[1],my_data$Latitude[1]))#/1000
      # BStx = 97.723722096#in watts
      # #sumMeasured = my_data$RSRPdbm[j] / my_data$SINR[j]
      # RXpower = 10^( my_data$RSRPdbm[j]/10)/1000
      # SINR =  10^( my_data$SINR[j]/10)/1000
      # SINRThreshold=  0.001995262315 #3db #0.0079432823472# 9 db
      # interferenceThreshold= 3.1622776602e-12# -85dbm0.001995# 3.1622776602e-10# that they used on the results 0.001995262315 
      # denominatorMin =12*(4*pi*3.6e9/3e8)^2* interferenceThreshold # Cjs(t)*(4pif/c)^2*I
      # denominatorMax =12*(4*pi*3.6e9/3e8)^2* SINRThreshold
      # distanceSquared = distance^2
      # #sumFadingMin = 0
      # #sumFadingMax = rrice(n=1,vee=1.59,sigma=1)#fading of the link
      # 
      # #for (eita in 2:3){#the sum of the measured power from other BSs (interference) an their fading. - I am using the 2 strongest BS because on the paper they mention "# of interferers L = 2"
      # #  sumFadingMin = sumFadingMin+rrice(n=1,vee=1.59,sigma=1)   #rician parameter K = 1.59 [35] #1 if we assume no interference
      # #}
      # raiz =  (RXpower /denominatorMin) - distanceSquared
      # if (raiz > 0){
      #   hmin = sqrt(raiz)#as.integer()
      # }else {hmin = 0}
      # print(hmin) 
      # raiz = (SINR/denominatorMax) - distanceSquared
      # if (raiz > 0){
      #   hmax = sqrt(raiz)#as.integer()
      # }else {hmax = 0}
      # print(hmax) 
      # #Here hmin will always be the maximum height possible, in our case 120m const - the hmin values are at least 908m, being more if we consider the riccian fading.
      # #######################################################################################

      foo = which(h==randUAVh)
      #if(randUAVh<maxH & randUAVh>mi
      # if (mcSINR[foo,j]>0){
        randReward= mcRate[foo,j]#log2(1+mcSINR[foo,j])#mcSINR[foo,j]
      # }else{randReward = 0}

      foo = which(h==genieUAVh)
    #  if (mcSINR[foo,j]>0){
        genieReward= mcRate[foo,j]#log2(1+mcSINR[foo,j])#mcSINR[foo,j]
     # }else{genieReward= 0}

      foo = which(h==as.integer(constUAVh))
     # if (mcSINR[foo,j]>0){
        constReward= mcRate[foo,j]#log2(1+mcSINR[foo,j])#mcSINR[foo,j]
      #}else{constReward= 0}


      foo = which(h==as.integer(const1UAVh))
     # if (mcSINR[foo,j]>0){
        const1Reward= mcRate[foo,j]#log2(1+mcSINR[foo,j])#mcSINR[foo,j]
    #  }else{const1Reward= 0}

      foo = which(h==as.integer(const2UAVh))
     # if (mcSINR[foo,j]>0){
        const2Reward = mcRate[foo,j]#log2(1+mcSINR[foo,j])#mcSINR[foo,j]
    #  }else{const2Reward = 0}


      rand_episode_reward= randReward+rand_episode_reward
      genie_episode_reward = genieReward+ genie_episode_reward
      const_episode_reward = constReward+ const_episode_reward
      const1_episode_reward = const1Reward+ const1_episode_reward
      const2_episode_reward = const2Reward+ const2_episode_reward

      randAchieveableRate[j] = randReward
      maxAchieveableRate[j] = genieReward
      constAchieveableRate[j] = constReward
      const1AchieveableRate[j] = const1Reward
      const2AchieveableRate[j] = const2Reward

      randheight[j] = randUAVh
      maxheight[j] = genieUAVh
      constheight[j] = constUAVh
      const1height[j] = const1UAVh
      const2height[j] = const2UAVh

    }#other solutions
    if(j%%5==0){
      #once every 5 episodes update the targetModel to have the same weights as the model
      targetModel  %>% set_weights(model %>% get_weights())
    }
    
  
  }#steps
 
  print(mean(achieveableRate))
  if (length(height)== (steps-1)){
    # if (k==1){  save(k,history1,maxAchieveableRate,randAchieveableRate,constAchieveableRate,const1AchieveableRate,const2AchieveableRate,optmAchieveableRate,achieveableRate,maxheight,randheight,constheight,const1height,const2height,optmheight,height,file="Real_canal_tudo1.Rdata")
    
    #} else {
    save_append(k,history1,maxAchieveableRate,randAchieveableRate,constAchieveableRate,const1AchieveableRate,const2AchieveableRate,optmAchieveableRate,achieveableRate,maxheight,randheight,constheight,const1height,const2height,optmheight,height,file="Real_canal_tudo1.Rdata")
    #}
  }
}#episodes



#load("Real_grancanal_result.Rdata")
# steps = 171
# #Calculating the const 2 again because I did wrong on the oficial reults
# const2UAVh = minH
# const2_episode_reward=0
# for(j in 1:(steps-1)){
#   if (const2UAVh < 11){const2UAVh = const2UAVh+1}
#   foo = which(h==as.integer(const2UAVh))
#   # if (mcSINR[foo,j]>0){
#   const2Reward = mcRate[foo,j]#log2(1+mcS
#   const2_episode_reward = const2Reward+ const2_episode_reward
#   const2AchieveableRate[j] = const2Reward
#   const2height[j] = const2UAVh
# }
# 
# 
# real_height = c(20,30,40,50,60,70,80,90,100,110,120)
# for (step in 1:(steps-1)){
#   for (alt in 1:length(real_height)){
#     if(height[step]==alt){
#       #print(alt)
#       height[step]=real_height[alt]}
#     if(optmheight[step]==alt){optmheight[step]=real_height[alt]}
#     if(maxheight[step]==alt){maxheight[step]=real_height[alt]}
#     if(randheight[step]==alt){randheight[step]=real_height[alt]}
#     if(constheight[step]==alt){constheight[step]=real_height[alt]}
#     if(const1height[step]==alt){const1height[step]=real_height[alt]}
#     if(const2height[step]==alt){const2height[step]=real_height[alt]}
#   }
# }
# optmheight =  head( optmheight, -1)
# 
print(mean(achieveableRate))
print(mean(optmAchieveableRate))
print(mean(maxAchieveableRate))
print(mean(randAchieveableRate))
print(mean(constAchieveableRate))
print(mean(const1AchieveableRate))
print(mean(const2AchieveableRate))
print(sd(achieveableRate))
print(sd(optmAchieveableRate))
print(sd(maxAchieveableRate))
print(sd(randAchieveableRate))
print(sd(constAchieveableRate))
print(sd(const1AchieveableRate))
print(sd(const2AchieveableRate))

print(sum(diff(height) != 0))
print(sum(diff(optmheight) != 0))
print(sum(diff(maxheight) != 0))
print(sum(diff(randheight) != 0))
print(sum(diff(constheight) != 0))
print(sum(diff(const1height) != 0))
print(sum(diff(const2height) != 0))
# 
# 
# # Give the chart file a name.
# png(file = "GranCanalPath.png")
# 
# # Plot the bar chart.
# # Plot the bar chart.
# plot(height[1:100],type = "S",col = "red", xlab = "Step", ylab = "Height (m)",  ylim=c(10, 160))
# lines(constheight[1:100], type = "S", col = "green")
# lines(maxheight[1:100], type = "S", col = "blue")
# #lines(const1height, type = "S", col = "darkgreen")
# #lines(const2height, type = "S", col = "lightgreen")
# lines(randheight[1:100], type = "S", col = "purple")
# lines(optmheight[1:100], type = "S", col = "orange")
# legend("topright", col=c("blue","green","red","purple","orange"),
#        legend=c("One-step-ahead","S. Zhang","RL","Random Walk","Optimal movement"), lty = 1:1, cex=0.8)
# 
# # Save the file.
# dev.off()















save(achieveableRate,height,history,action_replay,const_episode_reward,episode_reward,const1_episode_reward,const2_episode_reward,const1AchieveableRate,constAchieveableRate,const2AchieveableRate,maxAchieveableRate,
     randAchieveableRate,randheight,maxheight,genieReward,constheight,const1height,const2height,minibatch,reward_replay,file="Real_grancanal_result1.Rdata")
