import math
import pandas as pd
import os
#iN CASE OF RUNNING IN THE GRANCANNAL DATA
#path=os.path.dirname(os.getcwd())#os.getcwd()
#change = path+'/Grand Canal Docks Dedicated throughput'
#os.chdir(change)


#df_sheet_name = pd.read_excel('30m.xls', sheet_name='Metric Group 1')
files = [f for f in os.listdir('.') if 'mdedi.xls'in f ]#'mdedicated.xls'
dic_file={}
for file in files:
    #print(file)
    dic_file[file]=pd.read_excel(file,sheet_name='Metric Group 1')
    dic_file[file] = dic_file[file].drop(['Unnamed: 0','Time','Date','Serving Cell Identity','Cell Identity: Top #1'],1) #['Unnamed: 0','Time','Date'],1) #

#steps=160 # minimum value from the 60m file
base = '60mdedi.xls'#m.xls'
dic_final = {}
dic_final[base]=dic_file[base]
for file in dic_file:
    if file != base:
        temp=[]
        df_temp = pd.DataFrame(columns=dic_file[base].columns)
        for a in range(len(dic_file[base])):
            lign = 100#just to start with a big number

            for b in range(len(dic_file[file])):
                dist = math.sqrt((dic_file[base].Latitude[a] - dic_file[file].Latitude[b]) ** 2 + (
                            dic_file[base].Longitude[a] - dic_file[file].Longitude[b]) ** 2)
                if dist < lign:
                    lign = dist
                    smaller = b
            #order of the data

            df_temp = pd.concat([df_temp,dic_file[file][smaller:smaller+1]], ignore_index=True)
    df_temp=df_temp.rename(columns={'RS SINR Carrier 1 (dB)':'SINR','Serving Cell RSRP (dBm)':'RSRPdbm','Serving Cell RSRQ (dB)':'RSRQdb' })#'PDSCH Phy Throughput (kbps)':'Rate' })#
    dic_final[file]=df_temp
    name = file[0:-8:1]+'sinr.xls'#file[0:-14:1]+'n.xls'#'dedi.xls'#
    df_temp.to_excel(name)
    #df_temp.to_excel(name, sheet_name='sheet1', index=False)
        #dist = math.sqrt((dic_file['60m.xls'].Latitude[0] - dic_file['40m.xls'].Latitude[0]) ** 2 + (dic_file['60m.xls'].Longitude[0] - dic_file['40m.xls'].Longitude[0])** 2 )

        #dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
