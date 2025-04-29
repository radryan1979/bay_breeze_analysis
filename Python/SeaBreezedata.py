#John Lodise
#Sea breeze database
#YY  MM DD hh mm WDIR WSPD GST PRES  ATMP  WTMP
#yr  mo dy hr mn degT m/s  m/s  hPa  degC  degC 
import numpy as np              #imports python libraries
import matplotlib.pyplot as pl

brnd12 = open("brnd1h2012.txt","r").readlines() #Reads in txt file and counts lines



numline = len(brnd12)    #finds number of lines and saves in variable
Wtmp = np.zeros(numline)  # sets up array of 0's as long as txt file lines
Wdir = np.zeros(numline)
time = np.zeros(numline)
y=0
bJWdir = np.zeros(10*24) #sets up array of 0's for one day of data from brandywine in June
bJtime = np.zeros(10*24)
for k in range(numline):   #loops through everyline of file
    ob = brnd12[k]              #looks at each line individually 
    year = int(ob.split()[0])
    month = float(ob.split()[1])
    day = float(ob.split()[2])
    hour = float(ob.split()[3])                 #splits every line by spaces into variables
    minute = float(ob.split()[4])
    wdir = float(ob.split()[5])
    if wdir == 999.0:                           
        wdir = None        # assigns NaN to missing values
    Wdir[k] = wdir          # fills array for entire year
    wspd = ob.split()[6]
    pres = float(ob.split()[12])
    atmp = ob.split()[13]
    wtmp = float(ob.split()[14])
    if wtmp == 999.0:
        wtmp = None
    Wtmp[k] = wtmp  
    time[k] = year + (month-1)/12.0 + (day-1.0)/365.25 

    if month == 6 and day  == 8:        # if june 8th, puts wind direction and time into array
        bJWdir[y] = wdir
        bJtime[y] = (hour/24) + (minute/3600) #cumulative time....might be a better way
        y=y+1  #moves onto next indice


cpmy12 = open("capemay4h2012edit.txt","r").readlines()


nline = len(cpmy12)

cmWtmp = np.zeros(nline)
cmtime = np.zeros(nline)
cmWdir = np.zeros(nline)
cJWdir = np.zeros(10*24)
Jtime = np.zeros(10*24)
x=0

for k in range(nline):
    cm = cpmy12[k]
    cmyear = int(cm.split()[0])
    cmmonth = float(cm.split()[1])
    cmday = float(cm.split()[2])
    cmhour = float(cm.split()[3])
    cmminute = float(cm.split()[4])
    cmwdir = float(cm.split()[5])
    if cmwdir == 999.0:
       cmwdir = None
    cmWdir[k] = cmwdir
    cmwspd = cm.split()[6]
    cmpres = float(cm.split()[8])
    cmatmp = cm.split()[9]
    cmwtmp = cm.split()[10]
    if cmwtmp == 999.0:
        cmwtmp = None
    cmWtmp[k] = cmwtmp
    cmtime[k] = cmyear + (cmmonth-1)/12.0 + (cmday-1.0)/365.25
    
    if cmmonth == 6 and cmday  == 8:
        cJWdir[x] = cmwdir
        Jtime[x] = (cmhour/24) + (cmminute/3600)
        x=x+1

print cJWdir  #check for array accuracy 
print Jtime


pl.plot(Jtime,cJWdir,'b',bJtime,bJWdir,'g')    #combines arrays from both textfiles in one plot
pl.axis(xmin= 0,xmax=1,ymin=0,ymax=370)
pl.xlabel("Time")
pl.ylabel("Wind Direction (degrees)")
pl.title("July 3rd Seabreeze from Capemay and Brandywine")

#pl.plot(cmtime,cmWdir,'r',cmtime,cmWtmp,'b')
pl.show()
