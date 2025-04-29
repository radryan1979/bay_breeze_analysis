#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 23:02:34 2018

@author: allenea
"""


"""
From James:
https://github.com/jsimkins2/UD_SRS/blob/master/goesR/nexrad_kdox.py
"""

from siphon.radarserver import RadarServer
from datetime import datetime,timedelta
import cartopy
import matplotlib
import numpy as np
from siphon.cdmr import Dataset
import warnings
warnings.filterwarnings("ignore", category=matplotlib.cbook.MatplotlibDeprecationWarning)
from metpy.plots import ctables  # For NWS colortable
import matplotlib.pyplot as plt
import pyart
from dateutil import tz
import numpy.ma as ma
import time
from time import mktime
import pyproj
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl

#select the radar site
site = 'KDOX'
outdir = '/Users/allenea/Documents/Eric_Allen/Ferry_Data/RADAR/PyCase1/'
ref_norm, ref_cmap = ctables.registry.get_with_steps('NWSReflectivity', 5, 5)

rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')

query=rs.query()
dt=datetime(2014, 6, 4, 12)
query.stations('KDOX').time_range(dt,dt+timedelta(hours=15))
rs.validate_query(query)

cat = rs.get_catalog(query)
cat.datasets

#%% Plotting Set-up

#set up the plotting parameters (NWSReflectivity colormap, contour levels,
# and colorbar tick labels)

cmap = 'pyart_NWSRef'
levs = np.linspace(0,80,41,endpoint=True)
ticks = np.linspace(0,80,9,endpoint=True)
label = 'Radar Reflectivity Factor ($\mathsf{dBZ}$)'
#normalize the colormap based on the levels provided above
norm = mpl.colors.BoundaryNorm(levs,256)
same10Min = ""
#count = 0
for scan in range(0,len(cat.datasets),4):
    ds = list(cat.datasets.values())[scan]
    #break
    #print str(ds)[13:16],"  ",same10Min
    if str(ds)[13:16] != same10Min:
        same10Min = str(ds)[13:16]
        #print "USING THIS SCAN", ds
        #count +=1
        
        loc = pyart.io.nexrad_common.get_nexrad_location(site)
        lon0 = loc[1] ; lat0 = loc[0]
        
        radar = pyart.io.read_nexrad_cdm(ds.access_urls['OPENDAP'])
        
        # create timestamp
        timestamp = radar.time['units'].split(' ')[-1].split('T')
        timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
        timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        from_zone = tz.gettz('UTC')
        #to_zone = tz.gettz('America/New_York')
        utc = timestamp.replace(tzinfo=from_zone)
        dt= utc.strftime('%Y-%m-%d %H:%M:%S %Z')
        
        #local = utc.astimezone(to_zone)
        #lt = time.localtime()
        #dst = lt.tm_isdst
        #lt = time.localtime()
        #dst = lt.tm_isdst
        #if dst == 0:
        #    et = "EDT"
        #else:
        #    et = "EST"
            
        # Grid data into geographic coordinates
        display = pyart.graph.RadarMapDisplay(radar)
        x,y = display._get_x_y(0,True,None)
        
        
        #set up a 1x1 figure for plotting
        fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(9,9),dpi=100)
        #set up a basemap with a lambert conformal projection centered 
        # on the radar location, extending 1 degree in the meridional direction
        # and 1.5 degrees in the longitudinal in each direction away from the 
        # center point.
        m = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
                   llcrnrlat=lat0-1.25,llcrnrlon=lon0-1.5,
                   urcrnrlat=lat0+1.25,urcrnrlon=lon0+1.5,resolution='h')
        
        #get the plotting grid into lat/lon coordinates
        x0,y0 = m(lon0,lat0)
        glons,glats = m((x0+x*1000.), (y0+y*1000.),inverse=True)
        #read in the lowest scan angle reflectivity field in the NEXRAD file 
        refl = np.squeeze(radar.get_field(sweep=0,field_name='reflectivity'))
        #define the plot axis to the be axis defined above
        ax = axes
        #create a colormesh of the reflectivity using with the plot settings defined above
        cs = m.pcolormesh(glons,glats,refl,norm=norm,cmap=cmap,ax=ax,latlon=True)
        #add geographic boundaries and lat/lon labels
        m.drawparallels(np.arange(20,70,0.5),labels=[1,0,0,0],fontsize=12,
                        color='k',ax=ax,linewidth=0.001)
        m.drawmeridians(np.arange(-150,-50,1),labels=[0,0,1,0],fontsize=12,
                       color='k',ax=ax,linewidth=0.001)
        m.drawcounties(linewidth=0.5,color='gray',ax=ax)
        m.drawstates(linewidth=1.5,color='k',ax=ax)
        m.drawcoastlines(linewidth=1.5,color='k',ax=ax)
        #mark the radar location with a black dot
        m.scatter(lon0,lat0,marker='o',s=20,color='k',ax=ax,latlon=True)
        #add the colorbar axes and create the colorbar based on the settings above
        cax = fig.add_axes([0.075,0.075,0.85,0.025])
        cbar = plt.colorbar(cs,ticks=ticks,norm=norm,cax=cax,orientation='horizontal')
        cbar.set_label(label,fontsize=12)
        cbar.ax.tick_params(labelsize=11)
        #add a title to the figure
        fig.text(0.5,0.92, site + ' (0.5$^{\circ}$) Reflectivity\n ' + 
                 dt ,horizontalalignment='center',fontsize=16)
        #save and display the figure
        plt.savefig(outdir+'Radar'+site+'_'+dt+'.png')
        #plt.show()
        plt.close()
        #os.remove(localfile.name)
    else:
        same10Min = str(ds)[13:16]
