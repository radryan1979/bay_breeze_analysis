"""
Source: http://nbviewer.jupyter.org/gist/dopplershift/356f2e14832e9b676207
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



"""
def raw_to_masked_float(var, data):
    # Values come back signed. If the _Unsigned attribute is set, we need to convert
    # from the range [-127, 128] to [0, 255].
    if var._Unsigned:
        data = data & 255

    # Mask missing points
    data = np.ma.array(data, mask=data<0 )
    data = np.ma.array(data, mask=data>80)

    # Convert to float using the scale and offset
    return data * var.scale_factor + var.add_offset

def polar_to_cartesian(az, rng):
    az_rad = np.deg2rad(az)[:, None]
    x = rng * np.sin(az_rad)
    y = rng * np.cos(az_rad)
    return x, y

def new_map(fig, lon, lat):
    # Create projection centered on the radar. This allows us to use x
    # and y relative to the radar.
    proj = cartopy.crs.LambertConformal(central_longitude=lon, central_latitude=lat)

    # New axes with the specified projection
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    # Add coastlines
    ax.coastlines('50m', 'black', linewidth=2, zorder=2)

    # Grab state borders
    state_borders = cartopy.feature.NaturalEarthFeature(
        category='cultural', name='admin_1_states_provinces_lines',
        scale='50m', facecolor='none')
    ax.add_feature(state_borders, edgecolor='black', linewidth=1, zorder=3)
    
    return ax
"""
    
def togrid(polar, x, y, gridsize=1024, lim=460):
    src = np.column_stack((x.ravel(), y.ravel()))
    grid = np.linspace(-lim * 1000, lim * 1000, gridsize)
    mgrid = np.meshgrid(grid, grid[::-1])
    trg = np.column_stack((mgrid[0].ravel(), mgrid[1].ravel()))
    grid = ipol_nearest(src, trg, polar.ravel())
    return grid.reshape(gridsize, gridsize)



ref_norm, ref_cmap = ctables.registry.get_with_steps('NWSReflectivity', 5, 5)

rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')

query=rs.query()
dt=datetime(2017, 6, 1, 16)
query.stations('KDOX').time_range(dt,dt+timedelta(minutes=15))

cat = rs.get_catalog(query)
cat.datasets


"""
From James:
https://github.com/jsimkins2/UD_SRS/blob/master/goesR/nexrad_kdox.py
"""
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]

radar = pyart.io.read_nexrad_cdm(ds.access_urls['OPENDAP'])

# create timestamp
timestamp = radar.time['units'].split(' ')[-1].split('T')
timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
from_zone = tz.gettz('UTC')
to_zone = tz.gettz('America/New_York')
utc = timestamp.replace(tzinfo=from_zone)
local = utc.astimezone(to_zone)
lt = time.localtime()
dst = lt.tm_isdst
lt = time.localtime()
dst = lt.tm_isdst
if dst == 0:
    et = "EDT"
else:
    et = "EST"
    
# Handle data 
display = pyart.graph.RadarMapDisplay(radar)
x,y = display._get_x_y(0,True,None)
x0,y0 = kdoxH(lon0,lat0)
glons,glats = kdoxH((x0+x*1000.), (y0+y*1000.),inverse=True)

# Reflectivity
ref = radar.get_field(0, 'reflectivity')
x, y, _ = radar.get_gate_x_y_z(0)
grid = togrid(ref, x, y, gridsize=256, lim=250)
dBZ = np.ma.masked_where(clutter, ref)

#Define transect lines
translat=[[39.15,   39.15],
          [39.125,  39.125],
          [39.1,    39.1],
          [39.075,  39.075],
          [39.05,   39.05],
          [39.025,  39.025],
          [39.0,    39.0],
          [38.975,  38.975],
          [38.95,   38.95],
          [38.925,  38.925],
          [38.9,    38.9],
          [38.875,  38.875],
          [38.85,   38.85],
          [38.825,  38.825],
          [38.8,    38.8],
          [38.775,  38.775],
          [38.75,   38.75],#Begin ocean transects
          [38.725,  38.725],
          [38.7,    38.7],
          [38.675,  38.675],
          [38.65,   38.65],
          [38.625,  38.625],
          [38.6,    38.6],
          [38.575,  38.575],
          [38.55,   38.55],
          [38.525,  38.525],
          [38.5,    38.5],
          [38.475,  38.475],
          [38.45,   38.45],
          ]
translon=[[-75.25,      -76.0],
          [-75.2375,    -76.0],
          [-75.225,     -76.0],
          [-75.2125,    -76.0],
          [-75.2,       -76.0],
          [-75.1875,    -76.0],
          [-75.175,     -76.0],
          [-75.1625,    -76.0],
          [-75.15,      -76.0],
          [-75.1375,    -76.0],
          [-75.125,     -76.0],
          [-75.1125,    -76.0],
          [-75.1,       -76.0],
          [-75.0875,    -76.0],
          [-75.075,     -76.0],
          [-75.0625,    -76.0],
          [-75.05,      -76.0],
          [-74.7,       -76.0],#Begin ocean transects
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          [-74.7,       -76.0],
          ]

"""
Grab the first dataset so that we can get the longitude and 
latitude of the station and make a map for plotting. We'll go 
ahead and specify some longitude and latitude bounds for the map.
"""
"""Option of not using pyart:
ds = list(cat.datasets.values())[0]
data = Dataset(ds.access_urls['CdmRemote'])
fig = plt.figure(figsize=(10, 10))
ax = new_map(fig, data.StationLongitude, data.StationLatitude)



# Set limits in lat/lon space
ax.set_extent([-77.25, -73.75, 37.75, 40.5])

# Add ocean and land background
ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean', scale='50m',
                                            edgecolor='face',
                                            facecolor=cartopy.feature.COLORS['water'])
land = cartopy.feature.NaturalEarthFeature('physical', 'land', scale='50m',
                                          edgecolor='face',
                                          facecolor=cartopy.feature.COLORS['land'])

ax.add_feature(ocean, zorder=-1)
ax.add_feature(land, zorder=-1)


for i in range(len(translat)):
    ax.plot([translon[i][0], translon[i][1]], [translat[i][0], translat[i][1]],
         color='blue', linewidth=1, 
         transform=cartopy.crs.Geodetic(),
         )

#Convert to x,y grid
#transx, transy = geog_to_cartesian(translon, translat)




meshes = []
for item in sorted(cat.datasets.items()):
    # After looping over the list of sorted datasets, pull the actual Dataset object out
    # of our list of items and access over CDMRemote
    ds = item[1]
    data = Dataset(ds.access_urls['CdmRemote'])

    # Pull out the data of interest
    sweep = 0
    rng = data.variables['distanceR_HI'][:]
    az = data.variables['azimuthR_HI'][sweep]
    ref_var = data.variables['Reflectivity_HI']

    # Convert data to float and coordinates to Cartesian
    ref = raw_to_masked_float(ref_var, ref_var[sweep])
    x, y = polar_to_cartesian(az, rng)
    
    # Find SB in x,y grid
    # sb = find_sea_breeze(transx, transy, x, y, ref)

    # Plot the data and the timestamp
    mesh = ax.pcolormesh(x, y, ref, cmap=ref_cmap, norm=ref_norm, zorder=0)
    text = ax.text(0.65, 0.03, data.time_coverage_start, transform=ax.transAxes,
                  fontdict={'size':16})
    
    # Collect the things we've plotted so we can animate
    meshes.append((mesh, text))
"""