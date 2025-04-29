"""
http://nbviewer.jupyter.org/gist/dopplershift/356f2e14832e9b676207
"""


from siphon.radarserver import RadarServer
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
import cartopy
#from metpy.plots import ctables
import numpy as np
from siphon.cdmr import Dataset
import warnings
warnings.filterwarnings("ignore", category=matplotlib.cbook.MatplotlibDeprecationWarning)




def raw_to_masked_float(var, data):
    # Values come back signed. If the _Unsigned attribute is set, we need to convert
    # from the range [-127, 128] to [0, 255].
    if var._Unsigned:
        data = data & 255

    # Mask missing points
    data = np.ma.array(data, mask=data==0)

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



rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')

query=rs.query()
dt=datetime(2017, 5, 1)
query.stations('KDOX').time_range(dt,dt+timedelta(hours=1))

cat = rs.get_catalog(query)
cat.datasets

"""
Grab the first dataset so that we can get the longitude and 
latitude of the station and make a map for plotting. We'll go 
ahead and specify some longitude and latitude bounds for the map.
"""
ds = list(cat.datasets.values())[0]
data = Dataset(ds.access_urls['CdmRemote'])
fig = plt.figure(figsize=(10, 10))
ax = new_map(fig, data.StationLongitude, data.StationLatitude)

# Set limits in lat/lon space
ax.set_extent([-78, -73.5, 37, 40.5])

# Add ocean and land background
ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean', scale='50m',
                                            edgecolor='face',
                                            facecolor=cartopy.feature.COLORS['water'])
land = cartopy.feature.NaturalEarthFeature('physical', 'land', scale='50m',
                                           edgecolor='face',
                                           facecolor=cartopy.feature.COLORS['land'])

ax.add_feature(ocean, zorder=-1)
ax.add_feature(land, zorder=-1)