from siphon.radarserver import RadarServer
from datetime import datetime,timedelta
import pyart




rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')

query=rs.query()
dt=datetime(2017, 6, 1, 16)
query.stations('KDOX').time_range(dt,dt+timedelta(minutes=15))

cat = rs.get_catalog(query)
raw_list = list(cat.catalog_refs)
ds = list(cat.datasets.values())[0]

radar = pyart.io.read_nexrad_cdm(ds.access_urls['OPENDAP'])

ref = radar.get_field(0, 'reflectivity')

