# ancillary_time.py
#
# Short Python routine to find time stamp of GDAS or GFS file
#
# Virginia Sawyer
# Created 2019-08-08
# Updated 2019-08-08
 
import sys
from datetime import datetime,timedelta
 
# Check for correct number of arguments
if len(sys.argv) != 5:
    print('Usage: python ancillary_time.py YYYYMMDD HH MM dt | YYYYJJJ HH MM dt')
    sys.exit()
 
# Arguments to variables
ydate = sys.argv[1]
hh = sys.argv[2]
mm = sys.argv[3]
dt = float(sys.argv[4])
 
# Convert arguments to datetime
if len(ydate) == 8:
    tstamp = datetime.strptime(ydate+' '+hh+' '+mm,'%Y%m%d %H %M')
if len(ydate) == 7:
    tstamp = datetime.strptime(ydate+' '+hh+' '+mm,'%Y%j %H %M')
 
# Find 6-hourly interval for GDAS/GFS file
dhr = tstamp.hour + tstamp.minute/60.0
#ti = round(dhr/6.)*6 
ti = round(dhr/dt)*int(dt)
atime = tstamp+timedelta(hours=(ti-dhr))
# Write time stamp to ASCII file
#with open('ancillary_time','w') as outfile:
#    outfile.write(datetime.strftime(atime,'%Y %m %d %H'))
with open('gdas_basename','w') as outfile:
    outfile.write(datetime.strftime(atime,'%Y/%m/gdas1.PGrbF00.%y%m%d.%Hz'))
    