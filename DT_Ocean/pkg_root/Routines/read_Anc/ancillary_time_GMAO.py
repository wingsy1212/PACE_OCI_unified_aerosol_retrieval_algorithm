# ancillary_time_GMAO.py
#
# Short Python routine to find GMAO file with correct time stamp
#
# Virginia Sawyer
# Created 2019-08-08
# Updated 2022-06-28
 
import sys
from datetime import datetime,timedelta
from glob import glob

# Check for correct number of arguments
if len(sys.argv) != 6:
    print('Usage: python ancillary_time_GMAO.py input_path YYYYMMDD HH MM dt | YYYYJJJ HH MM dt')
    sys.exit()
 
# Arguments to variables
apath = sys.argv[1]
ydate = sys.argv[2]
hh = sys.argv[3]
mm = sys.argv[4]
dt = float(sys.argv[5])
 
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
print(atime)
# Write time stamp to ASCII file
#with open('ancillary_time','w') as outfile:
#    outfile.write(datetime.strftime(atime,'%Y %m %d %H'))
#with open('GMAO_basename','w') as outfile:
#    outfile.write(datetime.strftime(atime,'%Y/%m/gdas1.PGrbF00.%y%m%d.%Hz'))

aname = datetime.strftime(atime,'GEOS.fpit.asm.inst3_2d_asm_Nx.GEOS5124.%Y%m%d_%H%M.')
fname = glob(apath+datetime.strftime(atime,'/%Y/%j/')+aname+'*.nc4')[0]

with open('input_anc_file','w') as outfile:
    outfile.write(fname)
    

                                   
    
    
    
   
   
