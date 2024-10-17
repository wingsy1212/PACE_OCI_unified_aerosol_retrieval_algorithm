	# ancillary_time_2GMAO.py
#
# Short Python routine to find two GMAO files bracketing a given time stamp
#
# Virginia Sawyer
# Created 2019-08-08
# Updated 2024-04-29
 
import sys,math
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

# Find dt-hourly interval for ancillary file
# This version prints two lines, one with the GMAO file before the input timestamp
# and one with the GMAO file after the input timestamp
dhr = tstamp.hour + tstamp.minute/60.0
ti0 = math.floor(dhr/dt) * int(dt)
ti1 = math.ceil(dhr/dt) * int(dt)

atime0 = tstamp+timedelta(hours=(ti0-dhr))
atime1 = tstamp+timedelta(hours=(ti1-dhr))
print (atime0)
print (atime1)
 
#aname0 = datetime.strftime(atime0,'GEOS.fpit.asm.inst3_2d_asm_Nx.GEOS5124.%Y%m%d_%H%M.')
#aname1 = datetime.strftime(atime1,'GEOS.fpit.asm.inst3_2d_asm_Nx.GEOS5124.%Y%m%d_%H%M.') 
aname0 = datetime.strftime(atime0,'GMAO_MERRA2.%Y%m%dT%H%M') 
aname1 = datetime.strftime(atime1,'GMAO_MERRA2.%Y%m%dT%H%M')
print (aname0)
print (aname1) 
print ((apath+datetime.strftime(atime0,'/%Y/%m/%d/')+aname0+'*.nc'))
try:
    fname0 = glob(apath+datetime.strftime(atime0,'/%Y/%m/%d/')+aname0+'*.nc')[0]
    print(fname0)
    with open('input_anc_file1','w') as outfile:
        outfile.write(fname0+'\n')
except:
    print('Cannot find '+aname0)
 
try:
    fname1 = glob(apath+datetime.strftime(atime1,'/%Y/%m/%d/')+aname1+'*.nc')[0]
    print(fname1)
    with open('input_anc_file2','w') as outfile:
        outfile.write(fname1+'\n')
except:
    print('Cannot find '+aname1)
 
 