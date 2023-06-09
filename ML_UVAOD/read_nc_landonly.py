#import argparse
from sklearn.ensemble import RandomForestRegressor
import pandas as pd
import os
from datetime import datetime
#import plotting.boxPlot
#import plotting.mapPlot
#import stats.stats
import numpy as np
import pickle
import glob
import netCDF4 as nc4
import h5py 

configfile='../run/OCIUAAER_Config'
with open(configfile) as f:
    for line in f:
        if 'ml_model=' in line:
            modelpath=line[10:-1]

print (modelpath)
#modelpath='../ML_UVAOD/data/'
#landmodel=['landmodel_340.pkl']
#oceanmodel=['oceanmodel_340.pkl']
modelname='landmodel_both_nogeo_'

path='./'

ncfiles=glob.glob(path+'Interm_file.nc')
#ncfiles=glob.glob(path+'Interm_L2*.nc')



def loadModel(path):

    if os.path.exists(path):
        with open(path, 'rb') as f:
            rf = pickle.load(f)
            return rf
    else:
        print ('cannot find ML model: '+path)
    return None

for ncfile in ncfiles:
 #   print (ncfile)
#    date=datetime.strptime(ncfile[-15:-8],"%Y%j")
#    print (date.month)
    ifile=nc4.Dataset(ncfile, "r+", format='NETCDF4')
    geo=ifile['geolocation_data']
    month=geo.variables['month'][:]
    month=float(month[0])
#    print (month)
     
    slats = geo.variables['latitude'][:]
    slons = geo.variables['longitude'][:]
    slats = slats[:,:].astype(float).flatten()
    slons = slons[:,:].astype(float).flatten()
  #  print (slats.max(),slats.min())
    dimhdf=slats.shape
  #  print (slats.shape)
     
    gph=ifile['geophysical_data']
    landaod=gph.variables['Optical_Depth_Land'][:]
#    print (landaod.shape)
    dim1name=(ifile.dimensions['number_of_lines_8x8'].name)
    dim2name=(ifile.dimensions['number_of_pixels_8x8'].name)

    ifile.close()
    landaod=landaod[:,:,:].astype(float)
    dimhdf3=landaod.shape
#    print (dimhdf3)#5, 340,380,470,550,670 only take last 3wav
    dfland=pd.DataFrame()
    for i in range(0,dimhdf3[2]-2):
#        print ('AOD'+str(i))
        if i == 0:
            wavname='470nm'
        if i == 1:
            wavname='550nm'
        if i == 2:
            wavname='670nm'
        dfland['DB_'+wavname]=landaod[:,:,i+2].flatten()
    final = None
    df=pd.DataFrame()
#    print ('land',month)
#    features = ['month', 'Corrected_Optical_Depth_Land1', 'Corrected_Optical_Depth_Land2', 'Corrected_Optical_Depth_Land3', 'MOD-LAT', 'MOD-LON']
    df=pd.concat([df,dfland])
    df['month']=month
    cols=['month']+[col for col in dfland]
    df=df[cols]
    #df['lat']=slats
    #df['lon']=slons
    df=df.replace(-9999.0,np.nan)
    x= df.dropna(how='any')
#    print (x)
    out=x.copy(deep=True)
    for wavelength in [340,380]:
        print (modelpath+modelname + str(wavelength) + ".pkl")
        regr = loadModel(modelpath+modelname + str(wavelength) + ".pkl")
        out['predicted_land_'+str(wavelength)] = regr.predict(x)
    out=out.reindex(list(range(0,dfland.index.max()+1)),fill_value=-9999.)
#        print (out)
    cols=[col for col in out if  'predicted' not in col]
    #print (cols)
    if final is None: 
        final=out.drop(cols,axis=1)
    else:
        final=pd.concat([final,out.drop(cols,axis=1)],axis=1)
#    print (final)
    '''
    <class 'netCDF4._netCDF4.Variable'>
    int16 Corrected_Optical_Depth_Land(number_of_lines_8x8, number_of_pixels_8x8, Wavelength_Used_Land_2)
    valid_range: [ -50 5000]
    _FillValue: -9999
    long_name: Retrieved AOT at 0.48, 0.55, 0.67, 2.25 microns
    units: None
    scale_factor: 0.001
    add_offset: 0.0
    Parameter_Type: Output
    Geolocation_Pointer: Internal geolocation arrays
    coordinates: /geolocation_data/longitude /geolocation_data/latitude
    path = /geophysical_data
    unlimited dimensions: number_of_lines_8x8
    current shape = (404, 400, 4)
    filling on
    '''
    fid=nc4.Dataset(ncfile,'r+')
    for iw,wavelength in enumerate([340,380]):
        temp=final['predicted_land_'+str(wavelength)].to_numpy()
        ttp=np.reshape(temp,(dimhdf3[0],dimhdf3[1]))
#        print ('predicted_land_'+str(wavelength))
#            gph.createDimension(dim1name,(dim1))
#            gph.createDimension(dim2name,(dim2))
#            var=gph.createVariable('predicted_'+lo+'_'+str(wavelength),\
#                    'f',(dim1name,dim2name))
#            var[:]=ttp[:,:]
        gph=fid['geophysical_data']

        gph['Optical_Depth_Land'][:,:,iw]=ttp
        #fid.create_dataset('predicted_land_'+str(wavelength),data=ttp)
    fid.close()
    print ('Interm_file changed')
#    exit(1)
    os.rename(ncfile,path+'Interm_file.nc')


