'''
* NAME: files2nc.py
*
* DESCRIPTION: Generates UV netcdf LUT from list of input files.

* USAGE: files2nc.py -i [input file] -o [output file] 

* Created on December, 2023

* Author: Samuel Anderson
'''

import os
import sys
import argparse
import xarray as xr
import datatree as dtree
import configparser as cp
import numpy as np
from pyhdf.SD import SD, SDC
from scipy.io import FortranFile
from afrt import afrt_pace, afrt_viirs_ocean, afrt_viirs_land, afrt_misc
      
def s2num(s):
    try: return int(s)
    except ValueError:
        try: return float(s)
        except ValueError:
            return None            
    
def read_text(key, value):
    if key == 'sfc':
        ds = np.zeros((2160,1080), dtype='int32')
        f = open(value, 'r')
        lines = f.readlines()
        a = [list(map(s2num, s.split())) for s in lines]
        ds = np.asarray(a).T
        data = xr.DataArray(ds, dims=["dim_lon","dim_lat"])   
        return xr.Dataset(dict(data=data))

def read_binary(key, value):
    if key == 'sfc':
        ds = np.zeros((2160,1080), dtype='int32')
        f = np.fromfile(value, dtype='int8')
        p=0
        for i in range(1080):
            for j in range(2160):
                ds[j,i] = int.from_bytes(f[p:p+3],'big',signed=True)
                p+=3
            p+=1
        data = xr.DataArray(ds, dims=["dim_lon","dim_lat"])   
        return xr.Dataset(dict(data=data))
    elif key == 'viirs_ler':
        ds = xr.Dataset(
            data_vars=dict(        
                logi0 = xr.DataArray(np.zeros((20800)), dims=["dim_20800"]),
                z1i0 = xr.DataArray(np.zeros((20800)), dims=["dim_20800"]),
                z2i0 = xr.DataArray(np.zeros((20800)), dims=["dim_20800"]),
                ti0 = xr.DataArray(np.zeros((20800)), dims=["dim_20800"]),
                sb = xr.DataArray(np.zeros((260)), dims=["dim_260"]),
                logi0r = xr.DataArray(np.zeros((160)), dims=["dim_160"]),
                z1i0r = xr.DataArray(np.zeros((160)), dims=["dim_160"]),
                z2i0r = xr.DataArray(np.zeros((160)), dims=["dim_160"]),
                ti0r = xr.DataArray(np.zeros((160)), dims=["dim_160"]),
                sbr = xr.DataArray(np.zeros((2)), dims=["dim_2"]),
            )
        )
        f = np.fromfile(value, dtype='>f4')
        s = 20800
        r = 160
        ds['logi0'][:] = f[1:1+s]
        ds['z1i0'][:] = f[1+s:1+2*s]
        ds['z2i0'][:] = f[1+2*s:1+3*s]
        ds['ti0'][:] = f[1+3*s:1+4*s]
        ds['sb'][:] = f[1+4*s:1+4*s+260]
        ds['logi0r'][:] = f[3+4*s+260:3+4*s+260+r]
        ds['z1i0r'][:] = f[3+4*s+260+r:3+4*s+260+2*r]
        ds['z2i0r'][:] = f[3+4*s+260+2*r:3+4*s+260+3*r]
        ds['ti0r'][:] = f[3+4*s+260+3*r:3+4*s+260+4*r]
        ds['sbr'][:] = f[3+4*s+260+4*r:3+4*s+260+4*r+2]
        return ds
    elif key == 'airsco_clm':
        f = np.fromfile(value, dtype='<f4')
        ds = np.reshape(f[1:-1],(12,180,360))
        data = xr.DataArray(ds, dims=["dim_month","dim_lat", "dim_lon", ])   
        return xr.Dataset(dict(data=data))
    else:
        f = np.fromfile(value, dtype='float32')
        ds = np.reshape(f[1:-1],(180,360,12))
        data = xr.DataArray(ds, dims=["dim_lat", "dim_lon","dim_month" ])   
        return xr.Dataset(dict(data=data))
 
def read_hdf4(key, value):
    ds4 = SD(value, SDC.READ)
    fdict = ds4.datasets()
    info = ds4.info()
    attr = ds4.attributes()
    odict ={}
    for name in fdict.keys():
        data = ds4.select(name).get()
        dattr = ds4.select(name).attributes()
        dinfo = ds4.select(name).info()
        ddims=fdict[name][0]
        dxr = xr.DataArray(data, dims=ddims, attrs=dattr)
        odict.update({name: dxr})
    oxr = xr.Dataset( data_vars=odict, attrs= attr)
    return oxr
   
###################################################
#                  Main Function                  #
###################################################

def main():

    args = sys.argv
    
    description = '''Convert h5, h4, text and fortran binary LUTS to single NetCDF4 file.'''
    usage = "usage: files2nc.py -i [input par file] -o [output file]"
    version = "v1"

    parser = argparse.ArgumentParser()
    # Mandatory arguments
    parser.add_argument("-i", "--par", type=argparse.FileType('r'), 
                        help="config file", required=True)
    parser.add_argument("-o", "--ofile", type=argparse.FileType('w'), 
                        help="output file", required=False)

    # Parse the arguments from the command line
    args = parser.parse_args()

    config = cp.ConfigParser()
    config.read(args.par.name)    
    
    dict = {}
    for section in config.sections():
        if section == 'dbland_h4':
            for (key, val) in config.items(section):
                print (key + " : " + val)
                ds4 = read_hdf4(key, val)
                dict.update({key: ds4})
        if section == 'dbland_surf_coeffs':
            dblist = []
            for (key, val) in config.items(section):
                print (key + " : " + val)
                dblist.append(read_hdf4(key, val))
            dx = xr.concat(dblist, dim='season', coords='minimal')
            dict.update({'surface_coeffs': dx})
        if section == 'dbland_modis_surf':
            dblist = []
            for (key, val) in config.items(section):
                print (key + " : " + val)
                dblist.append(read_hdf4(key, val))
            dx = xr.concat(dblist, dim='season', coords='minimal')
            dict.update({'modis_surface': dx})
        if section == 'dbland_viirs_surf':
            dblist = []
            for (key, val) in config.items(section):
                print (key + " : " + val)
                dblist.append(read_hdf4(key, val))
            dx = xr.concat(dblist, dim='season', coords='minimal')
            dict.update({'viirs_surface': dx})
        if section == 'dbland_swir_vis_surf_coeffs':
            dblist = []
            for (key, val) in config.items(section):
                print (key + " : " + val)
                dblist.append(read_hdf4(key, val))
            dx = xr.concat(dblist, dim='season', coords='minimal')
            dict.update({'swir_vis_surface_coeffs': dx})
        if section == 'dbland_binary':
            for (key, val) in config.items(section):
                print (key + " : " + val)
                dsb = read_binary(key, val)
                dict.update({key: dsb})
        if section == 'dtland_aerosol':
            try: afl
            except NameError:
                afl = afrt_viirs_land()
            for (key, val) in config.items(section):
                print (key + " : " + val)
                afl.merge_file(key, val)
            dict.update({'land_aerosol': afl.ds})
        if section == 'dtocean_viirs_aerosol':
            try: afo
            except NameError:
                afo = afrt_viirs_ocean()
            for (key, val) in config.items(section):
                print (key + " : " + val)
                afo.merge_file(key, val)
            dict.update({'viirs_ocean_aerosol': afo.ds})
        if section == 'dtocean_pace_aerosol':
            try: afp
            except NameError:
                afp = afrt_pace()
            for (key, val) in config.items(section):
                print (key + " : " + val)
                afp.merge_file(key, val)
            dict.update({'pace_ocean_aerosol': afp.ds})
        if section == 'dt_misc':
            for (key, val) in config.items(section):
                aft = afrt_misc(key, val)
                print (key + " : " + val)
                dst = aft.read_file()
                dict.update({key: dst})
        if section == 'uv_h5':
            for (key, val) in config.items(section):
                print (key + " : " + val)
                dti = dtree.open_datatree(val, engine="h5netcdf", phony_dims='access')
                dict.update({key: dti[dti.groups[1]]})
        if section == 'uv_binary':
            for (key, val) in config.items(section):
                print (key + " : " + val)
                dsb = read_binary(key, val)
                dict.update({key: dsb})
        if section == 'uv_text':
            for (key, val) in config.items(section):
                print (key + " : " + val)
                dst = read_text(key, val)
                dict.update({key: dst})

    dto = dtree.DataTree.from_dict(dict) 
    dto.to_netcdf(args.ofile.name, engine="netcdf4")
    
    print("\n nc4 lut file created - " + args.ofile.name)   
    return 0

if __name__=='__main__':
    sys.exit(main())        
        
