#/usr/bin/env python
# encoding: utf-8

import os 
import sys
import math
import time
import argparse
import numpy as np
import xarray as xr 
from PyQt5.Qt import dec

# ocean
NWAV = 7
NWAV_UV = 2
NMOMENTS = 4
NSZA = 11
NVZA = 16
NRAA = 16
NAOT = 6
NBIG = 6
NSMALL = 4
NCASES = NBIG+NSMALL
NHEIGHT = 4
NOMEGA = 4
NWS = 4

#land
NLWAV=4
NLSZA=11
NLVZA=15
NLRAA=16
NLTAU=7
NLTABLE=5
NSEASONS=4
NLATS=180
NLONS=360

def s2num(s):
    try: return int(s)
    except ValueError:
        try: return float(s)
        except ValueError:
            return None                

class afrt_ocean(object):
    '''
        RGS                Radius
        SIGMA              Standard deviation
        EXT                Extinction coefficient
        MOMENTS            Moments
        EXTNOR             Extinction coefficient normalized
        BACKSCTT           Backscattering ratio
        ASSYM              Assymetry factor
        ALBEDO             Albedo single scattering
        ALBEDO_R           Reflected albedo
        ALBEDO_T           Transmitted albedo
        PHC                Azimuth angle
        THET               View angle.
        THET0              Solar zenith angle.
        AINT               Radiance(l/fo),fo=1 large mode
        TAUA               Aerosol optical thickness 
        JPHI               Azimuth
        ref_rayall         Radiance(l/fo),fo=1 FOR RAYLEIGH TAUA=0.0
    '''
   
    def __init__(self, nwav, ncase, mode):
        
        self.nwav = nwav
        self.ncase = ncase
        self.mode = mode
        self.tstep = 610
        self.astep = 11
        self.ds = xr.Dataset(
            data_vars=dict(        
                PHI = xr.DataArray(np.zeros((NRAA)), dims=["dim_raa"]),
                THET = xr.DataArray(np.zeros((NVZA)), dims=["dim_vza"]),
                THET0 = xr.DataArray(np.zeros((NSZA)), dims=["dim_sza"]),
                WAVE = xr.DataArray(np.zeros((nwav)), dims=["dim_wav"]),
                ALBEDO_R_RAY = xr.DataArray(np.zeros((NSZA,nwav)), 
                                    dims=["dim_sza","dim_wav"]),
                ALBEDO_T_RAY = xr.DataArray(np.zeros((NSZA,nwav)), 
                                    dims=["dim_sza","dim_wav"]),
                REF_RAYALL = xr.DataArray(np.zeros((nwav,NSZA,NVZA,NRAA)), 
                                    dims=["dim_wav","dim_sza","dim_vza","dim_raa"]),
                EXT = xr.DataArray(np.zeros((nwav,ncase)), 
                                    dims=["dim_wav","dim_case"]),
                RGS = xr.DataArray(np.zeros((ncase)), 
                                    dims=["dim_case"]),
                SIGMA = xr.DataArray(np.zeros((ncase)), 
                                    dims=["dim_case"]),
                MOMENTS = xr.DataArray(np.zeros((NMOMENTS,ncase)), 
                                    dims=["dim_mom","dim_case"]),
                CCN = xr.DataArray(np.zeros((ncase)), 
                                    dims=["dim_case"]),
                BACKSCTT = xr.DataArray(np.zeros((nwav,ncase)), 
                                    dims=["dim_wav","dim_case"]),
                ASSYM = xr.DataArray(np.zeros((nwav,ncase)), 
                                    dims=["dim_wav","dim_case"]),
                ALBEDO = xr.DataArray(np.zeros((nwav,ncase)), 
                                    dims=["dim_wav","dim_case"]),
                ALBEDO_R = xr.DataArray(np.zeros((ncase,nwav,NAOT,NSZA)), 
                                    dims=["dim_case","dim_wav","dim_aot","dim_sza"]),
                ALBEDO_T = xr.DataArray(np.zeros((ncase,nwav,NAOT,NSZA)), 
                                    dims=["dim_case","dim_wav","dim_aot","dim_sza"]),
                TAUA = xr.DataArray(np.zeros((NAOT,nwav,ncase)), 
                                    dims=["dim_aot","dim_wav","dim_case"]),
                AINT = xr.DataArray(np.zeros((ncase,nwav,NAOT,NSZA,NVZA,NRAA)), 
                                    dims=["dim_case","dim_wav","dim_aot","dim_sza","dim_vza","dim_raa"]),
                EFFRAD = xr.DataArray(np.zeros((ncase)), 
                                    dims=["dim_case"])
            )
        )

    def read_file(self,filepath):    
    
        f = open(filepath, 'r')
        lines = f.readlines()
        a = [list(map(s2num, s.split())) for s in lines]
        t = 0
        if self.mode == "small":
            for iwav in range(self.nwav):
                self.ds['WAVE'][iwav] = a[t+1][1]
                self.ds['THET0'][:] = a[t+12][1:] + a[t+13][:]
                self.ds['ALBEDO_R_RAY'][:,iwav] = a[t+22][1:] + a[t+23][:]
                self.ds['ALBEDO_T_RAY'][:,iwav] = a[t+24][1:] + a[t+25][:]
                self.ds['THET'][:] = a[t+28][2:] + a[t+29][:] + a[t+30][:]
                for itho in range(NSZA):
                    for iphi in range(NRAA):
                        p = 3*iphi;
                        self.ds['PHI'][iphi] = a[t+31+p][0]
                        self.ds['REF_RAYALL'][iwav][itho][:][iphi] = \
                            a[t+31+p][1:] + a[t+32+p][:] + a[t+33+p][:]
                t += (self.tstep - self.astep)
        for icase in range(self.ncase):
            for iwav in range(self.nwav):
                for itau in range(1,NAOT):
                    self.ds['WAVE'][iwav] = a[t+1][1]
                    self.ds['RGS'][icase] = a[t+5][2]
                    self.ds['SIGMA'][icase] = a[t+5][4]
                    self.ds['EFFRAD'][icase] = a[t+7][1]
                    self.ds['MOMENTS'][0][icase] = a[t+8][3]
                    self.ds['MOMENTS'][1][icase] = a[t+8][7]
                    self.ds['MOMENTS'][2][icase] = a[t+9][3]
                    self.ds['MOMENTS'][3][icase] = a[t+9][7]
                    self.ds['ALBEDO'][iwav][icase] = a[t+10][3]
                    self.ds['ASSYM'][iwav][icase] = a[t+10][7]
                    self.ds['BACKSCTT'][iwav][icase] = a[t+11][5]
                    self.ds['EXT'][iwav][icase] = a[t+12][5]
                    self.ds['TAUA'][itau][iwav][icase] = a[t+18][3]
                    self.ds['THET0'][:] = a[t+23][1:] + a[t+24][:]
                    self.ds['ALBEDO_R'][icase][iwav][itau][:] = a[t+33][1:] + a[t+34][:]
                    self.ds['ALBEDO_T'][icase][iwav][itau][:] = a[t+35][1:] + a[t+36][:]
                    self.ds['THET'][:] = a[t+39][2:] + a[t+40][:] + a[t+41][:]
                    for itho in range(NSZA):
                        for iphi in range(NRAA):
                            p = 3*iphi;
                            self.ds['PHI'][iphi] = a[t+42+p][0]
                            self.ds['AINT'][icase][iwav][itau][itho][:][iphi] = \
                            a[t+42+p][1:] + a[t+43+p][:] + a[t+44+p][:]
                    t += self.tstep
 
        return self.ds

    
class afrt_land(object):
 
    '''
      INT          Reflectance
      FD           Flux down
      T            Transmission factor
      OPTH         Optical depth
      SBAR         Sbar
      MASSCOEF     Mass Concentration coeficient
      EXTNORM      Normalized extinction coeficient
    '''
   
    def __init__(self):
        
        self.tstep = 1409
        self.bstep = 201
        self.sstep = 18
        self.ds = xr.Dataset(
            data_vars=dict(        
                PHI = xr.DataArray(np.zeros((NLRAA)), dims=["dim_raa"]),
                THET = xr.DataArray(np.zeros((NLVZA)), dims=["dim_vza"]),
                THET0 = xr.DataArray(np.zeros((NLSZA)), dims=["dim_sza"]),
                MU0 = xr.DataArray(np.zeros((NLSZA)), dims=["dim_sza"]),
                WAVE = xr.DataArray(np.zeros((1)), dims=["dim_wav"]),
                ROD = xr.DataArray(np.zeros((1)), dims=["dim_wav"]),
                GOD = xr.DataArray(np.zeros((1)), dims=["dim_wav"]),
                SSA = xr.DataArray(np.zeros((NLTABLE,NLTAU)), 
                            dims=["dim_table","dim_tau"]),
                QEXT = xr.DataArray(np.zeros((NLTABLE,NLTAU)), 
                            dims=["dim_table","dim_tau"]),
                BEXT = xr.DataArray(np.zeros((NLTABLE,NLTAU)), 
                            dims=["dim_table","dim_tau"]),
                VEXT = xr.DataArray(np.zeros((NLTABLE,NLTAU)), 
                            dims=["dim_table","dim_tau"]),
                MEXT = xr.DataArray(np.zeros((NLTABLE,NLTAU)), 
                            dims=["dim_table","dim_tau"]),
                MASSCOEF = xr.DataArray(np.zeros((NLTABLE,NLTAU)), 
                            dims=["dim_table","dim_tau"]),
                OPTH = xr.DataArray(np.zeros((NLTABLE,NLTAU)), 
                            dims=["dim_table","dim_tau"]),
                EXTNORM = xr.DataArray(np.zeros((NLTABLE,NLTAU)), 
                            dims=["dim_table","dim_tau"]),
                SBAR = xr.DataArray(np.zeros((NLTABLE,NLTAU,NLSZA)), 
                            dims=["dim_table","dim_tau","dim_sza"]),
                INT = xr.DataArray(np.zeros((NLTABLE,NLTAU,NLSZA,NLVZA,NLRAA)), 
                            dims=["dim_table","dim_tau","dim_sza","dim_vza","dim_raa"]),
                FD = xr.DataArray(np.zeros((NLTABLE,NLTAU,NLSZA)), 
                            dims=["dim_table","dim_tau","dim_sza"]),
                T = xr.DataArray(np.zeros((NLTABLE,NLTAU,NLSZA,NLVZA)), 
                            dims=["dim_table","dim_tau","dim_sza","dim_vza"]),
            )
        )

    def read_file(self,filepath):    
    
        f = open(filepath, 'r')
        lines = f.readlines()
        a = [list(map(s2num, s.split())) for s in lines]

        t = 0
        for itab in range(NLTABLE):
            self.ds['THET'][:] = a[t][1:] 
            self.ds['PHI'][:] = a[t+1][1:]
            b=0 
            for itau in range(NLTAU):
                self.ds['SSA'][itab][itau] = a[t+b+3][1]
                self.ds['QEXT'][itab][itau] = a[t+b+3][3]
                self.ds['BEXT'][itab][itau] = a[t+b+3][5]
                self.ds['VEXT'][itab][itau] = a[t+b+3][7]
                self.ds['MEXT'][itab][itau] = a[t+b+3][9]
                self.ds['MASSCOEF'][itab][itau] = a[t+b+3][11]
                self.ds['WAVE'] = a[t+b+4][1]
                self.ds['OPTH'][itab][itau] = a[t+b+4][3]
                self.ds['ROD'] = a[t+b+4][5]
                self.ds['GOD'] = a[t+b+4][7]
                s=0
                for ith0 in range(NLSZA):
                    self.ds['THET0'][ith0] = a[t+b+s+5][1]
                    self.ds['MU0'][ith0] = a[t+b+s+5][3]
                    self.ds['SBAR'][itab][itau][ith0] = a[t+b+s+5][5]
                    self.ds['FD'][itab][itau][ith0] = a[t+b+s+5][7]
                    self.ds['T'][itab][itau][ith0][:] = a[t+b+s+6][1:]
                    p=0 
                    for ithe in range(NLVZA):
                        try:
                            self.ds['INT'][itab][itau][ith0][ithe][:] = a[t+b+s+p+8]
                        except:
                            print( "check for negative values, minus sign occupying the blank spaces at line: ",str(t+b+s+p+9))
                        p += 1
                    s += self.sstep
                b += self.bstep
            t += self.tstep
 
        return self.ds
    

class afrt_pace(object):
   
    def __init__(self):
        
        self.ds = xr.Dataset(
            data_vars=dict(        
                PHI = xr.DataArray(np.zeros((NRAA)), dims=["dim_raa"]),
                THET = xr.DataArray(np.zeros((NVZA)), dims=["dim_vza"]),
                THET0 = xr.DataArray(np.zeros((NSZA)), dims=["dim_sza"]),
                WAVE = xr.DataArray(np.zeros((NWAV_UV)), dims=["dim_wav"]),
                TAUA = xr.DataArray(np.zeros((NAOT,NWAV_UV,NCASES)), 
                            dims=["dim_aot","dim_wav","dim_case"]),
                ALBEDO = xr.DataArray(np.zeros((NWAV_UV,NCASES,NHEIGHT,NOMEGA)), 
                            dims=["dim_wav","dim_case","dim_height","dim_omega"]),
                EXT = xr.DataArray(np.zeros((NWAV_UV,NCASES,NHEIGHT,NOMEGA)), 
                            dims=["dim_wav","dim_case","dim_height","dim_omega"]),
                REF_RAYALL = xr.DataArray(np.zeros((NHEIGHT,NWAV_UV,NSZA,NVZA,NRAA,NWS,NOMEGA)), 
                            dims=["dim_height","dim_wav","dim_sza","dim_vza","dim_raa","dim_wspd","dim_omega"]),
                AINT = xr.DataArray(np.zeros((NCASES,NWAV_UV,NAOT,NSZA,NVZA,NRAA,NWS,NHEIGHT,NOMEGA)), 
                            dims=["dim_case","dim_wav","dim_aot","dim_sza","dim_vza","dim_raa","dim_wspd","dim_height","dim_omega"])
            )
        )

    def merge_file(self, key, value):
         
        ht = int(key[9])-1
        om = int(key[10])-1
        ws = (int(key[11])-1)%4
        if int(key[11]) < 5:
            mode = "big"
            af = afrt_ocean(NWAV_UV,NBIG,mode)
        else:
            mode = "small"
            af = afrt_ocean(NWAV_UV,NSMALL,mode)
        
        dr = af.read_file(value)
        self.ds['PHI'] = dr.PHI
        self.ds['THET'] = dr.THET
        self.ds['THET0'] = dr.THET0
        self.ds['WAVE'] = dr.WAVE
        if mode == "small":
            self.ds['EXT'][:,:NSMALL,ht,om] = dr['EXT'] 
            self.ds['ALBEDO'][:,:NSMALL,ht,om] = dr['ALBEDO']
            self.ds['TAUA'][:,:,:NSMALL] = dr['TAUA']
            self.ds['AINT'][:NSMALL,:,:,:,:,:,ws,ht,om] = dr['AINT']
            self.ds['REF_RAYALL'][ht,:,:,:,:,ws,om] = dr['REF_RAYALL']
        else:
            self.ds['EXT'][:,NSMALL:,ht,om] = dr['EXT'] 
            self.ds['ALBEDO'][:,NSMALL:,ht,om] = dr['ALBEDO']
            self.ds['TAUA'][:,:,NSMALL:] = dr['TAUA']
            self.ds['AINT'][NSMALL:,:,:,:,:,:,ws,ht,om] = dr['AINT']
            self.ds['REF_RAYALL'][ht,:,:,:,:,ws,om] = dr['REF_RAYALL']

            
class afrt_viirs_ocean(object):
   
    def __init__(self):
        
        self.ds = xr.Dataset(
            data_vars=dict(        
                PHI = xr.DataArray(np.zeros((NRAA)), dims=["dim_raa"]),
                THET = xr.DataArray(np.zeros((NVZA)), dims=["dim_vza"]),
                THET0 = xr.DataArray(np.zeros((NSZA)), dims=["dim_sza"]),
                WAVE = xr.DataArray(np.zeros((NWAV)), dims=["dim_wav"]),
                ALBEDO_R_RAY = xr.DataArray(np.zeros((NSZA,NWAV)), 
                            dims=["dim_sza","dim_wav"]),
                ALBEDO_T_RAY = xr.DataArray(np.zeros((NSZA,NWAV)), 
                            dims=["dim_sza","dim_wav"]),
                REF_RAYALL = xr.DataArray(np.zeros((NWAV,NSZA,NVZA,NRAA,NWS)), 
                            dims=["dim_wav","dim_sza","dim_vza","dim_raa","dim_wspd",]),
                EXT = xr.DataArray(np.zeros((NWAV,NCASES,NWS)), 
                            dims=["dim_wav","dim_case","dim_wspd"]),
                RGS = xr.DataArray(np.zeros((NCASES)), 
                            dims=["dim_case"]),
                SIGMA = xr.DataArray(np.zeros((NCASES)), 
                            dims=["dim_case"]),
                MOMENTS = xr.DataArray(np.zeros((NMOMENTS,NCASES,NWS)), 
                            dims=["dim_mom","dim_case","dim_wspd"]),
                CCN = xr.DataArray(np.zeros((NCASES,NWS)), 
                            dims=["dim_case","dim_wspd"]),
                BACKSCTT = xr.DataArray(np.zeros((NWAV,NCASES,NWS)), 
                            dims=["dim_wav","dim_case","dim_wspd"]),
                ASSYM = xr.DataArray(np.zeros((NWAV,NCASES,NWS)), 
                            dims=["dim_wav","dim_case","dim_wspd"]),
                ALBEDO = xr.DataArray(np.zeros((NWAV,NCASES,NWS)), 
                            dims=["dim_wav","dim_case","dim_wspd"]),
                ALBEDO_R = xr.DataArray(np.zeros((NCASES,NWAV,NAOT,NSZA,NWS)), 
                            dims=["dim_case","dim_wav","dim_aot","dim_sza","dim_wspd"]),
                ALBEDO_T = xr.DataArray(np.zeros((NCASES,NWAV,NAOT,NSZA,NWS)), 
                            dims=["dim_case","dim_wav","dim_aot","dim_sza","dim_wspd"]),
                TAUA = xr.DataArray(np.zeros((NAOT,NWAV,NCASES)), 
                            dims=["dim_aot","dim_wav","dim_case"]),
                AINT = xr.DataArray(np.zeros((NCASES,NWAV,NAOT,NSZA,NVZA,NRAA,NWS)), 
                            dims=["dim_case","dim_wav","dim_aot","dim_sza","dim_vza","dim_raa","dim_wspd"]),
                EFFRAD = xr.DataArray(np.zeros((NCASES)), 
                            dims=["dim_case"])
            )
        )

    def merge_file(self, key, value):
         
        ws = (int(key[12])-1)%4
        if int(key[12]) < 5:
            mode = "big"
            af = afrt_ocean(NWAV,NBIG,mode)
        else:
            mode = "small"
            af = afrt_ocean(NWAV,NSMALL,mode)
        
        dr = af.read_file(value)
        
        self.ds['PHI'] = dr.PHI
        self.ds['THET'] = dr.THET
        self.ds['THET0'] = dr.THET0
        self.ds['WAVE'] = dr.WAVE
        if mode == "small":
            self.ds['EXT'][:,:NSMALL,ws] = dr['EXT'] 
            self.ds['RGS'][:NSMALL] = dr['RGS']
            self.ds['SIGMA'][:NSMALL] = dr['SIGMA']
            self.ds['MOMENTS'][:,:NSMALL,ws] = dr['MOMENTS']
            self.ds['CCN'][:NSMALL,ws] = dr['CCN']
            self.ds['BACKSCTT'][:,:NSMALL,ws] = dr['BACKSCTT']
            self.ds['ASSYM'][:,:NSMALL,ws] = dr['ASSYM']
            self.ds['ALBEDO'][:,:NSMALL,ws] = dr['ALBEDO']
            self.ds['ALBEDO_R'][:NSMALL,:,:,:,ws] = dr['ALBEDO_R']
            self.ds['ALBEDO_T'][:NSMALL,:,:,:,ws] = dr['ALBEDO_T']
            self.ds['TAUA'][:,:,:NSMALL] = dr['TAUA']
            self.ds['AINT'][:NSMALL,:,:,:,:,:,ws] = dr['AINT']
            self.ds['EFFRAD'][:NSMALL] = dr['EFFRAD']
        else:
            self.ds['EXT'][:,NSMALL:,ws] = dr['EXT'] 
            self.ds['RGS'][NSMALL:] = dr['RGS']
            self.ds['SIGMA'][NSMALL:] = dr['SIGMA']
            self.ds['MOMENTS'][:,NSMALL:,ws] = dr['MOMENTS']
            self.ds['CCN'][NSMALL:,ws] = dr['CCN']
            self.ds['BACKSCTT'][:,NSMALL:,ws] = dr['BACKSCTT']
            self.ds['ASSYM'][:,NSMALL:,ws] = dr['ASSYM']
            self.ds['ALBEDO'][:,NSMALL:,ws] = dr['ALBEDO']
            self.ds['ALBEDO_R'][NSMALL:,:,:,:,ws] = dr['ALBEDO_R']
            self.ds['ALBEDO_T'][NSMALL:,:,:,:,ws] = dr['ALBEDO_T']
            self.ds['TAUA'][:,:,NSMALL:] = dr['TAUA']
            self.ds['AINT'][NSMALL:,:,:,:,:,:,ws] = dr['AINT']
            self.ds['EFFRAD'][NSMALL:] = dr['EFFRAD']

           
class afrt_viirs_land(object):
   
    def __init__(self):
        
        self.ds = xr.Dataset(
            data_vars=dict(        
                PHI = xr.DataArray(np.zeros((NLRAA)), dims=["dim_raa"]),
                THET = xr.DataArray(np.zeros((NLVZA)), dims=["dim_vza"]),
                THET0 = xr.DataArray(np.zeros((NLSZA)), dims=["dim_sza"]),
                MU0 = xr.DataArray(np.zeros((NLSZA)), dims=["dim_sza"]),
                WAVE = xr.DataArray(np.zeros((NLWAV)), dims=["dim_wav"]),
                ROD = xr.DataArray(np.zeros((NLWAV)), dims=["dim_wav"]),
                GOD = xr.DataArray(np.zeros((NLWAV)), dims=["dim_wav"]),
                OPTH = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU)), 
                        dims=["dim_table","dim_wav","dim_tau"]),
                MASSCOEF = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU)), 
                        dims=["dim_table","dim_wav","dim_tau"]),
                EXTNORM = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU)), 
                        dims=["dim_table","dim_wav","dim_tau"]),
                SSA = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU)), 
                        dims=["dim_table","dim_wav","dim_tau"]),
                QEXT = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU)), 
                        dims=["dim_table","dim_wav","dim_tau"]),
                BEXT = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU)), 
                        dims=["dim_table","dim_wav","dim_tau"]),
                VEXT = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU)), 
                        dims=["dim_table","dim_wav","dim_tau"]),
                MEXT = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU)), 
                        dims=["dim_table","dim_wav","dim_tau"]),
                SBAR = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU,NLSZA)), 
                        dims=["dim_table","dim_wav","dim_tau","dim_sza"]),
                INT = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU,NLSZA,NLVZA,NLRAA)), 
                        dims=["dim_table","dim_wav","dim_tau","dim_sza","dim_vza","dim_raa"]),
                FD = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU,NLSZA)), 
                        dims=["dim_table","dim_wav","dim_tau","dim_sza"]),
                T = xr.DataArray(np.zeros((NLTABLE,NLWAV,NLTAU,NLSZA,NLVZA)), 
                        dims=["dim_table","dim_wav","dim_tau","dim_sza","dim_vza"]),
            )
        )

    def merge_file(self, key, value):
         
        wav = int(key[11])-1
        af = afrt_land()        
        dr = af.read_file(value)
        
        self.ds['PHI'] = dr['PHI']
        self.ds['THET'] = dr['THET']
        self.ds['THET0'] = dr['THET0']
        self.ds['MU0'] = dr['MU0']
        self.ds['WAVE'][wav] = dr['WAVE']
        self.ds['ROD'][wav] = dr['ROD']
        self.ds['GOD'][wav] = dr['GOD']
        self.ds['OPTH'][:,wav,:] = dr['OPTH'] 
        self.ds['MASSCOEF'][:,wav,:] = dr['MASSCOEF']
        self.ds['EXTNORM'][:,wav,:] = dr['EXTNORM']
        self.ds['SSA'][:,wav,:] = dr['SSA']
        self.ds['QEXT'][:,wav,:] = dr['QEXT']
        self.ds['BEXT'][:,wav,:] = dr['BEXT']
        self.ds['VEXT'][:,wav,:] = dr['VEXT']
        self.ds['MEXT'][:,wav,:] = dr['MEXT']
        self.ds['SBAR'][:,wav,:,:] = dr['SBAR']
        self.ds['INT'][:,wav,:,:,:,:] = dr['INT']
        self.ds['FD'][:,wav,:,:] = dr['FD']
        self.ds['T'][:,wav,:,:,:] = dr['T']


class afrt_misc(object):
   
    def __init__(self, name, filepath):
        self.name = name
        self.file = filepath
        
    def read_file(self):    
    
        f = open(self.file, 'r')
        lines = f.readlines()
        if self.name == "gas_coeffs":
            NBANDS=10
            a = [list(map(s2num, s.split())) for s in lines]
            dgc = xr.Dataset(
                data_vars=dict(        
                    MBAND = xr.DataArray(np.zeros((NBANDS),dtype=np.int32), dims=["dim_bands"]),
                    VBAND = xr.DataArray(np.zeros((NBANDS),dtype=np.int32), dims=["dim_bands"]),
                    WAVE = xr.DataArray(np.zeros((NBANDS)), dims=["dim_bands"]),
                    MOL = xr.DataArray(np.zeros((NBANDS)), dims=["dim_bands"]),
                    OPT_O3_CLIM = xr.DataArray(np.zeros((NBANDS)), dims=["dim_bands"]),
                    OPT_H2O_CLIM = xr.DataArray(np.zeros((NBANDS)), dims=["dim_bands"]),
                    OPT_CO2_CLIM = xr.DataArray(np.zeros((NBANDS)), dims=["dim_bands"]),
                    O3_COEF = xr.DataArray(np.zeros((2,NBANDS)), dims=["dim_coef2","dim_bands"]),
                    H2O_COEF = xr.DataArray(np.zeros((3,NBANDS)), dims=["dim_coef3","dim_bands"]),
                )
            )
            for ib in range(NBANDS):
                dgc['MBAND'][ib] = a[ib+1][0]
                dgc['VBAND'][ib] = a[ib+1][1]
                dgc['WAVE'][ib] = a[ib+1][2]
                dgc['MOL'][ib] = a[ib+1][3]
                dgc['OPT_O3_CLIM'][ib] = a[ib+1][4]
                dgc['OPT_H2O_CLIM'][ib] = a[ib+1][5]
                dgc['OPT_CO2_CLIM'][ib] = a[ib+1][6]
                dgc['O3_COEF'][0][ib] = a[ib+1][7]
                dgc['O3_COEF'][1][ib] = a[ib+1][8]
                dgc['H2O_COEF'][0][ib] = a[ib+1][9]
                dgc['H2O_COEF'][1][ib] = a[ib+1][10]
                dgc['H2O_COEF'][2][ib] = a[ib+1][11]
                
            return dgc
        
        if self.name == "ext_coeffs":
            a = [list(map(s2num, s.split())) for s in lines]
            dec = xr.Dataset(data_vars=dict(
                EXT_COEFFS = xr.DataArray(np.zeros((NCASES,NWS)), dims=["dim_case","dim_wspd"])))
            for iw in range(NWS):
                dec['EXT_COEFFS'][0][iw] = a[iw+1][0]
                dec['EXT_COEFFS'][1][iw] = a[iw+1][1]
                dec['EXT_COEFFS'][2][iw] = a[iw+1][2]
                dec['EXT_COEFFS'][3][iw] = a[iw+1][3]
                dec['EXT_COEFFS'][4][iw] = a[iw+5][0]
                dec['EXT_COEFFS'][5][iw] = a[iw+5][1]
                dec['EXT_COEFFS'][6][iw] = a[iw+5][2]
                dec['EXT_COEFFS'][7][iw] = a[iw+5][3]
                dec['EXT_COEFFS'][8][iw] = a[iw+5][4]
                dec['EXT_COEFFS'][9][iw] = a[iw+5][5]
            
            return dec
           
        if self.name == "aerosol_land_map":
            a = [list(map(s2num, s.split(','))) for s in lines]
            dlm = xr.Dataset(data_vars=dict(
                LAND_MAP = xr.DataArray(np.zeros((NSEASONS,NLATS,NLONS),dtype=np.int32), 
                        dims=["dim_seasons","dim_lat","dim_lon"])))
            
            for i in range(NSEASONS):
                for j in range(NLATS):
                    dlm["LAND_MAP"][i][j][:] = a[1+i*181+j]

            return dlm      