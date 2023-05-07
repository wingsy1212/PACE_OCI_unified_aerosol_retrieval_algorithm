 PRO CAL_JULIAN_DAY,YYYY,MM,DD,JDAY

SPAWN, 'cdate ' + YYYY + ' ' + MM + ' ' + DD, YYJDAY

year_jday = strsplit(YYJDAY, /regex, /extract)

jday = year_jday(1)

END
 PRO POSTPROCESS_GIANT,yyyy,mo,dd,hh,mm,day
  dayInMonth = [31,28,31,30,31,30,31,31,30,31,30,31]
  
; get ancillary data
    

    
;   yyyy = 2012
;   mo   = 04
;   dd   = 13

   IF ((yyyy MOD 4) EQ 0) THEN dayInMonth[1] = 29
; mmp is set to set so that it keeps the zero   
    mmp   = STRING(format = '(I2.2)', mm) 
    
    
   
   time = strtrim(string(hh),2)+':'+ strtrim(string(mmp),2) 
      

   timeparts = STR_SEP(time, ':')
  utc = FLOAT(timeparts[0]) + FLOAT(timeparts[1])/60.0
   
  hh  = ROUND(utc/6.)*6
 
  IF (hh EQ 24) THEN BEGIN

    hh = 0
    dd = dd+1

    IF (dd GT dayInMonth[mo-1]) THEN BEGIN

      dd = 1
      mo = mo + 1

      IF( mo EQ 13) THEN BEGIN

         mo  =1
         yyyy = yyyy+1
      ENDIF

   ENDIF

 ENDIF


  yyyy  = STRING(format = '(I4.4)', yyyy)
  hhz   = STRING(format = '(I2.2)', hh) + 'z'
  mo    = STRING(mo, format = '(I2.2)')
  yy    = STRMID(yyyy,2,2)
  dd    = STRING(format = '(i2.2)', dd)
  day   = STRING(format = '(i3.3)', day)
  
  yymmdd = yy+mo+dd 
  iyymmdd = LONG(yymmdd)
 
  ;------------------------------------------------------------------------
  ; Get GDAS data
  ;------------------------------------------------------------------------

  ;raid = '/products/GDAS1'
   gdas_dir = '/products/GDAS1/'+yyyy+'/'+mo
  ;gdas_dir = '/mnt/peate/ingest/ancillary/'+yyyy+'/'+day+'/GDAS_0ZF'
    CAL_JULIAN_DAY,YYYY,mo,DD,JDAY
  ;gdas_dir = '/mnt/peate/ingest/ancillary/'+yyyy+'/'+JDAY+'/GDAS_0ZF'
  ; gdas_dir = '/mnt/dawg/ancillary/GDAS_0ZF/'+yyyy+'/'+JDAY
   ;gdas_dir = '/f001/smattoo/Geo/ABI_AHI/ancillary/GDAS_0ZF/'+yyyy+'/'+JDAY
   gdas_file_root = gdas_dir +'/gdas1.PGrbF00.'+yymmdd+'.'+hhz 
                               
    
  gdas_files = findfile(gdas_file_root, count = gcount)
  
  print,'gdas_dir  ',gdas_file_root
 
  ;------------------------------------------------------------------------
  ; Get ozone
  ;------------------------------------------------------------------------

   
;  yearmonth = FIX(yyyy)*12 + FIX(mo) 
;  IF (yearmonth GE 24064) THEN BEGIN
;    oz_dir = '/products/TOAST'
;    oz_file_root = oz_dir + '/TOAST_'+yymmdd+'.GRB'  
; ENDIF
;  print,'oz_file_root',oz_file_root 
  OPENW, FLUN, 'gdas_file_name', /GET_LUN
  PRINTF, FLUN, gdas_file_root
  FREE_LUN, FLUN
;  OPENW, FLUN, 'oz_file_name', /GET_LUN
;  PRINTF, FLUN, oz_file_root
;  FREE_LUN, FLUN

   End
   
