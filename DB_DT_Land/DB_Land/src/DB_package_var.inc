  integer Year,day,month
  INTEGER IX1KM_B,IY1KM_B,nwave
  PARAMETER(IX1KM_B=1271,IY1KM_B=1722,nwave=14) 
	REAL Latitude_in(IX1KM_B,IY1KM_B),Longitude_in(IX1KM_B,IY1KM_B)
	REAL TerrainHeight(IX1KM_B,IY1KM_B),SolarZenithAngle(IX1KM_B,IY1KM_B)
	REAL ViewingZenithAngle(IX1KM_B,IY1KM_B),SolarAzimuthAngle(IX1KM_B,IY1KM_B)
	REAL ViewingAzimuthAngle(IX1KM_B,IY1KM_B)
	REAL UVtoSWIR_Radiances(nwave,IX1KM_B,IY1KM_B),UVtoSWIR_IRRadiances(nwave)
	REAL UVAI (IX1KM_B,IY1KM_B)
	Integer Land_Sea_Flag (IX1KM_B,IY1KM_B)
	