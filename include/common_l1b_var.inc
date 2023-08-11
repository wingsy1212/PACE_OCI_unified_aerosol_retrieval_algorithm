! 
! These are common variables from both proxy and synthetic 
! data at the original resolution. 
! Vinay Kayetha (SSAI)
!
Integer(Kind=4) :: l1b_nXtrack, l1b_nLines
Character(Len=255) :: l1b_date_time_str
Real(Kind=4), Dimension(:,:), Allocatable :: Latitude
Real(Kind=4), Dimension(:,:), Allocatable :: Longitude
Real(Kind=4), Dimension(:,:), Allocatable :: SolarAzimuthAngle
Real(Kind=4), Dimension(:,:), Allocatable :: SolarZenithAngle
Real(Kind=4), Dimension(:,:), Allocatable :: ViewingAzimuthAngle
Real(Kind=4), Dimension(:,:), Allocatable :: ViewingZenithAngle
Real(Kind=4), Dimension(:,:), Allocatable :: RelativeAzimuthAngle
Real(Kind=4), Dimension(:,:), Allocatable :: TerrainHeight
