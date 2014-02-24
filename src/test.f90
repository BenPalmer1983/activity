Program test

! Setup Modules
Use kinds				!data kinds
Use constants
Use stringfunctions		!string functions
Use maths

Real(kind=DoubleReal) :: a, b
Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChange
Real(kind=QuadrupoleReal), Dimension( : ), Allocatable :: coefficients	

print *,"Test"

Allocate(decayDataArray(1:3,1:4))
decayDataArray(1,1) = 1
decayDataArray(1,2) = 96477038085D0
decayDataArray(1,3) = 63108D0
decayDataArray(1,4) = 1
decayDataArray(2,1) = 2
decayDataArray(2,2) = 4.9789604658574400D+7
decayDataArray(2,3) = 125000
decayDataArray(2,4) = 1
decayDataArray(3,1) = 3
decayDataArray(3,2) = 4511243
decayDataArray(3,3) = 500000
decayDataArray(3,4) = 1
!decayDataArray(4,1) = 4
!decayDataArray(4,2) = 800000
!decayDataArray(4,3) = 10000000
!decayDataArray(4,4) = 1
!decayDataArray(5,1) = 5
!decayDataArray(5,2) = 70000
!decayDataArray(5,3) = 1D9
!decayDataArray(5,4) = 1
!decayDataArray(6,1) = 5
!decayDataArray(6,2) = 200000
!decayDataArray(6,3) = -1
!decayDataArray(6,4) = 1
isotopeChange = CalcIsotopeAmount(1.0D8,decayDataArray,0D0,300D0,250000D0)

!coefficients = GaverStehfestCoeffs(14)

End