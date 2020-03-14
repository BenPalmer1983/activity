Module constants
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
! Force declaration of all variables
  Implicit None
! Useful parameters for source code
  Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000




! Legacy - no naming standard
! Physical and Mathematical constants
  Real(kind=DoubleReal), Parameter :: avogadrosConstant = 6.02214129D23
  Real(kind=DoubleReal), Parameter :: elementaryCharge = 1.60217656D-19
  Real(kind=DoubleReal), Parameter :: pi = 3.1415926535898D0
  Real(kind=DoubleReal), Parameter :: lnTwo = 0.693147180559945D0
  Real(kind=QuadrupoleReal), Parameter :: lnTwoQ = 0.6931471805599453094172321214581765680755D0
  Real(kind=DoubleReal), Parameter :: eulersNumber = 2.7182818284590452353602874713527D0
  Real(kind=DoubleReal), Parameter :: sqrtTwo = 1.4142135623731D0
  Real(kind=DoubleReal), Parameter :: sqrtTwoPi = 2.506628274631D0
  Real(kind=DoubleReal), Parameter :: gravitationalConstant = 6.67408D-11                       ! http://physics.nist.gov/cgi-bin/cuu/Value?bg
! Privacy of functions/subroutines/variables
  Private
! Useful parameters for source code
  Public :: maxFileRows
! Physical and Mathematical constants
  Public :: avogadrosConstant
  Public :: elementaryCharge
  Public :: pi
  Public :: lnTwo, lnTwoQ
  Public :: eulersNumber
  Public :: sqrtTwo
  Public :: sqrtTwoPi
  Public :: gravitationalConstant
 End Module constants