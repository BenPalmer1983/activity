Module solveFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use calcFunctions
  Use rng
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: SolvePolynomial
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  
  Function SolvePolynomial(coefficients, lower, upper, convergenceThresholdIn) RESULT (output)
! Solves the polynomial p(x) = 0, in region close to p(x) = 0
! For 3rd order polynomials or above (quadratic equation for 2nd order)
! Newton's method
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal) :: upper, lower, output
    Real(kind=DoubleReal) :: xN,y,dydx
    Real(kind=DoubleReal) :: convergence, convergenceThreshold, convergenceTarget
    Real(kind=DoubleReal), Optional :: convergenceThresholdIn
    Integer(kind=StandardInteger) :: maxLoops
! Optional argument
    convergenceThreshold = 1.0D-7
    If(Present(convergenceThresholdIn))Then
      convergenceThreshold = convergenceThresholdIn
    End If
! Set values
    convergenceTarget = 0.0D0
    convergence = 1000
    maxLoops = 0
! set start value for x
    xN = lower + RandomLCG()*(upper-lower)
    Do while(convergence.gt.convergenceThreshold.and.maxLoops.le.10000)
      maxLoops = maxLoops + 1
      y = CalcPolynomial(coefficients, xN, 0)
      dydx = CalcPolynomial(coefficients, xN, 1)
      xN = xN - (y/dydx)
      convergence = abs(convergenceTarget - y)
    End Do
    output = xN
  End Function solvePolynomial
  
  
  

End Module solveFunctions