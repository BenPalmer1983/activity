Module calcFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: Gaussian
  Public :: MaxwellBoltzman
  Public :: CalcPolynomial
  Public :: CalcPolynomialExp
  Public :: ExpCalc
  Public :: BirchMurnCalc
  Public :: Zbl
  Public :: ZblFull
  Public :: EmbeddingA
  Public :: EmbeddingB
  Public :: EmbeddingC
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  
  Function Gaussian(x, sigma, mu) RESULT (y)
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: x,y,mu,sigma
! Calculation
    y = (1.0D0/(sigma*sqrtTwoPi))*exp(-1*((x-mu)**2/(2*sigma**2)))
  End Function Gaussian

  Function MaxwellBoltzman(x, a) RESULT (y)
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: x,y,a
! Calculation
    y = sqrtTwoPi*((x**2*exp(-1.0D0*((x*x)/(2.0D0*a*a))))/(a**3))
  End Function MaxwellBoltzman

  Function CalcPolynomial(polyCoefficientsIn, x, derivIn) RESULT (y)
! Calculates p(x) by default, p'(x) for derivativeIn = 1 etc
! p(x) = polyCoefficientsIn(1) x^0 +  polyCoefficientsIn(2) x^1  + ... + polyCoefficientsIn(n) x^(n-1)
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: polyCoefficientsIn
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger), Optional :: derivIn
! Out
    Real(kind=DoubleReal) :: y
! Private Variables
    Real(kind=DoubleReal), Dimension(1:size(polyCoefficientsIn,1)) :: polyCoefficients
    Integer(kind=StandardInteger) :: deriv
    Integer(kind=StandardInteger) :: i, j
! Optional
    deriv = 0
    If(Present(derivIn))Then
      deriv = derivIn
    End If
! Transfer coeffs to internal/private array
    polyCoefficients = polyCoefficientsIn
! Modify as necessary
    Do i=1,deriv
      Do j=1,size(polyCoefficients,1)
        polyCoefficients(j) = (j-i)*polyCoefficients(j)
      End Do
    End Do
! calculate output
    y = 0.0D0
    Do j=1,size(polyCoefficients,1)
      y = y + polyCoefficients(j)*x**(j-1-deriv)
    End Do
  End Function CalcPolynomial
  
  Function CalcPolynomialExp(polyCoefficientsIn, x, derivIn) RESULT (y)
! Calculates f(x) =exp(a+bx+...)
! just for f(x), f'(x) and f''(x)
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: polyCoefficientsIn
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger), Optional :: derivIn
! Out
    Real(kind=DoubleReal) :: y, yP, yP_coeff, yPP, yPP_coeff
! Private Variables
    Real(kind=DoubleReal), Dimension(1:size(polyCoefficientsIn,1)) :: polyCoefficients
    Integer(kind=StandardInteger) :: deriv
    Integer(kind=StandardInteger) :: i
! Optional
    deriv = 0
    If(Present(derivIn))Then
      deriv = derivIn
    End If
! Transfer coeffs to internal/private array
    polyCoefficients = polyCoefficientsIn
! Init vars
    y = 0.0D0
    yP = 0.0D0
    yP_coeff = 0.0D0
    yPP = 0.0D0
    yPP_coeff = 0.0D0
! f(x) = exp(p(x))
    If(deriv.ge.0)Then
      Do i=1,size(polyCoefficients,1)
        y = y + polyCoefficients(i)*(x**(i-1))
      End Do
      y = exp(y)
    End If
! f'(x) = p'(x) exp(p(x))
    If(deriv.ge.1)Then
! update coeffs
      Do i=1,size(polyCoefficients,1)
        polyCoefficients(i) = (i-1)*polyCoefficients(i)
      End Do
      Do i=1,size(polyCoefficients,1)
        yP_coeff = yP_coeff + polyCoefficients(i)*(x**(i-2))
      End Do
      yP = yP_coeff*y
    End If
! f''(x) = p'(x) exp(p(x))
    If(deriv.eq.2)Then
! update coeffs
      Do i=1,size(polyCoefficients,1)
        polyCoefficients(i) = (i-2)*polyCoefficients(i)
      End Do
      Do i=1,size(polyCoefficients,1)
        yPP_coeff = yPP_coeff + polyCoefficients(i)*(x**(i-3))
      End Do
      yPP = yP_coeff*yP+yPP_coeff*y
    End If
! Update output
    If(deriv.eq.1)Then
      y = yP
    End If
    If(deriv.eq.2)Then
      y = yPP
    End If
  End Function CalcPolynomialExp
  
  Function ExpCalc(x,coefficients) RESULT (y)  
! Single or multiple term exponential functions
    Implicit None  !Force declaration of all variables
! Declare variables  
! In:      Declare variables
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension(:) :: coefficients
! Out:     Declare variables
    Real(kind=DoubleReal) :: y
! Private: Declare variables
    Integer(kind=StandardInteger) :: i, j, terms
! Calculate
    terms = ceiling(size(coefficients,1)/2.0D0)     
    y = 0.0D0
    Do i=1,terms
      j = 2*(i-1)+1
      y = y + coefficients(j)*exp(coefficients(j+1)*x)
    End Do
  End Function ExpCalc

  Function BirchMurnCalc(volume,coefficients) RESULT (energy)
! Calculate energy from volume using Murnaghan EoS
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: volume, energy, eta
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients    ! E0, V0, B0, B'0
    Real(kind=DoubleReal) :: E, V, B, BP
! Shorten var names
    E = coefficients(1)
    V = coefficients(2)
    B = coefficients(3)
    BP = coefficients(4)
! Calculate energy
    eta = ((1.0D0*volume)/(1.0D0*V))**(1.0D0/3.0D0)
    energy = E+(9.0D0/16.0D0)*(B*V)*&
    ((eta**2-1.0D0)**2)*(6.0D0+BP*(eta**2-1.0D0)-4.0D0*eta**2)
! Rearranged:
! energy = E+(9.0D0/16.0D0)*(B)*(&
!  V**(-1.0D0)*volume**2.0D0*(BP-4.0D0)&
!  +V**(-1.0D0/3.0D0)*volume**(4.0D0/3.0D0)*(14.0D0-3.0D0*BP)&
!  +V**(1.0D0/3.0D0)*volume**(2.0D0/3.0D0)*(3.0D0*BP-16.0D0)&
!  +V*(6.0D0-BP))
  End Function BirchMurnCalc
  
  
  Function Zbl (x, qA, qB) RESULT (y)
! ZBL potential, separation x, charges qA and qB
    Implicit None  ! Force declaration of all variables
! declare variables
    Integer(kind=StandardInteger) :: qA, qB
    Real(kind=DoubleReal) :: xVal, x, y, xa, xs, exa
! Force none infinite result for 0
    If(x.eq.0.0D0)Then
      xVal = 0.00001D0
    Else
      xVal = x
    End If
! Calculate y
    xs = 0.4683766 * (qA**(2.0D0/3.0D0)+qB**(2.0D0/3.0D0))**0.5
    xa = 1.0D0*xVal/xs
    exa = 0.1818D0*exp(-3.2D0*xa)+0.5099D0*exp(-0.9423D0*xa)+&
    0.2802D0*exp(-0.4029*xa)+0.02817*exp(-0.2016D0*xa)
    y = ((1.0D0*qA*qB)/xVal)*exa
  End Function Zbl

  Function ZblFull (x, qA, qB) RESULT (yArray)
! y(x), y'(x), y''(x)
    Implicit None ! Force declaration of all variables
! declare variables
    Integer(kind=StandardInteger) :: qA, qB
    Real(kind=DoubleReal) :: xVal, x, y, dy, ddy, xs
    Real(kind=DoubleReal) :: termFa, termFb, termFc, termGa, termGb, termGc
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Force none infinite result for 0
    If(x.eq.0.0D0)Then
      xVal = 0.00001D0
    Else
      xVal = x
    End If
    xs = 0.4683766 * (qA**(2.0D0/3.0D0)+qB**(2.0D0/3.0D0))**0.5
! Calculate y
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                        !g(x)
    0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    0.02817*exp((-0.2016D0/xs)*xVal)
    y = termFa * termGa
    yArray(1) = y
! Calculate dy
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termFb = (1.0D0*qA*qB)*(xVal)**(-2.0D0)*(-1.0D0)                 !f'(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                        !g(x)
    0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    0.02817*exp((-0.2016D0/xs)*xVal)
    termGb = (-3.2D0/xs)*0.1818D0*exp((-3.2D0/xs)*xVal)+&            !g'(x)
    (-0.9423D0/xs)*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    (-0.4029D0/xs)*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    (-0.2016D0/xs)*0.02817*exp((-0.2016D0/xs)*xVal)
    dy = termFa*termGb+termFb*termGa
    yArray(2) = dy
! Calculate ddy
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termFb = (1.0D0*qA*qB)*(xVal)**(-2.0D0)*(-1.0D0)                        !f'(x)
    termFc = (1.0D0*qA*qB)*(xVal)**(-3.0D0)*(-1.0D0)*(-2.0D0)               !f''(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                             !g(x)
    0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    0.02817*exp((-0.2016D0/xs)*xVal)
    termGb = (-3.2D0/xs)*0.1818D0*exp((-3.2D0/xs)*xVal)+&                 !g'(x)
    (-0.9423D0/xs)*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    (-0.4029D0/xs)*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    (-0.2016D0/xs)*0.02817*exp((-0.2016D0/xs)*xVal)
    termGc = (-3.2D0/xs)**2*0.1818D0*exp((-3.2D0/xs)*xVal)+&                 !g''(x)
    (-0.9423D0/xs)**2*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    (-0.4029D0/xs)**2*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    (-0.2016D0/xs)**2*0.02817*exp((-0.2016D0/xs)*xVal)
    ddy = termFa*termGc+2*termFb*termGb+termFc*termGa
    yArray(3) = ddy
  End Function ZblFull
  
      
  Function EmbeddingA(x, parameters) RESULT (yArray)
! E(p) = a + bp^0.5 , E'(p) = 0.5bp^-0.5 , E''(p) = -0.25bp^-1.5
    Implicit None   ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension(1:2) :: parameters
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Private: Declare variables
! Init
    yArray = 0.0D0
! Calculation
    yArray(1) = parameters(1) + parameters(2)*x**0.5D0
    yArray(2) = 0.5D0 * parameters(2)*x**(-0.5D0)
    yArray(3) = -0.25D0 * parameters(2)*x**(-1.5D0)
  End Function EmbeddingA
        
  Function EmbeddingB(x, parameters) RESULT (yArray)
! E(p) = a + bp^0.5 + cp^2 , E'(p) = 0.5bp^-0.5 + 2cp , E''(p) = -0.25bp^-1.5 + 2c
    Implicit None   ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension(1:3) :: parameters
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Private: Declare variables
! Init
    yArray = 0.0D0
! Calculation
    yArray(1) = parameters(1) + parameters(2)*x**0.5D0 + parameters(3)*x**2.0D0
    yArray(2) = 0.5D0*parameters(2)*x**(-0.5D0) + 2.0D0*parameters(3)*x
    yArray(3) = -0.25D0*parameters(2)*x**(-1.5D0) + 2.0D0*parameters(3)
  End Function EmbeddingB
        
  Function EmbeddingC(x, parameters) RESULT (yArray)
! E(p) = a + bp^0.5 + cp^2 + dc^4 , E'(p) = 0.5bp^-0.5 + 2cp + 4dc^3 , E''(p) = -0.25bp^-1.5 + 2c + 12dc^2
    Implicit None   ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension(1:4) :: parameters
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Private: Declare variables
! Init
    yArray = 0.0D0
! Calculation
    yArray(1) = parameters(1) + parameters(2)*x**0.5D0 + parameters(3)*x**2.0D0 + parameters(4)*x**4.0D0
    yArray(2) = 0.5D0*parameters(2)*x**(-0.5D0) + 2.0D0*parameters(3)*x + 4.0D0*parameters(4)*x**3.0D0
    yArray(3) = -0.25D0*parameters(2)*x**(-1.5D0) + 2.0D0*parameters(3) + 12.0D0*parameters(4)*x**2.0D0
  End Function EmbeddingC
  
  
  
  
  
  
  
  
  
  
  
  
  
  

End Module calcFunctions