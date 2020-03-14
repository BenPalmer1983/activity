Module regression
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
  Use linearAlgebra
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! --variables--!
! --functions--!
  Public :: PolyFit
  Public :: LinearRegression
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

  Function PolyFit(points,order) RESULT (coefficients)
! Fits a polynomial of order to the points input
! Uses Vandermonde Matrix
    Implicit None  !Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Integer(kind=StandardInteger) :: order
! Out
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficients
! Private
    Integer(kind=StandardInteger) :: k,col,row,exponentValue
    Real(kind=DoubleReal), Dimension(1:(order+1),1:(order+1)) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: yMatrix
! Build Least Squares Fitting Vandermonde matrix
    Do row=1,(order+1)
      Do col=1,(order+1)
        exponentValue = row+col-2
        xMatrix(row,col) = 0.0D0
        Do k=1,size(points,1)
          xMatrix(row,col) = 1.0D0*xMatrix(row,col)+1.0D0*points(k,1)&
          **exponentValue
        End Do
      End Do
    End Do
    Do row=1,(order+1)
      exponentValue = row-1
      yMatrix(row) = 0.0D0
      Do k=1,size(points,1)
        yMatrix(row) = 1.0D0*yMatrix(row)+1.0D0*points(k,2)*&
        points(k,1)**exponentValue
      End Do
    End Do
! invert xMatrix
    xMatrix = InvertMatrix(xMatrix)
! multiply inverse by y to get coefficients
    coefficients = matMul(xMatrix,yMatrix)
  End Function PolyFit


  Function LinearRegression(X,Y) RESULT (parameters)
!  Linear regression
!  (XTX)B = XTY
!  f(x) = a_1 x_1 + a_2 x_2 +...+ a_n x_n
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: Y       ! rows        Y input points
    Real(kind=DoubleReal), Dimension(:,:) :: X     ! rows,cols   X input points
! Out
    Real(kind=DoubleReal), Dimension(1:size(X,2)) :: parameters
! Private
    Real(kind=DoubleReal), Dimension(1:size(X,2),1:size(X,1)) :: XT
    Real(kind=DoubleReal), Dimension(1:size(X,2),1:size(X,2)) :: XTX
    Real(kind=DoubleReal), Dimension(1:size(X,2)) :: XTY
! Init vars
    parameters = 0.0D0
! Check the input is OK
    If(size(Y,1).eq.size(X,1))Then
      XT = TransposeMatrix(X)
      XTX = matmul(XT,X)
      XTY = matmul(XT,Y)
      parameters = SolveLinearSet(XTX,XTY)
    End If
  End Function LinearRegression




!  Function LR(dataPoints,paramCount,calcFunction,startPointIn,endPointIn) RESULT (parameters)
!  Linear regression
!  B = (XTX)-1 XTY
!  f(x) = a_1 g_1(x) + a_2 g_2(x) +...+ a_n g_n(x)
!  Some linear function
!    Implicit None  ! Force declaration of all variables
! In:      Declare variables
!    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
!    Integer(kind=StandardInteger) :: paramCount
!    Real(kind=DoubleReal), External :: calcFunction
!    Integer(kind=StandardInteger), Optional :: startPointIn, endPointIn
! Out:     Declare variables
!    Real(kind=DoubleReal), Dimension(1:paramCount) :: parameters
! Private: Declare variables
!    Real(kind=DoubleReal), Dimension(1:Size(dataPoints,1),1:paramCount) :: X
!    Real(kind=DoubleReal), Dimension(1:Size(dataPoints,1)) :: Y
!    X = calcFunction(paramCount,dataPoints)
!  End Function LR


!---------------------------------------------------------------------------------------------------------------------------------------
! Linear Regression functions
!---------------------------------------------------------------------------------------------------------------------------------------

!  Function LR_Embed1(paramCount,dataPoints) RESULT (X)
!  Linear regression function for Embed1
!  f(p) = ap^0.5 + bp^2
!    Implicit None  ! Force declaration of all variables
! In:      Declare variables
!    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
!  End Function LR_Embed1



End Module regression
