Module basicMaths
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
  Public :: RoundDP
  Public :: BinomialCoefficient, BinomialCoefficientDP, BinomialCoefficientQ
  Public :: Odd
  Public :: Even
  Public :: RSSCalc
  Public :: RSSPoints
  Public :: Modulus
  Public :: CompareSign
! Primitive recursive
  Public :: Factorial, FactorialDP, FactorialQ
  Public :: Fib
! Recursive
  Public :: Ackermann
! Combinations and Permutations
  Public :: combinationsPrint

! Interfaces
  Interface Modulus
    Module Procedure Modulus_R, Modulus_I
  End Interface Modulus
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

  Function RoundDP(dpIn) RESULT (intOut)
! Round DP to nearest int
    Implicit None ! Force declaration of all variables
    Real(kind=DoubleReal) :: dpIn
    Integer(kind=StandardInteger) :: intOut
    intOut = Floor(dpIn+0.5D0)
  End Function RoundDP

  Function BinomialCoefficient(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: c,n,k
! calculate factorial
    c = Factorial(n)/(Factorial(n-k)*Factorial(k))
  End Function BinomialCoefficient

  Function BinomialCoefficientDP(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: n,k
    Real(kind=DoubleReal) :: c, nDP, nkDP, kDP
! calculate factorial
    nDP = FactorialDP(n)
    nkDP = Factorial(n-k)
    kDP = FactorialDP(k)
    c = 1.0D0*nDP/(nkDP*kDP)
  End Function BinomialCoefficientDP

  Function BinomialCoefficientQ(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: n,k
    Real(kind=QuadrupoleReal) :: c, nDP, nkDP, kDP
! calculate factorial
    nDP = FactorialDP(n)
    nkDP = Factorial(n-k)
    kDP = FactorialDP(k)
    c = 1.0D0*nDP/(nkDP*kDP)
  End Function BinomialCoefficientQ

  Function Odd(input) RESULT (output)
! Returns true if odd, false if even
    Implicit None  ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: input
    Real(kind=DoubleReal) :: dpA, dpB
    Logical :: output
    output = .false.
    dpA = 1.0D0*input
    dpB = 2.0D0*ceiling(dpA/2.0D0)
    If(dpB.gt.dpA)Then
      output = .true.
    End If
  End Function Odd

  Function Even(input) RESULT (output)
! Returns true if even, false if odd
    Implicit None  ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: input
    Real(kind=DoubleReal) :: dpA, dpB
    Logical :: output
    output = .true.
    dpA = 1.0D0*input
    dpB = 2.0D0*ceiling(dpA/2.0D0)
    If(dpB.gt.dpA)Then
      output = .false.
    End If
  End Function Even

  Function RSSCalc(inputA, inputB, factorIn) RESULT (output)
! Get value of function at x
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: inputA, inputB
    Real(kind=DoubleReal), optional :: factorIn
    Real(kind=DoubleReal) :: factor, output
    factor = 1.0D0
    If(Present(factorIn))Then
      factor = factorIn
    End If
    output = 0.0D0
    If(inputA.gt.-2.0D20.and.inputB.gt.-2.0D20)Then
      output = factor*(inputA-inputB)**2
    End If
  End Function RSSCalc

  Function RSSPoints(xArr, yArr) RESULT (output)
! Calculate RSS between values in arrays
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: xArr
    Real(kind=DoubleReal), Dimension(:,:) :: yArr
! Out:     Declare variables
    Real(kind=DoubleReal) :: output
! Private: Declare variables
    Integer(kind=StandardInteger) :: n
! Loop through points (if data sets same size)
    output = 0.0D0
    If(size(xArr,1).eq.size(yArr,1))Then
      Do n=1,size(xArr,1)
        output = output + (xArr(n,2)-yArr(n,2))**2
      End Do
    End If
  End Function RSSPoints

  Function Modulus_I(x, divisor) RESULT (y)
! calc modulus
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: x, y, divisor, factor
    If(x.lt.0)Then
      factor = ceiling(abs(x/(1.0D0*divisor)))
      y = x+factor*divisor
    Else
      factor = floor(x/(1.0D0*divisor))
      y = x-factor*divisor
    End If
  End Function Modulus_I

  Function Modulus_R(x, divisor) RESULT (y)
! calc modulus
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: x, y, divisor, factor
    If(x.lt.0.0D0)Then
      factor = ceiling(abs(x/(1.0D0*divisor)))
      y = x+factor*divisor
    Else
      factor = floor(x/(1.0D0*divisor))
      y = x-factor*divisor
    End If
  End Function Modulus_R

  Function CompareSign(x,y) Result (output)
! True is the same, false if + and -, "0" will be classed as positive
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: x, y
    Logical :: output
    output = .false.
    If(x.ge.0.0D0.and.y.lt.0.0D0)Then
      output = .true.
    Else
      If(x.lt.0.0D0.and.y.ge.0.0D0)Then
        output = .true.
      End If
    End If
  End Function CompareSign




Function Factorial(input) RESULT (output)
! force declaration of all variables
  Implicit None
! declare variables
  Integer(kind=StandardInteger) :: i,input
  Integer(kind=StandardInteger) :: output
! calculate factorial
  output = 1
  Do i=1,input
    output = i * output
  End Do
End Function Factorial

Function FactorialDP(input) RESULT (output)
! force declaration of all variables
  Implicit None
! declare variables
  Integer(kind=StandardInteger) :: i,input
  Integer(kind=VeryLongInteger) :: tempInt
  Real(kind=DoubleReal) :: output
! calculate factorial
  tempInt = 1
  Do i=1,input
    tempInt = i * tempInt
  End Do
  output = 1.0D0*tempInt
End Function FactorialDP

Function FactorialQ(input) RESULT (output)
! force declaration of all variables
  Implicit None
! declare variables
  Integer(kind=StandardInteger) :: i,input
  Integer(kind=VeryLongInteger) :: tempInt
  Real(kind=QuadrupoleReal) :: tempQ
  Real(kind=QuadrupoleReal) :: output
! calculate factorial
  tempInt = 1
  tempQ = 1
  Do i=1,input
    If(i.le.33)Then
      tempInt = i * tempInt
      tempQ = 1.0D0*tempInt
    End If
    If(i.eq.34)Then
      tempQ = 1.0D0*i*tempInt
    End If
    If(i.ge.35)Then
      tempQ = 1.0D0*i*tempQ
    End If
  End Do
  output = tempQ
End Function FactorialQ


  Function Fib(input) RESULT (output)
! Fibonacci sequence
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Integer(kind=StandardInteger) :: input
! Out:     Declare variables
    Integer(kind=StandardInteger) :: output
! Private: Declare variables
    Integer(kind=StandardInteger) :: i, nA, nB, nC
! calculate fib
    If(input.eq.1)Then
      output = 1
    Elseif(input.eq.2)Then
      output = 1
    Else
      nA = 1
      nB = 1
      Do i=3,input
        nC = nA + nB
        nA = nB
        nB = nC
      End Do
      output = nC
    End If
  End Function Fib

  Recursive Function Ackermann(x,y) RESULT (z)
! Ackermann recursive function
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Integer(kind=StandardInteger) :: x, y
! Out:     Declare variables
    Integer(kind=StandardInteger) :: z

    If(x.eq.0)Then
      z = y + 1
    Elseif(y.eq.0)Then
      z = Ackermann(x-1,1)
    Else
      z = Ackermann(x-1,Ackermann(x,y-1))
    End If
  End Function Ackermann


! ----------------------------------------------
! Combinations and Permutations
! ----------------------------------------------
  Subroutine combinationsPrint(setSize, numbersUsed)
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Integer(kind=StandardInteger) :: setSize, numbersUsed
! Vars Private
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger), Dimension(1:numbersUsed) :: combinationSet
    Character(Len=1) :: charTemp
    Character(Len=16) :: outTemp
    Logical :: loopCombinations

    Do i=1,numbersUsed
      combinationSet(i) = i
    End Do

    loopCombinations = .true.
    k = 0
    Do while(loopCombinations)
      k = k + 1
      If(k.gt.1)Then
        j = numbersUsed
        Do i=1,numbersUsed
          loopCombinations = .false.
          If(combinationSet(j).lt.(setSize+1-i))Then
            combinationSet(j) = combinationSet(j) + 1
            loopCombinations = .true.
            Exit
          End If
          j = j - 1
        End Do
      End If
      If(loopCombinations)Then
        outTemp = "                "
        Do i=1,numbersUsed
          Write(charTemp,"(I1)") combinationSet(i)
          outTemp(i:i) = charTemp
        End Do
        print *,outTemp
      End If
    End Do
  End Subroutine combinationsPrint





End Module basicMaths
