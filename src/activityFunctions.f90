! --------------------------------------------------------------!
! Plot
! plotKinds, plotTypes
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Uses Matplotlib python library to build charts
! Requires matplotlib to be installed
!
! Two modules must be used - plotTypes to load the types, and plot to use the plot functions
! and subroutines.  It must be compiled as a part of this library of functions as there are
! quite a lot of modules it depends on.
!
! There is an option to perform one or multiple fits on the input data.
!
! ----------------------------------------
! Updated: 21st September 2015
! ----------------------------------------

Module activityFunctionsTypes
! Setup Modules
  Use kinds

  Type :: decayChainObj
    Real(kind=DoubleReal) :: time = 0.0D0
    !Real(kind=DoubleReal) :: productionRate = 0.0D0
    Integer(kind=StandardInteger) :: isotopes
    Character(Len=16), Dimension(1:100) :: label
    Real(kind=DoubleReal), Dimension(1:100) :: productionRate = 0.0D0
    Real(kind=DoubleReal), Dimension(1:100) :: branchFactor = 1.0D0    ! from isotope parent
    Real(kind=DoubleReal), Dimension(1:100) :: decayConstant = -1.0D0  ! negative for stable
    Real(kind=DoubleReal), Dimension(1:100) :: halfLife = -1.0D0       ! negative for stable
    Real(kind=DoubleReal), Dimension(1:100) :: amountStart = 0.0D0     ! N(0)
    Real(kind=DoubleReal), Dimension(1:100) :: amountEnd = 0.0D0       ! N(t)
    Real(kind=DoubleReal), Dimension(1:100) :: activity = 0.0D0        ! lambda N(t)
  End Type

  Type :: activityTimeObj
    Real(kind=DoubleReal) :: time = 0.0D0
    Real(kind=DoubleReal), Dimension(1:100) :: amount = 0.0D0
    Real(kind=DoubleReal), Dimension(1:100) :: activity = 0.0D0
    Real(kind=DoubleReal), Dimension(1:100) :: amount_gs = 0.0D0
    Real(kind=DoubleReal), Dimension(1:100) :: activity_gs = 0.0D0
  End Type

End Module activityFunctionsTypes


Module activityFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use strings
  Use calcFunctions
  Use solveFunctions
  Use laplaceTransforms
  Use activityFunctionsTypes
  Use printModTypes
  Use printMod
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! --functions--!
  Public :: CalcIsotopeAmount
  Public :: CalcIsotopeAmountGS
  Public :: MaxTrajDepth
  Public :: ActivityCompareGS
  Public :: decayChainComplete
  Public :: decayChainPrint
  Public :: CalcIsotopeChainGS
  Public :: CalcIsotopeChain
  Public :: CalcActivities
  Public :: CalcActivitiesPrint
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------

  Function CalcIsotopeAmount(w,decayDataArray,t,calcOptionIn) RESULT (isotopeChange)
! Force declaration of all variables
    Implicit None
! Declare variables
! Vars:  in
    Real(kind=DoubleReal) :: w, t
    Integer(kind=StandardInteger), optional :: calcOptionIn
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
! Vars:  out
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChange
    Integer(kind=StandardInteger) :: i,decaySteps
    Integer(kind=StandardInteger) :: calcOption
    Type(decayChainObj) :: decayChain
! -------------------------------------------------
! decaySteps really means decay isotopes in chain (steps = decaySteps-1)
! -------------------------------------------------
! Input decay chain array:
! decayDataArray(i,1) !Tally key
! decayDataArray(i,2) !No. Atoms
! decayDataArray(i,3) !Half life
! decayDataArray(i,4) !branching factor
! decayDataArray(i,5) !isotope Z
! decayDataArray(i,6) !isotope A
! -------------------------------------------------
! Output decay chain array:
! isotopeChange(i,1)    !Tally key
! isotopeChange(i,2)    !Change in isotope amount
! isotopeChange(i,3)    !Start amount
! isotopeChange(i,4)    !End amount
! isotopeChange(i,5)    !Isotope Z
! isotopeChange(i,6)    !Isotope A
! isotopeChange(i,7)    !T1/2
! isotopeChange(i,8)    !Decay constant
! isotopeChange(i,9)    !Branching factor
! isotopeChange(i,10)   !Parent production rate
! isotopeChange(i,11)   !Time
! isotopeChange(i,12)   !GS End
! -------------------------------------------------
! Optional arguments
    calcOption = 1  ! 1 = Analytic, 2 = Numeric
    If(Present(calcOptionIn))Then
      calcOption = calcOptionIn
    End If
! set decay steps/isotopes
    decaySteps = size(decayDataArray,1)
! allocate isotopeChange array
    Allocate(isotopeChange(1:decaySteps,1:12))
! Fill with starting data
    Do i=1,decaySteps
      isotopeChange(i,1) = decayDataArray(i,1)
      isotopeChange(i,2) = 0.0D0          !default no change
      isotopeChange(i,3) = decayDataArray(i,2)
      isotopeChange(i,4) = decayDataArray(i,2)    !default no change
      isotopeChange(i,5) = decayDataArray(i,5)
      isotopeChange(i,6) = decayDataArray(i,6)
      isotopeChange(i,7) = decayDataArray(i,3)
      isotopeChange(i,8) = log(2.0D0)/decayDataArray(i,3)
      isotopeChange(i,9) = decayDataArray(i,4)
      isotopeChange(i,10) = w
      isotopeChange(i,11) = t
      isotopeChange(i,12) = 0.0D0          !default no change
    End Do
! Init decay chain
    decayChain%time = t
    Do i=1,decaySteps
      decayChain%label(i) = " "
      decayChain%halfLife(i) = decayDataArray(i,3)
      decayChain%amountStart(i) = decayDataArray(i,2)
      If(i.eq.1)Then
        decayChain%branchFactor(i) = 1.0D0  ! Not used
        decayChain%productionRate(i) = w
      Else
        decayChain%branchFactor(i) = decayDataArray(i,4)
        decayChain%productionRate(i) = 0.0D0
      End If
    End Do
! Calculate
    Call CalcIsotopeChain(decayChain)
! Store changes in isotope amounts
    Do i=1,size(isotopeChange,1)
      isotopeChange(i,4) = decayChain%amountEnd(i)
      isotopeChange(i,2) = decayChain%amountEnd(i) - decayChain%amountStart(i)
    End Do
  End Function CalcIsotopeAmount

  Function CalcIsotopeAmountGS(t,isotopeStep,isotopeChangeIn) RESULT (output)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i, isotopeStep, M, k
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChangeIn
    Real(kind=QuadrupoleReal), Dimension(1:50) :: weightingQ
    Real(kind=QuadrupoleReal), Dimension(1:20) :: L ! Lambda
    Real(kind=QuadrupoleReal), Dimension(1:20) :: N ! Starting number of atoms
    Real(kind=QuadrupoleReal), Dimension(1:20) :: B ! Starting number of atoms
    Real(kind=QuadrupoleReal) :: kQ, w, t, ft, s, FS, output
! -------------------------------------------------
! Output decay chain array:
! isotopeChange(i,1)    !Tally key
! isotopeChange(i,2)    !Change in isotope amount
! isotopeChange(i,3)    !Start amount
! isotopeChange(i,4)    !End amount
! isotopeChange(i,5)    !Isotope Z
! isotopeChange(i,6)    !Isotope A
! isotopeChange(i,7)    !T1/2
! isotopeChange(i,8)    !Decay constant
! isotopeChange(i,9)    !Branching factor
! isotopeChange(i,10)   !Parent production rate
! -------------------------------------------------
! Init variables
    M = 8
    weightingQ = GaverStehfestWeightingQ(M,weightingQ)
    w = isotopeChangeIn(1,10)
    output = 0.0D0
! Adjust the isotope values
    Do i=1,isotopeStep
      If(isotopeChangeIn(i,4).lt.0.0D0)Then
        isotopeChangeIn(i,4) = 0.0D0
      End If
    End Do
! Store lambda starting atom number data
    Do i=1,isotopeStep
      L(i) = lnTwoQ/isotopeChangeIn(i,7)
      N(i) = isotopeChangeIn(i,3)
      If(i.eq.1)Then
        B(i) = 1.0D0
      Else
        B(i) = isotopeChangeIn(i,9)
      End If
    End Do
! Perform calculation
    ft = 0.0D0
    Do k=1,2*M
      kQ = 1.0D0 * k
      s = (kQ*lnTwoQ)/t
! -----------------------
      FS = (1.0D0/(s+L(1)))*(w/s+N(1))
      Do i=2,isotopeStep
        FS = (1.0D0/(s+L(i)))*(B(i)*L(i-1)*FS+N(2))
      End Do
! FS = (1.0D0/(s+L(1)))*(w/s+N(1))
! -----------------------
      ft = ft + weightingQ(k)*FS
    End Do
    ft = (lnTwoQ/t)*ft
    output = Dble(ft)
! isotopeChangeOut(isotopeStep,4) = Dble(ft)
  End Function CalcIsotopeAmountGS

  Function MaxTrajDepth(coefficients, maxDepthIn) RESULT (maxDepth)
! Calc max depth (E = 0) for ion trajectory, described by polynomial
    Implicit None  ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, j, nMax
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal) :: c, x, y, maxDepth, maxDepthL
    Real(kind=DoubleReal), Optional :: maxDepthIn
! Set optional argument
    maxDepth = 1.0D10     ! 1m in ang
    If(Present(maxDepthIn))Then
      maxDepth = maxDepthIn
    End If
! Do three refinement loops
    Do i=1,4
      nMax = 20+(i*10)
      c = 10**(log10(maxDepth)/(1.0D0*nMax))
      Do j=1,50
        If(i.lt.3)Then
          x = 1.0D0*c**j
        Else
          x = (maxDepth+1.0D0)-c**(nMax-j)
        End If
        y = CalcPolynomial(coefficients, x)
        If(y.lt.0.0D0)Then
          maxDepth = x
          Exit
        End If
        If(i.eq.4)Then
          maxDepthL = x
        End If
      End Do
    End Do
    maxDepth = SolvePolynomial (coefficients, maxDepthL, maxDepth, 1.0D-6)
  End Function MaxTrajDepth

  Function ArraySize1DDouble (inputArray,arraySize) RESULT (outputArray)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: arraySize
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: outputArray
! Allocate output array
    Allocate(outputArray(1:arraySize))
! transfer data
    Do i=1,arraySize
      If(i.le.size(inputArray))Then
        outputArray(i) = inputArray(i)
      Else
        outputArray(i) = 0.0D0
      End If
    End Do
  End Function ArraySize1DDouble

  Function ArraySize2DDouble (inputArray,arraySizeHeight,arraySizeWidthIn) &
    RESULT (outputArray)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i, j
    Integer(kind=StandardInteger) :: arraySizeHeight
    Integer(kind=StandardInteger), optional :: arraySizeWidthIn
    Integer(kind=StandardInteger) :: arraySizeWidth
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: outputArray
! catch optional width
    If(present(arraySizeWidthIn))Then
      arraySizeWidth = arraySizeWidthIn
    Else
      arraySizeWidth = size(inputArray,2)
    End If
! Allocate output array
    Allocate(outputArray(1:arraySizeHeight,1:arraySizeWidth))
! transfer data
    Do i=1,arraySizeHeight
      Do j=1,arraySizeWidth
        If(i.le.size(inputArray,1).and.j.le.size(inputArray,2))Then
          outputArray(i,j) = inputArray(i,j)
        Else
          outputArray(i,j) = 0.0D0
        End If
      End Do
    End Do
  End Function ArraySize2DDouble



! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------

  Subroutine ActivityCompareGS(lambda, dataPointsCount, endTime, nAtoms, productionRate)
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal) :: lambda
    Integer(kind=StandardInteger) :: dataPointsCount
    Real(kind=DoubleReal) :: endTime, nAtoms, productionRate
! Private: Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: analyticResult, numericResult, numericResultB, t
    Real(kind=DoubleReal), Dimension(1:10) :: p
    p(1) = productionRate
    p(2) = lambda
    p(3) = nAtoms
!
    Do i=1,dataPointsCount
      t = (i-1)*(endTime/(dataPointsCount-1))
      analyticResult = productionRate*((1/lambda)-(exp(-1.0D0*lambda*t)/lambda))+nAtoms*exp(-1.0D0*lambda*t)
      If(t.eq.0.0D0)Then
        numericResult = nAtoms
      Else
        !numericResult = GaverStehfest(ltDecay,t,p,7)
        !numericResultB = GaverStehfest(ltDecay,t,p,8)
      End If
      print *,t,analyticResult,numericResult,numericResultB
    End Do
  End Subroutine ActivityCompareGS

  Subroutine decayChainComplete(decayChain)
! Completes the decay chain object
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Type(decayChainObj) :: decayChain
! Vars Private
    Integer(kind=StandardInteger) :: i, n
    n = 0
    Do i=1,100
      n = n + 1
      If(decayChain%decayConstant(i).eq.-1.0D0.and.decayChain%halfLife(i).gt.0.0D0)Then
! complete decay constant from half life
        decayChain%decayConstant(i) = lnTwo/decayChain%halfLife(i)
      End If
      If(decayChain%halfLife(i).le.0.0D0.and.decayChain%decayConstant(i).gt.0.0D0)Then
! complete decay constant from half life
        decayChain%halfLife(i) = lnTwo/decayChain%decayConstant(i)
      End If
! Adjust for stable isotope
      If(decayChain%decayConstant(i).lt.0.0D0)Then
        decayChain%halfLife(i) = -1.0D0
        decayChain%decayConstant(i) = 0.0D0
      End If
      If(decayChain%halfLife(i).lt.0.0D0)Then
        decayChain%halfLife(i) = -1.0D0
        decayChain%decayConstant(i) = 0.0D0
      End If
! Break out if stable
      If(decayChain%decayConstant(i).eq.0.0D0)Then
        Exit
      End If
    End Do
    decayChain%isotopes = n
  End Subroutine decayChainComplete


  Subroutine decayChainPrint(decayChain)
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Type(decayChainObj) :: decayChain
    Integer(kind=StandardInteger) :: i
! Output
    Print "(A96)","================================================================================================"
    Print "(A32,F16.8)","Time:                      ",decayChain%time
    Print "(A32,I8)","Isotopes in chain: ",decayChain%isotopes
    Print "(A8,A20,A8,A20,A20,A20,A20)","n","w_n","B_n-1,n","L_n","T1/2_n","nStart","nEnd"
    Do i=1,decayChain%isotopes
      Print "(I8,E20.8,F8.4,E20.8,E20.8,E20.8,E20.8)",&
      i,decayChain%productionRate(i),decayChain%branchFactor(i),decayChain%decayConstant(i),decayChain%halfLife(i),&
      decayChain%amountStart(i),decayChain%amountEnd(i)
    End Do
    Print "(A96)","================================================================================================"
    Print *,""
  End Subroutine decayChainPrint


  Subroutine CalcIsotopeChainGS(decayChain)
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Type(decayChainObj) :: decayChain
! Vars Private
    Integer(kind=StandardInteger) :: n, k
    Real(kind=DoubleReal) :: nEnd
    Real(kind=DoubleReal) :: t
    Real(kind=DoubleReal), Dimension(1:100) :: p
! Complete decay chain data
    Call decayChainComplete(decayChain)
! Store time
    t = decayChain%time
! Prepare parameters array
    p(1) = 0   ! Chain to calculate
    p(2) = decayChain%isotopes
    Do n=1,decayChain%isotopes
      k = 2+4*(n-1)+1
      p(k) = decayChain%productionRate(n)
      k = 2+4*(n-1)+2
      If(decayChain%decayConstant(n).le.0.0D0)Then  ! Set very small decay constant if stable
        p(k) = 1.0D-30
      Else
        p(k) = decayChain%decayConstant(n)
      End If
      k = 2+4*(n-1)+3
      p(k) = decayChain%branchFactor(n)
      k = 2+4*(n-1)+4
      p(k) = decayChain%amountStart(n)
    End Do
! Calculate isotope amounts
    Do n=1,decayChain%isotopes
      p(1) = n
      nEnd = GaverStehfest(ltDecay, t, p, 8)
      decayChain%amountEnd(n) = nEnd
    End Do
  End Subroutine CalcIsotopeChainGS







  Subroutine CalcIsotopeChain(decayChain)
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Type(decayChainObj) :: decayChain
! Vars Private
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: t, nEnd
    Real(kind=DoubleReal), Dimension(1:100) :: W ! Production Rate
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: N ! Starting number of atoms
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Exp
! Complete decay chain data
    Call decayChainComplete(decayChain)
! store input in shortned name arrays to make equations clearer
    t = decayChain%time
    Do i=1,decayChain%isotopes
      W(i) = decayChain%productionRate(i)
      L(i) = decayChain%decayConstant(i)
      B(i) = decayChain%branchFactor(i)
      N(i) = decayChain%amountStart(i)
    End Do
! Break infinities
    Call DecayBreakInfinities(L,decayChain%isotopes)
! Run analytic calculations
! Loop through isotopes
! Using L-1(1/(q+ps) = (1/p)*exp(-1*(q/p)*t) and partial fractions
    Do i=1,decayChain%isotopes
      nEnd = CalcIsotopeChainCalc(t,i,W,L,N,B)
      decayChain%amountEnd(i) = nEnd
      decayChain%activity(i) = nEnd*decayChain%decayConstant(i)
    End Do
  End Subroutine CalcIsotopeChain

  Subroutine DecayBreakInfinities(L,n)
!
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Real(kind=DoubleReal), Dimension(:) :: L ! Lambda
    Integer(kind=StandardInteger) :: n
! Vars Private
    Integer(kind=StandardInteger) :: i,j
    Logical :: breaking
! Loop and alter matching decay constants slightly
    breaking = .true.
    Do While(breaking)
      breaking = .false.
      Do i=1,n-1
        Do j=i+1,n
          If(L(i).eq.L(j))Then
            breaking = .true.
            L(i) = L(i)*1.0000001D0  ! Vary by 0.00001%
          End If
        End Do
      End Do
    End Do
  End Subroutine DecayBreakInfinities

  Function CalcIsotopeChainCalc(t,m,W,L,N,B) Result (nEnd)
    Implicit None ! Force declaration of all variables
! Vars In
    Real(kind=DoubleReal) :: t
    Integer(kind=StandardInteger) :: m
    Real(kind=DoubleReal), Dimension(1:100) :: W ! Production Rate
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: N ! Starting number of atoms
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Branch factor
! Vars Out
    Real(kind=DoubleReal) :: nEnd
! Vars Private
    Integer(kind=StandardInteger) :: k
    Real(kind=DoubleReal) :: mult, multR
    Real(kind=DoubleReal) :: nChange
! Init
    nEnd = 0.0D0
! ---------------------------------------------
! UNSTABLE Isotopes
! ---------------------------------------------
    If(L(m).gt.0.0D0)Then
! Loop through terms
      Do k=1,m
        multR = CalcIsotopeChainMultR(k,m,L,B)
        mult = multR * N(k)
! decay of starting matter
        nChange = CalcIsotopeChainF_Unstable(k,m,t,mult,L)
        nEnd = nEnd + nChange
! production term
        mult = multR * W(k)
        nChange = CalcIsotopeChainG_Unstable(k,m,t,mult,L)
        nEnd = nEnd + nChange
      End Do
    Else
! ---------------------------------------------
! STABLE Isotopes
! ---------------------------------------------
! Loop through terms
      nEnd = N(m)+t*W(m)
      Do k=1,m-1
        multR = CalcIsotopeChainMultR(k,m,L,B)
        mult = multR * N(k)
! decay of starting matter
        nChange = CalcIsotopeChainF_Stable(k,m,t,mult,L)
        nEnd = nEnd + nChange
! production term
        mult = multR * W(k)
        nChange = CalcIsotopeChainG_Stable(k,m,t,mult,L)
        nEnd = nEnd + nChange
      End Do
    End If
  End Function CalcIsotopeChainCalc

  Function CalcIsotopeChainMultR(k,m,L,B) Result (multR)
! Vars In
    Integer(kind=StandardInteger) :: k, m
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Branching Factor
! Vars Out
    Real(kind=DoubleReal) :: multR
! Private
    Integer(kind=StandardInteger) :: i
! Result
    multR = 1.0D0
    If(k.lt.m)Then
      Do i=k,m-1
        multR = multR * B(i+1) * L(i)
      End Do
    End If
  End Function CalcIsotopeChainMultR

! --------------------------
! Unstable
! --------------------------

  Function CalcIsotopeChainF_Unstable(k,m,t,mult,L) Result (nChange)
! Vars In
    Integer(kind=StandardInteger) :: k, m
    Real(kind=DoubleReal) :: t, mult
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
! Vars Out
    Real(kind=DoubleReal) :: nChange
! Private
    Real(kind=DoubleReal) :: multP
    Real(kind=DoubleReal) :: r
    Integer(kind=StandardInteger) :: i, j
! Calculate isotope amount change
    nChange = 0.0D0
    multP = (-1.0D0)**(m-k)
! Loop through pfrac
    Do i=k,m
      r = 1.0D0
      Do j=k,m
        If(j.ne.i)Then
          r = r * (L(i)-L(j))
        End If
      End Do
      nChange = nChange + (1.0D0/r)*exp(-1.0D0*L(i)*t)*multP*mult
    End Do
  End Function CalcIsotopeChainF_Unstable

  Function CalcIsotopeChainG_Unstable(k,m,t,mult,L) Result (nChange)
! Vars In
    Integer(kind=StandardInteger) :: k, m
    Real(kind=DoubleReal) :: t, mult
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
! Vars Out
    Real(kind=DoubleReal) :: nChange
! Private
    Real(kind=DoubleReal) :: multP
    Real(kind=DoubleReal) :: r
    Integer(kind=StandardInteger) :: i, j
! Calculate isotope amount change
    nChange = 0.0D0
    multP = (-1.0D0)**(m-k+1)
  ! term A
    r = 1.0D0
    Do i=k,m
      r = r * L(i)
    End Do
    nChange = nChange + (1.0D0/r)*mult
! Loop through pfrac
    Do i=k,m
      r = 1.0D0*L(i)
      Do j=k,m
        If(j.ne.i)Then
          r = r * (L(i)-L(j))
        End If
      End Do
      nChange = nChange + (1.0D0/r)*exp(-1.0D0*L(i)*t)*multP*mult
    End Do
  End Function CalcIsotopeChainG_Unstable

! --------------------------
! Stable
! --------------------------

  Function CalcIsotopeChainF_Stable(k,mIn,t,mult,L) Result (nChange)
! Vars In
    Integer(kind=StandardInteger) :: k, mIn
    Real(kind=DoubleReal) :: t, mult
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
! Vars Out
    Real(kind=DoubleReal) :: nChange
! Private
    Integer(kind=StandardInteger) :: m
    Real(kind=DoubleReal) :: multP
    Real(kind=DoubleReal) :: r
    Integer(kind=StandardInteger) :: i, j
! In
    m = mIn-1
! Calculate isotope amount change
    nChange = 0.0D0
    multP = (-1.0D0)**(m-k+1)
! term A
    r = 1.0D0
    Do i=k,m
      r = r * L(i)
    End Do
    nChange = nChange + (1.0D0/r)*mult
! Loop through pfrac
    Do i=k,m
      r = 1.0D0*L(i)
      Do j=k,m
        If(j.ne.i)Then
          r = r * (L(i)-L(j))
        End If
      End Do
      nChange = nChange + (1.0D0/r)*exp(-1.0D0*L(i)*t)*multP*mult
    End Do
  End Function CalcIsotopeChainF_Stable


  Function CalcIsotopeChainG_Stable(k,mIn,t,mult,L) Result (nChange)
  ! Vars In
      Integer(kind=StandardInteger) :: k, mIn
      Real(kind=DoubleReal) :: t, mult
      Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
  ! Vars Out
      Real(kind=DoubleReal) :: nChange
  ! Private
      Integer(kind=StandardInteger) :: m
      Real(kind=DoubleReal) :: multP
      Real(kind=DoubleReal) :: p, q, r
      Integer(kind=StandardInteger) :: i, j
  ! In
      m = mIn-1
  ! Calculate isotope amount change
      nChange = 0.0D0
      multP = (-1.0D0)**(m-k)
  ! term A
      r = 1.0D0
      Do i=k,m
        r = r * L(i)
      End Do
      nChange = nChange + (1.0D0/r)*t*mult
  ! term B
      p = CalcIsotopeChainC(L,k,m)
      q = 1.0D0
      Do i=k,m
        q = q*L(i)*L(i)
      End Do
      r = (-1.0D0)*(p/q)
      nChange = nChange + r*mult
  ! Loop through pfrac
      Do i=k,m
        r = 1.0D0*L(i)*L(i)
        Do j=k,m
          If(j.ne.i)Then
            r = r * (L(i)-L(j))
          End If
        End Do
        nChange = nChange + (1.0D0/r)*exp(-1.0D0*L(i)*t)*multP*mult
      End Do
    End Function CalcIsotopeChainG_Stable

  Function CalcIsotopeChainC(L,k,m) Result (numerator)
! Calculates numerator in isotope activity function
    Implicit None ! Force declaration of all variables
! Vars In
    Real(kind=DoubleReal), Dimension(:) :: L
    Integer(kind=StandardInteger) :: k, m
! Vars Out
    Real(kind=DoubleReal) :: numerator
! Vars Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: tempMult
    numerator = 0.0D0
    Do i=k,m
      tempMult = 1.0D0
      Do j=k,m
        If(j.ne.i)Then
          tempMult = tempMult * L(j)
        End If
      End Do
      numerator = numerator + tempMult
    End Do
  End Function CalcIsotopeChainC

! --------------------------
! Calc Activities
! --------------------------

  Subroutine CalcActivities(decayChain,activityTime,endTime,zeroProductionTimeIn)
! calculate activities over time
! Vars In
    Type(decayChainObj) :: decayChain
    Type(activityTimeObj), Dimension(:) :: activityTime
    Real(kind=DoubleReal) :: endTime
    Real(kind=DoubleReal), Optional :: zeroProductionTimeIn
! Vars Private
    Integer(kind=StandardInteger) :: i, j, timeSteps
    Real(kind=DoubleReal) :: t
    Real(kind=DoubleReal) :: zeroProductionTime
    Logical :: productionOff
    Real(kind=DoubleReal), Dimension(1:100) :: n0w0
    Real(kind=DoubleReal), Dimension(1:100) :: n0w0_gs
! Optional arguments
    productionOff = .false.
    zeroProductionTime = 2.0D0*endTime
    If(Present(zeroProductionTimeIn))Then
      print *,endTime,zeroProductionTimeIn
      zeroProductionTime = zeroProductionTimeIn
    End If
! Time steps
    timeSteps = Size(activityTime,1)
! Loop through steps and calculate activity
    Do i=1,timeSteps
      t = (1.0D0*(i/(1.0D0*timeSteps)))*endTime
      activityTime(i)%time = t
      decayChain%time = t
      If((productionOff.eqv..false.).and.(t.ge.zeroProductionTime))Then
        productionOff = .true.
! Store post irradiation/source starting amounts
        Call CalcIsotopeChain(decayChain)
        Do j=1,decayChain%isotopes
          n0w0(j) = decayChain%amountEnd(j)
        End Do
! Store post irradiation/source starting amounts - gs
        Call CalcIsotopeChainGS(decayChain)
        Do j=1,decayChain%isotopes
          n0w0_gs(j) = decayChain%amountEnd(j)
          If(n0w0_gs(j).lt.0.0D0)Then
            n0w0_gs(j) = 0.0D0
          End If
        End Do
! Turn off source
        Do j=1,decayChain%isotopes
          decayChain%productionRate(j) = 0.0D0
        End Do
      End If
! Analytic calculation
      If(productionOff)Then
        decayChain%time = t-zeroProductionTime
        Do j=1,decayChain%isotopes
          decayChain%amountStart(j) = n0w0(j)
        End Do
      End If
      Call CalcIsotopeChain(decayChain)
      Do j=1,decayChain%isotopes
        activityTime(i)%amount(j) = decayChain%amountEnd(j)
        activityTime(i)%activity(j) = decayChain%activity(j)
      End Do
! Numeric calculation
      If(productionOff)Then
        decayChain%time = t-zeroProductionTime
        Do j=1,decayChain%isotopes
          decayChain%amountStart(j) = n0w0_gs(j)
        End Do
      End If
      Call CalcIsotopeChainGS(decayChain)
      If(t.eq.zeroProductionTime)Then
        Do j=1,decayChain%isotopes
          activityTime(i)%amount_gs(j) = n0w0_gs(j)
          activityTime(i)%activity_gs(j) = n0w0_gs(j)*decayChain%decayConstant(j)
        End Do
      Else
        Do j=1,decayChain%isotopes
          activityTime(i)%amount_gs(j) = decayChain%amountEnd(j)
          activityTime(i)%activity_gs(j) = decayChain%activity(j)
        End Do
      End If
    End Do
  End Subroutine CalcActivities



  Subroutine CalcActivitiesPrint(decayChain,activityTime,tableIn)
! calculate activities over time
! Vars In
    Type(decayChainObj) :: decayChain
    Type(activityTimeObj), Dimension(:) :: activityTime
    Logical, Optional :: tableIn
! Vars Private
    Integer(kind=StandardInteger) :: i, j, timeSteps
    Character(Len=16) :: tempLine
    Character(Len=1024) :: headerLine, dataLine
    Type(tableObj) :: activityTable
    Logical :: table
! Optional
    table = .true.
    If(Present(tableIn))Then
      table = tableIn
    End If
! Time steps
    timeSteps = Size(activityTime,1)
    If(table)Then
! Init table
      Call printTableInit(activityTable)
! Table settings
      activityTable%colAutoWidth = .false.
      activityTable%printHeaderRow = .true.
      activityTable%printHeaderColumn = .true.
! Headers
      activityTable%headerRowColumn = "Time/s"
      Do j=1,decayChain%isotopes
        activityTable%headerRow(j) = decayChain%label(j)
      End Do
      Do j=1,decayChain%isotopes
        tempLine = trim(decayChain%label(j))//" (GS)"
        activityTable%headerRow(decayChain%isotopes+j) = tempLine
      End Do
! Loop through steps and calculate activity
      Do i=1,timeSteps
        Write(tempLine,"(E14.7)") activityTime(i)%time
        activityTable%headerColumn(i) = tempLine
        Do j=1,decayChain%isotopes
          Write(tempLine,"(E14.7)") activityTime(i)%amount(j)
          activityTable%tableData(i,j) = tempLine
          Write(tempLine,"(E14.7)") activityTime(i)%amount_gs(j)
          activityTable%tableData(i,j+decayChain%isotopes) = tempLine
        End Do
      End Do
      activityTable%columns = 2*decayChain%isotopes
      activityTable%rows = timeSteps
      Call printTableMake(activityTable)
    Else
      headerLine = BlankString(headerLine)
      headerLine = "Time/s"
      Do j=1,decayChain%isotopes
        tempLine = trim(decayChain%label(j))
        headerLine = trim(headerLine)//","//trim(tempLine)
      End Do
      Do j=1,decayChain%isotopes
        tempLine = trim(decayChain%label(j))//" (GS)"
        headerLine = trim(headerLine)//","//trim(tempLine)
      End Do
      print *,trim(headerLine)
! Loop through steps and calculate activity
      Do i=1,timeSteps
        dataLine = BlankString(dataLine)
        Write(tempLine,"(E14.7)") activityTime(i)%time
        dataLine = trim(tempLine)
        Do j=1,decayChain%isotopes
          Write(tempLine,"(E14.7)") activityTime(i)%amount(j)
          activityTable%tableData(i,j) = tempLine
          dataLine = trim(dataLine)//","//trim(tempLine)
        End Do
        Do j=1,decayChain%isotopes
          Write(tempLine,"(E14.7)") activityTime(i)%amount_gs(j)
          activityTable%tableData(i,j) = tempLine
          dataLine = trim(dataLine)//","//trim(tempLine)
        End Do
        print *,trim(dataLine)
      End Do


    End If
  End Subroutine CalcActivitiesPrint






















End Module activityFunctions
