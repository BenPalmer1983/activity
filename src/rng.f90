Module rng
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use arrayFunctions
  Use constants
! Force declaration of all variables
  Implicit None
! Public variables  
  Integer(kind=LongInteger) :: randomLCG_n=0
  Integer(kind=LongInteger) :: randomLCG_xn
  Integer(kind=LongInteger) :: randomLCG_R_n=0
  Integer(kind=LongInteger) :: randomLCG_R_xn
! Make private
  Private
! Public
! ---- Variables
  Public :: randomLCG_n
  Public :: randomLCG_xn
  Public :: randomLCG_R_n
  Public :: randomLCG_R_xn
! ---- Functions
  Public :: RandomLCG
  Public :: RandomLCG_R
  Public :: RandomInteger
  Public :: RandomFloat
  Public :: IntegerList
! ---- Subroutines    
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  
  Function RandomLCG(seedIn) RESULT (output)
! Random number - linear congruential generator
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=LongInteger) :: m, a, c, clockTime
    Integer(kind=LongInteger) :: seed
    Integer(kind=StandardInteger), Optional :: seedIn
    Real(kind=DoubleReal) :: output
! init
    seed = 0
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! Read input, reset counter
    If(Present(seedIn))Then
      seed = seedIn
      If(seed.lt.0)Then ! If seed -1 (for example) get random seed
        Call SYSTEM_CLOCK(clockTime) ! "nano seconds" - well, an estimate
        seed = mod(clockTime,m) ! fit in m
      End If
      randomLCG_n = 0
    End If
! If first iteration
    If(randomLCG_n.eq.0)Then
      If(seed.eq.0)Then
        seed = 12791244+45778951 ! Use default seed
      End If
      randomLCG_n = 0
      randomLCG_xn = seed
    End If
! Increment counter
    randomLCG_n = randomLCG_n + 1
! calculate
    randomLCG_xn = mod((a*randomLCG_xn+c),m)
    output = (1.0D0*randomLCG_xn)/(1.0D0*m)
  End Function RandomLCG  
  
  Function RandomLCG_R() RESULT (output)
! Random number - linear congruential generator
! This function starts with a random seed
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=LongInteger) :: m, a, c, clockTime
    Integer(kind=LongInteger) :: seed
    Real(kind=DoubleReal) :: output
! init
    seed = 0
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! If first iteration
    If(randomLCG_R_n.eq.0)Then
! Make "random" seed
      Call SYSTEM_CLOCK(clockTime) ! "nano seconds" - well, an estimate
      seed = mod(clockTime,m)      
      Do i=1,10
        seed = mod((a*seed+c),m)
      End Do
      randomLCG_R_n = 0
      randomLCG_R_xn = seed
    End If
! Increment counter
    randomLCG_R_n = randomLCG_R_n + 1
! calculate
    randomLCG_R_xn = mod((a*randomLCG_R_xn+c),m)
    output = (1.0D0*randomLCG_R_xn)/(1.0D0*m)
  End Function RandomLCG_R
  
  Function RandomInteger(lower,upper) RESULT (randInt)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: lower, upper
    Integer(kind=StandardInteger) :: randInt
    Integer(kind=StandardInteger) :: diff, tempInt
    Real(kind=DoubleReal) :: randDouble
! Make random integer
    If(lower.gt.upper)Then
      tempInt = lower
      lower = upper
      upper = tempInt
    End If
    diff = (upper - lower) + 1
! Call RANDOM_NUMBER(randDouble)
    randDouble = RandomLCG()
    randInt = lower+floor(1.0D0*diff*randDouble)
  End Function RandomInteger

  Function RandomFloat(lower,upper) RESULT (randFloat)
! force declaration of all variables
    Implicit None
! declare variables
    Real(kind=DoubleReal) :: lower, upper
    Real(kind=DoubleReal) :: randFloat
    Real(kind=DoubleReal) :: diff, tempFloat
    Real(kind=DoubleReal) :: randDouble
! Make random float
    If(lower.eq.upper)Then
      lower = 0.0D0
    Else If(lower.gt.upper)Then
      tempFloat = lower
      lower = upper
      upper = tempFloat
    End If
    diff = 1.0D0*(upper - lower)
    randDouble = RandomLCG()
    randFloat = lower+1.0D0*diff*randDouble
  End Function RandomFloat

  Function IntegerList(listStart,listEnd,shuffles) RESULT (list)
! Array filled with integers, possibly shuffled
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: listStart, listEnd, listSize, shuffles, rowA, rowB, shuffleCount
    Integer(kind=StandardInteger), Dimension(1:(listEnd-listStart+1)) :: list
! Initialise variables
    shuffleCount = 0
    listSize = listEnd-listStart+1
! Make list
    Do i=1,listSize
      list(i) = listStart+i-1
    End Do
! Shuffle list
    If(shuffles.gt.0)Then
      Do While(shuffleCount.lt.shuffles)
        rowA = RandomInteger(1,listSize)
        rowB = RandomInteger(1,listSize)
        If(rowA.ne.rowB)Then
          Call swapRows(list,rowA,rowB)
          shuffleCount = shuffleCount + 1
        End If
      End Do
    End If
  End Function IntegerList

End Module rng