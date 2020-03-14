Module units
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! 
! 
! 
  Use kinds
  Use strings
! Force declaration of all variables
  Implicit None
! Define type
  Type :: convType
    Character(len=8) :: unitName="ZZZZZZZZ"
    Real(kind=DoubleReal) :: factor=0.0D0
  End Type
! Variables  
  Integer(kind=StandardInteger), Parameter :: convTableSize = 100
  Type(convType), Dimension(1:convTableSize) :: convTable
  Integer(kind=StandardInteger) :: convTableN = 0
! Make private
  Private
! Public variables  
  Public :: convTableSize, convTable, convTableN
! Public
  Public :: UnitConvert
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  Function UnitConvert(inputValue, inputUnitArg, outputUnitArg) RESULT (outputValue)
! Converts from one unit to another
    Implicit None  ! Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal) :: inputValue
    Character(*) :: inputUnitArg
    Character(*) :: outputUnitArg
! Out
    Real(kind=DoubleReal) :: outputValue
! Private
    Character(len=8) :: inputUnit
    Character(len=8) :: outputUnit    
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: factor
! Init variables    
    inputUnit = "        "
    inputUnit = StrToUpper(inputUnitArg)
    outputUnit = "        "
    outputUnit = StrToUpper(outputUnitArg)
! Load units into table    
    If(convTableN.eq.0)Then
      Call LoadUnits()
    End If
! Convert
    factor = 1.0D0
    Do i=1,convTableSize 
      If(convTable(i)%unitName.eq."ZZZZZZZZ")Then
        Exit
      End If
      If(convTable(i)%unitName.eq.inputUnitArg)Then
        factor = factor * convTable(i)%factor    
      End If      
      If(convTable(i)%unitName.eq.outputUnitArg)Then
        factor = factor * (1.0D0/convTable(i)%factor)
      End If  
    End Do
    outputValue = factor * inputValue
  End Function UnitConvert
  
  Subroutine LoadUnits()
! Returns factorial of input
    Implicit None  ! Force declaration of all variables
! ------------------------------    
! Add units
!
! ------------------------------  
! Lengths
    Call AddUnit("M", 1.0D0)
    Call AddUnit("CM", 1.0D-2)
    Call AddUnit("MM", 1.0D-3)
    Call AddUnit("UM", 1.0D-6)
    Call AddUnit("NM", 1.0D-9)
    Call AddUnit("ANG", 1.0D-10,"ANGS")
    Call AddUnit("BOHR", 5.2917721D-11)
! Volumes
    Call AddUnit("M3", 1.0D0)
    Call AddUnit("CM3", 1.0D-6)
    Call AddUnit("MM3", 1.0D-9)
    Call AddUnit("UM3", 1.0D-18)
    Call AddUnit("NM3", 1.0D-27)
    Call AddUnit("ANG3", 1.0D-30,"ANGS3")
! Times
    Call AddUnit("S", 1.0D0)
    Call AddUnit("MIN", 60.0D0)
    Call AddUnit("HR", 3600.0D0)
! Energies
    Call AddUnit("J", 1.0D0)                ! J, kg m3 s-2, Nm
    Call AddUnit("RY", 2.179872172697D-18)
    Call AddUnit("EV", 1.602176568D-19)
! Stresses
    Call AddUnit("PA", 1.0D0)               ! Pa, N/m2, J/m3
    Call AddUnit("RYBOH3", 1.471050658D13,"RYBOHR3")
    Call AddUnit("GPA", 1.0D9)
    Call AddUnit("EVAN3", 1.602176568D11,"EVANG3")   ! eV/ang3
! Forces
    Call AddUnit("N", 1.0D0)                     ! N, J/m
    Call AddUnit("RYBOHR", 4.11936139295036D-8)
    Call AddUnit("EVANG", 1.602176568D-9)       
  End Subroutine LoadUnits
  
  Subroutine AddUnit(unitNameArg, siFactor, unitNameArgB_in, unitNameArgC_in)
! Returns factorial of input
    Implicit None  ! Force declaration of all variables
! Declare variables  
! In
    Character(*) :: unitNameArg
    Real(kind=DoubleReal) :: siFactor
    Character(*), Optional :: unitNameArgB_in
    Character(*), Optional :: unitNameArgC_in
! Private
    Character(Len=8) :: unitName
    Character(Len=8) :: unitNameB
    Character(Len=8) :: unitNameC
! Prepare name    
    unitName = "        "
    unitName = StrToUpper(unitNameArg)
! Optional arguments
    unitNameB = "        "
    If(Present(unitNameArgB_in))Then
      unitNameB = StrToUpper(unitNameArgB_in)
    End If
    unitNameC = "        "
    If(Present(unitNameArgC_in))Then
      unitNameC = StrToUpper(unitNameArgC_in)
    End If
! Increment counter
    convTableN = convTableN + 1
! Add to table    
    convTable(convTableN)%unitName = unitName
    convTable(convTableN)%factor = siFactor
! Add other names to table
    If(unitNameB.ne."        ")Then
! Increment counter
      convTableN = convTableN + 1
! Add to table    
      convTable(convTableN)%unitName = unitNameB
      convTable(convTableN)%factor = siFactor
    End If
    If(unitNameC.ne."        ")Then
! Increment counter
      convTableN = convTableN + 1
! Add to table    
      convTable(convTableN)%unitName = unitNameC
      convTable(convTableN)%factor = siFactor
    End If
  End Subroutine AddUnit
  
End Module units
