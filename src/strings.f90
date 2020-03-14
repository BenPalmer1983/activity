Module strings
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use logicalMod
! Force declaration of all variables
  Implicit None
! Public variables
  Integer(kind=LongInteger) :: randomLCG_n_strings=0
  Integer(kind=LongInteger) :: randomLCG_xn_strings
  Integer(kind=LongInteger) :: randomLCG_R_n_strings=0
  Integer(kind=LongInteger) :: randomLCG_R_xn_strings
! Make private
  Private
! Public
! ---- Variables
  Public :: randomLCG_n_strings
  Public :: randomLCG_xn_strings
  Public :: randomLCG_R_n_strings
  Public :: randomLCG_R_xn_strings
! ---- Functions
  Public :: StrToUpper      ! All chars to upper
  Public :: StrToLower      ! All chars to lower
  Public :: StrInStr        !
  Public :: StrReplace
  Public :: NumericOnly
  Public :: SingleSpaces
  Public :: RemoveSpaces
  Public :: RemoveSpacesQ
  Public :: TrimSpaces
  Public :: TrimStr
  Public :: TrimmedLength
  Public :: BlankString
  Public :: BlankStringArray
  Public :: WipeString
  Public :: WipeStringArray
  Public :: Spaces
  Public :: SpacesRight
  Public :: RemoveComments
  Public :: RemoveQuotes
  Public :: RemoveTrailing
  Public :: IntToStr
  Public :: DpToStr
  Public :: StrToInt
  Public :: StrToDp
  Public :: StrToBool
  Public :: RandName
  Public :: TempFileName
  Public :: CleanString
  Public :: TimeToHuman
  Public :: IfStringEmpty
  Public :: ConcatStr
  Public :: StrMatch
  Public :: TestBoolStr
! ---- Subroutines
  Public :: explode
  Public :: randCharacter
  Public :: TrimString
  Public :: StrCenter
  Public :: StrAlign
! Interfaces
  Interface BlankStringArray
    Module Procedure BlankString1DArray, BlankString2DArray
  End Interface BlankStringArray
  Interface WipeStringArray
    Module Procedure WipeString1DArray, BlankString2DArray
  End Interface WipeStringArray
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------

! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------

  Function StrToUpper (input) RESULT (output)
! Returns factorial of input
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(IN) :: input
    Character(LEN(input)) :: output
! Private
    Integer(kind=StandardInteger) :: i, n
! Loop through characters
    output = "   "
    Do i=1,Len(input)
      n = Iachar(input(i:i))
      If(n.ge.97.and.n.le.122)Then
        n = n - 32
      End If
      output(i:i) = char(n)
    End Do
  End Function StrToUpper
!---------------------------------------------------------------------------------------------------------------------------------------
  Function StrToLower (input) RESULT (output)
! Returns factorial of input
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(IN) :: input
    Character(Len(input)) :: output
! Private
    Integer(kind=StandardInteger) :: i, n
! Loop through characters
    Do i=1,Len(input)
      n = Iachar(input(i:i))
      If(n.ge.65.and.n.le.90)Then
        n = n + 32
      End If
      output(i:i) = char(n)
    End Do
  End Function StrToLower
!---------------------------------------------------------------------------------------------------------------------------------------
    Function StrInStr (haystack, needle, caseSensitiveIn) RESULT (inString)
! Replace string
      Implicit None  !Force declaration of all variables
! Vars:  In
      Character(*) :: haystack
      Character(*) :: needle
      Logical, Optional :: caseSensitiveIn
! Vars:  Out
      Logical :: inString
! Vars:  Private
      Integer(kind=StandardInteger) :: i
      Logical :: caseSensitive
      Character(Len(haystack)) :: haystackWorking
      Character(Len(needle)) :: needleWorking
! Init vars
      inString = .false.
      caseSensitive = .false.
! Optional
      If(Present(caseSensitiveIn))Then
        caseSensitive = caseSensitiveIn
      End If
! If needle longer than haystack
      If(Len(needle).gt.Len(haystack))Then
        inString = .false.
      Else
! Both to uppercase
        If(caseSensitive)Then
          haystackWorking = haystack
          needleWorking = needle
        Else
          haystackWorking = StrToUpper(haystack)
          needleWorking = StrToUpper(needle)
        End If
! Loop through characters
        Do i=1,Len(haystackWorking)-len(needleWorking)+1
          If(haystackWorking(i:(i+len(needleWorking)-1)).eq.needleWorking)Then
            inString = .true.
          End If
          If(inString)Then
            Exit
          End If
        End Do
      End If
    End Function StrInStr
!---------------------------------------------------------------------------------------------------------------------------------------
  Function StrReplace (input, needle, replace) RESULT (output)
! Replace string
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*) :: input
    Character(*) :: needle
    Character(*) :: replace
! Out
    Character(LEN(input)) :: output
! Private
    Integer(kind=StandardInteger) :: n, i, j, k
! Loop through characters
    j = 0
    n = 0
    Do i=1,Len(input)-len(needle)+1
      n = n + 1
      If((n+len(needle)-1).gt.Len(input))Then
        Exit
      End If
      If(input(n:(n+len(needle)-1)).eq.needle)Then
        n = n + len(needle)-1
        Do k=1,len(replace)
          j = j + 1
          output(j:j) = replace(k:k)
        End Do
      Else
        j = j + 1
        output(j:j) = input(n:n)
      End If
    End Do
  End Function StrReplace
!---------------------------------------------------------------------------------------------------------------------------------------
  Function NumericOnly (input) RESULT (output)
! Returns factorial of input
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(IN) :: input
    Character(Len(input)) :: outputTemp
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i, n
! Remove characters
    outputTemp = input
    Do i = 1, Len( outputTemp )
      output( i:i ) = " "
    End Do
    n = 0
    Do i=1,Len(outputTemp)
      If(outputTemp( i:i ).eq.".".or.(iachar(outputTemp( i:i )).ge.48.and.iachar(outputTemp( i:i )).le.57))Then
        n = n + 1
        output( n:n ) = outputTemp( i:i )
      End If
    End Do
  End Function NumericOnly
!---------------------------------------------------------------------------------------------------------------------------------------
  Function SingleSpaces (input) RESULT (output)
! Remove spaces
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i, j
    Logical :: writeChar
! Blank output
    output = BlankString(output)
! transfer outputtemp to output without spaces
    j = 0
    Do i=1,Len(input)
      writeChar = .true.
      If(i.lt.Len(input))Then
        If(input(i:i).eq." ".and.input(i+1:i+1).eq." ")Then
          writeChar = .false.
        End If
      End If
      If(writeChar)Then
        j = j + 1
        output(j:j) = input(i:i)
      End If
    End Do
  End Function SingleSpaces
!---------------------------------------------------------------------------------------------------------------------------------------
  Function RemoveSpaces (input) RESULT (output)
! Remove spaces
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(IN) :: input
    Character(Len(input)) :: outputTemp
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i, j
! Remove spaces
    outputTemp = input
! Blank output
    Do i=1,Len(outputTemp)
      output(i:i) = " "
    End Do
! transfer outputtemp to output without spaces
    j = 0
    Do i=1,Len(outputTemp)
      If(outputTemp( i:i ).ne." ")Then
        j = j + 1
        output(j:j) = outputTemp(i:i)
      End If
    End Do
  End Function RemoveSpaces
!---------------------------------------------------------------------------------------------------------------------------------------
  Function RemoveSpacesQ (input) RESULT (output)
! Removes spaces, not within quotation marks
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(IN) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i, j, inLen
    Logical :: inQuotes
! Remove spaces
    output = BlankString(output)
    inQuotes = .false.
    inLen = len(input)
    j = 0
    Do i=1,inLen
      If(ichar(input(i:i)).eq.34)Then
        j = j + 1
        inQuotes = FlipLogical(inQuotes)
        output(j:j) = input(i:i)
      Else
        If(input(i:i).eq." ")Then
          If(inQuotes)Then
            j = j + 1
            output(j:j) = input(i:i)
          Else
            ! skip
          End If
        Else
          j = j + 1
          output(j:j) = input(i:i)
        End If
      End If
    End Do
  End Function RemoveSpacesQ
!---------------------------------------------------------------------------------------------------------------------------------------
  Function TrimSpaces(trimStr, padCharIn) Result (workingStr)
! Trims
    Implicit None  ! Force declaration of all variables
! In/Out:      Declare variables
    Character(*) :: trimStr
    Character(Len=1), Optional :: padCharIn
! Private:     Declare variables
    Character(Len(trimStr)) :: workingStr
    Integer(kind=StandardInteger) :: i, j, k, inputLen
    Logical :: store
    Character(Len=1) :: padChar
! Optional Argument
    padChar(1:1) = ""
    If(Present(padCharIn))Then
      padChar(1:1) = padCharIn(1:1)
    End If
! trim
    store = .false.
    j = 0
    inputLen = len(trimStr)
    Do i=1,inputLen
      If(trimStr(i:i).ne." ")Then
        store = .true.
      End If
      If(store)Then
        j = j + 1
        workingStr(j:j) = trimStr(i:i)
      End If
    End Do
    j = j + 1
    If(j.lt.inputLen)Then
      Do k=j,inputLen
        workingStr(k:k) = " "
      End Do
    End If
    Do i=1,inputLen
      j = i-1
      k = inputLen-j
      If(workingStr(k:k).eq." ")Then
        workingStr(k:k) = padChar(1:1)
      Else
        Exit
      End If
    End Do
  End Function TrimSpaces
!---------------------------------------------------------------------------------------------------------------------------------------
  Function TrimStr(strIn) Result (strOut)
! Trims both ends of the string
    Implicit None  ! Force declaration of all variables
! Vars:  In
    Character(*) :: strIn
! Vars:  Out
    Character(Len(strIn)) :: strOut
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, flag
! Init
    flag = 0
    j = 0
    strOut = WipeString(strOut)
! Trim leading spaces
    Do i=1,len(strIn)
      If(flag.eq.0)Then
        If(strIn(i:i).eq.char(0))Then
! Do nothing
        ElseIf(strIn(i:i).eq.char(32))Then
! Do nothing
        Else
          flag = flag + 1
        End If
      End If
      If(flag.eq.1)Then
        j= j + 1
        strOut(j:j) = strIn(i:i)
      End If
    End Do
! Trim trailing spaces
    j = len(strOut)+1
    Do i=1,len(strOut)
      j = j - 1
      If((iChar(strOut(j:j)).eq.0))Then
        strOut(j:j) = char(0)
      ElseIf((iChar(strOut(j:j)).eq.32))Then
        strOut(j:j) = char(0)
      Else
        Exit
      End If
    End Do
  End Function TrimStr
!---------------------------------------------------------------------------------------------------------------------------------------
  Function TrimmedLength(strIn) Result (strLen)
! Trims both ends of the string
    Implicit None  ! Force declaration of all variables
! Vars:  In
    Character(*) :: strIn
! Vars:  Out
    Integer(kind=StandardInteger) :: strLen
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
! Init
    strLen = Len(strIn)
! Subtract leading spaces/blanks
    Do i=1,len(strIn)
      If(strIn(i:i).eq.char(0))Then
        strLen = strLen - 1
      ElseIf(strIn(i:i).eq.char(32))Then
        strLen = strLen - 1
      Else
        Exit
      End If
    End Do
! Subtract trailing spaces
    j = len(strIn)+1
    Do i=1,len(strIn)
      j = j - 1
      If((iChar(strIn(j:j)).eq.0))Then
        strLen = strLen - 1
      ElseIf((iChar(strIn(j:j)).eq.32))Then
        strLen = strLen - 1
      Else
        Exit
      End If
    End Do
  End Function TrimmedLength
!---------------------------------------------------------------------------------------------------------------------------------------
  Function BlankString (input) RESULT (output)
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i
    Do i=1,Len(input)
      output(i:i) = " "
    End Do
  End Function BlankString
!---------------------------------------------------------------------------------------------------------------------------------------
  Function BlankString1DArray (input) RESULT (output)
    Character(*), Dimension(:), INTENT(IN) :: input
    Character(Len(input)) :: line
    Character(Len(input)), Dimension(1:size(input,1)) :: output
    Integer(kind=StandardInteger) :: i
    Do i=1,Len(input)
      line(i:i) = " "
    End Do
    Do i=1,size(input,1)
      output(i) = line
    End Do
  End Function BlankString1DArray
!---------------------------------------------------------------------------------------------------------------------------------------
  Function BlankString2DArray (input) RESULT (output)
    Character(*), Dimension(:,:), INTENT(IN) :: input
    Character(Len(input)) :: line
    Character(Len(input)), Dimension(1:size(input,1),1:size(input,2)) :: output
    Integer(kind=StandardInteger) :: i, j
    Do i=1,Len(input)
      line(i:i) = " "
    End Do
    Do i=1,size(input,1)
      Do j=1,size(input,2)
        output(i,j) = line
      End Do
    End Do
  End Function BlankString2DArray

  Function WipeString (input) RESULT (output)
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i
    Do i=1,Len(input)
      output(i:i) = char(0)
    End Do
  End Function WipeString

  Function WipeString1DArray (input) RESULT (output)
    Character(*), Dimension(:), INTENT(IN) :: input
    Character(Len(input)) :: line
    Character(Len(input)), Dimension(1:size(input,1)) :: output
    Integer(kind=StandardInteger) :: i
    Do i=1,Len(input)
      line(i:i) = char(0)
    End Do
    Do i=1,size(input,1)
      output(i) = line
    End Do
  End Function WipeString1DArray

  Function WipeString2DArray (input) RESULT (output)
    Character(*), Dimension(:,:), INTENT(IN) :: input
    Character(Len(input)) :: line
    Character(Len(input)), Dimension(1:size(input,1),1:size(input,2)) :: output
    Integer(kind=StandardInteger) :: i, j
    Do i=1,Len(input)
      line(i:i) = char(0)
    End Do
    Do i=1,size(input,1)
      Do j=1,size(input,2)
        output(i,j) = line
      End Do
    End Do
  End Function WipeString2DArray





  Function Spaces (length) RESULT (output)
    Integer(kind=StandardInteger), INTENT(IN) :: length
    Character(Len=length) :: output
    Integer(kind=StandardInteger) :: i
    Do i=1,length
      output(i:i) = " "
    End Do
  End Function Spaces

  Function SpacesRight (input) RESULT (output)
! Adds spaces to right of string
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: tempStr, output
    Integer(kind=StandardInteger) :: i
    tempStr = trim(adjustl(input))
    Do i=1,Len(tempStr)
      output(i:i) = " "
    End Do
    Do i=1,Len(trim(tempStr))
      output(i:i) = tempStr(i:i)
    End Do
  End Function SpacesRight

  Function RemoveComments (input) RESULT (output)
! Removes comments from
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i
    output = BlankString(output)
! Comment character !
    Do i=1,Len(input)
      If(input(i:i).eq."!")Then
        Exit
      Else
        output(i:i) = input(i:i)
      End If
    End Do
  End Function RemoveComments

  Function RemoveQuotes (input) RESULT (output)
! Removes comments from
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i, j
    output = BlankString(output)
! Comment character !
    j = 0
    Do i=1,Len(input)
      If(ichar(input(i:i)).eq.34.or.ichar(input(i:i)).eq.39)Then
! Do nothing
      Else
        j = j + 1
        output(j:j) = input(i:i)
      End If
    End Do
  End Function RemoveQuotes

  Function RemoveTrailing (input, trailingIn) RESULT (output)
! Returns TRUE if A larger than B
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Character(*), INTENT(IN) :: input
    Character(Len=1), Optional :: trailingIn
! Out:     Declare variables
    Character(Len(input)) :: output
! Private
    Character(Len=1) :: trailing
    Integer(kind=StandardInteger) :: i, j, inLen, trailingFlag
! Optional
    trailing = "."
    If(Present(trailingIn))Then
      trailing = trailingIn
    End If
! Read backwards
    inLen = Len(input)
    j = inLen
    trailingFlag = 1
    Do i=1,inLen
      If(trailingFlag.eq.1.and.ichar(input(j:j)).ne.32)Then
        trailingFlag = 2
      End If
      If(trailingFlag.eq.2)Then
        If(input(j:j).eq.trailing(1:1))Then
          output(j:j) = " "
        Else
          output(j:j) = input(j:j)
        End If
        trailingFlag = 3
        j = j - 1
      Else
        output(j:j) = input(j:j)
        j = j - 1
      End If
    End Do
  End Function RemoveTrailing

  Function IntToStr (input) RESULT (output)
! Apply style to last dataset added (unless otherwise specified)
    Implicit None  ! Force declaration of all variables
  ! In:      Declare variables
    Integer(kind=StandardInteger) :: input
  ! Out:     Declare variables
    Character(Len=16) :: output
    Write(output,"(I16)") input
    output = trim(adjustl(output))
  End Function IntToStr

  Function DpToStr (input, numFormatIn) RESULT (output)
! Apply style to last dataset added (unless otherwise specified)
    Implicit None  ! Force declaration of all variables
  ! In:      Declare variables
    Real(kind=DoubleReal) :: input
    Character(*), Optional :: numFormatIn
  ! Out:     Declare variables
    Character(Len=16) :: output
  ! Private: Declare variables
    Character(Len=12) :: numFormat
    numFormat = "(E16.8)"
    If(Present(numFormatIn))Then
      numFormat = numFormatIn
    End If
    Write(output,numFormat) input
    output = trim(adjustl(output))
  End Function DpToStr

  Function StrToInt (input) RESULT (output)
! Apply style to last dataset added (unless otherwise specified)
    Implicit None  ! Force declaration of all variables
  ! In:      Declare variables
    Character(*) :: input
  ! Out:     Declare variables
    Integer(kind=StandardInteger) :: output
    output = 0
    Read(input,*) output
  End Function StrToInt

  Function StrToDp (input) RESULT (outputDouble)
! Apply style to last dataset added (unless otherwise specified)
    Implicit None  ! Force declaration of all variables
! Vars:  In
    Character(*) :: input
! Vars:  Out
    Real(kind=DoubleReal) :: outputDouble
! Vars:  Private
    Character(Len=32) :: inputRewritten
    Integer(kind=StandardInteger) :: i, j
! Check format
    input = StrToUpper(input)
    input = Trim(adjustl(input))
    input = StrReplace(input, "E", "D")
    input = StrReplace(input, "D+", "D")
! rewrite
    inputRewritten = BlankString(inputRewritten)
    j = 0
    Do i=1,len(input)
      If(CharCheckReal(input(i:i)))Then
        j = j + 1
        inputRewritten(j:j) = input(i:i)
      End If
    End Do
! Read into var
    Read(inputRewritten,*) outputDouble
  End Function StrToDp

  Function StrToBool (inputIn) RESULT (output)
! converts a user input into true/false bool
    Implicit None  ! Force declaration of all variables
  ! In:      Declare variables
    Character(*) :: inputIn
  ! Out:     Declare variables
    Logical :: output
  ! Private: Declare variables
    Character(len=6) :: input
    input = inputIn
    output = .false.
    input = StrToUpper(Trim(Adjustl(input)))
    If(input(1:1).eq."1")Then
      output = .true.
    End If
    If(input(1:1).eq."+")Then
      output = .true.
    End If
    If(input(1:1).eq."Y")Then
      output = .true.
    End If
    If(input(1:4).eq."TRUE")Then
      output = .true.
    End If
    If(input(1:6).eq.".TRUE.")Then
      output = .true.
    End If
  End Function StrToBool

  Function RandName(randSwitchIn, prefixIn) Result (randNameOut)
! Make a random 8 character "name"
    Implicit None  ! Force declaration of all variables
! In:      Declare variables
    Logical, Optional :: randSwitchIn
    Character(*), Optional :: prefixIn
! Out:     Declare variables
    Character(len=8) :: randNameOut
! Private: Declare variables
    Logical :: randSwitch
    Character(len=8) :: prefix
    Integer(kind=StandardInteger) :: i
    Character(len=1) :: randChar
    Logical :: writePrefix = .true.
! Optional
    randSwitch = .true.
    prefix = "        "
    If(Present(randSwitchIn))Then
      randSwitch = randSwitchIn
    End If
    If(Present(prefixIn))Then
      prefix = prefixIn
    End If
! Init output
    randNameOut = "        "
! Loop through letters
    Do i=1,8
      If(prefix(i:i).eq." ")Then
        writePrefix = .false.
      End If
      If(writePrefix)Then
        randNameOut(i:i) = prefix(i:i)
      Else
        Call randCharacter(randChar, randSwitchIn, 2)
        randNameOut(i:i) = randChar
      End If
    End Do
  End Function RandName

  Function TempFileName(randSwitchIn) Result (fileNameOut)
! Make random name for temp file
    Implicit None  ! Force declaration of all variables
! In:      Declare variables
    Logical, Optional :: randSwitchIn
! Out:     Declare variables
    Character(len=8) :: fileNameOut
! Private: Declare variables
    Logical :: randSwitch
! Optional
    randSwitch = .true.
    If(Present(randSwitchIn))Then
      randSwitch = randSwitchIn
    End If
! Init
    fileNameOut = "        "
! Make name
    fileNameOut = RandName(randSwitch, "tmp")
  End Function TempFileName

  Function CleanString(stringIn) Result (stringOut)
! Clean a string - a-z A-Z 0-9 space
    Implicit None  ! Force declaration of all variables
! In:      Declare variables
    Character(*) :: stringIn
! Out:     Declare variables
    Character(Len(stringIn)) :: stringOut
! Private: Declare variables
    Integer(kind=StandardInteger) :: i, n, charVal
! Blank output string
    stringOut = BlankString(stringOut)
    n = 0
    Do i=1,len(stringIn)
      charVal = ichar(stringIn(i:i))
      If(charVal.eq.32.or.&
      (charVal.ge.48.and.charVal.le.57).or.&
      (charVal.ge.65.and.charVal.le.90).or.&
      (charVal.ge.97.and.charVal.le.122))Then
        n = n + 1
        stringOut(n:n) = stringIn(i:i)
      End If
    End Do
  End Function CleanString


  Function TimeToHuman(timeIn) Result (stringOut)
! Converst time in seconds to more human friendly version
    Implicit None  ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal) :: timeIn
! Out:     Declare variables
    Character(Len=128) :: stringOut
! Private: Declare variables
    Real(kind=DoubleReal) :: timeSeconds
    Integer(kind=StandardInteger) :: years=0
    Integer(kind=StandardInteger) :: days=0
    Integer(kind=StandardInteger) :: hours=0
    Integer(kind=StandardInteger) :: minutes=0
    Real(kind=DoubleReal) :: seconds=0.0D0
!
    timeSeconds = timeIn
    If(timeSeconds.gt.31557600.0D0)Then
      years = Floor(timeSeconds/31557600.0D0)
      timeSeconds = timeSeconds-31557600.0D0*years
    End If
    If(timeSeconds.gt.86400.0D0)Then
      days = Floor(timeSeconds/86400.0D0)
      timeSeconds = timeSeconds-86400.0D0*days
    End If
    If(timeSeconds.gt.3600.0D0)Then
      hours = Floor(timeSeconds/3600.0D0)
      timeSeconds = timeSeconds-3600.0D0*hours
    End If
    If(timeSeconds.gt.60.0D0)Then
      minutes = Floor(timeSeconds/60.0D0)
      timeSeconds = timeSeconds-60.0D0*minutes
    End If
    seconds = timeSeconds
    stringOut = BlankString(stringOut)
    If(years.gt.0)Then
      stringOut = trim(stringOut)//" "//adjustl(trim(IntToStr(years)))//" yr "
    End If
    If(days.gt.0)Then
      stringOut = trim(stringOut)//" "//adjustl(trim(IntToStr(days)))//" d "
    End If
    If(hours.gt.0)Then
      stringOut = trim(stringOut)//" "//adjustl(trim(IntToStr(hours)))//" hr "
    End If
    If(minutes.gt.0)Then
      stringOut = trim(stringOut)//" "//adjustl(trim(IntToStr(minutes)))//" min "
    End If
    stringOut = trim(stringOut)//" "//adjustl(trim(DpToStr(seconds,"(F10.3)")))//" s "
    stringOut = trim(adjustl(stringOut))
  End Function TimeToHuman

  Function IfStringEmpty(stringIn) Result (result)
! Check for an empty string
! Spaces (iChar = 32) and iChar = 0 are empty
    Implicit None  ! Force declaration of all variables
! In:      Declare variables
    Character(*) :: stringIn
! Out:     Declare variables
    Logical :: result
    Integer(kind=StandardInteger) :: i, charAscii
! Assume empty
    result = .true.
! Test string
    Do i=1,Len(stringIn)
      charAscii = iChar(stringIn(i:i))
      If(charAscii.ne.0)Then
        If(charAscii.ne.32)Then
          result = .false.
          Exit
        End If
      End If
    End Do
  End Function IfStringEmpty

  Function ConcatStr(stringA, stringB, forceTrimIn) Result (stringC)
! Join two strings, returning in first string length
    Implicit None  ! Force declaration of all variables
! Vars:  In
    Character(*) :: stringA
    Character(*) :: stringB
    Logical, Optional :: forceTrimIn
! Vars:  Out
    Character(Len(stringA)) :: stringC
! Vars:  Private
    Logical :: forceTrim
    Integer(kind=StandardInteger) :: i, j, flag
! Init
    forceTrim = .true.
    stringC = WipeString(stringC)
! Optional arguments
    If(Present(forceTrimIn))Then
      forceTrim = forceTrimIn
    End If
! Trim strings if required
    If(forceTrim)Then
      stringA = TrimStr(stringA)
      stringB = TrimStr(stringB)
    End If
    flag = 1
    j = 0
    Do i=1,Len(stringA)
      j = j + 1
      If(flag.eq.1)Then
        If(ichar(stringA(j:j)).eq.0)Then
          flag = 2
          j = 1
        End If
      End If
      If(flag.eq.2)Then
        If(j.gt.len(stringB))Then
          Exit
        End If
        If(ichar(stringB(j:j)).eq.0)Then
          Exit
        End If
      End If
      If(flag.eq.1)Then
        stringC(i:i) = stringA(j:j)
      End If
      If(flag.eq.2)Then
        stringC(i:i) = stringB(j:j)
      End If
    End Do
  End Function ConcatStr

  Function StrMatch(strA_In,strB_In,caseSensitiveIn) Result (matchResult)
! Trims spaces and zeros
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Character(*) :: strA_In
    Character(*) :: strB_In
    Logical, Optional :: caseSensitiveIn
! Vars:  Out
    Logical :: matchResult
! Vars:  Private
    Logical :: caseSensitive
    Integer(kind=StandardInteger) :: i
    Character(Len(strA_In)) :: strA
    Character(Len(strB_In)) :: strB
! Optional Arguments
    caseSensitive = .true.
    If(Present(caseSensitiveIn))Then
      caseSensitive = caseSensitiveIn
    End if
! Init
    matchResult = .true.
! Adjust strings
    If(caseSensitive)Then
      strA = strA_In
      strB = strB_In
    Else
      strA = StrToUpper(strA_In)
      strB = StrToUpper(strB_In)
    End If
! Trim
    strA = TrimStr(strA)
    strB = TrimStr(strB)
! Check for match
    Do i=1,len(strA)
      If(i.gt.len(strB))Then
        If(iChar(strA(i:i)).ne.0)Then
          matchResult = .false.
        End If
        Exit
      End If
      If(strA(i:i).ne.strB(i:i))Then
        matchResult = .false.
        Exit
      End If
      If(ichar(strA(i:i)).eq.0)Then
        If(iChar(strB(i:i)).eq.0)Then
          Exit
        Else
          matchResult = .false.
          Exit
        End If
      End If
    End Do
  End Function StrMatch

  Function TestBoolStr(str) Result (result)
! Trims spaces and zeros
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Character(*) :: str
! Vars:  Out
    Logical :: result
! Trim
    str = TrimStr(str)
! Convert to UC
    str = StrToUpper(str)
! Init
    result = .false.
! Test
    If(str(1:1).eq."Y")Then
      result = .true.
    End If
    If(str(1:1).eq."T")Then
      result = .true.
    End If
    If(str(1:2).eq.".T")Then
      result = .true.
    End If
  End Function TestBoolStr

! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------

  Subroutine explode(inputString, fieldSplit, outputArray, outputCount)
! In/Out:  Declare variables
    Character(*) :: inputString
    Character(*) :: fieldSplit
    Character(*), Dimension(:) :: outputArray
    Integer(kind=StandardInteger) :: outputCount
! Private: Declare variables
    Character(Len(fieldSplit)) :: trialSegment
    Integer(kind=StandardInteger) :: fieldCount
    Integer(kind=StandardInteger) :: lenInput
    Integer(kind=StandardInteger) :: lenSplit
    Integer(kind=StandardInteger) :: i, charI, n, k
! Init
    outputCount = 0
    Call TrimString(inputString,lenInput," ")
    !lenInput = len(inputString)
    lenSplit = len(fieldSplit)
    If(lenInput.gt.lenSplit)Then
      n = 0
      fieldCount = 1
      charI = 0
      Do i=1,lenInput
        charI = charI + 1
        trialSegment = inputString(charI:(charI+lenSplit-1))
        If(trialSegment.eq.fieldSplit)Then
          Do k=n+1,len(outputArray)
            outputArray(fieldCount)(k:k) = " "
          End Do
          fieldCount = fieldCount + 1
          n = 0
          charI = charI+lenSplit-1
        Else
          n = n + 1
          outputArray(fieldCount)(n:n) = inputString(charI:charI)
        End If
        If(charI.ge.lenInput)Then
          Exit
        End If
      End Do
! process last field
      Do k=n+1,len(outputArray)
        outputArray(fieldCount)(k:k) = " "
      End Do
! store field count
      outputCount = fieldCount
    End If
  End Subroutine explode

! Character Functions

  Subroutine randCharacter(letter, randSwitchIn, setIn)
! In/Out:      Declare variables
    Character(len=1) :: letter
    Logical, Optional :: randSwitchIn
    Integer(kind=StandardInteger), Optional :: setIn
! Private: Declare variables
    Logical :: randSwitch
    Integer(kind=StandardInteger) :: set
    Integer(kind=StandardInteger) :: characterNum
    Real(kind=DoubleReal) :: randNumber
    Character(len=52), Parameter :: alpha = &
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(len=62), Parameter :: alphaNum = &
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890'
! Optional arguments
    randSwitch = .false.   ! false use repeatable rng, true use random based off cpu clock
    set = 1                ! 1 = alpha, 2 = alpha numeric
    If(Present(randSwitchIn))Then
      randSwitch = randSwitchIn
    End If
    If(Present(setIn))Then
      set = setIn
    End If
    letter = " "
! Make character
    If(randSwitch)Then
      randNumber = RandomLCG_R_strings()
    Else
      randNumber = RandomLCG_strings()
    End If
    If(set.eq.1)Then
      characterNum = Ceiling(62.0D0*randNumber+1.0D0)
      If(characterNum.lt.1)Then
        characterNum = 1
      End If
      If(characterNum.gt.52)Then
        characterNum = 52
      End If
      letter = alpha(characterNum:characterNum)
    End If
    If(set.eq.2)Then
      characterNum = Ceiling(52.0D0*randNumber+1.0D0)
      If(characterNum.lt.1)Then
        characterNum = 1
      End If
      If(characterNum.gt.62)Then
        characterNum = 62
      End If
      letter = alphaNum(characterNum:characterNum)
    End If
  End Subroutine randCharacter

  Function CharCheckAlpha(charIn) Result (boolResult)
! Check if character is a letter
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Character(Len=1), Optional :: charIn
! Vars:  Out
    Logical :: boolResult
! Init
    boolResult = .false.
! Check if alpha Upper Case
    If(iChar(charIn).ge.65.and.iChar(charIn).le.90)Then
      boolResult = .true.
    End If
! Check if alpha Lower Case
    If(iChar(charIn).ge.97.and.iChar(charIn).le.122)Then
      boolResult = .true.
    End If
  End Function CharCheckAlpha

  Function CharCheckNumeric(charIn) Result (boolResult)
! Check if character is a letter
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Character(Len=1), Optional :: charIn
! Vars:  Out
    Logical :: boolResult
! Init
    boolResult = .false.
! Check if numeric
    If(iChar(charIn).ge.48.and.iChar(charIn).le.57)Then
      boolResult = .true.
    End If
  End Function CharCheckNumeric

  Function CharCheckAlphaNumeric(charIn) Result (boolResult)
! Check if character is a letter
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Character(Len=1), Optional :: charIn
! Vars:  Out
    Logical :: boolResult
! Init
    boolResult = .false.
    boolResult = CharCheckAlpha(charIn)
    If(.not.boolResult)Then
      boolResult = CharCheckNumeric(charIn)
    End If
  End Function CharCheckAlphaNumeric

  Function CharCheckReal(charIn) Result (boolResult)
  ! Check if character is a letter
    Implicit None   ! Force declaration of all variables
  ! Vars:  In
    Character(Len=1), Optional :: charIn
  ! Vars:  Out
    Logical :: boolResult
  ! Init
    boolResult = .false.
    If(charIn.eq.".")Then
      boolResult = .true.
    End If
    If(charIn.eq."-")Then
      boolResult = .true.
    End If
    If(charIn.eq."+")Then
      boolResult = .true.
    End If
    If(charIn.eq."D")Then
      boolResult = .true.
    End If
    If(.not.boolResult)Then
      boolResult = CharCheckNumeric(charIn)
    End If
  End Function CharCheckReal


  Subroutine TrimString(trimStr, outputLength, padCharIn)
! Trims
    Implicit None  ! Force declaration of all variables
! In/Out:      Declare variables
    Character(*) :: trimStr
    Integer(kind=StandardInteger) :: outputLength
    Character(Len=1), Optional :: padCharIn
! Private:     Declare variables
    Character(Len(trimStr)) :: workingStr
    Integer(kind=StandardInteger) :: i, j, k, inputLen
    Logical :: store
    Character(Len=1) :: padChar
! Optional Argument
    padChar(1:1) = ""
    If(Present(padCharIn))Then
      padChar(1:1) = padCharIn(1:1)
    End If
! trim
    store = .false.
    j = 0
    inputLen = len(trimStr)
    Do i=1,inputLen
      If(trimStr(i:i).ne." ")Then
        store = .true.
      End If
      If(store)Then
        j = j + 1
        workingStr(j:j) = trimStr(i:i)
      End If
    End Do
    j = j + 1
    If(j.lt.inputLen)Then
      Do k=j,inputLen
        workingStr(k:k) = " "
      End Do
    End If
    outputLength = 0
    Do i=1,inputLen
      j = i-1
      k = inputLen-j
      If(workingStr(k:k).eq." ")Then
        workingStr(k:k) = padChar(1:1)
      Else
        Exit
      End If
    End Do
    outputLength = inputLen-j
    Do i=1,inputLen
      trimStr(i:i) = workingStr(i:i)
    End Do
  End Subroutine TrimString

!----------

  Subroutine StrCenter(line, lineLengthIn)
! Trims
    Implicit None  ! Force declaration of all variables
! Vars:  In/Out
    Character(*) :: line
    Integer(kind=StandardInteger), Optional :: lineLengthIn
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Integer(kind=StandardInteger) :: lineLength
    Character(Len(line)) :: lineTemp
    Integer(kind=StandardInteger) :: trimStrLength, paddingL, paddingR
! Optional arguments
    lineLength = Len(line)
    If(Present(lineLengthIn))Then
      lineLength = lineLengthIn
    End If
! Init vars
    lineTemp = wipeString(lineTemp)
! Trim string
    line = TrimStr(line)
! Padding values
    trimStrLength = TrimmedLength(line)
    If(trimStrLength.lt.lineLength)Then
      paddingL = floor((lineLength-trimStrLength)/2.0D0)
      paddingR = (lineLength-trimStrLength)-paddingL
      j=0
      Do i=1,paddingL
        j = j + 1
        lineTemp(j:j) = " "
      End Do
      Do i=1,trimStrLength
        j = j + 1
        lineTemp(j:j) = line(i:i)
      End Do
      Do i=1,paddingR
        j = j + 1
        lineTemp(j:j) = " "
      End Do
    Else
      lineTemp(1:lineLength) = line(1:lineLength)
    End If
    line = lineTemp
  End Subroutine StrCenter

  Subroutine StrAlign(line, align, lineLengthIn)
! Trims
    Implicit None  ! Force declaration of all variables
! Vars:  In/Out
    Character(*) :: line
    Character(Len=1) :: align
    Integer(kind=StandardInteger), Optional :: lineLengthIn
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Integer(kind=StandardInteger) :: lineLength
    Character(Len(line)) :: lineTemp
    Integer(kind=StandardInteger) :: trimStrLength, paddingL, paddingR
! Optional arguments
    lineLength = Len(line)
    If(Present(lineLengthIn))Then
      lineLength = lineLengthIn
    End If
! Init vars
    lineTemp = wipeString(lineTemp)
! Trim string
    line = TrimStr(line)
! Padding values
    trimStrLength = TrimmedLength(line)
    If(trimStrLength.lt.lineLength)Then
      If(align.eq."C")Then
        paddingL = floor((lineLength-trimStrLength)/2.0D0)
        paddingR = (lineLength-trimStrLength)-paddingL
      ElseIf(align.eq."R")Then
        paddingL = (lineLength-trimStrLength)
        paddingR = 0
      Else
        paddingL = 0
        paddingR = (lineLength-trimStrLength)
      End If
      j=0
      Do i=1,paddingL
        j = j + 1
        lineTemp(j:j) = " "
      End Do
      Do i=1,trimStrLength
        j = j + 1
        lineTemp(j:j) = line(i:i)
      End Do
      Do i=1,paddingR
        j = j + 1
        lineTemp(j:j) = " "
      End Do
    Else
      lineTemp(1:lineLength) = line(1:lineLength)
    End If
    line = lineTemp
  End Subroutine StrAlign







!---------------------------------------------------------------------------------------------------------------------------------------
! String random number functions (as these are loaded AFTER strings MOD with the rng MOD
!---------------------------------------------------------------------------------------------------------------------------------------

  Function RandomLCG_strings(seedIn) RESULT (output)
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
      randomLCG_n_strings = 0
    End If
! If first iteration
    If(randomLCG_n_strings.eq.0)Then
      If(seed.eq.0)Then
        seed = 12791244+45778951 ! Use default seed
      End If
      randomLCG_n_strings = 0
      randomLCG_xn_strings = seed
    End If
! Increment counter
    randomLCG_n_strings = randomLCG_n_strings + 1
! calculate
    randomLCG_xn_strings = mod((a*randomLCG_xn_strings+c),m)
    output = (1.0D0*randomLCG_xn_strings)/(1.0D0*m)
  End Function RandomLCG_strings

  Function RandomLCG_R_strings() RESULT (output)
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
    If(randomLCG_R_n_strings.eq.0)Then
! Make "random" seed
      Call SYSTEM_CLOCK(clockTime) ! "nano seconds" - well, an estimate
      seed = mod(clockTime,m)
      Do i=1,10
        seed = mod((a*seed+c),m)
      End Do
      randomLCG_R_n_strings = 0
      randomLCG_R_xn_strings = seed
    End If
! Increment counter
    randomLCG_R_n_strings = randomLCG_R_n_strings + 1
! calculate
    randomLCG_R_xn_strings = mod((a*randomLCG_R_xn_strings+c),m)
    output = (1.0D0*randomLCG_R_xn_strings)/(1.0D0*m)
  End Function RandomLCG_R_strings

End Module strings