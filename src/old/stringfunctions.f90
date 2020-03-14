Module stringfunctions

! Setup Modules
  Use kinds

! force declaration of all variables
  Implicit None
  Private
  Public :: StrToUpper
  Public :: NumericOnly
  Public :: RemoveSpaces
  Public :: CorrectFilePath

  character( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  character( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  Contains
  Function StrToUpper (input) RESULT (output)
! -- Argument and result
    CHARACTER(*), INTENT(IN) :: input
    CHARACTER(LEN(input)) :: output
! -- Local variables
    INTEGER :: i, n
! -- Copy input string
    output = input
! -- Loop over string elements
    Do i = 1, LEN( output )
! -- Find location of letter in lower case constant string
      n = INDEX( LOWER_CASE, output( i:i ) )
! -- If current substring is a lower case letter, make it upper case
      If( n /= 0 ) output( i:i ) = UPPER_CASE( n:n )
      End Do
    End Function StrToUpper

    Function NumericOnly (input) RESULT (output)
! -- Argument and result
      CHARACTER(*), INTENT(IN) :: input
      CHARACTER(LEN(input)) :: outputTemp
      CHARACTER(LEN(input)) :: output
! -- Local variables
      INTEGER :: i,n
! -- Copy input string
      outputTemp = input
      Do i = 1, LEN( outputTemp )
        output( i:i ) = " "
      End Do
      n = 0
      Do i = 1, LEN( outputTemp )
        If(outputTemp( i:i ).eq.".".or.(iachar(outputTemp( i:i )).ge.48.and.iachar(outputTemp( i:i )).le.57))Then
          n = n + 1
          output( n:n ) = outputTemp( i:i )
        Else
          output( i:i ) = " "
        End If
      End Do
    End Function NumericOnly

    Function RemoveSpaces (input) RESULT (output)
      CHARACTER(*), INTENT(IN) :: input
      CHARACTER(LEN(input)) :: outputTemp
      CHARACTER(LEN(input)) :: output
! -- Local variables
      INTEGER :: i,j
! -- Copy input string
      outputTemp = input
! Blank output
      Do i = 1, LEN( outputTemp )
        output( i:i ) = " "
      End Do
! transfer outputtemp to output without spaces
      j = 0
      Do i = 1, LEN( outputTemp )
        If(outputTemp( i:i ).ne." ")Then
          j = j + 1
          output( j:j ) = outputTemp( i:i )
        End If
      End Do
    End Function RemoveSpaces

    Function CorrectFilePath (input) RESULT (output)
      CHARACTER(*), INTENT(IN) :: input
      CHARACTER(LEN(input)) :: outputTemp
      CHARACTER(LEN(input)) :: output
! -- Local variables
      INTEGER :: i,n
! -- Copy input string
      outputTemp = input
      Do i = 1, LEN( outputTemp )
        output( i:i ) = " "
      End Do
      n = 0
      Do i = 1, LEN( outputTemp )
        If((iachar(outputTemp( i:i )).ge.32.and.iachar(outputTemp( i:i )).le.126))Then
          If(outputTemp( i:i ).ne."?".and.outputTemp( i:i ).ne."*".and.&
            outputTemp( i:i ).ne."%".and.outputTemp( i:i ).ne."+")Then
            n = n + 1
            output( n:n ) = outputTemp( i:i )
          End If
        End If
      End Do
      output = trim (output)
    End Function CorrectFilePath

  End Module stringfunctions
