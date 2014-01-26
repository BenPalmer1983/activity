Module initialise

! Setup Modules
  Use kinds
  Use constants
  Use stringfunctions		!string functions
  Use maths


!force declaration of all variables
  Implicit None
  
!declare global variables  
  Real :: programStartTime
  Character(len=255) :: currentWorkingDirectory
  Character(len=255) :: outputFile

!Privacy of functions/subroutines/variables
  Private
  Public :: programStartTime		!Variable
  Public :: outputFile      		!Variable
  Public :: currentWorkingDirectory	!Variable
  Public :: runInitialise		    !Subroutine
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines

  Subroutine runInitialise()
	
	!Internal subroutine variables
	Integer :: i, j, k
	Integer, Dimension(1:3) :: theTime, theDate
	Character(len=255) :: outputFilePath
	
	!store start time
	CALL cpu_time(programStartTime)
	call idate(theDate)   ! theDate(1)=day, (2)=month, (3)=year
    call itime(theTime)   ! theDate(1)=hour, (2)=minute, (3)=second
	
	!get the working directory
	CALL getcwd(currentWorkingDirectory)
	
	!Create output file
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	!open(unit=999,file=trim(outputFile),status="new",action="write")
	open(unit=999,file=trim(outputFile))
	write(999,"(A38)") "======================================"
	write(999,"(A38)") "ACTIVITY code University of Birmingham"
	write(999,"(A38)") "======================================"
	write(999,"(A1)") " "
	write(999,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
	"Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)	
	close(999)

  End Subroutine runInitialise
  
  
  
  

End Module initialise