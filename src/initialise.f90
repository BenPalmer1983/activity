Module initialise

! Setup Modules
  Use kinds
  Use msubs
  Use general
  Use globals

! force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Privacy of functions/subroutines/variables
  Private
! Subroutines
  Public :: runInitialise        !Subroutine
! Functions
  Public :: ProgramTime             !Function

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

  contains

! Run all the input subroutines

  Subroutine runInitialise()
! Internal subroutine variables
    Integer, Dimension(1:3) :: theTime, theDate
    Integer(kind=StandardInteger) :: error
! store start time
    CALL cpu_time(programStartTime)
    call idate(theDate)   ! theDate(1)=day, (2)=month, (3)=year
    call itime(theTime)   ! theDate(1)=hour, (2)=minute, (3)=second
! get the working directory
    CALL getcwd(currentWorkingDirectory)
! Set output and temp/scratch directories
    outputDirectory = Trim(currentWorkingDirectory)//"/output"
    tempDirectory = Trim(currentWorkingDirectory)//"/temp"
! Make directories if not there
    Call makeDir(outputDirectory)
    Call makeDir(tempDirectory)
! MPI variables (public)
    Call MPI_Comm_size( MPI_COMM_WORLD ,mpiProcessCount,error)
    Call MPI_Comm_rank(MPI_COMM_WORLD,mpiProcessID,error)
! Set file paths
    fileOutputData = Trim(outputDirectory)//"/"//"output.dat"
    fileActivityHistory = Trim(outputDirectory)//"/"//"activityHistory.dat"
    fileInOut = Trim(outputDirectory)//"/"//"inOutFile.dat"
    fileIsotopeActivity = Trim(outputDirectory)//"/"//"isotopeActivityFile.dat"
    fileIsotopeActivityG = Trim(outputDirectory)//"/"//"isotopeActivityFileG.dat"
! Output program details
    If(mpiProcessID.eq.0)Then
      call idate(theDate)   ! theDate(1)=day, (2)=month, (3)=year
      call itime(theTime)   ! theDate(1)=hour, (2)=minute, (3)=second
      print *,"----------------------------------------------------------------------"
      print "(A47)","    Activity Code University of Birmingham 2014"
      print "(A14,A24)","    Compiled: ",compileLine
      print "(A19,I4)","    MPI Processes: ",mpiProcessCount
      print "(A16,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)",&
      "    Started at: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      print *,"----------------------------------------------------------------------"
    End If
! Initialise files
    If(mpiProcessID.eq.0)Then
! Create output file
      open(unit=999,file=trim(fileOutputData))
      write(999,"(A38)") "======================================"
      write(999,"(A38)") "ACTIVITY code University of Birmingham"
      write(999,"(A38)") "======================================"
      write(999,"(A1)") " "
      write(999,"(A10,A24)") "Compiled: ",compileLine
      write(999,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(999,"(A1)") " "
      write(999,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      close(999)
! Create activity history file
      open(unit=998,file=trim(fileActivityHistory))
      write(998,"(A38)") "======================================"
      write(998,"(A38)") "        ACTIVITY HISTORY FILE         "
      write(998,"(A38)") "======================================"
      write(998,"(A1)") " "
      write(998,"(A10,A24)") "Compiled: ",compileLine
      write(998,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(998,"(A1)") " "
      write(998,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      write(998,"(A1)") " "
      write(998,"(A21,A25)") "Time/s               ","Total isotope activity/Bq"
      close(998)
! Create in-out file
      open(unit=997,file=trim(fileInOut))
      write(997,"(A38)") "======================================"
      write(997,"(A38)") "              In-Out FILE             "
      write(997,"(A38)") "======================================"
      write(997,"(A1)") " "
      write(997,"(A10,A24)") "Compiled: ",compileLine
      write(997,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(997,"(A1)") " "
      write(997,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      close(997)
! Create isotopeActivityFile file
      open(unit=996,file=trim(fileIsotopeActivity))
      write(996,"(A38)") "======================================"
      write(996,"(A38)") "        ISOTOPE ACTIVITY FILE         "
      write(996,"(A38)") "======================================"
      write(996,"(A1)") " "
      write(996,"(A10,A24)") "Compiled: ",compileLine
      write(996,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(996,"(A1)") " "
      write(996,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      close(996)
! Create isotopeActivityFile file
      open(unit=996,file=trim(fileIsotopeActivityG))
      write(996,"(A38)") "======================================"
      write(996,"(A38)") "        ISOTOPE ACTIVITY FILE         "
      write(996,"(A38)") "======================================"
      write(996,"(A1)") " "
      write(996,"(A10,A24)") "Compiled: ",compileLine
      write(996,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(996,"(A1)") " "
      write(996,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      close(996)
! Create isotopeActivityFile file
      open(unit=995,file=Trim(outputDirectory)//"/"//"decayChains.log")
      write(995,"(A38)") "======================================"
      write(995,"(A38)") "              DECAY CHAINS            "
      write(995,"(A38)") "======================================"
      write(995,"(A1)") " "
      write(995,"(A10,A24)") "Compiled: ",compileLine
      write(995,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(995,"(A1)") " "
      write(995,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      close(995)
! Create isotopeActivityFile file
      open(unit=995,file=Trim(outputDirectory)//"/"//"ionTraj.dat")
      write(995,"(A38)") "======================================"
      write(995,"(A38)") "       AVERAGE ION TRAJECTORY         "
      write(995,"(A38)") "======================================"
      write(995,"(A1)") " "
      write(995,"(A10,A24)") "Compiled: ",compileLine
      write(995,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(995,"(A1)") " "
      write(995,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      close(995)
! Create cross section data file
      open(unit=995,file=Trim(outputDirectory)//"/"//"xs.dat")
      write(995,"(A38)") "======================================"
      write(995,"(A38)") "       CROSS SECTION DATA USED        "
      write(995,"(A38)") "======================================"
      write(995,"(A1)") " "
      write(995,"(A10,A24)") "Compiled: ",compileLine
      write(995,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(995,"(A1)") " "
      write(995,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      close(995)
! Create total reaction xs vs energy file
      open(unit=995,file=Trim(outputDirectory)//"/"//"xsTotal.dat")
      write(995,"(A38)") "======================================"
      write(995,"(A38)") "       CROSS SECTION DATA USED        "
      write(995,"(A38)") "======================================"
      write(995,"(A1)") " "
      write(995,"(A10,A24)") "Compiled: ",compileLine
      write(995,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(995,"(A1)") " "
      write(995,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      close(995)
! Create total reaction xs vs energy file
      open(unit=995,file=Trim(outputDirectory)//"/"//"gammaLines.dat")
      write(995,"(A38)") "======================================"
      write(995,"(A38)") "             GAMMA LINES              "
      write(995,"(A38)") "======================================"
      write(995,"(A1)") " "
      write(995,"(A10,A24)") "Compiled: ",compileLine
      write(995,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(995,"(A1)") " "
      write(995,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
      close(995)
    End If
  End Subroutine runInitialise

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

  Function ProgramTime () RESULT (outputTime)
! -- Argument and result
    Real(kind=DoubleReal) :: inputTime, outputTime
    Call cpu_time(inputTime)
    outputTime = inputTime - programStartTime
  End Function ProgramTime

End Module initialise
