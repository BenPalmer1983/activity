Module input

! Setup Modules
  Use kinds
  Use strings
  Use msubs
  Use globals
  Use initialise

! force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! declare global variables

! Privacy of functions/subroutines/variables
  Private
  Public :: runInput        !Subroutine
  Public :: makeIsotopeKey      !function

!  Output data arrays
!  elements(n) = ELEMENT
!  materialComposition(n) = fraction of composition e.g. 0.7
!  exyzKey(n) = ion number
!  exyzData(n,1) = distance travelled
!  exyzData(n,2) = energy
!
!
!
!

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

  contains

! Run all the input subroutines

  Subroutine runInput()
! Internal subroutine variables
! Read in data on all processes
    Call readUserInput()
    Call readIsotopeData()
    Call readActivityIn()
    Call readTraj()
    Call readXS()
    Call readDecayModes()
    Call readGammaLines()
! Return data on just master process
    If(mpiProcessId.eq.0)Then
      Call returnInputData()
    End If
! Wait for all processes to reach the same point
    Call M_synchProcesses()
  End Subroutine runInput

! read in user input data, input from the command line
  Subroutine readUserInput()
! force declaration of all variables
    Implicit None
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
      write(999,"(A34,F8.4)") "Read user input                   ",ProgramTime()
    End If
! Read in name of input file
    call get_command_argument(1,activityFile)
! close output file
    If(mpiProcessID.eq.0)Then
      close(999)
    End If
  End Subroutine readUserInput

! read in isotope data
  Subroutine readIsotopeData()
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: ios,i,isotopeCounter,fileRows
    Integer(kind=StandardInteger) :: maxZ
    Character(len=25) :: buffera,bufferb,bufferc,bufferd,buffere
    Character(len=255) :: bufferLong
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
      write(999,"(A36,F8.4)") "Read Isotope Data                   ",ProgramTime()
    End If
! read isotope file path from input file
    Open(UNIT=1,FILE=activityFile)
    Do i=1,10000000
! Read in line
      Read(1,*,IOSTAT=ios) buffera
! break out
      If(ios /= 0)Then
        EXIT
      End If
      If(buffera(1:9).eq."#isotopes")Then
        Read(1,*,IOSTAT=ios) bufferLong
        isotopeFile = bufferLong
      End If
    End Do
    CLOSE(1) !close file
! count isotopes
    isotopeCounter = 0
    fileRows = 0
    Open(UNIT=1,FILE=isotopeFile)
    Do i=1,10000000
! Read in line
      Read(1,*,IOSTAT=ios) buffera, bufferb
! break out
      If(ios /= 0)Then
        EXIT
      End If
      fileRows = fileRows + 1
      If(buffera(1:6).eq."Atomic".and.bufferb(1:6).eq."Number")Then
        isotopeCounter = isotopeCounter + 1
      End If
    End Do
    CLOSE(1) !close file
! print *,isotopeCounter
! Allocate Arrays
    Allocate(isotopesChar(1:isotopeCounter))
    Allocate(isotopesInt(1:isotopeCounter,1:10))
    Allocate(isotopesReal(1:isotopeCounter,1:10))
! Read in data
    isotopeCounter = 0
    Open(UNIT=1,FILE=isotopeFile)
    maxZ = 0
    Do i=1,fileRows
! Read in line
      Read(1,*,IOSTAT=ios) buffera, bufferb
! break out
      If(ios /= 0)Then
        EXIT
      End If
      If(buffera(1:6).eq."Atomic".and.bufferb(1:6).eq."Number")Then
        isotopeCounter = isotopeCounter + 1
! Re-read file line
        BACKSPACE(1)
        Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
        read(bufferd,*) isotopesInt(isotopeCounter,1)
        If(isotopesInt(isotopeCounter,1).gt.maxZ)Then
          maxZ = isotopesInt(isotopeCounter,1)
        End If
      End If
      If(buffera(1:6).eq."Atomic".and.bufferb(1:6).eq."Symbol")Then
! Re-read file line
        BACKSPACE(1)
        Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
        isotopesChar(isotopeCounter) = StrToUpper(bufferd(1:2))
      End If
      If(buffera(1:4).eq."Mass".and.bufferb(1:6).eq."Number")Then
! Re-read file line
        BACKSPACE(1)
        Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
        read(bufferd,*) isotopesInt(isotopeCounter,2)
      End If
      If(buffera(1:8).eq."Relative".and.bufferb(1:6).eq."Atomic")Then
! Re-read file line
        BACKSPACE(1)
        Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, buffere
        buffere = NumericOnly(buffere)
        If(buffere(1:5).eq."     ")Then
          isotopesReal(isotopeCounter,1) = 1.0 * isotopesInt(isotopeCounter,2)
        Else
          read(buffere,*) isotopesReal(isotopeCounter,1)
        End If
      End If
      If(buffera(1:8).eq."Isotopic".and.bufferb(1:11).eq."Composition")Then
! Re-read file line
        BACKSPACE(1)
        Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
        bufferd = NumericOnly(bufferd)
        If(bufferd(1:5).eq."     ")Then
          isotopesReal(isotopeCounter,2) = 0.0
        Else
          read(bufferd,*) isotopesReal(isotopeCounter,2)
        End If
      End If
    End Do
    CLOSE(1) !close file
! make element array
    Allocate(elementSymbol(0:maxZ))
! store Z and element symbol
    elementSymbol(0) = "NN" !Neutron
    Do i=1,size(isotopesInt,1)
      elementSymbol(isotopesInt(i,1)) = isotopesChar(i)
    End Do
! close output file
    If(mpiProcessID.eq.0)Then
      close(999)
    End If
  End Subroutine readIsotopeData

! read in activity input file
  Subroutine readActivityIn()
! force declaration of all variables
    Implicit None
! declare variables
    Integer :: ios,i,j,elementCounter,headerRow
    Real(kind=DoubleReal) :: gammaResolutionTemp
    Character(len=30) :: buffera,bufferb
    Character(len=255) :: bufferLong
    Character(len=1) :: storeType
    real :: sumMaterialComposition
    Integer(kind=StandardInteger) :: fileRows
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
      write(999,"(A45,F8.4)") "Read Activity User Input File                ",ProgramTime()
    End If
! Set any defaults
    verboseTerminal = .true.
    vpi = 0.0D0
    targetDPA = 0.0D0
! count file rows
    fileRows = 0
    Open(UNIT=1,FILE=activityFile)
    Do i=1,10000000
! Read in line
      fileRows = fileRows + 1
      Read(1,*,IOSTAT=ios) buffera
! break out
      If(ios /= 0)Then
        EXIT
      End If
    End Do
    CLOSE(1) !close file
! count to allocate array rows
    elementCounter = 0
    storeType = "U"
    Open(UNIT=1,FILE=activityFile)
    Do i=1,fileRows
! Read in line
      Read(1,*,IOSTAT=ios) buffera
! break out
! If(ios /= 0)Then
!  EXIT
! End If
      headerRow  = 0
      If(buffera(1:1).eq."#")Then
        storeType = "U"
        headerRow = 1
      End If
      If(buffera(1:9).eq."#elements")Then
        storeType = "E"
      End If
      If(headerRow.eq.0)Then
        If(storeType.eq."E")Then
          elementCounter = elementCounter + 1
        End If
      End If
    End Do
    CLOSE(1) !close file
! Allocate arrays
    Allocate(elements(1:elementCounter))
    Allocate(materialComposition(1:elementCounter))
! store data
    elementCounter = 0
    Open(UNIT=1,FILE=activityFile,STATUS='OLD',ACTION='READ')
    Do i=1,fileRows
! Read in line
      Read(1,*,IOSTAT=ios) buffera
! break out if end of file
! If(ios /= 0)Then
!  EXIT
! End If
! read in files + directories
      If(buffera(1:9).eq."#trajfile")Then
        Read(1,*,IOSTAT=ios) bufferLong
        trajFile = trim(bufferLong)
      End If
      If(buffera(1:11).eq."#decaymodes")Then
        Read(1,*,IOSTAT=ios) bufferLong
        decayModesFile = trim(bufferLong)
      End If
      If(buffera(1:14).eq."#gammaenergies")Then
        Read(1,*,IOSTAT=ios) bufferLong
        gammaEnergiesFile = trim(bufferLong)
      End If
      If(buffera(1:8).eq."#xsfiles")Then
        Read(1,*,IOSTAT=ios) bufferLong
        xsFiles = trim(bufferLong)
      End If
! Read in other parameters
      If(buffera(1:13).eq."#polyfitorder")Then
        Read(1,*,IOSTAT=ios) buffera
        read(buffera,*) polyFitOrder
      End If
      If(buffera(1:23).eq."#integrationgranularity")Then
        Read(1,*,IOSTAT=ios) buffera
        read(buffera,*) integrationGranularity
      End If
      If(buffera(1:9).eq."#beamflux")Then
        Read(1,*,IOSTAT=ios) buffera
        read(buffera,*) beamFlux
      End If
      If(buffera(1:11).eq."#beamenergy")Then
        Read(1,*,IOSTAT=ios) buffera
        read(buffera,*) beamEnergy
      End If
      If(buffera(1:13).eq."#beamduration")Then
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) beamDuration
        beamDuration = NormaliseTime (1.0D0*beamDuration, bufferb)
      End If
      If(buffera(1:9).eq."#beamarea")Then
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) beamArea
      End If
      If(buffera(1:7).eq."#amtime")Then   ! Activity Measurement Time
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) amTime
        amTime = NormaliseTime (1.0D0*amTime, bufferb)
      End If
      If(buffera(1:9).eq."#timestep")Then
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) timeStep
        timeStep = NormaliseTime (1.0D0*timeStep, bufferb)
      End If
      If(buffera(1:11).eq."#projectile")Then
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) projectileZ
        read(bufferb,*) projectileA
      End If
      If(buffera(1:16).eq."#targetthickness")Then
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) targetThickness
        targetThickness = NormaliseLength (1.0D0*targetThickness, bufferb)
      End If
      If(buffera(1:16).eq."#materialdensity")Then
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) materialDensity
        materialDensity = NormaliseDensity (1.0D0*materialDensity, bufferb)
      End If
      If(buffera(1:4).eq."#vpi")Then
        Read(1,*,IOSTAT=ios) buffera
        read(buffera,*) vpi
      End If
      If(buffera(1:17).eq."#activitytimestep")Then
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) activityTimeStep
        activityTimeStep = NormaliseTime (1.0D0*activityTimeStep, bufferb)
      End If
      If(buffera(1:26).eq."#individualisotopeactivity")Then
        Read(1,*,IOSTAT=ios) buffera
        individualIsotopeActivity = StrToUpper(buffera(1:1))
      End If
      If(buffera(1:16).eq."#verboseterminal")Then
        Read(1,*,IOSTAT=ios) buffera
        If(StrToUpper(buffera(1:1)).eq."Y")Then
          verboseTerminal = .true.
        Else
          verboseTerminal = .false.
        End If
      End If
      If(buffera(1:10).eq."#targetdpa")Then
        Read(1,*,IOSTAT=ios) buffera
        read(buffera,*) targetDPA
      End If
      If(buffera(1:21).eq."#gammachartresolution")Then
        Read(1,*,IOSTAT=ios) buffera
        read(buffera,*) gammaResolutionTemp
        gammaResolution = ceiling(1.0D0*gammaResolutionTemp)
      End If
! Read in elements
      If(buffera(1:9).eq."#elements")Then
        Do j=1,1000
! check not a new header
          Read(1,*,IOSTAT=ios) buffera
          If(buffera(1:1).eq." ".or.buffera(1:1).eq."#")Then
            Backspace(1)
            exit
          Else
            Backspace(1)
            Read(1,*,IOSTAT=ios) buffera, bufferb
            elementCounter = elementCounter + 1
            elements(elementCounter) = StrToUpper(buffera(1:2))
            read(bufferb,*) materialComposition(elementCounter)
          End If
        End Do
      End If
    End Do
    CLOSE(1) !close file
    If(mpiProcessID.gt.0)Then
      verboseTerminal = .false.
    End If
! convert to fraction of 1
    sumMaterialComposition = 0
    Do i=1,size(materialComposition)
      sumMaterialComposition = sumMaterialComposition + Int(materialComposition(i))
    End Do
    Do i=1,size(materialComposition)
      materialComposition(i) = materialComposition(i) / sumMaterialComposition
    End Do
! close output file
    If(mpiProcessID.eq.0)Then
      close(999)
    End If
  End Subroutine readActivityIn

  Subroutine readTraj()
! force declaration of all variables
    Implicit None
! declare variables
    Integer :: ios,i,startCounter,dataCounter,fileRows
    Integer :: ionNumber, ionNumberMax, lastIonNumber
    Integer :: startDataPoint, dataLength
    Real :: x, y, z
    Character(len=256) :: fileRow
    Character(len=20) :: buffera,bufferb,bufferc,bufferd,buffere
    If(projectileZ.gt.0)Then    ! Do not read trajectory for Neutrons
      If(mpiProcessID.eq.0)Then
! open output file
        open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
        write(999,"(A36,F8.4)") "Read Trajectory File                ",ProgramTime()
      End If
! print if verbose on
      If(verboseTerminal.eqv..true.)Then
        print "(A22,F8.4)","    Reading exyz file ",ProgramTime()
      End If
! count data lines
      startCounter = 0
      dataCounter = 0
      fileRows = 0
      ionNumberMax = 0
      Open(UNIT=1,FILE=trajFile)
      Do i=1,10000000
        fileRows = fileRows + 1
        Read(1,"(A255)",IOSTAT=ios) fileRow
        If(ios /= 0)Then
          EXIT !break out
        End If
! count data row
        If(fileRow(1:1).eq."#".or.fileRow(1:1).eq." ".or.fileRow(1:1).eq."!")Then
! skip blank row or comment
        Else
          If(startCounter.eq.1)Then
! increment exyz data point counter
            dataCounter = dataCounter + 1
! read the ion number, to find number of ions
            Read(fileRow,*) buffera
            Read(buffera,*) ionNumber
            If(ionNumber.gt.ionNumberMax)Then
              ionNumberMax = ionNumber
            End If
          End If
        End If
! start counting as in data
        If(fileRow(1:7).eq."-------")Then
          startCounter = 1
        End If
      End Do
      CLOSE(1) !close file
      dataCounter = dataCounter + 1
! Allocate arrays
      Allocate(exyzKey(1:dataCounter))
      Allocate(exyzIonKey(1:ionNumberMax,1:2))
      Allocate(exyzData(1:dataCounter,1:2))
! Read in trajectory data
      startCounter = 0
      dataCounter = 0
      ionNumber = 1
      lastIonNumber = 1
      startDataPoint = 1
      dataLength = 0
      Open(UNIT=1,FILE=trajFile)
      Do i=1,fileRows
        Read(1,"(A255)",IOSTAT=ios) fileRow
        If(fileRow(1:4).ne."    ")Then
          If(startCounter.eq.1)Then
            dataCounter = dataCounter + 1
            dataLength = dataLength + 1
! Read in line
            Read(fileRow,*) buffera, bufferb, bufferc, bufferd, buffere
! Read ion number
            Read(buffera,*) ionNumber
            If(ionNumber.ne.lastIonNumber)Then
              dataLength = dataLength - 1 !don't count this row for last ion
! print *,startDataPoint,dataLength,ionNumber,lastIonNumber
! Store to key
              exyzIonKey(lastIonNumber,1) = startDataPoint
              exyzIonKey(lastIonNumber,2) = dataLength
! Update start counter and length
              startDataPoint = dataCounter
              dataLength = 1  !count this row for next ion
            End If
! Store point to data array
! read(buffera,*) exyzKey(dataCounter)
            read(bufferb,*) exyzData(dataCounter,2)
            read(bufferc,*) x
            read(bufferd,*) y
            read(buffere,*) z
            exyzData(dataCounter,1) = (x**2+y**2+z**2)**0.5
! Last ion number
            lastIonNumber = ionNumber
          End If
! Identify when in data section of exyz file
          If(fileRow(1:7).eq."-------")Then
            startCounter = 1
          End If
        End If
      End Do
! Store last ion
! print *,startDataPoint,(dataLength-1),ionNumber,lastIonNumber
      exyzIonKey(lastIonNumber,1) = startDataPoint
      exyzIonKey(lastIonNumber,2) = (dataLength-1)
! write to output file
      If(mpiProcessID.eq.0)Then
        write(999,"(A18,I8,A8,F8.4)") "Data point count: ",dataCounter,"        ",ProgramTime()
      End If
! close output file
      If(mpiProcessID.eq.0)Then
        close(999)
      End If
    End If
  End Subroutine readTraj

  Subroutine readXS()
! read in xs data, and calculate number density
! xsKey(tarZ, tarA, tarM, prodZ, prodA, prodM, start, length)
! force declaration of all variables
    Implicit None
! declare variables
    Integer :: ios,i,j,k,l,xsCounter
    Integer :: reactionCounter, xsDataCounter, xsDataCounterTemp
    Integer :: xsCounterStart, xsCounterLength, xsCounterEnd, switch
    Character(len=24) :: xsFile
    Character(len=3) :: tempA, tempB, tempC, tempD
    Character(len=20) :: buffera,bufferb,bufferc,bufferd
    Double Precision xsDataTempA, xsDataTempB
    Real(kind=DoubleReal) :: energyMin, energyMax, energyVal
    Real(kind=DoubleReal) :: energyValStore, xsValStore
    Real(kind=DoubleReal), Dimension(0:200,1:2) :: xsTotal
! Init vars
    energyMin = 0.0D0
    energyMax = 0.0D0
! open output file
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
      write(999,"(A1)") " "
      write(999,"(A21)") "XS Data Files Loaded:"
      open(unit=998,file=trim(Trim(outputDirectory)//"/"//"xs.dat"),&
      status="old",position="append",action="write")
    End If
! print if verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A20,F8.4)","    Reading xs Data ",ProgramTime()
    End If
! count xs data
    reactionCounter = 0
    xsDataCounter = 0
! loop through elements
    Do i=1,size(elements)
! loop through stable isotopes
      Do j=1,size(isotopesChar)
        If(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.0)Then
! construct xs file name
          write(tempA,'(i3)') projectileZ
          write(tempB,'(i3)') projectileA
          write(tempC,'(i3)') isotopesInt(j,1)
          write(tempD,'(i3)') isotopesInt(j,2)
          xsFile = tempA//"_"//tempB//"-"//tempC//"_"//tempD//"_0.dat"
          xsFile = RemoveSpaces(xsFile)
          If(mpiProcessID.eq.0)Then
            write(999,"(A16,A24)") "Reading xs file ",xsFile
          End If
! count data in file
          Open(UNIT=1,FILE=trim(xsFiles)//'/'//xsFile)
          Do k=1,10000000
            Read(1,*,IOSTAT=ios) buffera
            If(ios /= 0)Then
              EXIT !break out
            End If
            If(buffera(1:7).eq."#Header")Then
              reactionCounter = reactionCounter + 1
              Read(1,*,IOSTAT=ios) buffera
              Read(1,*,IOSTAT=ios) buffera
              Read(1,*,IOSTAT=ios) buffera, bufferb
              read(bufferb,*) xsDataCounterTemp
              xsDataCounter = xsDataCounter + xsDataCounterTemp
            End If
          End Do
          CLOSE(1) !close file
        End If
      End Do
    End Do
! Allocate arrays
    Allocate(xsKey(1:reactionCounter,1:8))
    Allocate(xsMCData(1:reactionCounter,1:1))
    Allocate(xsData(1:xsDataCounter,1:2))
! start xs counter
    xsCounterStart = 1
    xsCounterLength = 0
    reactionCounter = 0
    xsCounter = 0
! loop through elements
    Do i=1,size(elements)
! loop through stable isotopes
      Do j=1,size(isotopesChar)
        If(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.0)Then
! construct xs file name
          write(tempA,'(i3)') projectileZ
          write(tempB,'(i3)') projectileA
          write(tempC,'(i3)') isotopesInt(j,1)
          write(tempD,'(i3)') isotopesInt(j,2)
          xsFile = tempA//"_"//tempB//"-"//tempC//"_"//tempD//"_0.dat"
          xsFile = RemoveSpaces(xsFile)
! load data from file
          Open(UNIT=1,FILE=trim(xsFiles)//'/'//xsFile)
          Do k=1,10000000
            Read(1,*,IOSTAT=ios) buffera
            If(ios /= 0)Then
              EXIT !break out
            End If
            If(buffera(1:7).eq."#Header")Then
              reactionCounter = reactionCounter + 1
              Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
              read(bufferb,*) xsKey(reactionCounter,1)
              read(bufferc,*) xsKey(reactionCounter,2)
              read(bufferd,*) xsKey(reactionCounter,3)
              Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
              read(bufferb,*) xsKey(reactionCounter,4)
              read(bufferc,*) xsKey(reactionCounter,5)
              read(bufferd,*) xsKey(reactionCounter,6)
              Read(1,*,IOSTAT=ios) buffera, bufferb
              read(bufferb,*) xsCounterLength
              xsKey(reactionCounter,7) = xsCounterStart
              xsKey(reactionCounter,8) = xsCounterLength
              xsCounterStart = xsCounterStart + xsCounterLength
! loop through data rows and store to array
              Do l=1,xsCounterLength
                xsCounter = xsCounter + 1
                Read(1,*,IOSTAT=ios) buffera, bufferb
                read(buffera,*) xsDataTempA
                read(bufferb,*) xsDataTempB
                xsData(xsCounter,1) = xsDataTempA / 1000    ! Save as KeV
                xsData(xsCounter,2) = xsDataTempB
                If(mpiProcessID.eq.0)Then
                  write(998,"(I4,I4,I4,I4,I4,E14.6,E14.6)") &
                  projectileZ, projectileA, isotopesInt(j,1), isotopesInt(j,2),&
                  xsCounter,xsData(xsCounter,1),xsData(xsCounter,2)
                End If
                If(xsData(xsCounter,1).gt.energyMax)Then
                  energyMax = xsData(xsCounter,1)
                End If
              End Do
            End If
          End Do
          CLOSE(1) !close file
        End If
      End Do
    End Do
    If(mpiProcessID.eq.0)Then
      xsTotal = 0.0D0
! Output total reaction cross section vs energy
      open(unit=995,file=trim(Trim(outputDirectory)//"/"//"xsTotal.dat"),&
      status="old",position="append",action="write")
      Do k=1,size(xsKey,1)
        xsCounterStart = xsKey(k,7)
        xsCounterLength = xsKey(k,8)
        xsCounterEnd = xsCounterStart + xsCounterLength - 1
        j = xsCounterStart
        Do i=0,200
          energyVal = i*(energyMax/(1.0D0*200))
          If(energyVal.gt.xsData(xsCounterEnd,1))Then    ! If above available xs data
            exit
          End If
          If(energyVal.lt.xsData(xsCounterStart,1))Then  ! If below available cross section data
            energyValStore = energyVal
            xsValStore = 0.0D0
          Else
            If(j.lt.xsCounterEnd)Then
              switch = 0
              Do while(switch.eq.0)
                If(energyVal.ge.xsData(j+1,1))Then
                  j = j + 1
                Else
                  switch = 1
                End If
                If(j.eq.xsCounterEnd)Then
                  switch = 1
                End If
              End Do
            End If
            energyValStore = xsData(j,1)
            xsValStore = xsData(j,2)
          End If
          xsTotal(i,1) = energyValStore
          xsTotal(i,2) = xsTotal(i,2) + xsValStore
        End Do
      End Do
      Do i=0,200
        write(995,"(E14.6,E14.6)") xsTotal(i,1),xsTotal(i,2)
      End Do
      close(995)
    End If
    If(mpiProcessID.eq.0)Then
      write(999,"(A1)") " "
      Close(999)
      Close(998)
    End If
  End Subroutine readXS

  Subroutine readDecayModes()
! force declaration of all variables
    Implicit None
! declare variables
    Integer :: ios,i,decayCounter,rowCounter
    Character(len=25) :: buffera, bufferb, bufferc, bufferd
    Character(len=25) :: buffere, bufferf, bufferg, bufferh
! count to allocate array rows
    decayCounter = 0
    rowCounter = 0
    Open(UNIT=1,FILE=decayModesFile)
    Do i=1,10000000
! Read in line
      Read(1,*,IOSTAT=ios) buffera
! break out
      If(ios /= 0)Then
        EXIT
      End If
      rowCounter = rowCounter + 1
      decayCounter = decayCounter + 1
    End Do
    CLOSE(1) !close file
! Allocate arrays
    Allocate(decayChar(1:decayCounter))
    Allocate(decayInt(1:decayCounter,1:8))
    Allocate(decayDouble(1:decayCounter,1:2))
!  decayChar(n) = Element e.g. FE, CR
!  decayInt(n,1) = parent A       !A,Z wrong way round
!  decayInt(n,2) = parent Z
!  decayInt(n,3) = parent meta
!  decayInt(n,4) = child A
!  decayInt(n,5) = child Z
!  decayInt(n,6) = child meta
!  decayDouble(n,1) = branching factor
!  decayDouble(n,2) = half life
! count to allocate array rows
    decayCounter = 0
    Open(UNIT=1,FILE=decayModesFile)
    Do i=1,rowCounter
! Read in line
      Read(1,*,IOSTAT=ios) buffera
! break out
      If(ios /= 0)Then
        EXIT
      End If
      rowCounter = rowCounter + 1
      decayCounter = decayCounter + 1
      BACKSPACE(1)
      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, &
      buffere, bufferf, bufferg, bufferh
      decayChar(decayCounter) = buffera(1:2)
      read(bufferb,*) decayInt(decayCounter,1)
      read(bufferc,*) decayInt(decayCounter,2)
      read(bufferd,*) decayInt(decayCounter,3)
      read(buffere,*) decayInt(decayCounter,4)
      read(bufferf,*) decayInt(decayCounter,5)
      decayInt(decayCounter,6) = 0
      decayInt(decayCounter,7) = &
      makeIsotopeKey(decayInt(decayCounter,2),decayInt(decayCounter,1),&
      decayInt(decayCounter,3))
      decayInt(decayCounter,8) = &
      makeIsotopeKey(decayInt(decayCounter,5),decayInt(decayCounter,4),&
      decayInt(decayCounter,6))
      read(bufferg,*) decayDouble(decayCounter,1)
      read(bufferh,*) decayDouble(decayCounter,2)
    End Do
    CLOSE(1) !close file
  End Subroutine readDecayModes

  Subroutine readGammaLines()
! force declaration of all variables
    Implicit None
! declare variables
    Integer :: ios,i,j,k,dataCounter,fileRows
    Integer :: dataCounterTemp
    Integer :: z,a,m
    Character(len=25) :: buffera, bufferb, bufferc, bufferd
! count file rows and data points
    dataCounter = 0
    fileRows = 0
    Open(UNIT=1,FILE=gammaEnergiesFile)
    Do i=1,10000000
! Read file row
      Read(1,*,IOSTAT=ios) buffera
      If(ios /= 0)Then
        EXIT !break out
      End If
      fileRows = fileRows + 1
      If(buffera(1:11).eq."#Datapoints")Then
! re-read file line
        BACKSPACE(1)
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(bufferb,*) dataCounterTemp
        dataCounter = dataCounter + dataCounterTemp
      End If
    End Do
    CLOSE(1) !close file
! Allocate array for gamma line data
    Allocate(gammaLines(1:dataCounter,1:2))
    Allocate(gammaLinesKey(1:dataCounter,1:3))
! store data
    k = 0
    Open(UNIT=1,FILE=gammaEnergiesFile)
    Do i=1,10000000
! Read file row
      Read(1,*,IOSTAT=ios) buffera
      If(ios /= 0)Then
        EXIT !break out
      End If
      If(buffera(1:7).eq."#Parent")Then
! re-read file line
        BACKSPACE(1)
        Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
        read(bufferb,*) z
        read(bufferc,*) a
        read(bufferd,*) m
      End If
      If(buffera(1:11).eq."#Datapoints")Then
! re-read file line
        BACKSPACE(1)
        Read(1,*,IOSTAT=ios) buffera, bufferb
        read(bufferb,*) dataCounterTemp
        If(dataCounterTemp.ge.1)Then
          Do j=1,dataCounterTemp
! store data
            k = k + 1  !increment store key
! read file line
            Read(1,*,IOSTAT=ios) buffera, bufferb
            read(buffera,*) gammaLines(k,1) !energy
            read(bufferb,*) gammaLines(k,2) !intensity
            gammaLinesKey(k,1) = z
            gammaLinesKey(k,2) = a
            gammaLinesKey(k,3) = m
! print *,z,a,m,buffera,bufferb
          End Do
        End If
      End If
    End Do
    CLOSE(1) !close file
  End Subroutine readGammaLines

  Subroutine returnInputData()
! force declaration of all variables
    Implicit None
! declare variables
! Print to file
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
      write(999,"(A20)") "Calculation Settings"
      write(999,"(A1)") " "
      write(999,"(A32,d17.10)") "Beam Duration/s:                ",beamDuration
      write(999,"(A32,d17.10)") "Beam Current/uA:                ",beamFlux
      write(999,"(A32,d17.10)") "Beam Energy/MeV:                ",beamEnergy
      write(999,"(A32,d17.10)") "Beam Area/mm2:                  ",beamArea
      write(999,"(A32,d17.10)") "Target Thickness/Angstrom:      ",targetThickness
      write(999,"(A32,d17.10)") "Target Density/kgm-3:           ",materialDensity
      write(999,"(A32,d17.10)") "Activity Measurement Time/s:    ",amTime
      write(999,"(A32,A60)") "Input File:                     ",activityFile
      write(999,"(A32,A60)") "Isotope File:                   ",isotopeFile
      write(999,"(A32,A60)") "Decay Modes File:               ",decayModesFile
      write(999,"(A32,A60)") "Gamma Energies File:            ",gammaEnergiesFile
      write(999,"(A32,A60)") "XS File Directory:              ",xsFiles
      write(999,"(A1)") " "
      write(999,"(A1)") " "
      close(999)
    End If
! Print to screen
    If(verboseTerminal.eqv..true.)Then
      print "(A27)","    User Input Sim Details:"
      print "(A38,d17.10)","      Beam Duration/s:                ",beamDuration
      print "(A38,d17.10)","      Beam Current/uA:                ",beamFlux
      print "(A38,d17.10)","      Beam Energy/MeV:                ",beamEnergy
      print "(A38,d17.10)","      Beam Area/mm2:                  ",beamArea
      print "(A38,d17.10)","      Target Thickness/Angstrom:      ",targetThickness
      print "(A38,d17.10)","      Target Density/kgm-3:           ",materialDensity
      print "(A38,d17.10)","      Activity Measurement Time/s:    ",amTime
    End If
  End Subroutine returnInputData

  Subroutine sendInput()
! force declaration of all variables
    Implicit None
! declare variables
! Data to send: exyzIonKey, exyzData
! send
! Call MPI_sendData2DDP_A(exyzData,size(exyzData,1),size(exyzData,2))
    Call M_distDouble2D(exyzData)
  End Subroutine sendInput

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

  Function NormaliseTime (inputTime, timeUnit) RESULT (outputTime)
! -- Argument and result
    Real(kind=DoubleReal) :: inputTime, factor
    CHARACTER(*) :: timeUnit
    Real(kind=DoubleReal) :: outputTime
! -- Local variables
    factor = 1.0D0
    timeUnit = StrToUpper(timeUnit)
    If(timeUnit(1:2).eq."HR")Then
      factor = 3600.0D0
    ElseIf(timeUnit(1:2).eq."MS")Then
      factor = 0.001D0
    ElseIf(timeUnit(1:1).eq."M")Then
      factor = 60.0D0
    ElseIf(timeUnit(1:1).eq."S")Then
      factor = 1.0D0
    ElseIf(timeUnit(1:1).eq."D")Then
      factor = 86400.0D0
    End If
    outputTime = factor * inputTime
  End Function NormaliseTime

  Function NormaliseLength (inputLength, lengthUnit) RESULT (outputLength)
! -- Argument and result
    Real(kind=DoubleReal) :: inputLength, factor
    CHARACTER(*) :: lengthUnit
    Real(kind=DoubleReal) :: outputLength
! -- Local variables
    factor = 1.0D0
    lengthUnit = StrToUpper(lengthUnit)
    If(lengthUnit(1:2).eq."MM")Then
      factor = 1.0D7
    ElseIf(lengthUnit(1:2).eq."CM")Then
      factor = 1.0D8
    ElseIf(lengthUnit(1:1).eq."M")Then
      factor = 1.0D10
    ElseIf(lengthUnit(1:1).eq."A")Then
      factor = 1.0D0
    End If
    outputLength = factor * inputLength
  End Function NormaliseLength

  Function NormaliseDensity (inputDensity, densityUnit) RESULT (outputDensity)
! -- Argument and result
    Real(kind=DoubleReal) :: inputDensity, factor
    CHARACTER(*) :: densityUnit
    Real(kind=DoubleReal) :: outputDensity
! -- Local variables
    factor = 1.0D0
    densityUnit = StrToUpper(densityUnit)
    If(densityUnit(1:4).eq."KGM3")Then
      factor = 1.0D0
    ElseIf(densityUnit(1:4).eq."GCM3")Then
      factor = 1000.0D0
    End If
    outputDensity = factor * inputDensity
  End Function NormaliseDensity

  Function makeIsotopeKey (z, a, m) RESULT (key)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: z,a,m,key
    key = 590 * z + 290 * m + a
! returns y
  End Function makeIsotopeKey

End Module input
