Module input

! Setup Modules
  Use kinds
  Use constants
  Use stringfunctions		!string functions
  Use maths
  Use initialise


!force declaration of all variables
  Implicit None
  
!declare global variables  
  Character(len=2), Dimension( : ), Allocatable :: isotopesChar
  Integer, Dimension( : , : ), Allocatable :: isotopesInt
  Real, Dimension( : , : ), Allocatable :: isotopesReal
  Character(len=2), Dimension( : ), Allocatable :: decayChar
  Integer, Dimension( : , : ), Allocatable :: decayInt
  Double Precision, Dimension( : , : ), Allocatable :: decayDouble
  Character(len=2), Dimension( : ), Allocatable :: elements
  Real, Dimension( : ), Allocatable :: materialComposition
  Integer, Dimension( : ), Allocatable :: exyzKey
  Real, Dimension( : , : ), Allocatable :: exyzData
  Integer :: polyFitOrder
  Integer, Dimension( : , : ), Allocatable :: xsKey
  Real, Dimension( : , : ), Allocatable :: xsData
  Real, Dimension( : , : ), Allocatable :: xsMCData
  Real, Dimension( : , : ), Allocatable :: gammaLines
  Integer, Dimension( : , : ), Allocatable :: gammaLinesKey
  Integer :: integrationGranularity
  Real :: beamEnergy
  Real :: beamFlux
  Real :: beamDuration
  Real :: beamArea
  Real :: amTime
  Real :: timeStep
  Real :: targetThickness
  Real :: numberDensity
  Real :: materialDensity
  Real :: vpi
  Integer :: projectileZ
  Integer :: projectileA
!Program Files
  Character(len=255) :: isotopeFile
  Character(len=255) :: activityFile
  Character(len=255) :: decayModesFile
  Character(len=255) :: gammaEnergiesFile
  Character(len=255) :: xsFiles
  Character(len=255) :: trajFile
  
!Privacy of functions/subroutines/variables
  Private
  Public :: runInput				!Subroutine
  Public :: isotopesChar			!Variable
  Public :: isotopesInt  			!Variable
  Public :: isotopesReal			!Variable
  Public :: decayChar				!Variable
  Public :: decayInt				!Variable
  Public :: decayDouble				!Variable
  Public :: elements				!Variable
  Public :: materialComposition	    !Variable
  Public :: exyzKey			    	!Variable
  Public :: exyzData				!Variable
  Public :: polyFitOrder			!Variable
  Public :: xsKey			    	!Variable
  Public :: xsData					!Variable
  Public :: xsMCData				!Variable
  Public :: gammaLines  			!Variable
  Public :: gammaLinesKey  			!Variable
  Public :: integrationGranularity  !Variable
  Public :: beamEnergy  			!Variable
  Public :: beamFlux  				!Variable
  Public :: beamDuration  			!Variable
  Public :: beamArea     			!Variable
  Public :: amTime  				!Variable
  Public :: timeStep  				!Variable
  Public :: targetThickness  		!Variable
  Public :: numberDensity  		    !Variable
  Public :: materialDensity  		!Variable
  Public :: vpi  		!Variable
  Public :: projectileZ  			!Variable
  Public :: projectileA  			!Variable
  
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
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines

  Subroutine runInput()
	
	!Internal subroutine variables
	Integer :: i, j, k
		
	Call readUserInput()
	Call readIsotopeData()
	Call readActivityIn()
	Call readTraj()
	Call readNutab()
	Call readXS()
	Call readDecayModes()
	Call readGammaLines()
	Call returnInputData()

  End Subroutine runInput
  
  
!read in user input data, input from the command line
  Subroutine readUserInput()
  !force declaration of all variables
	Implicit None
  
!Read in name of input file
    call get_command_argument(1,activityFile)
	
  
  End Subroutine readUserInput
  
  
  
!read in isotope data  
  Subroutine readIsotopeData()
!force declaration of all variables
	Implicit None
	
!declare variables
	Integer, Parameter :: maxFileRows = 1E8 
	Integer :: ios, i, j, k, isotopeCounter, fileRows
	Character(len=25) :: buffera, bufferb, bufferc, bufferd, buffere, bufferf
	Character(len=255) :: bufferLong
	 
!read isotope file path from input file
    Open(UNIT=1,FILE=activityFile) 
    do i=1,maxFileRows 
	  !Read in line
	  Read(1,*,IOSTAT=ios) buffera
	  !break out
	  if (ios /= 0) then
	    EXIT 
	  end if 
	  if(buffera(1:9).eq."#isotopes")then
	    Read(1,*,IOSTAT=ios) bufferLong
	    isotopeFile = bufferLong
	  endif	  
    enddo
	CLOSE(1) !close file	
	
	
!count isotopes
	isotopeCounter = 0
	fileRows = 0
	Open(UNIT=1,FILE=isotopeFile) 
    do i=1,maxFileRows 
	  !Read in line
	  Read(1,*,IOSTAT=ios) buffera, bufferb
	  !break out
	  if (ios /= 0) then
	    EXIT 
	  end if
	  fileRows = fileRows + 1
	  if(buffera(1:6).eq."Atomic".and.bufferb(1:6).eq."Number")then
		isotopeCounter = isotopeCounter + 1
	  endif	  
    enddo
	CLOSE(1) !close file	
	!print *,isotopeCounter
	
!Allocate Arrays
    Allocate(isotopesChar(1:isotopeCounter))
    Allocate(isotopesInt(1:isotopeCounter,1:10))
    Allocate(isotopesReal(1:isotopeCounter,1:10))

!Read in data
	isotopeCounter = 0
	Open(UNIT=1,FILE=isotopeFile) 
    do i=1,fileRows 
	  !Read in line
	  Read(1,*,IOSTAT=ios) buffera, bufferb
	  !break out
	  if (ios /= 0) then
	    EXIT 
	  end if
	  if(buffera(1:6).eq."Atomic".and.bufferb(1:6).eq."Number")then
		isotopeCounter = isotopeCounter + 1
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		read(bufferd,*) isotopesInt(isotopeCounter,1) 		
	  endif	  
	  if(buffera(1:6).eq."Atomic".and.bufferb(1:6).eq."Symbol")then
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		isotopesChar(isotopeCounter) = StrToUpper(bufferd(1:2))		
	  endif	 
	  if(buffera(1:4).eq."Mass".and.bufferb(1:6).eq."Number")then
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		read(bufferd,*) isotopesInt(isotopeCounter,2) 		
	  endif	  
	  if(buffera(1:8).eq."Relative".and.bufferb(1:6).eq."Atomic")then
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, buffere
		buffere = NumericOnly(buffere)
		if(buffere(1:5).eq."     ")then
		  isotopesReal(isotopeCounter,1) = 1.0 * isotopesInt(isotopeCounter,2)
		else
		  read(buffere,*) isotopesReal(isotopeCounter,1) 		
		endif
	  endif	  
	  if(buffera(1:8).eq."Isotopic".and.bufferb(1:11).eq."Composition")then
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		bufferd = NumericOnly(bufferd)
		if(bufferd(1:5).eq."     ")then
		  isotopesReal(isotopeCounter,2) = 0.0
		else
		  read(bufferd,*) isotopesReal(isotopeCounter,2) 		
		endif
	  endif	
    enddo
	CLOSE(1) !close file	
	
	!must convert from isotopic composition to percentage by number
	
	!do i=1,isotopeCounter
	!  print *,isotopesInt(i,1),isotopesChar(i),isotopesInt(i,2),&
	!  isotopesReal(i,1),isotopesReal(i,2)
	!enddo
	

  End Subroutine readIsotopeData
  
  
  
  
!read in activity input file    
  Subroutine readActivityIn()
  
!force declaration of all variables
	Implicit None
	
!declare variables
	Integer, Parameter :: maxFileRows = 1E8 
	Integer :: ios, i, j, k, elementCounter, headerRow
	Character(len=25) :: buffera, bufferb, bufferc, bufferd
	Character(len=255) :: bufferLong
	Character(len=1) :: storeType
	real :: sumMaterialComposition
  	
 !count to allocate array rows 
    elementCounter = 0
	storeType = "U"
    Open(UNIT=1,FILE=activityFile) 
    do i=1,maxFileRows 
	  !Read in line
	  Read(1,*,IOSTAT=ios) buffera
	  !break out
	  if (ios /= 0) then
	    EXIT 
	  end if
      headerRow	= 0  
	  if(buffera(1:1).eq."#")then
	    storeType = "U"
		headerRow = 1
	  endif	  
	  if(buffera(1:9).eq."#elements")then
	    storeType = "E"
	  endif 
	  if(headerRow.eq.0)then
	    if(storeType.eq."E")then
          elementCounter = elementCounter + 1
		endif
	  endif
    enddo
	CLOSE(1) !close file	
	
!Allocate arrays
	Allocate(elements(1:elementCounter))
	Allocate(materialComposition(1:elementCounter))
	
	
!store data 
	elementCounter = 0
    Open(UNIT=1,FILE=activityFile,STATUS='OLD',ACTION='READ') 
    do i=1,maxFileRows 
	  !Read in line
	  Read(1,*,IOSTAT=ios) buffera
	  !break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  endif	  
	  
!read in files + directories
	  if(buffera(1:9).eq."#trajfile")then
	    Read(1,*,IOSTAT=ios) bufferLong
        trajFile = trim(bufferLong)
	  endif
	  if(buffera(1:11).eq."#decaymodes")then
	    Read(1,*,IOSTAT=ios) bufferLong
		decayModesFile = trim(bufferLong)
	  endif
	  if(buffera(1:14).eq."#gammaenergies")then
	    Read(1,*,IOSTAT=ios) bufferLong
		gammaEnergiesFile = trim(bufferLong)
	  endif
	  if(buffera(1:8).eq."#xsfiles")then
	    Read(1,*,IOSTAT=ios) bufferLong
		xsFiles = trim(bufferLong)
	  endif
	  
!Read in other parameters	    
	  if(buffera(1:13).eq."#polyfitorder")then
	    Read(1,*,IOSTAT=ios) buffera
	    read(buffera,*) polyFitOrder
	  endif
	  if(buffera(1:23).eq."#integrationgranularity")then
	    Read(1,*,IOSTAT=ios) buffera
	    read(buffera,*) integrationGranularity
	  endif
	  if(buffera(1:9).eq."#beamflux")then
	    Read(1,*,IOSTAT=ios) buffera
		read(buffera,*) beamFlux
	  endif
	  if(buffera(1:11).eq."#beamenergy")then
	    Read(1,*,IOSTAT=ios) buffera
        read(buffera,*) beamEnergy
	  endif
	  if(buffera(1:13).eq."#beamduration")then
		Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) beamDuration
		beamDuration = NormaliseTime (beamDuration, bufferb)
	  endif
	  if(buffera(1:9).eq."#beamarea")then
		Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) beamArea
	  endif
	  if(buffera(1:7).eq."#amtime")then
		Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) amTime
		amTime = NormaliseTime (amTime, bufferb)
	  endif
	  if(buffera(1:9).eq."#timestep")then
		Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) timeStep
		timeStep = NormaliseTime (timeStep, bufferb)
	  endif
	  if(buffera(1:11).eq."#projectile")then
	    Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) projectileZ
        read(bufferb,*) projectileA
	  endif
	  if(buffera(1:16).eq."#targetthickness")then
	    Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) targetThickness
		targetThickness = NormaliseLength (targetThickness, bufferb)
	  endif
	  if(buffera(1:16).eq."#materialdensity")then
	    Read(1,*,IOSTAT=ios) buffera, bufferb
        read(buffera,*) materialDensity
		materialDensity = NormaliseDensity (materialDensity, bufferb)
	  endif
	  if(buffera(1:4).eq."#vpi")then
	    Read(1,*,IOSTAT=ios) buffera
        read(buffera,*) vpi
	  endif
	  
!Read in elements
	  if(buffera(1:9).eq."#elements")then
	    do j=1,1000
		  !check not a new header
		  Read(1,*,IOSTAT=ios) buffera
		  if(buffera(1:1).eq." ".or.buffera(1:1).eq."#")then
		    Backspace(1)
			exit
		  else
		    Backspace(1)
			Read(1,*,IOSTAT=ios) buffera, bufferb
            elementCounter = elementCounter + 1
            elements(elementCounter) = StrToUpper(buffera(1:2))
		    read(bufferb,*) materialComposition(elementCounter)
		  endif
		enddo
      endif	  
	  
    enddo
	CLOSE(1) !close file	
	  
	!convert to fraction of 1  
	sumMaterialComposition = 0
	do i=1,size(materialComposition)
	  sumMaterialComposition = sumMaterialComposition + materialComposition(i)
	enddo
	do i=1,size(materialComposition)
	  materialComposition(i) = materialComposition(i) / sumMaterialComposition
	enddo

	
	
	
  
  End Subroutine readActivityIn  
  
  
  
  Subroutine readTraj()
  
  !force declaration of all variables
	Implicit None
	
!declare variables
	Integer, Parameter :: maxFileRows = 1E8 
	Integer :: ios, i, j, k, startCounter, dataCounter, fileRows
	Real :: x, y, z
	Character(len=20) :: buffera, bufferb, bufferc, bufferd, buffere, bufferf
  
!count data lines
    startCounter = 0
    dataCounter = 0
	fileRows = 0
    Open(UNIT=1,FILE=trajFile) 
    do i=1,maxFileRows 
	  fileRows = fileRows + 1
	  Read(1,*,IOSTAT=ios) buffera
	  if (ios /= 0) then
	    EXIT !break out
	  end if
	  !count data row
	  if(buffera(1:1).eq."#".or.buffera(1:1).eq." ".or.buffera(1:1).eq."!")then
	    !skip blank row or comment
	  else
	    if(startCounter.eq.1)then
		  dataCounter = dataCounter + 1
		endif
	  endif
	  !start counting as in data
	  if(buffera(1:7).eq."-------")then
	    startCounter = 1  
	  endif
    enddo
	CLOSE(1) !close file
    dataCounter = dataCounter + 1
  
 !Allocate arrays 
    Allocate(exyzKey(1:dataCounter))
    Allocate(exyzData(1:dataCounter,1:2))
	
!read in data
    startCounter = 0
    dataCounter = 0
    Open(UNIT=1,FILE=trajFile) 
    do i=1,fileRows 
	  if(startCounter.eq.0)then
	    Read(1,*,IOSTAT=ios) buffera
	  endif	
	  if(startCounter.eq.1)then
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, buffere
		dataCounter = dataCounter + 1
		read(buffera,*) exyzKey(dataCounter)
		read(bufferb,*) exyzData(dataCounter,2)
		read(bufferc,*) x
		read(bufferd,*) y
		read(buffere,*) z
		exyzData(dataCounter,1) = (x**2+y**2+z**2)**0.5
	  endif
	  !start counting as in data
	  if(buffera(1:7).eq."-------")then
	    startCounter = 1  
	  endif
    enddo
	CLOSE(1) !close file
  
  
    !do i=1,size(exyzKey)
	!  print *,exyzKey(i),exyzData(i,1),exyzData(i,2)
	!enddo
    
  
  End Subroutine readTraj  
  
  
  
    
  
  Subroutine readNutab()
  
	!Read in atom/isotope data
  
  
  
  
  
  End Subroutine readNutab  
  
  
  
  
  Subroutine readXS()
  
	!  
	! read in xs data, and calculate number density
	!
	! xsKey(tarZ, tarA, tarM, prodZ, prodA, prodM, start, length)
	!
    
	  
  !force declaration of all variables
	Implicit None
	
!declare variables
	Integer, Parameter :: maxFileRows = 1E8 
	Integer :: ios, i, j, k, l, dataCounter, fileRows, xsCounter, rowStart
	Integer :: reactionCounter, xsDataCounter, xsDataCounterTemp
	Integer :: xsCounterStart, xsCounterLength
	Character(len=24) :: xsFile
	Character(len=3) :: tempA, tempB, tempC, tempD
	Character(len=20) :: buffera, bufferb, bufferc, bufferd, buffere, bufferf
	Double Precision xsDataTempA, xsDataTempB
			
	!xsMCData		
			
!count xs data		
    reactionCounter = 0
	xsDataCounter = 0
!loop through elements
    do i=1,size(elements)
!loop through stable isotopes
	  do j=1,size(isotopesChar)
	    if(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.0)then
!construct xs file name
		  write(tempA,'(i3)') projectileZ
		  write(tempB,'(i3)') projectileA
		  write(tempC,'(i3)') isotopesInt(j,1)
		  write(tempD,'(i3)') isotopesInt(j,2)
		  xsFile = tempA//"_"//tempB//"-"//tempC//"_"//tempD//"_0.dat"
		  xsFile = RemoveSpaces(xsFile)	
!count data in file
          Open(UNIT=1,FILE=trim(xsFiles)//'/'//xsFile) 
          do k=1,maxFileRows 
	        Read(1,*,IOSTAT=ios) buffera
			if (ios /= 0) then
	          EXIT !break out
	        end if
			if(buffera(1:7).eq."#Header")then
		      reactionCounter = reactionCounter + 1
			  Read(1,*,IOSTAT=ios) buffera
			  Read(1,*,IOSTAT=ios) buffera
			  Read(1,*,IOSTAT=ios) buffera, bufferb
			  read(bufferb,*) xsDataCounterTemp
			  xsDataCounter = xsDataCounter + xsDataCounterTemp
			endif
          enddo
	      CLOSE(1) !close file		  
		endif  
	  enddo
	enddo
!Allocate arrays
    Allocate(xsKey(1:reactionCounter,1:8))
	Allocate(xsMCData(1:reactionCounter,1:1))	
    Allocate(xsData(1:xsDataCounter,1:2))		
!start xs counter
    xsCounterStart = 1
	xsCounterLength = 0
	reactionCounter = 0
	xsCounter = 0
!loop through elements
    do i=1,size(elements)
!loop through stable isotopes
	  do j=1,size(isotopesChar)
	    if(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.0)then
!construct xs file name
		  write(tempA,'(i3)') projectileZ
		  write(tempB,'(i3)') projectileA
		  write(tempC,'(i3)') isotopesInt(j,1)
		  write(tempD,'(i3)') isotopesInt(j,2)
		  xsFile = tempA//"_"//tempB//"-"//tempC//"_"//tempD//"_0.dat"
		  xsFile = RemoveSpaces(xsFile)
		  
!load data from file
          Open(UNIT=1,FILE=trim(xsFiles)//'/'//xsFile) 
          do k=1,maxFileRows 
	        Read(1,*,IOSTAT=ios) buffera
			if (ios /= 0) then
	          EXIT !break out
	        end if
			if(buffera(1:7).eq."#Header")then
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
!loop through data rows and store to array
              do l=1,xsCounterLength
			    xsCounter = xsCounter + 1
				Read(1,*,IOSTAT=ios) buffera, bufferb
				read(buffera,*) xsDataTempA
				read(bufferb,*) xsDataTempB
				xsData(xsCounter,1) = xsDataTempA / 1000
				xsData(xsCounter,2) = xsDataTempB
			  enddo
              
			endif
          enddo
	      CLOSE(1) !close file		  
		endif  
	  enddo
	enddo
	
	!do i=1,size(xsKey)/8
    !  print *,xsKey(i,1),xsKey(i,2),xsKey(i,3),xsKey(i,4),&
	!  xsKey(i,5),xsKey(i,6),xsKey(i,7),xsKey(i,8)
	!enddo
	
	!do i=1,size(xsData)/2
    !  print *,xsData(i,1),xsData(i,2)
	!enddo
 


	
  
  End Subroutine readXS  
  
  
  
  
  
  
  
  Subroutine readDecayModes()
    	  
!force declaration of all variables
	Implicit None
	
!declare variables
	Integer, Parameter :: maxFileRows = 1E8 
	Integer :: ios, i, j, k, decayCounter, rowCounter
	Character(len=25) :: buffera, bufferb, bufferc, bufferd
	Character(len=25) :: buffere, bufferf, bufferg, bufferh
	
 !count to allocate array rows 
    decayCounter = 0
	rowCounter = 0
    Open(UNIT=1,FILE=decayModesFile) 
    do i=1,maxFileRows 
	  !Read in line
	  Read(1,*,IOSTAT=ios) buffera
	  !break out
	  if (ios /= 0) then
	    EXIT 
	  end if
	  rowCounter = rowCounter + 1
	  decayCounter = decayCounter + 1
    enddo
	CLOSE(1) !close file	
	
!Allocate arrays
	Allocate(decayChar(1:decayCounter))
	Allocate(decayInt(1:decayCounter,1:6))
	Allocate(decayDouble(1:decayCounter,1:2))
  
!  decayChar(n) = Element e.g. FE, CR  
!  decayInt(n,1) = parent Z     
!  decayInt(n,2) = parent A    
!  decayInt(n,3) = parent meta      
!  decayInt(n,4) = child Z     
!  decayInt(n,5) = child A    
!  decayInt(n,6) = child meta    
!  decayDouble(n,1) = branching factor   
!  decayDouble(n,2) = half life   
  
 !count to allocate array rows 
    decayCounter = 0
    Open(UNIT=1,FILE=decayModesFile) 
    do i=1,rowCounter 
	  !Read in line
	  Read(1,*,IOSTAT=ios) buffera
	  !break out
	  if (ios /= 0) then
	    EXIT 
	  end if
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
      read(bufferg,*) decayDouble(decayCounter,1)
	  read(bufferh,*) decayDouble(decayCounter,2)
    enddo
	CLOSE(1) !close file
	
	!do i=1,size(decayChar)
    !  print *,i,decayChar(i),decayInt(i,1),decayInt(i,2),decayInt(i,3),decayInt(i,4),decayInt(i,5),&
	!  decayDouble(i,1),decayDouble(i,2)
	!enddo
	

  End Subroutine readDecayModes
  
    
  
  Subroutine readGammaLines()
    	  
!force declaration of all variables
	Implicit None
	
!declare variables
	Integer, Parameter :: maxFileRows = 1E8 
	Integer :: ios, i, j, k, lineCounter, rowCounter, dataCounter, fileRows
	Integer :: dataCounterTemp
	Integer :: z,a,m
	Character(len=25) :: buffera, bufferb, bufferc, bufferd
	Character(len=25) :: buffere, bufferf, bufferg, bufferh
	
!count file rows and data points
	dataCounter = 0
	fileRows = 0
    Open(UNIT=1,FILE=gammaEnergiesFile) 
    do i=1,maxFileRows 
!Read file row	  
	  Read(1,*,IOSTAT=ios) buffera
	  if (ios /= 0) then
	    EXIT !break out
	  end if
	  fileRows = fileRows + 1
      if(buffera(1:11).eq."#Datapoints")then
!re-read file line
	  BACKSPACE(1)
	    Read(1,*,IOSTAT=ios) buffera, bufferb
		read(bufferb,*) dataCounterTemp
		dataCounter = dataCounter + dataCounterTemp
	  endif		  
    enddo
	CLOSE(1) !close file
	
!Allocate array for gamma line data	
	Allocate(gammaLines(1:dataCounter,1:2))
	Allocate(gammaLinesKey(1:dataCounter,1:3))
	
!store data	
	k = 0
    Open(UNIT=1,FILE=gammaEnergiesFile) 
    do i=1,maxFileRows 
!Read file row	  
	  Read(1,*,IOSTAT=ios) buffera
	  if (ios /= 0) then
	    EXIT !break out
	  end if
      if(buffera(1:7).eq."#Parent")then
	    !re-read file line
	    BACKSPACE(1)
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		read(bufferb,*) z
		read(bufferc,*) a
		read(bufferd,*) m
	  endif
      if(buffera(1:11).eq."#Datapoints")then
!re-read file line
	    BACKSPACE(1)
	    Read(1,*,IOSTAT=ios) buffera, bufferb
		read(bufferb,*) dataCounterTemp
		if(dataCounterTemp.ge.1)then
		  do j=1,dataCounterTemp
!store data
			k = k + 1	!increment store key
!read file line
	        Read(1,*,IOSTAT=ios) buffera, bufferb
            read(buffera,*) gammaLines(k,1)
			read(bufferb,*) gammaLines(k,2)
			gammaLinesKey(k,1) = z
			gammaLinesKey(k,2) = a
			gammaLinesKey(k,3) = m
			!print *,z,a,m,buffera,bufferb
		  enddo
		endif
	  endif		  
    enddo
	CLOSE(1) !close file
	
  End Subroutine readGammaLines
  
  
  
  
  Subroutine returnInputData()
    	  
!force declaration of all variables
	Implicit None
	
!declare variables
	Integer :: ios, i, j, k
  
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
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
  
  
  End Subroutine returnInputData
  
  
  
  
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
  
function NormaliseTime (inputTime, timeUnit) RESULT (outputTime)
    ! -- Argument and result
	Real :: inputTime, factor
    CHARACTER(*) :: timeUnit
	Double Precision :: outputTime
    ! -- Local variables
    INTEGER :: i, n
  
	timeUnit = StrToUpper(timeUnit)
	
	if(timeUnit(1:2).eq."HR")then
	  factor = 3600
	elseif(timeUnit(1:2).eq."MS")then
	  factor = 0.001
	elseif(timeUnit(1:1).eq."M")then
	  factor = 60
	elseif(timeUnit(1:1).eq."S")then
	  factor = 1
	elseif(timeUnit(1:1).eq."D")then
	  factor = 86400
	endif
  
    outputTime = factor * inputTime
  
End Function NormaliseTime  
  
  
function NormaliseLength (inputLength, lengthUnit) RESULT (outputLength)
    ! -- Argument and result
	Real :: inputLength, factor
    CHARACTER(*) :: lengthUnit
	Double Precision :: outputLength
    ! -- Local variables
    INTEGER :: i, n
  
	lengthUnit = StrToUpper(lengthUnit)
	
	if(lengthUnit(1:2).eq."MM")then
	  factor = 1E7
	elseif(lengthUnit(1:2).eq."CM")then
	  factor = 1E8
	elseif(lengthUnit(1:1).eq."M")then
	  factor = 1E10
	elseif(lengthUnit(1:1).eq."A")then
	  factor = 1
	endif
  
    outputLength = factor * inputLength
  
End Function NormaliseLength    
  
  
  
function NormaliseDensity (inputDensity, densityUnit) RESULT (outputDensity)
    ! -- Argument and result
	Real :: inputDensity, factor
    CHARACTER(*) :: densityUnit
	Double Precision :: outputDensity
    ! -- Local variables
    INTEGER :: i, n
  
	densityUnit = StrToUpper(densityUnit)
	
	if(densityUnit(1:4).eq."KGM3")then
	  factor = 1
	elseif(densityUnit(1:4).eq."GCM3")then
	  factor = 1000
	endif
  
    outputDensity = factor * inputDensity
  
End Function NormaliseDensity  
  
  
  
  
  
  
  
  
  
  

End Module input