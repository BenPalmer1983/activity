Module productionLoss

! Setup Modules
  Use kinds
  Use constants
  Use stringfunctions		!string functions
  Use maths
  Use input
  Use prep
  Use initialise				! input


!force declaration of all variables
  Implicit None
  
!declare global variables  
  Double Precision, Dimension( : ), Allocatable :: reactionRate
  Real :: maxBeamDepth
  Real :: simulationScalingFactor
  Real :: affectedMaterialDepth
  Real :: affectedMaterialVolume
  Real :: dpa
  Double Precision, Dimension( : , : ), Allocatable :: gammaLineSpectra
  Double Precision :: totalActivity, totalActivityPowerOutput, affectedAtoms
  
  
!Privacy of functions/subroutines/variables
  Private
  Public :: runProductionLoss			!Subroutine
  Public :: reactionRate		        !Variable
  Public :: maxBeamDepth		        !Variable
  Public :: simulationScalingFactor     !Variable
  Public :: totalActivity               !Variable
  Public :: totalActivityPowerOutput    !Variable
  Public :: affectedMaterialDepth       !Variable
  Public :: affectedMaterialVolume      !Variable
  Public :: dpa                         !Variable
  Public :: affectedAtoms               !Variable
  
  
  
!Module Subroutines  
contains 

!------------------------------------------------------------------------!
! SUBROUTINE runProductionLoss                                           !
!------------------------------------------------------------------------!
  Subroutine runProductionLoss()
  
	!Calculate the production of all isotopes and loss of radioactive isotopes
	!Production is assumed to be time independent
	!Loss is assumed to vary with time/amount of radioactive isotope

!force declaration of all variables
	Implicit None	
!declare variables


!Run subroutines
    Call isotopeProductionRates()		!calculate the production rate
	Call runInOut()
	Call calculateDpa()
	Call calculateActivity()
  
  End Subroutine runProductionLoss

  
!------------------------------------------------------------------------!
! SUBROUTINE isotopeProductionRates                                      !
! calculate the isotope production rates                                 !
!------------------------------------------------------------------------!
  Subroutine isotopeProductionRates()
  
	!Calculate the production of all isotopes and loss of radioactive isotopes
	!Production is assumed to be time independent
	!Loss is assumed to vary with time/amount of radioactive isotope
	
	!force declaration of all variables
	Implicit None	
    !declare variables
	Integer :: i,j,k,z,a,m,key
	Integer :: lowerX, upperX
	Real :: solution
	Real :: x,y
	Real :: ionsPerSecond
	Real :: mgOfMaterial
	Real :: contentFactor
	Double precision :: maxDepth, trajectoryDepth, energyAtDepth, segmentLength, xs
	Double precision :: reactionProbability, affectedDepth,tempDoubleA,tempDoubleB
	Double precision :: averagedXS
	!Integer, Dimension( : , : ) :: reactionRateKey
	!Double Precision, Dimension( : , : ) :: reactionRateData
	!Integer, Dimension( : , : ) :: reactionKey
	!Double Precision, Dimension( : , : ) :: reactionData
	
	!  Overview
	!
	!  1. calculate average depth of ion
	!  2. calculate reaction rate for each target-product combination
	!  3. 
	
	
	!Allocate Arrays
	Allocate(reactionRate(1:size(xsKey)/8))
		
	!Production rate
	!solve equation = 0 for the trajectory
	x = 0
	y = 0
	do i=0,size(fitCoefficients)-1
	  y = y + x**(i) * fitCoefficients(i)
	enddo
	!find upper-lower boundaries
	do while(y.gt.0)
	  lowerX = x
	  x = x + 50
	  upperX = x
	  y = 0
	  do i=0,size(fitCoefficients)-1
	    y = y + x**(i) * fitCoefficients(i)	
	  enddo
	enddo
	
	!get max depth
	maxDepth = SolvePolynomial(fitCoefficients, Dble(lowerX), Dble(upperX))
	
	!affected depth
	if(maxDepth.gt.targetThickness)then
	  affectedDepth = targetThickness
	else  
	  affectedDepth = maxDepth
	endif
!set global variable
    affectedMaterialDepth = affectedDepth
	
	!maxDepth A
	!beamArea mm2  beamArea
	!density kgm3  materialDensity
	!simulationScalingFactor scales from large size to a tally of 1milligram of material
	mgOfMaterial = (beamArea * 1E-6 * affectedDepth * 1E-10) * (materialDensity * 1E6)
	simulationScalingFactor = 1 / mgOfMaterial
	
	!ions per second
	ionsPerSecond = (beamFlux * 1E-6) / elementaryCharge	!convert flux in uA to ions per second
	
!Loop over target-product combinations	
	do i=1,size(xsKey)/8		
!Get key and data for cross section	
	  z = xsKey(i,1)
	  a = xsKey(i,2)
	  m = xsKey(i,3)
	  key = 590 * z + 290 * m + a
	  tempDoubleA = 1.0 * isotopeTallyAtoms(key,1)
	  tempDoubleB = 1.0 * totalSimulationAtoms
	  contentFactor = tempDoubleA / tempDoubleB
	
!split up trajectory function	
      segmentLength = affectedDepth / integrationGranularity
	  reactionRate(i) = 0.0
	  averagedXS = 0.0
	  do j=1,integrationGranularity
	    !calculate trajectory depth
		trajectoryDepth = (affectedDepth * (j - 0.5)) / integrationGranularity
		!energy at depth		
	    energyAtDepth = 0
		do k=0,size(fitCoefficients)-1
	      energyAtDepth = energyAtDepth + trajectoryDepth**(k) * fitCoefficients(k)
	    enddo
! find reaction cross section for energy - sets xs 
		xs = 0 !make cross section equal to zero by default
		Call searchXS(energyAtDepth, xs, i)
!add to averaged cross section
		averagedXS = averagedXS + (1.0/(1.0*integrationGranularity)) * xs	
	  enddo
	  
!convert averaged cross section to m-2
      averagedXS = averagedXS * 1E-28
!store reaction rate (for 1mg of material in the beam area)
	  reactionRate(i) = ionsPerSecond * averagedXS * numberDensity * &
	  affectedDepth * 1E-10 * contentFactor * simulationScalingFactor
!store the reaction rate
	  z = xsKey(i,4)
	  a = xsKey(i,5)
	  m = xsKey(i,6)
	  key = 590 * z + 290 * m + a	  
	  isotopeTallyActive(key,3) = isotopeTallyActive(key,3) + 1.0 * reactionRate(i)	  
	enddo
  
  End Subroutine isotopeProductionRates


!------------------------------------------------------------------------!
! SUBROUTINE searchXS                                                    !
! find the xs from the loaded data                                       !
!------------------------------------------------------------------------!
  
  Subroutine searchXS(energy, xs, i)
    !force declaration of all variables
	Implicit None	
    !declare variables
	Integer :: i,j,k,startKey,endKey,key, difference
    
	Double precision :: xs, energy, tempEnergyL, tempEnergyU, factor
	Double precision :: eA,oA,eB,oB,eM,oM
	
	!xsKey(i,7)   data row start
	!xsKey(i,8)   data row length
	startKey = xsKey(i,7)
	endKey = xsKey(i,7) + xsKey(i,8) - 1
	factor = 0.5
	difference = ceiling(factor*(endKey - startKey))
	
	if(energy.lt.xsData(startKey,1))then
	  key = startKey	  
	  eA = 0
	  oA = 0
	  eB = xsData(key,1)
	  oB = xsData(key,2)
	  eM = energy
	  oM = oA + ((eM-eA)/(eB-eA))*(oB-oA)
	  xs = oM	  
    elseif(energy.gt.xsData(endKey,1))then
	  key = endKey	  
	  eA = xsData(key,1)
	  oA = xsData(key,2)
	  eB = 2*energy
	  oB = 0
	  eM = energy
	  oM = oA + ((eM-eA)/(eB-eA))*(oB-oA)
	  xs = oM	
	else
	  !start point		
	  key = startKey + difference
	  do j=1,10
	    !adjust key if too high or too low
	    if(key.ge.endKey)then
	      key = endKey-1
	    endif
	    if(key.lt.startKey)then
	      key = startKey
	    endif
	    tempEnergyL = xsData(key,1)
	    tempEnergyU = xsData(key+1,1)
	    if(energy.ge.tempEnergyL.and.energy.le.tempEnergyU)then
          !energy bound found
		  !Linear interpolate between start and end energy
	      xs = xsData(key,2) + &
		  ((energy-xsData(key,1))/(xsData(key+1,1)-xsData(key,1))) * &
		  (xsData(key+1,2)-xsData(key,2))
	      exit !break out of loop
	    else	  
	      if(energy.lt.tempEnergyL)then
	        difference = ceiling(factor * difference)
		    key = key - difference
		    !print *,"Decrease"
	      endif
	      if(energy.gt.tempEnergyU)then
	        difference = ceiling(factor * difference)
		    key = key + difference
		    !print *,"Decrease"
		  endif
	    endif
	  enddo
	endif
	
  End Subroutine searchXS
  
  

!------------------------------------------------------------------------!
! SUBROUTINE runInOut                                                    !
! Run the in-out production calculation time steps                       !
!------------------------------------------------------------------------!
  
  Subroutine runInOut()

!force declaration of all variables
	Implicit None	
!declare variables
 	Integer :: i,j,k
	Real :: simTime

!save output data table	to file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
	!open(unit=999,file=trim(outputFile))
	write(999,"(A22)") "Starting Isotope Tally"
	write(999,"(A8,A4,A4,A11,A20,A14)") "Element ","Z   ","A   ","Metastable ",& 
	"Atoms/mg            ","Activity      "
	do i=1,size(isotopeTallyChar)
      if(isotopeTallyChar(i).ne."ZZ".and.isotopeTallyAtoms(i,1).gt.0)then
	    write(999,"(A7,A1,I3.3,A1,I3.3,A1,I1.1,A10,I19.19,A1,d17.10)") &
		isotopeTallyChar(i)," ",isotopeTallyInt(i,1)," ",isotopeTallyInt(i,2)," ", &
		isotopeTallyInt(i,3),"          ",isotopeTallyAtoms(i,1)," ",&
		isotopeTallyActive(i,2)
	  endif
    enddo
	write(999,"(A1)") " "
	close(999)
	
	simTime = 0
	do i=1,1000000000
	  if(simTime.gt.amTime)then
        exit
	  endif
	  Call tallyInOut(simTime)
	  simTime = simTime + timeStep
	enddo


	
	
  End Subroutine runInOut  
  
    
  

!------------------------------------------------------------------------!
! SUBROUTINE tallyInOut                                                  !
! Run the in-out production calculation time steps                       !
!------------------------------------------------------------------------!
  
  Subroutine tallyInOut(simTime)
	
!force declaration of all variables
	Implicit None	
!declare variables
 	Integer :: i,j,k	
 	Integer :: z,m,a
    Integer :: isotopeCounter, childIsotopeCounter	
	Integer :: key
	Double Precision :: change 
	Integer(kind=ib) :: changeInteger 
	Integer, Dimension( : , : ), Allocatable :: parentIsotopes
	Integer, Dimension( : , : ), Allocatable :: childIsotopes
	Double Precision, Dimension( : , : ), Allocatable :: decayData
	Real :: simTime
	Double Precision :: tempDoubleA, tempDoubleB

!-----------------------------------	
!Notes
!-----------------------------------	
! Affected volume 1 milligram
! Number created by beam dependent only on beam, not volume affected by beam
! 
! 
	
	
!-----------------------------------	
!subtract isotopes lost due to decay
!-----------------------------------

    change = 0	
	
!Deallocate any allocated arrays	
	if(Allocated(parentIsotopes))then
	  Deallocate(parentIsotopes)
	endif
	if(Allocated(decayData))then
	  Deallocate(decayData)
	endif
	if(Allocated(childIsotopes))then
	  Deallocate(childIsotopes)
	endif

!count isotopes
    isotopeCounter = 0	
	do i=1,size(isotopeTallyChar)
      if(isotopeTallyChar(i).ne."ZZ".and.isotopeTallyAtoms(i,1).gt.0)then
	    isotopeCounter = isotopeCounter + 1
	  endif
	enddo
!allocate parent isotope array	
	Allocate(parentIsotopes(1:isotopeCounter,1:3))
!store parent isotopes
    isotopeCounter = 0	
	do i=1,size(isotopeTallyChar)
      if(isotopeTallyChar(i).ne."ZZ".and.isotopeTallyAtoms(i,1).gt.0)then
	    isotopeCounter = isotopeCounter + 1
	    parentIsotopes(isotopeCounter,1) = isotopeTallyInt(i,1)
	    parentIsotopes(isotopeCounter,2) = isotopeTallyInt(i,2)
	    parentIsotopes(isotopeCounter,3) = isotopeTallyInt(i,3)
	  endif
	enddo	
!allocate parent isotope array	
	Allocate(decayData(1:isotopeCounter,1:2))
	Allocate(childIsotopes(1:isotopeCounter,1:3))
!load decay data
    do j=1,size(parentIsotopes)/3
	  !assume isotope stable
	  decayData(j,1) = -1.0
	  decayData(j,2) = -1.0
	  childIsotopes(j,1) = 0
	  childIsotopes(j,2) = 0
	  childIsotopes(j,3) = 0
!find half-life and branching factor if unstable
      do i=1,size(decayChar)	
		if (decayInt(i,2).eq.parentIsotopes(j,1).and.decayInt(i,1).eq.parentIsotopes(j,2).and.&
		decayInt(i,3).eq.parentIsotopes(j,3)) then
		  decayData(j,1) = 1.0 * decayDouble(i,1)
		  decayData(j,2) = 1.0 * decayDouble(i,2)
	      childIsotopes(j,1) = decayInt(i,5)
	      childIsotopes(j,2) = decayInt(i,4)
	      childIsotopes(j,3) = decayInt(i,6)
		endif
	  enddo
	enddo
!calculate loss and gain from decay
    do i=1,size(parentIsotopes)/3
      if(decayData(i,2).gt.0)then
	    z = parentIsotopes(i,1)
	    a = parentIsotopes(i,2)
	    m = parentIsotopes(i,3)
        key = 590 * z + 290 * m + a
		!calculate change
		change = ceiling(isotopeTallyAtoms(key,1) * (1 - exp(-1*timeStep*(log(2.0)/decayData(i,2)))))
		!subtract from parent isotope
		isotopeTallyAtoms(key,1) = isotopeTallyAtoms(key,1) - change
		!add to child isotope
	    z = childIsotopes(i,1)
	    a = childIsotopes(i,2)
	    m = childIsotopes(i,3)
		key = 590 * z + 290 * m + a
		isotopeTallyAtoms(key,1) = isotopeTallyAtoms(key,1) + change	!tally the change in atoms
	  endif
	enddo
	
!process remaining steps in decay chain (loop)
    do k=1,10
!count child isotopes
      isotopeCounter = 0
	  do i=1,size(childIsotopes)/3
        if(decayData(i,2).gt.0)then
	      isotopeCounter = isotopeCounter + 1
	    endif
	  enddo
!clear and reallocate parent isotope array
	  Deallocate(parentIsotopes)
	  Allocate(parentIsotopes(1:isotopeCounter,1:3))
!store input child isotopes as next parent isotopes
	  isotopeCounter = 0
	  do i=1,size(childIsotopes)/3
        if(decayData(i,2).gt.0)then
	      isotopeCounter = isotopeCounter + 1
		  parentIsotopes(isotopeCounter,1) = childIsotopes(i,1)
		  parentIsotopes(isotopeCounter,2) = childIsotopes(i,2)
		  parentIsotopes(isotopeCounter,3) = childIsotopes(i,3)
	    endif
	  enddo
!clear and reallocate parent isotope array
	  Deallocate(childIsotopes)
	  Allocate(childIsotopes(1:isotopeCounter,1:3))
	  Deallocate(decayData)	
	  Allocate(decayData(1:isotopeCounter,1:2))		
!load decay data
      do j=1,size(parentIsotopes)/3
	    !assume isotope stable
	    decayData(j,1) = -1.0
	    decayData(j,2) = -1.0
	    childIsotopes(j,1) = 0
	    childIsotopes(j,2) = 0
	    childIsotopes(j,3) = 0
!find half-life and branching factor if unstable
        do i=1,size(decayChar)
		  if (decayInt(i,2).eq.parentIsotopes(j,1).and.decayInt(i,1).eq.parentIsotopes(j,2).and.&
		  decayInt(i,3).eq.parentIsotopes(j,3)) then
		    decayData(j,1) = 1.0 * decayDouble(i,1)
		    decayData(j,2) = 1.0 * decayDouble(i,2)
	        childIsotopes(j,1) = decayInt(i,5)
	        childIsotopes(j,2) = decayInt(i,4)
	        childIsotopes(j,3) = decayInt(i,6)
		  endif
	    enddo
	  enddo	
!calculate loss and gain from decay
      childIsotopeCounter = 0
      do i=1,size(parentIsotopes)/3
        if(decayData(i,2).gt.0)then
		  childIsotopeCounter = childIsotopeCounter + 1
	      z = parentIsotopes(i,1)
	      a = parentIsotopes(i,2)
	      m = parentIsotopes(i,3)
          key = 590 * z + 290 * m + a
		  !calculate change
		  change = ceiling(isotopeTallyAtoms(key,1) * (1 - exp(-1*timeStep*(log(2.0)/decayData(i,2)))))
		  !subtract from parent isotope
		  isotopeTallyAtoms(key,1) = isotopeTallyAtoms(key,1) - change
		  !add to child isotope
	      z = childIsotopes(i,1)
	      a = childIsotopes(i,2)
	      m = childIsotopes(i,3)
		  key = 590 * z + 290 * m + a
		  isotopeTallyAtoms(key,1) = isotopeTallyAtoms(key,1) + change	!tally the change in atoms
	    endif
	  enddo
!break out if no more decay steps to go
	  if(childIsotopeCounter.eq.0)then
	    exit
	  endif
!end decay loop
    enddo	  
	

!----------------------------------------	
!add/subtract isotopes due to irradiation
!----------------------------------------	

    change = 0
	if(simTime.le.beamDuration)then
      do i=1,size(xsKey)/8
	    !only do if reaction rate greater than zero
	    if(reactionRate(i).gt.0.and.simTime.lt.beamDuration)then
	      z = xsKey(i,1)
	      a = xsKey(i,2)
	      m = xsKey(i,3)
          key = 590 * z + 290 * m + a
	      change = reactionRate(i) * timeStep
		  changeInteger = change				
	      isotopeTallyAtoms(key,1) = isotopeTallyAtoms(key,1) - changeInteger
	      z = xsKey(i,4)
	      a = xsKey(i,5)
	      m = xsKey(i,6)
          key = 590 * z + 290 * m + a
	      isotopeTallyAtoms(key,1) = isotopeTallyAtoms(key,1) + changeInteger
	    endif
	    !
      enddo
	endif
	!reactionRate(i)
	!timeStep
	


  End Subroutine tallyInOut  
  
  
!------------------------------------------------------------------------!
! SUBROUTINE calculateDpa                                                !
! Displacements per Atom                                                 !
!------------------------------------------------------------------------!
  
   Subroutine calculateDpa()
	
!force declaration of all variables
	Implicit None	
!declare variables
 	Integer :: i,j,k
	Real :: time
	Double Precision :: ionsPerSecond, totalVacancies
	
  
!set dpa to 0
    dpa = 0.0
	
	affectedMaterialVolume = (affectedMaterialDepth * 1E-10) * (beamArea * 1E-6)
	affectedAtoms = affectedMaterialVolume * numberDensity	

    if(vpi.gt.0)then   
	  	  
	  if(amTime.lt.beamDuration)then
	    time = amTime
	  else
	    time = beamDuration
	  endif
	  ionsPerSecond = (beamFlux * 1E-6) / elementaryCharge	!convert flux in uA to ions per second
	  totalVacancies = ionsPerSecond * time * vpi
	  dpa = totalVacancies / affectedAtoms
	  
	  !numberDensity
	  
	  
    endif	
  
    
  End Subroutine calculateDpa  
  
  
  

!------------------------------------------------------------------------!
! SUBROUTINE calculateActivity                                           !
! Run the in-out production calculation time steps                       !
!------------------------------------------------------------------------!
  
  Subroutine calculateActivity()
	
!force declaration of all variables
	Implicit None	
!declare variables
 	Integer :: i,j,k	
	Integer :: key, finalIsotopeCount, gammaSpectraCount
	Integer, Dimension( : , : ), Allocatable :: finalIsotopeList
	Double Precision, Dimension( : , : ), Allocatable :: finalIsotopeActivity
	Real :: time
	
!Tally details
!isotopeTallyChar(i) = "ZZ"
!isotopeTallyInt(i,1) = -1		!Z
!isotopeTallyInt(i,2) = -1		!A
!isotopeTallyInt(i,3) = -1		!M
!isotopeTallyInt(i,4) = -1		!Decay calc marker
!isotopeTallyActivity(i,1) = -1.0		!half life of isotope (-1 if stable)
!isotopeTallyActivity(i,2) = 0.0		!activity of isotope in bq
!isotopeTallyActive(i,1) = -1.0		!half life of isotope (-1 if stable)
!isotopeTallyActive(i,2) = 0.0		!activity of isotope in bq
!isotopeTallyActive(i,3) = 0.0		!beam creation rate
!	
	
	
!calculate activity and count final isotopes
	finalIsotopeCount = 0
!loop through all isotopes with gt 0 atoms
	do i=1,size(isotopeTallyChar)
      if(isotopeTallyChar(i).ne."ZZ".and.isotopeTallyAtoms(i,1).gt.0&
	  .and.isotopeTallyActive(i,1).gt.0)then
!calculate activity for entire affected volume of material
	    isotopeTallyActive(i,2) = (log(2.0)/isotopeTallyActive(i,1)) * &
		isotopeTallyAtoms(i,1) * (1.0/simulationScalingFactor)
		finalIsotopeCount = finalIsotopeCount + 1
	  endif
    enddo
	

!save output data table	to file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
	write(999,"(A19)") "Final Isotope Tally"
	write(999,"(A8,A4,A4,A11,A20,A14,A14)") "Element ","Z   ","A   ","Metastable ",& 
	"Atoms/mg            ","Activity      ","T1/2          "
	totalActivity = 0.0D0
	do i=1,size(isotopeTallyChar)
      if(isotopeTallyChar(i).ne."ZZ".and.isotopeTallyAtoms(i,1).gt.0)then
	    write(999,"(A7,A1,I3.3,A1,I3.3,A1,I1.1,A10,I19.19,A1,d17.10,A1,d17.10)") &
		isotopeTallyChar(i)," ",isotopeTallyInt(i,1)," ",isotopeTallyInt(i,2)," ", &
		isotopeTallyInt(i,3),"          ",isotopeTallyAtoms(i,1)," ",&
		isotopeTallyActive(i,2)," ",isotopeTallyActive(i,1)
		!print *,isotopeTallyInt(i,1),isotopeTallyInt(i,2),isotopeTallyActive(i,2)
!calculate the total activity
		totalActivity = totalActivity + isotopeTallyActive(i,2)
	  endif
    enddo
	write(999,"(A1)") " "
	close(999)
	
	
!allocate array
    Allocate(finalIsotopeList(1:finalIsotopeCount,1:3))
    Allocate(finalIsotopeActivity(1:finalIsotopeCount,1:1))	
	
!save final radioactive isotope list
	k = 0
!loop through all isotopes with gt 0 atoms
	do i=1,size(isotopeTallyChar)
      if(isotopeTallyChar(i).ne."ZZ".and.isotopeTallyAtoms(i,1).gt.0&
	  .and.isotopeTallyActive(i,1).gt.0)then
	    k = k + 1
		finalIsotopeList(k,1) = isotopeTallyInt(i,1)
		finalIsotopeList(k,2) = isotopeTallyInt(i,2)
		finalIsotopeList(k,3) = isotopeTallyInt(i,3)
		finalIsotopeActivity(k,1) = isotopeTallyActive(i,2)
	  endif
    enddo
		
!count gamma line spectra data
    gammaSpectraCount = 0
	do i=1,(size(gammaLinesKey)/3)
		do j=1,(size(finalIsotopeList)/3)
			if(gammaLinesKey(i,1).eq.finalIsotopeList(j,1).and.&
			gammaLinesKey(i,2).eq.finalIsotopeList(j,2).and.&
			gammaLinesKey(i,3).eq.finalIsotopeList(j,3))then
			  gammaSpectraCount = gammaSpectraCount + 1
			endif
		enddo
	enddo
	
!allocate array
    Allocate(gammaLineSpectra(1:gammaSpectraCount,1:2))
	
!store data
	
	!save output data table	to file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
	write(999,"(A14)") "Gamma Spectra "
	write(999,"(A4,A4,A11,A16,A14)") "Z   ","A   ","Metastable ",& 
	"Gamma Energy/eV ","Count         "
	k = 0
	totalActivityPowerOutput = 0.0D0
	do i=1,(size(gammaLinesKey)/3)
		do j=1,(size(finalIsotopeList)/3)
			if(gammaLinesKey(i,1).eq.finalIsotopeList(j,1).and.&
			gammaLinesKey(i,2).eq.finalIsotopeList(j,2).and.&
			gammaLinesKey(i,3).eq.finalIsotopeList(j,3))then
			  k = k + 1
			  gammaLineSpectra(k,1) = gammaLines(i,1)
			  gammaLineSpectra(k,2) = gammaLines(i,2) * finalIsotopeActivity(j,1)
			  
			  write(999,"(I3.3,A1,I3.3,A1,I1.1,A10,d17.10,A1,d17.10)") &
		      finalIsotopeList(j,1)," ",finalIsotopeList(j,2),&
			  " ",finalIsotopeList(j,3),"          ", &
		      gammaLineSpectra(k,1)," ",gammaLineSpectra(k,2)
			  
			  !sum total power output
			  totalActivityPowerOutput = totalActivityPowerOutput + gammaLines(i,1) &
			  * gammaLines(i,2) * finalIsotopeActivity(j,1)
			endif
		enddo
	enddo
	write(999,"(A1)") " "
	close(999)
	
!Output totals
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
	write(999,"(A19)") "Activity Summary   "
	write(999,"(A32)") "--------------------------------"
	write(999,"(A19,d17.10)") "Total Activity/Bq: ",totalActivity
	write(999,"(A19,d17.10)") "Total power/eVs-1: ",totalActivityPowerOutput
	write(999,"(A1,d17.10)") " "
	
	write(999,"(A19)") "Calculated Values  "
	write(999,"(A32)") "--------------------------------"
	write(999,"(A19,d17.10)") "Number Density:    ",numberDensity
	write(999,"(A19,d17.10)") "Beam Depth/A:      ",maxBeamDepth
	write(999,"(A19,d17.10)") "Affected Depth/A:  ",affectedMaterialDepth
	write(999,"(A19,d17.10)") "Affected Volume/m3:",affectedMaterialVolume
	write(999,"(A19,d17.10)") "Affected Atoms     :",affectedAtoms
	if(vpi.gt.0)then
	  write(999,"(A19,d17.10)") "VPI:               ",vpi
	  write(999,"(A19,d17.10)") "DPA:               ",dpa
	endif
	write(999,"(A1)") "  "
	
	call cpu_time(time)
	time = time - programStartTime
	write(999,"(A20,d17.10)")  "Calculation Time/s: ",time
  
  End Subroutine calculateActivity    
  !isotopeTallyActivity(i)
  

  
End Module productionLoss