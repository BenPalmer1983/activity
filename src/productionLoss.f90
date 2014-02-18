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
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: reactionRate
  Integer(kind=StandardInteger), Dimension( : , :), Allocatable :: targetReactionRatesInt
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: targetReactionRates
  Integer(kind=StandardInteger), Dimension( : , :), Allocatable :: productReactionRatesInt
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: productReactionRates
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: averageXS
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: numberDensityArray
  Integer(kind=StandardInteger), Dimension( : , :), &
  Allocatable :: decayChainIsotopesArray
  Integer(kind=StandardInteger) :: isotopeCounter, publicX, publicY, publicZ 
  Integer(kind=StandardInteger), Dimension(:,:), Allocatable :: public2DTempArray 
  Integer(kind=StandardInteger) :: simIsotopes
  Real(kind=DoubleReal), Dimension(:,:,:), Allocatable :: public3DTempArray 
  Real :: maxBeamDepth
  Real :: simulationScalingFactor
  Real :: affectedMaterialDepth
  Real :: affectedMaterialVolume
  Real :: dpa
  Real :: totalActivityAtTime
  Double Precision, Dimension( : , : ), Allocatable :: gammaLineSpectra
  Double Precision :: totalActivity, totalActivityPowerOutput, affectedAtoms
  !Integer, Dimension( : , : ), Allocatable :: decayIntReduced
  !Double Precision, Dimension( : , : ), Allocatable :: decayDoubleReduced
  Double Precision, Dimension( : ), Allocatable :: decayTempTally
  !Integer(kind=StandardInteger), Dimension( : ), Allocatable :: simIsotopeKeys
  
  
!Privacy of functions/subroutines/variables
  Private
  Public :: runProductionLoss			!Subroutine
  Public :: reactionRate		        !Variable
  Public :: maxBeamDepth		        !Variable
  Public :: simulationScalingFactor     !Variable
  Public :: totalActivityAtTime
  Public :: totalActivity               !Variable
  Public :: totalActivityPowerOutput    !Variable
  Public :: affectedMaterialDepth       !Variable
  Public :: affectedMaterialVolume      !Variable
  Public :: dpa                         !Variable
  Public :: affectedAtoms               !Variable
  Public :: decayChainIsotopesArray
  Public :: isotopeCounter
  Public :: publicX
  Public :: publicY 
  Public :: publicZ
  Public :: simIsotopes
  Public :: productReactionRatesInt
  Public :: productReactionRates
  Public :: targetReactionRatesInt
  Public :: targetReactionRates
  
  
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
    Call preCalcOutput()		!calculate the production rate
    Call isotopeProductionRates()		!calculate the production rate
	Call prepareIsotopeTally()
	!Call makeReducedDecayList()
	Call runInOut()
  
  End Subroutine runProductionLoss

  
  
!------------------------------------------------------------------------!
! SUBROUTINE 
!                            
!------------------------------------------------------------------------!
  Subroutine preCalcOutput()
!force declaration of all variables
	Implicit None	
!declare variables
	Integer(kind=StandardInteger) :: i,j,k
  !open output file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	write(999,"(A1)") " "
    write(999,"(A70)") "----------------------------------------------------------------------"
	write(999,"(A22)") "Pre-Simulation Summary"
	write(999,"(F8.4)") ProgramTime()
	write(999,"(A70)") "----------------------------------------------------------------------"
	write(999,"(A1)") " "
	write(999,"(A30,E20.10)") "Beam flux ions/s:             ",((beamFlux*1D-6)/elementaryCharge)
  
!close output file
	write(999,"(A1)") " "
    close(999)
    
  
  End Subroutine preCalcOutput
  
  
  
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
	Integer(kind=StandardInteger) :: i,j,k,z,a,m,key
	Integer(kind=StandardInteger) :: lowerX, upperX
	Integer(kind=StandardInteger) :: zA,aA,mA,zB,aB,mB
	Integer(kind=StandardInteger) :: store, storeCounter
    Real(kind=DoubleReal) :: volumeAffected
	Real :: solution
	Real :: x,y
	Real :: ionsPerSecond
	Real :: mgOfMaterial
	Real :: contentFactor
	Real(kind=DoubleReal) :: maxDepth, trajectoryDepth, energyAtDepth, segmentLength, xs
	Real(kind=DoubleReal) :: reactionProbability, affectedDepth,tempDoubleA,tempDoubleB
	Real(kind=DoubleReal) :: averagedXS
	Character(len=255) :: fittingPolynomial	
	Character(len=255) :: tempChar	
	Integer(kind=StandardInteger), Dimension( : , :), Allocatable :: targetReactionRatesIntTemp
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: targetReactionRatesTemp
	Integer(kind=StandardInteger), Dimension( : , :), Allocatable :: productReactionRatesIntTemp
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: productReactionRatesTemp
	!  Overview
	!
	!  1. calculate average depth of ion
	!  2. calculate reaction rate for each target-product combination	
!open output file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	write(999,"(A1)") " "
	write(999,"(A51,F8.4)") "Calculate projectile-target-product reaction rates ",ProgramTime()	
!Allocate Arrays
	Allocate(reactionRate(1:size(xsKey,1)))
	Allocate(averageXS(1:size(xsKey,1)))
	Allocate(numberDensityArray(1:size(xsKey,1)))			
!Production rate
!solve equation = 0 for the trajectory
	x = 0
	y = 0
	fittingPolynomial = "y = "
	write(999,"(A25)") "Trajectory fit equation: "
	do i=0,size(fitCoefficients)-1
	  y = y + x**(i) * fitCoefficients(i)
	  write(999,"(I8,E20.10)") i,fitCoefficients(i)
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
	write(999,"(A25,E20.10)") "Affected material depth: ",affectedMaterialDepth
!ions per second
	ionsPerSecond = (beamFlux * 1.0D-6) / elementaryCharge	!convert flux in uA to ions per second
!volume of material affected by beam
    volumeAffected = 1.0D0 * beamArea * 1E-6 * affectedDepth * 1E-10
!update tally to show atoms in affected volume only
	totalSimulationAtoms = 0.0D0
	Do i=1,size(simIsotopeKeys) 
	  key = simIsotopeKeys(i)
	  totalSimulationAtoms = totalSimulationAtoms + isotopeTallyActive(key,4)
	End Do
	
!Loop over target-product combinations	
	do i=1,size(xsKey,1)		
!Get key and data for cross section	
	  z = xsKey(i,1)
	  a = xsKey(i,2)
	  m = xsKey(i,3)
	  key = 590 * z + 290 * m + a
	  tempDoubleA = 1.0D0 * isotopeTallyActive(key,4)
	  tempDoubleB = 1.0D0 * totalSimulationAtoms
	  contentFactor = tempDoubleA / tempDoubleB
	
!split up trajectory function	
      segmentLength = affectedDepth / integrationGranularity
	  reactionRate(i) = 0.0
	  averagedXS = 0.0D0
	  do j=1,integrationGranularity
	    !calculate trajectory depth
		trajectoryDepth = (affectedDepth * (j - 0.5)) / integrationGranularity
		!energy at depth		
	    energyAtDepth = 0
		do k=0,size(fitCoefficients)-1
	      energyAtDepth = energyAtDepth + trajectoryDepth**(k) * fitCoefficients(k)
	    enddo
! find reaction cross section for energy - sets xs 
		xs = 0.0D0 !make cross section equal to zero by default
		Call searchXS(energyAtDepth, xs, i)
!add to averaged cross section
		averagedXS = averagedXS + (1.0D0/(1.0D0*integrationGranularity)) * xs	
		!Testing
		!If(xsKey(i,1).eq.26.and.xsKey(i,2).eq.56.and.&
		!xsKey(i,4).eq.25.and.xsKey(i,5).eq.48)Then
		!  xs = searchXSa(energyAtDepth,i)
		!  print *,j,energyAtDepth,xs
		!End If
	  enddo
	  
!convert averaged cross section to m-2
      averagedXS = averagedXS * 1.0D-28
	  averageXS(i) = averagedXS
!store reaction rate (for 1mg of material in the beam area)
	  reactionRate(i) = ionsPerSecond * averagedXS * numberDensity * &
	  affectedDepth * 1.0D-10 * contentFactor
	  numberDensityArray(i) = 1.0D0 * numberDensity * contentFactor
	enddo	
!output variables to file	
    write(999,"(A1)") " "
    write(999,"(A70)") "----------------------------------------------------------------------"
    write(999,"(A20)") "Calculated variables"
    write(999,"(A70)") "----------------------------------------------------------------------"
	write(999,"(A30,E20.10)") "Beam flux ions/s:             ",ionsPerSecond
	write(999,"(A30,E20.10)") "Affected depth/m:             ",(affectedDepth*1.0D-10)
	write(999,"(A30,E20.10)") "Affected depth/ang:           ",(affectedDepth)
	write(999,"(A30,E20.10)") "Affected volume/m3:           ",volumeAffected
	write(999,"(A30,E20.10)") "Number Density:               ",numberDensity
    write(999,"(A1)") " "
!output reaction rates to output data file
    write(999,"(A1)") " "
    write(999,"(A70)") "----------------------------------------------------------------------"
	write(999,"(A38,F8.4)") "Reaction rates in the target material ",ProgramTime()
	write(999,"(A59,A40)") "       Sym Z   A   M   Sym Z   A   M   Reaction Rate       ",&
	                       "Average XS          Number Density      "
	write(999,"(A70)") "----------------------------------------------------------------------"
	do i=1,size(xsKey,1)
!only do if reaction rate greater than zero
	  if(reactionRate(i).gt.0)then
	    zA = xsKey(i,1)
	    aA = xsKey(i,2)
	    mA = xsKey(i,3)
	    zB = xsKey(i,4)
	    aB = xsKey(i,5)
	    mB = xsKey(i,6)
		write(999,"(I4,A3,A4,I4,I4,I4,A4,I4,I4,I4,E20.10,E20.10,E20.10)") i,":  ",&
		elementSymbol(zA),zA,aA,mA,elementSymbol(zB),zB,aB,mB,&
		(1.0D0*reactionRate(i)),averageXS(i),numberDensityArray(i)
	  endif
    enddo
!store target rates
	Allocate(targetReactionRatesIntTemp(1:size(xsKey,1),1:3))
	Allocate(targetReactionRatesTemp(1:size(xsKey,1)))	
	storeCounter = 0
	Do i=1,size(xsKey,1)
	  !Only store if greater than zero
	  If(reactionRate(i).ne.0.0D0)Then
	    store = 0 !store as new = yes 
	    Do j=1,size(xsKey,1)
	      If(targetReactionRatesIntTemp(j,1).eq.xsKey(i,1).and.&
		    targetReactionRatesIntTemp(j,2).eq.xsKey(i,2).and.&
		    targetReactionRatesIntTemp(j,3).eq.xsKey(i,3))Then
	          store = j !store as new = no 
	     End If
	    End Do
        If(store.eq.0)Then
          storeCounter = storeCounter + 1
		  targetReactionRatesIntTemp(storeCounter,1) = xsKey(i,1)
		  targetReactionRatesIntTemp(storeCounter,2) = xsKey(i,2)
		  targetReactionRatesIntTemp(storeCounter,3) = xsKey(i,3)
		  targetReactionRatesTemp(storeCounter) = -1.0D0*reactionRate(i)
	    Else	
		  targetReactionRatesTemp(store) = &
		    targetReactionRatesTemp(store) - 1.0D0*reactionRate(i)
        End If
	  End If
	End Do	  
	Allocate(targetReactionRatesInt(1:storeCounter,1:3))
	Allocate(targetReactionRates(1:storeCounter))	
	Do i=1,storeCounter
	  targetReactionRatesInt(i,1) = targetReactionRatesIntTemp(i,1)
	  targetReactionRatesInt(i,2) = targetReactionRatesIntTemp(i,2)
	  targetReactionRatesInt(i,3) = targetReactionRatesIntTemp(i,3)
	  targetReactionRates(i) = targetReactionRatesTemp(i)
	  !print *,i,targetReactionRatesInt(i,1),targetReactionRatesInt(i,2),&
	  !targetReactionRatesInt(i,3),targetReactionRatesTemp(i)
	End Do		
!store production rates
	Allocate(productReactionRatesIntTemp(1:size(xsKey,1),1:3))
	Allocate(productReactionRatesTemp(1:size(xsKey,1)))	
	storeCounter = 0
	!Print *,"Target"
	!Do i=1,size(xsKey,1)
	!  print *,i,xsKey(i,1),xsKey(i,2),&
	!  xsKey(i,3),xsKey(i,4),xsKey(i,5),&
	!  xsKey(i,6),reactionRate(i)
	!End Do
	!Print *,"Product"
	!Do i=1,size(xsKey,1)
	!  print *,i,xsKey(i,4),xsKey(i,5),&
	!  xsKey(i,6),reactionRate(i)
	!End Do
	!Print *," "
	!Print *," "
	!Print *," "
	Do i=1,size(xsKey,1)
	  !Only store if greater than zero
	  If(reactionRate(i).ne.0.0D0)Then
	    store = 0 !store as new = yes 
	    Do j=1,size(xsKey,1)
	      If(productReactionRatesIntTemp(j,1).eq.xsKey(i,4).and.&
		    productReactionRatesIntTemp(j,2).eq.xsKey(i,5).and.&
		    productReactionRatesIntTemp(j,3).eq.xsKey(i,6))Then
	          store = j !store as new = no 
	      End If
	    End Do
        If(store.eq.0)Then
          storeCounter = storeCounter + 1
		  productReactionRatesIntTemp(storeCounter,1) = xsKey(i,4)
		  productReactionRatesIntTemp(storeCounter,2) = xsKey(i,5)
		  productReactionRatesIntTemp(storeCounter,3) = xsKey(i,6)
		  productReactionRatesTemp(storeCounter) = 1.0D0*reactionRate(i)
	    Else	
		  productReactionRatesTemp(store) = &
		    productReactionRatesTemp(store) + 1.0D0*reactionRate(i)
        End If
	  End If
	End Do	  
	Allocate(productReactionRatesInt(1:storeCounter,1:3))
	Allocate(productReactionRates(1:storeCounter))	
	Do i=1,storeCounter
	  productReactionRatesInt(i,1) = productReactionRatesIntTemp(i,1)
	  productReactionRatesInt(i,2) = productReactionRatesIntTemp(i,2)
	  productReactionRatesInt(i,3) = productReactionRatesIntTemp(i,3)
	  productReactionRates(i) = productReactionRatesTemp(i)
	  !print *,i,productReactionRatesInt(i,1),productReactionRatesInt(i,2),&
	  !productReactionRatesInt(i,3),productReactionRatesTemp(i)
	End Do	
	
	!productReactionRatesInt
	!productReactionRates
	!do i=1,size(xsKey,1)
	  

	
!close output file
	write(999,"(A1)") " "
    close(999)	
  End Subroutine isotopeProductionRates


!------------------------------------------------------------------------!
! SUBROUTINE searchXS                                                    !
! find the xs from the loaded data                                       !
!------------------------------------------------------------------------!
  
  Subroutine searchXS(energy, xs, i)
!force declaration of all variables
	Implicit None	
!declare variables
	Integer(kind=StandardInteger) :: i,j,k,startKey,endKey,key, difference    
	Real(kind=DoubleReal) :: xs, energy, tempEnergyL, tempEnergyU, factor
	Real(kind=DoubleReal) :: eA,oA,eB,oB,eM,oM
!set start and end points
	startKey = xsKey(i,7)                    !xsKey(i,7)   data row start
	endKey = xsKey(i,7) + xsKey(i,8) - 1     !xsKey(i,8)   data row length
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
	      xs = 1.0D0 * xsData(key,2) + &
		    1.0D0 * ((energy-xsData(key,1))/(xsData(key+1,1)-xsData(key,1))) * &
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
! SUBROUTINE prepareIsotopeTally                                         
!                     
!
!
!------------------------------------------------------------------------!

  Subroutine prepareIsotopeTally()
  
!force declaration of all variables
	Implicit None	
!declare variables
 	Integer(kind=StandardInteger) :: i,j,k,key
	Integer(kind=StandardInteger) :: zA,aA,mA,zB,aB,mB
	Integer(kind=StandardInteger) :: z,a,m
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: isotopeList
!open output file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	write(999,"(A1)") " "
	write(999,"(A40,F8.4)") "Prepare full isotope tally              ",ProgramTime()	
!Fill in reaction products on isotope tally
    do i=1,size(xsKey,1)
!only do if reaction rate greater than zero
	  if(reactionRate(i).gt.0)then
	    zA = xsKey(i,1)
	    aA = xsKey(i,2)
	    mA = xsKey(i,3)
		key = 590 * zA + 290 * mA + aA	
!subtract reaction rate from target isotope
		isotopeTallyActive(key,3) = isotopeTallyActive(key,3) - 1.0D0*reactionRate(i)
	    zB = xsKey(i,4)
	    aB = xsKey(i,5)
	    mB = xsKey(i,6)
		key = 590 * zB + 290 * mB + aB	
!add reaction rate to the product isotope
		isotopeTallyActive(key,3) = isotopeTallyActive(key,3) + 1.0D0*reactionRate(i)
	  endif
    enddo	
	
!save starting isotope tally to output file
	write(999,"(A1)") " "
	write(999,"(A140)") "-----------------------------------------------------------------&
	---------------------------------------------------------------------------"
	write(999,"(A39)") "Full Starting Isotope Tally"
	write(999,"(A6,A8,A4,A4,A2,A4,A4,A18,A18,A18,A21)") &
	"Key   ","Element ",&
	"Z   ","A   ","M ",& 
	"mk  ","sim ",& 
	"Half life         ","Activity          ","Reaction Rate     ",&
	"Atoms               "
	write(999,"(A140)") "-----------------------------------------------------------------&
	---------------------------------------------------------------------------"
	simIsotopes = 0
	Do j=1,size(simIsotopeKeys) 
	  key = simIsotopeKeys(j)
	  isotopeTallyActive(key,5) = isotopeTallyActive(key,4)
	  write(999,"(I6.6,A7,A1,&
	  I3.3,A1,I3.3,A1,I1.1,A1,&
	  I3.3,A1,I3.3,A1,&
	  d17.10,A1,d17.10,A1,d17.10,A1,d17.10)") &
	  key,isotopeTallyChar(key)," ",&
	  isotopeTallyInt(key,1)," ",isotopeTallyInt(key,2)," ",isotopeTallyInt(key,3)," ",&
	  isotopeTallyInt(key,4)," ",isotopeTallyInt(key,5)," ",&
	  isotopeTallyActive(key,1)," ",isotopeTallyActive(key,2)," ",&
	  isotopeTallyActive(key,3)," ",isotopeTallyActive(key,4)		
	  simIsotopes = simIsotopes + 1
    enddo
!close output file
	write(999,"(A1)") " "
    close(999)
  End Subroutine prepareIsotopeTally
  
  
  
  

  
  
  
  
!------------------------------------------------------------------------!
! SUBROUTINE runInOut                                         
! Run the in-out production calculation time steps                       
!
!
!------------------------------------------------------------------------!   
  
  Subroutine runInOut()  
!force declaration of all variables
	Implicit None	
!declare variables
 	Integer(kind=StandardInteger) :: i,j,k
	Real(kind=DoubleReal) :: parentProductionRate,tStart,tEnd,tBeamEnd
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
	Real(kind=DoubleReal) :: atoms
	Real(kind=DoubleReal) :: rrTallyChange
	Real(kind=DoubleReal) :: w
	Real(kind=SingleReal) :: simTime, simTimeLast, workingTimeStep, simTimeNext
 	Integer(kind=StandardInteger) :: key,keyT,keyP
 	Integer(kind=StandardInteger) :: zT,aT,mT,zP,aP,mP
	Integer(kind=StandardInteger) :: beamOnOff, beamTransition, resetTimeStep
!Set start sim time 0
	simTime = 0
	workingTimeStep = timeStep	
!write start of simulation to output file
      open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A65,A75)") "=================================================================",&
	  "==========================================================================="
	  write(999,"(A1)") " " 
	  write(999,"(A19)") "START of SIMULATION"
	  write(999,"(A1)") " " 
	  write(999,"(A21,E20.10,A4,E20.10)") "Simulation time from ",simTime," to ",amTime
	  write(999,"(A10,E20.10,A4,E20.10)") "Beam from ",simTime," to ",beamDuration
	  write(999,"(A1)") " " 
	  write(999,"(A65,A75)") "=================================================================",&
	  "==========================================================================="
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
	  write(999,"(A1)") " " 
!close file
      close(999)	
!set beam on or off
	If(beamDuration.gt.0)Then
	  beamOnOff = 1 !assume beam on
	Else
	  beamOnOff = 0
	End If
!flags
	beamTransition = 0
	resetTimeStep = 0
!set up temporary tally array	  
    If(Allocated(decayTempTally))Then
	  Deallocate(decayTempTally)
	End If
	Allocate(decayTempTally(1:70800))		
!loop through time steps
	do i=1,1000000000
	  if(simTime.ge.amTime)then
        exit
	  endif
!store last sim time
      simTimeLast = simTime
!Check if sim time step exceeds activity measurement time
      If(simTime.lt.amTime.and.(simTime+workingTimeStep).gt.amTime)Then
	    workingTimeStep = amTime - simTime
	  EndIf
!check if beam ends in this time step      
      If (beamOnOff.eq.1.and.beamDuration.gt.simTime.and.&
	    beamDuration.lt.(simTime+workingTimeStep).and.beamTransition.eq.0) Then
	      beamTransition = 1
		  simTimeNext = simTime+workingTimeStep !store next simTime
		  workingTimeStep = beamDuration - simTime !adjust time step
	  ElseIf(beamTransition.eq.1)Then
	    beamTransition = 0
	    beamOnOff = 0  !switch off beam
		workingTimeStep = simTimeNext - simTime
		resetTimeStep = 1
!zero out reaction rates (beam off)
		Do j=1,size(simIsotopeKeys) 
		  key = simIsotopeKeys(j)
		  isotopeTallyActive(key,3) = 0.0D0
		End Do
	  End If
!store sim info to file - open output file
      open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	  write(999,"(A70)") "----------------------------------------------------------------------" 
	  write(999,"(A21,E20.10,A4,E20.10)") "Simulation time from ",simTime," to ",(simTime+workingTimeStep)
	  write(999,"(A70)") "----------------------------------------------------------------------"
!close file
      close(999)
!fill sim isotope slots with zeros
	  Do j=1,size(simIsotopeKeys,1) 
		key = simIsotopeKeys(j)
	    decayTempTally(key) = 0.0D0
	  End Do
!target atoms lost due to transmutation
	  If(beamOnOff.eq.1)Then
	    Do j=1,size(targetReactionRatesInt,1) 
		  zT = targetReactionRatesInt(j,1)
	      aT = targetReactionRatesInt(j,2)
	      mT = targetReactionRatesInt(j,3)
		  keyT=makeIsotopeKey(zT,aT,mT)
		  decayTempTally(keyT) = decayTempTally(keyT) + targetReactionRates(j)
		End Do  
	  End If
!Product atoms lost
	  Do j=1,size(productReactionRatesInt,1) 
		zP = productReactionRatesInt(j,1)
	    aP = productReactionRatesInt(j,2)
	    mP = productReactionRatesInt(j,3)
		w = 1.0D0*productReactionRates(j)
        If(beamOnOff.eq.1)Then
	      Call tallyDecay&
		    (zP,aP,mP,w,1.0D0*simTime,&
		    1.0D0*(simTime+workingTimeStep),1.0D0*beamDuration)	  
	    Else
	  	  Call tallyDecay&
		    (zP,aP,mP,0.0D0,1.0D0*simTime,&
		    1.0D0*(simTime+workingTimeStep),1.0D0*beamDuration)
	    End If
	  End Do
!Apply tally changes
      Do j=1,size(simIsotopeKeys) 
		key = simIsotopeKeys(j)
	    isotopeTallyActive(key,4) = isotopeTallyActive(key,4) + 1.0D0 * decayTempTally(key)
!ensure no negative tally counts
        If(isotopeTallyActive(key,4).lt.0.0D0)Then
		  isotopeTallyActive(key,4) = 0.0D0
		End If
	  End Do
!Calculate activity
      Call calcTallyActivity()
!Save activity/time to file - open output file
      open(unit=998,file=trim(activityHistoryFile),status="old",position="append",action="write")	  
	  write(998,"(E20.10,A1,E20.10,A1,E20.10)") simTime," ",&
	  (simTime+workingTimeStep)," ",totalActivityAtTime
	  close(998)
!output tally to output file
	  Call outputTally()
!increment time step
	  simTime = simTimeLast + workingTimeStep
!reset time step
      If (resetTimeStep.eq.1) Then
	    workingTimeStep = timeStep
		resetTimeStep = 0
	  End If
	enddo
!close file
    close(997)
		
		
	

  End Subroutine runInOut
  
   
!------------------------------------------------------------------------!
! SUBROUTINE tallyDecay                                      
!                      
!
!
!------------------------------------------------------------------------!    
  
  Subroutine tallyDecay(zP,aP,mP,parentProductionRate,tStart,tEnd,tBeamEnd)
!force declaration of all variables
	Implicit None	
!declare variables
 	Integer(kind=StandardInteger) :: i,j,k
	Real(kind=DoubleReal) :: parentProductionRate,tStart,tEnd,tBeamEnd
 	Integer(kind=StandardInteger) :: keyP, key
 	Integer(kind=StandardInteger) :: zP,aP,mP
	Real(kind=DoubleReal),Dimension(:,:,:),Allocatable::decayChainArray	
	Real(kind=DoubleReal),Dimension(:,:),Allocatable::decayDataArray
	!Real(kind=DoubleReal),Dimension(:,:),Allocatable::newTally		
	Real(kind=DoubleReal),Dimension(:,:),Allocatable::isotopeChange
	Integer(kind=StandardInteger) :: arrayHeight, arrayWidth, chainLength
		
!Parent/product key
	keyP = makeIsotopeKey (zP, aP, mP)			
!make isotope decay chain array
	decayChainArray = IsotopesDecayChains(zP,aP,mP)	
!for each chain, run decay and tally
    arrayHeight = size(decayChainArray,1)	!i, max chain length
    arrayWidth = size(decayChainArray,2)    !j, number of decay chains
!loop each chain
    Do j=1,arrayWidth
!Allocate array
      If(Allocated(decayDataArray))Then
	    Deallocate(decayDataArray)
	  End If
!loop through each isotope
      chainLength = 0
	  Do i=1,arrayHeight
	    If(decayChainArray(i,j,1).ne.-1)Then
		  chainLength = chainLength + 1
		Else
		  Exit
		End If
	  End Do
	  If(chainLength.gt.0)Then
!Allocate data array
	    Allocate(decayDataArray(1:chainLength,1:6))
!Store data to pass to activity calc function - build decay data array
        Do i=1,chainLength		 
		  key = decayChainArray(i,j,1) 
		  decayDataArray(i,1) = key 						!Tally key
		  decayDataArray(i,2) = isotopeTallyActive(key,4)   !No. Atoms
		  decayDataArray(i,3) = decayChainArray(i,j,2)      !Half life
		  decayDataArray(i,4) = decayChainArray(i,j,3)      !branching factor
		  decayDataArray(i,5) = isotopeTallyInt(key,1)		!isotope Z
		  decayDataArray(i,6) = isotopeTallyInt(key,2)		!isotope A	
          !print *,isotopeTallyInt(key,1),isotopeTallyInt(key,2),key,isotopeTallyActive(key,4)	  
		End Do
!calculate change in isotope amounts
		isotopeChange=CalcIsotopeAmount(parentProductionRate,decayDataArray,tStart,tEnd,tBeamEnd)
!apply changes to temp tally
	    Do i=1,size(isotopeChange,1)
	      key = isotopeChange(i,1)
		  If(key.gt.0.and.key.le.70800)Then
		    !print *,key,isotopeChange(i,2)
		    decayTempTally(key) = decayTempTally(key) + isotopeChange(i,2)
		  End If
	    End Do
	  End If
	  
	End Do

	
  End Subroutine tallyDecay

  
!------------------------------------------------------------------------!
! SUBROUTINE calcTallyActivity                                      
!                      
!
!
!------------------------------------------------------------------------! 
	
  Subroutine calcTallyActivity
!force declaration of all variables
	Implicit None	
!declare variables
 	Integer(kind=StandardInteger) :: i,key
!calculate activities
    totalActivityAtTime = 0.0D0
	Do i=1,size(simIsotopeKeys)
	  !isotopeTallyActive(i,6) decay constant
	  !isotopeTallyActive(i,4) number of atoms
	  !isotopeTallyActive(i,2) activity
	  key = simIsotopeKeys(i)
	  isotopeTallyActive(key,2) = isotopeTallyActive(key,6) *&
        isotopeTallyActive(key,4) 
	  totalActivityAtTime = totalActivityAtTime + isotopeTallyActive(key,2)	
    End Do
  End Subroutine calcTallyActivity
  
  
  
  
	
!------------------------------------------------------------------------!
! SUBROUTINE outputTally                                      
!                      
!
!
!------------------------------------------------------------------------! 
	
  Subroutine outputTally
!force declaration of all variables
	Implicit None	
!declare variables
 	Integer(kind=StandardInteger) :: i,j,k,key
!open output file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	
  !save starting isotope tally to output file
	write(999,"(A1)") " "
	write(999,"(A140)") "-----------------------------------------------------------------&
	---------------------------------------------------------------------------"
	write(999,"(A15)") "Isotope Tally  "
	write(999,"(A6,A8,A4,A4,A2,A4,A4,A18,A18,A18,A18,A18,A18)") &
	"Key   ","Element ",&
	"Z   ","A   ","M ",& 
	"mk  ","sim ",& 
	"Half life         ","Decay Constant    ","Activity          ","Reaction Rate     ",&
	"Start Atoms       ","End Atoms         "
	write(999,"(A140)") "-----------------------------------------------------------------&
	---------------------------------------------------------------------------"
	Do j=1,size(simIsotopeKeys)
	  key = simIsotopeKeys(j)
	  If(isotopeTallyActive(key,4).gt.0.0D0)Then
		i = key
	    write(999,"(I6.6,A7,A1,&
		I3.3,A1,I3.3,A1,I1.1,A1,&
		I3.3,A1,I3.3,A1,&
		d17.10,A1,d17.10,A1,d17.10,A1,d17.10,A1,d17.10,A1,d17.10)") &
		key,isotopeTallyChar(i)," ",&
		isotopeTallyInt(i,1)," ",isotopeTallyInt(i,2)," ",isotopeTallyInt(i,3)," ",&
		isotopeTallyInt(i,4)," ",isotopeTallyInt(i,5)," ",&
		isotopeTallyActive(i,1)," ",isotopeTallyActive(i,6)," ",isotopeTallyActive(i,2)," ",&
		isotopeTallyActive(i,3)," ",isotopeTallyActive(i,5)," ",isotopeTallyActive(i,4)		
	  End If
    End Do
  !close output file
	write(999,"(A1)") " "
    close(999)	
	
  End Subroutine outputTally
	
	
	
	
	
	
	
	
	
!------------------------------------------------------------------------!
! Decay Functions
!  
!------------------------------------------------------------------------!	
  
  Function IsotopesInDecayTree(z,a,m) RESULT (isotopeArray)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,z,a,m,complete
	Integer(kind=StandardInteger) :: zP, aP, mP
	Integer(kind=StandardInteger) :: zC, aC, mC
	Real(kind=DoubleReal), Dimension(:,:), &
    Allocatable :: isotopeArray
!Allocate an array to store isotope list	
    If(Allocated(decayChainIsotopesArray))Then
	  Deallocate(decayChainIsotopesArray)
	End If
	Allocate(decayChainIsotopesArray(1:300,1:3))
!fill with blank data
    do i=1,300
	  do j=1,3
	    decayChainIsotopesArray(i,j) = 0
	  enddo	
	enddo
!set counter = 1 and store parent isotope
    decayChainIsotopesArray(1,1) = z
    decayChainIsotopesArray(1,2) = a
    decayChainIsotopesArray(1,3) = m
    isotopeCounter = 1
!call function recursively
    complete = IsotopesInDecayTreeR(z,a,m)	
!transfer data to new array	
	Allocate(isotopeArray(1:isotopeCounter,1:3))
    do i=1,isotopeCounter
	  do j=1,3
	    isotopeArray(i,j) = decayChainIsotopesArray(i,j)
	  enddo
	enddo
!Deallocate array
	If(Allocated(decayChainIsotopesArray))Then
	  Deallocate(decayChainIsotopesArray)
	End If
  End Function IsotopesInDecayTree
!---------------------------------------------------------------------
  Recursive Function IsotopesInDecayTreeR(z,a,m) RESULT (lastRow)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,z,a,m,lastRow,nextRow
	Integer(kind=StandardInteger) :: zP, aP, mP
	Integer(kind=StandardInteger) :: zC, aC, mC
	Logical :: store
!loop through decay chains
	do i=1,size(decayChar)
	  zP = decayInt(i,2)
	  aP = decayInt(i,1)
	  mP = decayInt(i,3)
	  if(z.eq.zP.and.a.eq.aP.and.m.eq.mP)then
	    zC = decayInt(i,5)
	    aC = decayInt(i,4)
	    mC = decayInt(i,6)
!check not in array already
        store = .true.
		do j=1,size(decayChainIsotopesArray,1)
		  if(decayChainIsotopesArray(j,1).eq.zC.and.&
		  decayChainIsotopesArray(j,2).eq.aC.and.&
		  decayChainIsotopesArray(j,3).eq.mC)then
		    store = .false.
		  endif
		enddo
		if(store.eqv..true.)then
		  isotopeCounter = isotopeCounter + 1
		  decayChainIsotopesArray(isotopeCounter,1) = zC
		  decayChainIsotopesArray(isotopeCounter,2) = aC
		  decayChainIsotopesArray(isotopeCounter,3) = mC
	      !print *,isotopeCounter,zP,aP,mP,zC,aC,mC
		endif
		lastRow = IsotopesInDecayTreeR(zC,aC,mC)
	  endif	  
	enddo  
  End Function IsotopesInDecayTreeR
  
  
  
  
  
  Function IsotopesDecayChains(z,a,m) RESULT (decayChains)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: z, a, m
	Integer(kind=StandardInteger) :: nodes
	Integer(kind=StandardInteger) :: tempArraySize
	Real(kind=DoubleReal),Dimension(:,:,:),Allocatable::decayChains
	Integer(kind=StandardInteger) :: arrayHeight, arrayWidth
	Integer(kind=StandardInteger) :: keyP, keyC, key
	Logical :: arrayEnd
!make temporary square array 50x50
    tempArraySize = 25
    If(Allocated(public3DTempArray))Then
	  Deallocate(public3DTempArray)
	EndIf
    Allocate(public3DTempArray(1:tempArraySize,1:tempArraySize,1:3))
!fill with blank markers (-1)
    do i=1,tempArraySize
      do j=1,tempArraySize
	    public3DTempArray(i,j,1) = -1
	  enddo
	enddo	
!run recursive function to fill public3DTempArray
	isotopeCounter = 0
	publicX = 0
	publicY = 1
    nodes = IsotopesDecayChainsR(z,a,m)	
!measure public3DTempArray
    arrayHeight = 0
    do i=1,tempArraySize
	  arrayEnd = .true.
      do j=1,tempArraySize
	    If(public3DTempArray(i,j,1).ne.-1)Then
		  arrayEnd = .false.
		End If		
	  enddo
	  If(arrayEnd.eqv..true.)Then
	    arrayHeight = i - 1
		Exit
	  End If
	enddo
	arrayWidth = 0
    do j=1,tempArraySize
	  arrayEnd = .true.
      do i=1,tempArraySize
	    If(public3DTempArray(i,j,1).ne.-1)Then
		  arrayEnd = .false.
		End If		
	  enddo
	  If(arrayEnd.eqv..true.)Then
	    arrayWidth = j - 1
		Exit
	  End If
	enddo
!Allocate decay chains array
    Allocate(decayChains(1:arrayHeight,1:arrayWidth,1:3))	
!Transfer into new array of correct size
    do i=1,arrayHeight
	  arrayEnd = .true.
      do j=1,arrayWidth
	    decayChains(i,j,1) = public3DTempArray(i,j,1)	
	  enddo
	enddo
!Fill in blanks
    Do j=1,arrayWidth
	  Do i=1,arrayHeight
	    If(decayChains(i,j,1).ne.-1)Then
		  Exit
		Else  
		  If(j.gt.1)Then
		    decayChains(i,j,1) = decayChains(i,j-1,1)
		  End If		  
		End If
	  End Do
	End Do    
!Fill in half lives
    Do j=1,arrayWidth
	  Do i=1,arrayHeight
	    !Do k=1,size(decayInt,1)
		!  If(decayChains(i,j,1).eq.decayInt(k,7))Then
		!    decayChains(i,j,2) = decayDouble(k,2)
		!	Exit
		!  End If
	    !End Do		
		key = decayChains(i,j,1)
		decayChains(i,j,2) = isotopeTallyActive(key,1)
	  End Do
	End Do
!Fill in branching factors
    Do j=1,arrayWidth
	  Do i=1,arrayHeight
	    If(i.eq.1)Then
		  decayChains(i,j,3) = 0
		Else  
	      Do k=1,size(decayInt,1)
		    keyP = decayChains(i-1,j,1)
		    keyC = decayChains(i,j,1)
		    If(keyP.eq.decayInt(k,7).and.&
			keyC.eq.decayInt(k,8))Then
		      decayChains(i,j,3) = decayDouble(k,1)
			  Exit
		    End If
	      End Do		  
		End If
	  End Do
	End Do
	
  End Function IsotopesDecayChains
!-----------------------------------------------------------------
  Recursive Function IsotopesDecayChainsR(z,a,m) RESULT (nodes)
!force declaration of all variables
	Implicit None
!declare variables  
	Integer(kind=StandardInteger) :: nodes, lastIsotope
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: z, a, m
	Integer(kind=StandardInteger) :: zP, aP, mP
	Integer(kind=StandardInteger) :: zC, aC, mC
	Integer(kind=StandardInteger) :: key
	Logical :: unstable
!loop through decay chains
	unstable=.false.
	do i=1,size(decayChar)
	  zP = decayInt(i,2)
	  aP = decayInt(i,1)
	  mP = decayInt(i,3)
	  if(z.eq.zP.and.a.eq.aP.and.m.eq.mP)then
	    unstable=.true.
	    zC = decayInt(i,5)
	    aC = decayInt(i,4)
	    mC = decayInt(i,6)
        isotopeCounter = isotopeCounter + 1
		publicX = publicX + 1
		key = makeIsotopeKey(zP,aP,mP)
		public3DTempArray(publicX,publicY,1) = key
!store half life
        !public3DTempArray(publicX,publicY,2) = decayDouble(i,2)		
		publicZ = publicX + 1
		nodes = IsotopesDecayChainsR(zC,aC,mC)
		publicX = publicX - 1
	  endif	  
	enddo  
	if(unstable.eqv..false.)then
	  key = makeIsotopeKey(z,a,m)
	  public3DTempArray(publicZ,publicY,1) = key
	  publicY = publicY + 1
	endif  
  End Function IsotopesDecayChainsR
  
  
    
  Function BranchProductionLoss&
  (parentProductionRate,decayDataArray,tStart,tEnd,tBeamEnd) RESULT (branchTally)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,j,k,decaySteps
	Real(kind=DoubleReal) :: parentProductionRate,tStart,tEnd,tBeamEnd
	Real(kind=DoubleReal) :: timeBPL,timeStepBPL
	Real(kind=DoubleReal) :: lastLoss, numberOfAtoms, sourceAtoms, lossAtoms
	Real(kind=DoubleReal) :: lnTwo
	Real(kind=DoubleReal) :: atoms
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: branchTally  
	Integer(kind=StandardInteger) :: iterations    
    !decayDataArray(1,1) = Key			!Key
  	!decayDataArray(1,2) = 100			!Np(tStart)
	!decayDataArray(1,3) = 60     	    !T1/2p
	!decayDataArray(1,4) = 0	    		!
    !decayDataArray(2,1) = Key			!Key
	!decayDataArray(2,2) = 500			!Na(tStart)
	!decayDataArray(2,3) = 120	        !T1/2a
	!decayDataArray(2,4) = 0.8	    	!BF P to A
!set variables
    lnTwo = 0.693147180559945
!set iterations  
    iterations = 100
!decay steps
    decaySteps = size(decayDataArray,1)
!allocate branch tally array
    Allocate(branchTally(1:decaySteps,1:2))
!fill tally array
    do i=1,decaySteps
	  branchTally(i,1) = decayDataArray(i,1)	!Key
	  branchTally(i,2) = decayDataArray(i,2)	!No. atoms
	  !print *,branchTally(i,1),branchTally(i,2)
	enddo
!set time step
    timeStepBPL = (tEnd-tStart)/iterations
!set time	
	timeBPL = tStart
!iterate from start to end time
	do i=1,iterations
	  do j=1,decaySteps
		sourceAtoms = 0.0D0
		lossAtoms = 0.0D0 
!source parent
	    If(j.eq.1)Then  !parent
		  If(timeBPL.lt.tBeamEnd)Then
	        sourceAtoms = timeStepBPL * parentProductionRate
		  End If	
		End If
!source from previous isotope in decay chain
		if(j.gt.1)then
		  sourceAtoms = sourceAtoms + decayDataArray(j,4)*lastLoss !add bf * parent loss to source
		endif
        lastLoss = 0.0D0
!loss from decay
	    if(decayDataArray(1,3).gt.0)then
		  lossAtoms = branchTally(j,2)*timeStepBPL*(lnTwo/decayDataArray(j,3))
		endif
!save last loss amount
		lastLoss = lossAtoms
		if(lastLoss.gt.branchTally(j,2))then
		  lastLoss = branchTally(j,2)
		endif
!tally with source and loss atoms
        branchTally(j,2) = branchTally(j,2) + sourceAtoms - lossAtoms		
!check greater than or equal to 0
        if(branchTally(j,2).lt.0)then
		  branchTally(j,2) = 0.0D0
		endif
		
	  enddo
	  timeBPL = timeBPL + timeStepBPL
	enddo		
  !print
    print *,"Prod",parentProductionRate
    print *,"tStart",tStart
    print *,"tEnd",tEnd
    do i=1,decaySteps
	  print *,"S",i,branchTally(i,1),branchTally(i,2),decayDataArray(i,3),decayDataArray(i,4)
	enddo
  
  End Function BranchProductionLoss  
  
  
  
  
!------------------------------------------------------------------------!
! SUBROUTINE searchXSa                             
! Might need to develop, add 3 point interpolation
!------------------------------------------------------------------------!
  
  Function searchXSa(energy, i) RESULT (xs)
!force declaration of all variables
	Implicit None	
!declare variables
	Integer(kind=StandardInteger) :: i,j,k,startKey,endKey,key, difference    
	Real(kind=DoubleReal) :: xs, energy, tempEnergyL, tempEnergyU, factor
	Real(kind=DoubleReal) :: eA,oA,eB,oB,eM,oM
!set start and end points
	startKey = xsKey(i,7)                    !xsKey(i,7)   data row start
	endKey = xsKey(i,7) + xsKey(i,8) - 1     !xsKey(i,8)   data row length
	factor = 0.5
	difference = ceiling(factor*(endKey - startKey))
	print *,xsData(startKey,1),xsData(endKey,1)
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
	      xs = 1.0D0 * xsData(key,2) + &
		    1.0D0 * ((energy-xsData(key,1))/(xsData(key+1,1)-xsData(key,1))) * &
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
	
  End Function searchXSa
  
  
End Module productionLoss