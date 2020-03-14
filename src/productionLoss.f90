Module productionLoss

! Setup Modules
  Use kinds
  Use constants
  Use activityFunctions
  Use strings
  Use msubs
  Use input
  Use prep
  Use globals
  Use initialise        ! input

! force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! declare global variables

! Privacy of functions/subroutines/variables
  Private
  Public :: runProductionLoss      !Subroutine

! Module Subroutines
  contains

! ------------------------------------------------------------------------!
! SUBROUTINE runProductionLoss                                           !
! ------------------------------------------------------------------------!
  Subroutine runProductionLoss()
! Calculate the production of all isotopes and loss of radioactive isotopes
! Production is assumed to be time independent
! Loss is assumed to vary with time/amount of radioactive isotope
! force declaration of all variables
    Implicit None
! declare variables
    Call M_synchProcesses()
! Run subroutines
    If(mpiProcessID.eq.0)Then
      Call preCalcOutput()    !calculate the production rate
      Call isotopeProductionRates()    !calculate the production rate
      Call prepareIsotopeTally()
! Call makeReducedDecayList()
      Call runInOut()
      Call finishProductionLoss()
    End If
    Call M_synchProcesses()
  End Subroutine runProductionLoss

! ------------------------------------------------------------------------!
! SUBROUTINE
!
! ------------------------------------------------------------------------!
  Subroutine preCalcOutput()
! force declaration of all variables
    Implicit None
! declare variables
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
      write(999,"(A1)") " "
      write(999,"(A70)") "----------------------------------------------------------------------"
      write(999,"(A22)") "Pre-Simulation Summary"
      write(999,"(F8.4)") ProgramTime()
      write(999,"(A70)") "----------------------------------------------------------------------"
      write(999,"(A1)") " "
      write(999,"(A30,E20.10)") "Beam flux ions/s:             ",((beamFlux*1D-6)/elementaryCharge)
! close output file
      write(999,"(A1)") " "
      close(999)
    End If
  End Subroutine preCalcOutput

! ------------------------------------------------------------------------!
! SUBROUTINE isotopeProductionRates                                      !
! calculate the isotope production rates                                 !
! ------------------------------------------------------------------------!
  Subroutine isotopeProductionRates()
! Calculate the production of all isotopes and loss of radioactive isotopes
! Production is assumed to be time independent
! Loss is assumed to vary with time/amount of radioactive isotope
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,k,n,z,a,m,key
    Integer(kind=StandardInteger) :: zA,aA,mA,zB,aB,mB
    Integer(kind=StandardInteger) :: store, storeCounter
    Real(kind=DoubleReal) :: volumeAffected
    Real(kind=DoubleReal) :: x,y
    Real(kind=DoubleReal) :: ionsPerSecond
    Real(kind=DoubleReal) :: contentFactor
    Real(kind=DoubleReal) :: trajectoryDepth, energyAtDepth, segmentLength, xs
    Real(kind=DoubleReal) :: affectedDepth,tempDoubleA,tempDoubleB
    Real(kind=DoubleReal) :: averagedXS
    Real(kind=DoubleReal) :: DPA
    Real(kind=DoubleReal) :: reactionRateTemp
    Character(len=255) :: fittingPolynomial
    Integer(kind=StandardInteger), Dimension( : , :), Allocatable :: targetReactionRatesIntTemp
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: targetReactionRatesTemp
    Integer(kind=StandardInteger), Dimension( : , :), Allocatable :: productReactionRatesIntTemp
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: productReactionRatesTemp
    Real(kind=DoubleReal), Dimension(1:polyFitOrder+1) :: coefficients
!  Overview
!
!  1. calculate average depth of ion
!  2. calculate reaction rate for each target-product combination
! print if verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A34,F8.4)","    Calc isotope production rates ",ProgramTime()
    End If
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
      write(999,"(A1)") " "
      write(999,"(A51,F8.4)") "Calculate projectile-target-product reaction rates ",ProgramTime()
      write(999,"(A29,I6)") "Reaction rates to calculate: ",size(xsKey,1)
    End If
! Allocate Arrays
    Allocate(reactionRate(1:size(xsKey,1)))
    Allocate(averageXS(1:size(xsKey,1)))
    Allocate(numberDensityArray(1:size(xsKey,1)))
! Transfer array
    Do i=1,polyFitOrder+1
      coefficients(i) = fitCoefficients(i-1)
    End Do
! Production rate
! solve equation = 0 for the trajectory
    x = 0
    y = 0
    If(verboseTerminal.eqv..true.)Then
      print "(A26,F8.4)","    Trajectory equation:  ",ProgramTime()
      print "(A12)","    E(r) = "
      Do i=0,size(fitCoefficients)-1
        If(i.eq.0)Then
          print "(A4,E20.10,A2,I2)","    ",fitCoefficients(i),"r^",i
        Else
          print "(A4,E20.10,A3,I2)","    ",fitCoefficients(i),"+r^",i
        End If
      End Do
    End If
    fittingPolynomial = "y = "
    If(mpiProcessID.eq.0)Then
      write(999,"(A25)") "Trajectory fit equation: "
    End If
    Do i=0,size(fitCoefficients)-1
      y = y + x**(i) * fitCoefficients(i)
      If(mpiProcessID.eq.0)Then
        write(999,"(I8,E20.10)") i,fitCoefficients(i)
      End If
    End Do
! affected depth
    If(trajDepth.gt.targetThickness)Then
      affectedDepth = targetThickness
    Else
      affectedDepth = trajDepth
    End If
! set global variable
    affectedMaterialDepth = affectedDepth
    If(mpiProcessID.eq.0)Then
      write(999,"(A25,E20.10)") "Affected material depth: ",affectedMaterialDepth
    End If
    If(verboseTerminal.eqv..true.)Then
      print "(A4,A25,E20.10,F8.4)","    ","Affected material depth: ",&
      affectedMaterialDepth,ProgramTime()
    End If
! ions per second
    ionsPerSecond = (beamFlux * 1.0D-6) / elementaryCharge  !convert flux in uA to ions per second
! volume of material affected by beam
    volumeAffected = 1.0D0 * beamArea * 1E-6 * affectedDepth * 1E-10
! update tally to show atoms in affected volume only
    totalSimulationAtoms = 0.0D0
    Do i=1,size(simIsotopeKeys)
      key = simIsotopeKeys(i)
      If(key.gt.0.and.key.le.70800)Then
        totalSimulationAtoms = totalSimulationAtoms + isotopeTallyActive(key,4)
      End If
    End Do
    If(verboseTerminal.eqv..true.)Then
      print "(A4,A18,E20.10,F8.4)","    ","Total atom tally: ",&
      totalSimulationAtoms,ProgramTime()
    End If
! Loop over target-product combinations
    If(verboseTerminal.eqv..true.)Then
      print "(A4,A24,F8.4)","    ","Calculate reaction rates",ProgramTime()
    End If
    i = 0
    Do n=1,size(xsKey,1)
! Get key and data for cross section
      z = xsKey(n,1)
      a = xsKey(n,2)
      m = xsKey(n,3)
! key = 590 * z + 290 * m + a
      key=makeIsotopeKey(z,a,m)
      If(key.gt.0.and.key.le.70800)Then
        tempDoubleA = 1.0D0 * isotopeTallyActive(key,4)
        tempDoubleB = 1.0D0 * totalSimulationAtoms
        contentFactor = tempDoubleA / (tempDoubleB*1.0D0)
! split up trajectory function
        segmentLength = affectedDepth / integrationGranularity
        averagedXS = 0.0D0
        Do j=1,integrationGranularity
! calculate trajectory depth
          trajectoryDepth = (affectedDepth * (j - 0.5)) / integrationGranularity
! energy at depth
          energyAtDepth = 0
          Do k=0,size(fitCoefficients)-1
            energyAtDepth = energyAtDepth + trajectoryDepth**(k) * fitCoefficients(k)
          End Do
! find reaction cross section for energy - sets xs
          xs = 0.0D0 !make cross section equal to zero by default
          Call searchXS(energyAtDepth, xs, n)
! print *,energyAtDepth,xs,n
! add to averaged cross section
          averagedXS = averagedXS + (1.0D0/(1.0D0*integrationGranularity)) * xs
        End Do
! convert averaged cross section to m-2
        averagedXS = averagedXS * 1.0D-28
! calculate reaction rate
        reactionRateTemp = ionsPerSecond * averagedXS * numberDensity * &
        affectedDepth * 1.0D-10 * contentFactor
! print *,averagedXS,reactionRateTemp
! Store values
        i = i + 1
        reactionRate(i) = 1.0D0*reactionRateTemp
        averageXS(i) = 1.0D0*averagedXS
        numberDensityArray(i) = 1.0D0 * numberDensity * contentFactor
      End If
    End Do
! DPA
    If(vpi.gt.0.0D0)Then
      DPA = 1.0D0*(vpi*ionsPerSecond*beamDuration)/(volumeAffected*numberDensity)
    End If
! change beamDuration and amtime to meet target DPA
    If(vpi.gt.0.0D0.and.targetDPA.gt.0.0D0)Then
      beamDuration = 1.0D0*beamDuration*(targetDPA/DPA)
      If(amTime.lt.beamDuration)Then
        amTime = beamDuration
      End If
! print if verbose on
      If(verboseTerminal.eqv..true.)Then
        print "(A18,F20.10,A6,F8.4)","    Target DPA    ",targetDPA,"      ",ProgramTime()
        print "(A18,F20.10,A6,F8.4)","    Beam Duration ",beamDuration,"      ",ProgramTime()
      End If
! recalculate DPA
      DPA = 1.0D0*(vpi*ionsPerSecond*beamDuration)/(volumeAffected*numberDensity)
    End If
! DPA
    If(vpi.gt.0.0D0)Then
! print if verbose on
      If(verboseTerminal.eqv..true.)Then
        print "(A8,F20.10,A6,F8.4)","    VPI ",vpi,"      ",ProgramTime()
        print "(A8,F20.10,A6,F8.4)","    DPA ",DPA,"      ",ProgramTime()
      End If
    End If
! output variables to file
    If(mpiProcessID.eq.0)Then
      write(999,"(A1)") " "
      write(999,"(A70)") "----------------------------------------------------------------------"
      write(999,"(A20)") "Calculated variables"
      write(999,"(A70)") "----------------------------------------------------------------------"
      write(999,"(A30,E20.10)") "Beam flux ions/s:             ",ionsPerSecond
      write(999,"(A30,E20.10)") "Affected depth/m:             ",(affectedDepth*1.0D-10)
      write(999,"(A30,E20.10)") "Affected depth/ang:           ",(affectedDepth)
      write(999,"(A30,E20.10)") "Affected volume/m3:           ",volumeAffected
      write(999,"(A30,E20.10)") "Number Density:               ",numberDensity
      If(vpi.gt.0)Then
        write(999,"(A30,E20.10)") "VPI:                          ",vpi
        write(999,"(A30,E20.10)") "DPA:                          ",DPA
      End If
      If(vpi.gt.0.0D0.and.targetDPA.gt.0.0D0)Then
        write(999,"(A1)") " "
        write(999,"(A70)") "----------------------------------------------------------------------"
        write(999,"(A20)") "Adjusted variables"
        write(999,"(A70)") "----------------------------------------------------------------------"
        write(999,"(A30)") "Beam Duration:                ",beamDuration
        write(999,"(A30)") "Activity Measurement Time:    ",amTime
        write(999,"(A1)") " "
      End If
! output reaction rates to output data file
      write(999,"(A1)") " "
      write(999,"(A70)") "----------------------------------------------------------------------"
      write(999,"(A38,F8.4)") "Reaction rates in the target material ",ProgramTime()
      write(999,"(A59,A40)") "       Sym Z   A   M   Sym Z   A   M   Reaction Rate       ",&
      "Average XS          Number Density      "
      write(999,"(A70)") "----------------------------------------------------------------------"
    End If
    Do i=1,size(xsKey,1)
! only do if reaction rate greater than zero
      If(reactionRate(i).gt.0)Then
        zA = xsKey(i,1)
        aA = xsKey(i,2)
        mA = xsKey(i,3)
        zB = xsKey(i,4)
        aB = xsKey(i,5)
        mB = xsKey(i,6)
        If(mpiProcessID.eq.0)Then
          write(999,"(I4,A3,A4,I4,I4,I4,A4,I4,I4,I4,E20.10,E20.10,E20.10)") i,":  ",&
          elementSymbol(zA),zA,aA,mA,elementSymbol(zB),zB,aB,mB,&
          (1.0D0*reactionRate(i)),averageXS(i),numberDensityArray(i)
        End If
      End If
    End Do
! store target rates
    Allocate(targetReactionRatesIntTemp(1:size(xsKey,1),1:3))
    Allocate(targetReactionRatesTemp(1:size(xsKey,1)))
    storeCounter = 0
    Do i=1,size(xsKey,1)
! Only store if greater than zero
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
    End Do
! store production rates
    Allocate(productReactionRatesIntTemp(1:size(xsKey,1),1:3))
    Allocate(productReactionRatesTemp(1:size(xsKey,1)))
    storeCounter = 0
    Do i=1,size(xsKey,1)
! Only store if greater than zero
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
    End Do
! close output file
    If(mpiProcessID.eq.0)Then
      write(999,"(A1)") " "
      close(999)
    End If
  End Subroutine isotopeProductionRates

! ------------------------------------------------------------------------!
! SUBROUTINE searchXS                                                    !
! find the xs from the loaded data                                       !
! ------------------------------------------------------------------------!

  Subroutine searchXS(energy, xs, i)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,startKey,endKey,key, difference
    Real(kind=DoubleReal) :: xs, energy, tempEnergyL, tempEnergyU, factor
    Real(kind=DoubleReal) :: eA,oA,eB,oB,eM,oM
! set start and end points
    startKey = xsKey(i,7)                    !xsKey(i,7)   data row start
    endKey = xsKey(i,7) + xsKey(i,8) - 1     !xsKey(i,8)   data row length
    factor = 0.5
    difference = ceiling(factor*(endKey - startKey))
    If(energy.lt.xsData(startKey,1))Then
      key = startKey
      eA = 0
      oA = 0
      eB = xsData(key,1)
      oB = xsData(key,2)
      eM = energy
      oM = oA + ((eM-eA)/(eB-eA))*(oB-oA)
      xs = oM
    ElseIf(energy.gt.xsData(endKey,1))Then
      key = endKey
      eA = xsData(key,1)
      oA = xsData(key,2)
      eB = 2*energy
      oB = 0
      eM = energy
      oM = oA + ((eM-eA)/(eB-eA))*(oB-oA)
      xs = oM
    Else
! start point
      key = startKey + difference
      Do j=1,10
! adjust key if too high or too low
        If(key.ge.endKey)Then
          key = endKey-1
        End If
        If(key.lt.startKey)Then
          key = startKey
        End If
        tempEnergyL = xsData(key,1)
        tempEnergyU = xsData(key+1,1)
        If(energy.ge.tempEnergyL.and.energy.le.tempEnergyU)Then
! energy bound found
! Linear interpolate between start and end energy
          xs = 1.0D0 * xsData(key,2) + &
          1.0D0 * ((energy-xsData(key,1))/(xsData(key+1,1)-xsData(key,1))) * &
          (xsData(key+1,2)-xsData(key,2))
          exit !break out of loop
        Else
          If(energy.lt.tempEnergyL)Then
            difference = ceiling(factor * difference)
            key = key - difference
! print *,"Decrease"
          End If
          If(energy.gt.tempEnergyU)Then
            difference = ceiling(factor * difference)
            key = key + difference
! print *,"Decrease"
          End If
        End If
      End Do
    End If
  End Subroutine searchXS

! ------------------------------------------------------------------------!
! SUBROUTINE prepareIsotopeTally
!
!
!
! ------------------------------------------------------------------------!

  Subroutine prepareIsotopeTally()
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,key
    Integer(kind=StandardInteger) :: zA,aA,mA,zB,aB,mB
    Integer(kind=StandardInteger) :: z,a,m
! print if verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A41,F8.4)","    Prepare isotope tally for simulation ",ProgramTime()
    End If
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
      write(999,"(A1)") " "
      write(999,"(A40,F8.4)") "Prepare full isotope tally              ",ProgramTime()
    End If
! Fill in reaction products on isotope tally
    Do i=1,size(xsKey,1)
! only do if reaction rate greater than zero
      If(reactionRate(i).gt.0)Then
        zA = xsKey(i,1)
        aA = xsKey(i,2)
        mA = xsKey(i,3)
        key = makeIsotopeKey(z,a,m)
! subtract reaction rate from target isotope
        If(key.gt.0.and.key.le.70800)Then
          isotopeTallyActive(key,3) = isotopeTallyActive(key,3) - 1.0D0*reactionRate(i)
        End If
        zB = xsKey(i,4)
        aB = xsKey(i,5)
        mB = xsKey(i,6)
        key = makeIsotopeKey(z,a,m)
! add reaction rate to the product isotope
        If(key.gt.0.and.key.le.70800)Then
          isotopeTallyActive(key,3) = isotopeTallyActive(key,3) + 1.0D0*reactionRate(i)
        End If
      End If
    End Do
    If(mpiProcessID.eq.0)Then
! save starting isotope tally to output file
      write(999,"(A1)") " "
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
      write(999,"(A39)") "Full Starting Isotope Tally"
      write(999,"(A6,A8,A4,A4,A2,A4,A4,A18,A18,A18,A21)") &
      "Key   ","Element ",&
      "Z   ","A   ","M ",&
      "mk  ","sim ",&
      "Half life         ","Activity          ","Reaction Rate     ",&
      "Atoms               "
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
    End If
    simIsotopes = 0
    Do j=1,size(simIsotopeKeys)
      key = simIsotopeKeys(j)
      If(key.gt.0.and.key.le.70800)Then
        isotopeTallyActive(key,5) = isotopeTallyActive(key,4)
        If(mpiProcessID.eq.0)Then
          write(999,"(I6.6,A7,A1,I3.3,A1,I3.3,A1,I1.1,A1,I3.3,A1,I3.3,&
A1,D17.10,A1,D17.10,A1,d17.10,A1,D17.10)") &
          key,isotopeTallyChar(key)," ",&
          isotopeTallyInt(key,1)," ",isotopeTallyInt(key,2)," ",isotopeTallyInt(key,3)," ",&
          isotopeTallyInt(key,4)," ",isotopeTallyInt(key,5)," ",&
          isotopeTallyActive(key,1)," ",isotopeTallyActive(key,2)," ",&
          isotopeTallyActive(key,3)," ",isotopeTallyActive(key,4)
        End If
        simIsotopes = simIsotopes + 1
      End If
    End Do
    If(mpiProcessID.eq.0)Then
! close output file
      write(999,"(A1)") " "
      close(999)
    End If
  End Subroutine prepareIsotopeTally

! ------------------------------------------------------------------------!
! SUBROUTINE runInOut
! Run the in-out production calculation time steps
!
!
! ------------------------------------------------------------------------!

  Subroutine runInOut()
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,k,keyTime,saveIsotope,topFive
    Real(kind=DoubleReal) :: w
    Real(kind=DoubleReal) :: simTime, simTimeLast, workingTimeStep, simTimeNext
    Real(kind=DoubleReal) :: startTime, endTime, changeTime
    Real(kind=DoubleReal) :: startTimeDisplay, endTimeDisplay
    Real(kind=DoubleReal) :: gammaActivity
    Real(kind=DoubleReal) :: totalIsotopeActivity
    Integer(kind=StandardInteger) :: key,keyT
    Integer(kind=StandardInteger) :: zT,aT,mT,zP,aP,mP
    Integer(kind=StandardInteger) :: z,a,m
    Integer(kind=StandardInteger) :: isotopeActivityHeight, isotopeActivityWidth
    Integer(kind=StandardInteger) :: beamOnOff, beamTransition, resetTimeStep
    Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: isotopeActivityKey
    Integer(kind=StandardInteger) :: printCounter
    Integer(kind=StandardInteger), Dimension(1:5) :: topFiveKeys
    Real(kind=DoubleReal), Dimension( : , : , : ), Allocatable :: isotopeActivity
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: orderedActivity
    Real(kind=DoubleReal) :: tempAA, tempAB, tempBA, tempBB
    Integer(kind=StandardInteger) :: orderedActivityCount
    Character(len=7) :: tempStrA, tempStrB
    Character(len=14) :: tempStrC, tempStrD, tempStrE
    Character(len=24) :: dbleStr
    Character(len=32) :: timeStr
    Character(len=2048) :: fileRow
    Character(len=16384) :: plotA, plotB, plotC
    Character(len=32768) :: outputLine
    Character(len=65536) :: outputLineA
    Character(len=65536) :: outputLineB
    Character(len=65536) :: outputLineC
    Logical :: sorted
    Real(kind=DoubleReal) :: gammaEnergy, maxGammaEnergy, maxGammaCounts
    Integer(kind=StandardInteger) :: gammaBin
! cmd message
    Integer(kind=StandardInteger) :: termExitStat
! Init variables
    topFiveKeys = 0
    isotopeActivityHeight = 0
! Set start sim time 0
    simTime = 0
    workingTimeStep = timeStep
! starting simulation - print if verbose on
    If(verboseTerminal.eqv..true.)Then
      Print"(A15)","    Sim started"
    End If
! print if verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A23,F8.4)","    Simulation Started ",ProgramTime()
    End If
    If(mpiProcessID.eq.0)Then
! write start of simulation to output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
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
! close file
      close(999)
    End If
! set beam on or off
    If(beamDuration.gt.0)Then
      beamOnOff = 1 !assume beam on
    Else
      beamOnOff = 0
    End If
! flags
    beamTransition = 0
    resetTimeStep = 0
! set up temporary tally array
    If(Allocated(decayTempTally))Then
      Deallocate(decayTempTally)
    End If
    Allocate(decayTempTally(1:70800))
! calculate amounts at beamDuration (needs start atoms, and atoms when beam turned off)
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
      write(999,"(A70)") "----------------------------------------------------------------------"
      write(999,"(A20,E20.10)") "Simulation tally at ",1.0D0*beamDuration
      write(999,"(A70)") "----------------------------------------------------------------------"
! close file
      close(999)
    End If
! fill sim isotope slots with zeros
    Do j=1,size(simIsotopeKeys,1)
      key = simIsotopeKeys(j)
      If(key.gt.0.and.key.le.70800)Then
        decayTempTally(key) = 0.0D0
      End If
    End Do
! target atoms lost due to transmutation
    Do j=1,size(targetReactionRatesInt,1)
      zT = targetReactionRatesInt(j,1)
      aT = targetReactionRatesInt(j,2)
      mT = targetReactionRatesInt(j,3)
      w = 1.0D0*targetReactionRates(j)
      keyT=makeIsotopeKey(zT,aT,mT)
      If(keyT.gt.0.and.keyT.le.70800)Then
        decayTempTally(keyT) = decayTempTally(keyT) + 1.0D0*w*(beamDuration-simTime)
      End If
    End Do
! Product atoms lost
    Do j=1,size(productReactionRatesInt,1)
      zP = productReactionRatesInt(j,1)
      aP = productReactionRatesInt(j,2)
      mP = productReactionRatesInt(j,3)
      w = 1.0D0*productReactionRates(j)
      Call tallyDecay&
      (zP,aP,mP,w,1.0D0*simTime,&
      1.0D0*beamDuration,1.0D0*beamDuration,4)
    End Do
! Apply tally changes
    Do j=1,size(simIsotopeKeys)
      key = simIsotopeKeys(j)
      If(key.gt.0.and.key.le.70800)Then
        isotopeTallyActive(key,4) = isotopeTallyActive(key,5) + 1.0D0 * decayTempTally(key)
! ensure no negative tally counts
        If(isotopeTallyActive(key,4).lt.0.0D0)Then
          isotopeTallyActive(key,4) = 0.0D0
        End If
! store in the beamDuration tally
        isotopeTallyActive(key,7) = isotopeTallyActive(key,4)
      End If
    End Do
! Calculate activity
    Call calcTallyActivity()
! output tally to output file
    Call outputTally()
! make empty isotopeActivity arrays
    If(individualIsotopeActivity(1:1).eq."Y")Then
      isotopeActivityHeight = Int((amTime/timeStep)+10)
      isotopeActivityWidth = size(simIsotopeKeys)
      Allocate(isotopeActivityKey(1:isotopeActivityHeight,1:isotopeActivityWidth))
      Allocate(isotopeActivity(1:isotopeActivityHeight,1:isotopeActivityWidth,1:3))
      !isotopeActivityKey = 0
      !isotopeActivity = 0.0D0
! fill with -1s
      Do i=1,isotopeActivityHeight
        Do j=1,isotopeActivityWidth
          isotopeActivityKey(i,j) = -1
          Do k=1,3
            isotopeActivity(i,j,k) = -1.0D0
          End Do
        End Do
      End Do
    End If
    If(mpiProcessID.eq.0)Then
! loop through time steps
      printCounter = 0
      Do i=1,1000000000
        If(simTime.ge.amTime)Then
          exit
        End If
! print if verbose on
        If(verboseTerminal.eqv..true.)Then
          If(printCounter.eq.0)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 0% complete, sim time   ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
          If(simTime.ge.(0.1*amTime).and.printCounter.eq.1)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 10% complete, sim time  ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
          If(simTime.ge.(0.2*amTime).and.printCounter.eq.2)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 20% complete, sim time  ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
          If(simTime.ge.(0.3*amTime).and.printCounter.eq.3)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 30% complete, sim time  ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
          If(simTime.ge.(0.4*amTime).and.printCounter.eq.4)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 40% complete, sim time  ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
          If(simTime.ge.(0.5*amTime).and.printCounter.eq.5)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 50% complete, sim time  ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
          If(simTime.ge.(0.6*amTime).and.printCounter.eq.6)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 60% complete, sim time  ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
          If(simTime.ge.(0.7*amTime).and.printCounter.eq.7)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 70% complete, sim time  ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
          If(simTime.ge.(0.8*amTime).and.printCounter.eq.8)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 80% complete, sim time  ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
          If(simTime.ge.(0.9*amTime).and.printCounter.eq.9)Then
            print "(A39,F14.6,A6,F8.4)","    Simulation 90% complete, sim time  ",&
            simTime,"      ",ProgramTime()
            printCounter = printCounter + 1
          End If
        End If
! store last sim time
        simTimeLast = simTime
! Check if sim time step exceeds activity measurement time
        If(simTime.lt.amTime.and.(simTime+workingTimeStep).gt.amTime)Then
          workingTimeStep = amTime - simTime
        End If
! check if beam ends in this time step
        If(simTime.ge.beamDuration)Then
          beamOnOff = 0
        End If
        If(beamOnOff.eq.1.and.beamDuration.gt.simTime.and.&
          beamDuration.lt.(simTime+workingTimeStep).and.beamTransition.eq.0)Then
          beamTransition = 1
          simTimeNext = simTime+workingTimeStep !store next simTime
          workingTimeStep = beamDuration - simTime !adjust time step
        ElseIf(beamTransition.eq.1)Then
          beamTransition = 0
          beamOnOff = 0  !switch off beam
          workingTimeStep = simTimeNext - simTime
          resetTimeStep = 1
! zero out reaction rates (beam off)
          Do j=1,size(simIsotopeKeys)
            key = simIsotopeKeys(j)
            If(key.gt.0.and.key.le.70800)Then
              isotopeTallyActive(key,3) = 0.0D0
            End If
          End Do
        End If
! start-end times
        If(beamOnOff.eq.1)Then
          startTime = 0.0D0
          endTime = simTime + workingTimeStep
          startTimeDisplay = startTime
          endTimeDisplay = endTime
        Else
          startTime = 0.0D0
          endTime = simTime + workingTimeStep - beamDuration
          startTimeDisplay = startTime + beamDuration
          endTimeDisplay = endTime + beamDuration
        End If
        changeTime = endTime - startTime
! fill sim isotope slots with zeros
        Do j=1,size(simIsotopeKeys,1)
          key = simIsotopeKeys(j)
          If(key.gt.0.and.key.le.70800)Then
            decayTempTally(key) = 0.0D0
          End If
        End Do
! target atoms lost due to transmutation
        If(beamOnOff.eq.1)Then
          Do j=1,size(targetReactionRatesInt,1)
            zT = targetReactionRatesInt(j,1)
            aT = targetReactionRatesInt(j,2)
            mT = targetReactionRatesInt(j,3)
            w = 1.0D0*targetReactionRates(j)
            keyT=makeIsotopeKey(zT,aT,mT)
            If(keyT.gt.0.and.keyT.le.70800)Then
              decayTempTally(keyT) = decayTempTally(keyT) + w * changeTime
            End If
          End Do
        End If
! print *,simTime,beamOnOff,startTime,endTime
! Product atoms lost
        Do j=1,size(productReactionRatesInt,1)
          zP = productReactionRatesInt(j,1)
          aP = productReactionRatesInt(j,2)
          mP = productReactionRatesInt(j,3)
          w = 1.0D0*productReactionRates(j)
          If(beamOnOff.eq.1)Then
            Call tallyDecay&
            (zP,aP,mP,w,1.0D0*startTime,&
            1.0D0*endTime,1.0D0*beamDuration,5)
          Else
            Call tallyDecay&
            (zP,aP,mP,0.0D0,1.0D0*startTime,&
            1.0D0*endTime,1.0D0*beamDuration,7)
! Call tallyDecay&
!  (zP,aP,mP,0.0D0,0.0D0,&
!  0.0001D0,1.0D0*beamDuration,7)
          End If
        End Do
! Apply tally changes
        Do j=1,size(simIsotopeKeys)
          key = simIsotopeKeys(j)
          If(key.gt.0.and.key.le.70800)Then
            If(beamOnOff.eq.1)Then
              isotopeTallyActive(key,4) = isotopeTallyActive(key,5) + 1.0D0 * decayTempTally(key)
            Else
              isotopeTallyActive(key,4) = isotopeTallyActive(key,7) + 1.0D0 * decayTempTally(key)
            End If
! ensure no negative tally counts
            If(isotopeTallyActive(key,4).lt.0.0D0)Then
              isotopeTallyActive(key,4) = 0.0D0
            End If
! store individual isotope activity tally
            If(individualIsotopeActivity(1:1).eq."Y")Then
              If(isotopeTallyActive(key,1).gt.0.and.isotopeTallyActive(key,4).gt.0)Then
                isotopeActivityKey(i,j) = key
                isotopeActivity(i,j,1) = endTimeDisplay
                isotopeActivity(i,j,2) = isotopeTallyActive(key,4)
                isotopeActivity(i,j,3) = isotopeTallyActive(key,6)*isotopeTallyActive(key,4)
              End If
            End If
          End If
        End Do
! Calculate activity
        Call calcTallyActivity()
! Total activity
        totalIsotopeActivity = 0.0D0
        Do j=1,size(simIsotopeKeys)
          key = simIsotopeKeys(j)
          If(key.gt.0.and.key.le.70800)Then
            If(isotopeTallyActive(key,1).gt.0.and.isotopeTallyActive(key,4).gt.0)Then
              totalIsotopeActivity=totalIsotopeActivity+1.0D0*&
              isotopeTallyActive(key,6)*isotopeTallyActive(key,4)
            End If
          End If
        End Do
! Save activity/time to file - open output file
        If(mpiProcessID.eq.0)Then
          open(unit=998,file=trim(fileActivityHistory),status="old",position="append",action="write")
          write(998,"(E20.10,A1,E20.10)") &
          (simTime+workingTimeStep)," ",totalIsotopeActivity
          close(998)
        End If
! increment time step
        simTime = simTimeLast + workingTimeStep
! reset time step
        If(resetTimeStep.eq.1)Then
          workingTimeStep = timeStep
          resetTimeStep = 0
        End If
      End Do
    End If
! output final tally
    If(mpiProcessID.eq.0)Then
! write to output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
      write(999,"(A70)") "----------------------------------------------------------------------"
      write(999,"(A20,E20.10)") "Simulation tally at ",1.0D0*amTime
      write(999,"(A70)") "----------------------------------------------------------------------"
      close(999)
    End If
    Call outputTally()
! If print verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A39,F14.6,A6,F8.4)","    Simulation 100% complete, sim time  ",&
      amTime,"      ",ProgramTime()
      printCounter = printCounter + 1
    End If
    If(individualIsotopeActivity(1:1).eq."Y")Then
! save isotope activities to file - open file
      If(mpiProcessID.eq.0)Then
        open(unit=996,file=trim(fileIsotopeActivity),status="old",position="append",action="write")
      End If
! loop through isotopes
      Do j=1,isotopeActivityWidth
        saveIsotope = 1
        If(isotopeActivityKey(1,j).ne.-1)Then
          Do i=1,isotopeActivityHeight
            If(isotopeActivityKey(i,j).ne.-1.and.&
              isotopeActivity(i,j,2).gt.0.0D0.and.&
              isotopeActivity(i,j,3).gt.0.0D0)Then
              saveIsotope = 1
            End If
          End Do
        End If
        If(saveIsotope.eq.1)Then
          If(mpiProcessID.eq.0)Then
            write(996,"(A65)") "-----------------------------------------------------------------"
          End If
          key = isotopeActivityKey(1,j)
          If(key.gt.0.and.key.le.70800)Then
            If(mpiProcessID.eq.0)Then
              tempStrA = "       "
              tempStrB = "       "
              write(tempStrA,"(I2)") isotopeTallyInt(key,1)
              tempStrB = trim(tempStrA)//elementSymbol(isotopeTallyInt(key,1))
              tempStrA = "       "
              write(tempStrA,"(I2)") isotopeTallyInt(key,2)
              tempStrB = trim(tempStrB)//trim(RemoveSpaces(tempStrA))
              If(isotopeTallyInt(key,3).eq.1)Then
                tempStrB = trim(tempStrB)//"m"
              End If
              write(996,"(A7)") trim(tempStrB)
              write(996,"(I8,A1,I8,A1,I8,A1,I8)") key," ",isotopeTallyInt(key,1)," ",&
              isotopeTallyInt(key,2)," ",isotopeTallyInt(key,3)
              write(996,"(A65)") "-----------------------------------------------------------------"
            End If
! loop through times
            Do i=1,isotopeActivityHeight
              If(isotopeActivityKey(i,j).ne.-1)Then
                If(mpiProcessID.eq.0)Then
                  write(996,"(E20.10,A1,E20.10,A1,E20.10)") &
                  isotopeActivity(i,j,1)," ",isotopeActivity(i,j,2)," ",isotopeActivity(i,j,3)
                End If
              End If
            End Do
            If(mpiProcessID.eq.0)Then
              write(996,"(A1)") " "
              write(996,"(A1)") " "
            End If
          End If
        End If
      End Do
      If(mpiProcessID.eq.0)Then
        close(996)
      End If
    End If
    If(mpiProcessID.eq.0)Then
      open(unit=996,file=trim(fileIsotopeActivityG),status="old",position="append",action="write")
! Activity only tally
! print *,isotopeActivityWidth
! loop through times
      fileRow = BlankString(fileRow)
      fileRow = "Time/s"
! loop through isotopes to make heading
      k = 0
      Do j=1,isotopeActivityWidth     ! Loop through isotopes
        key = isotopeActivityKey(1,j)
        If(key.gt.0.and.key.le.70800)Then
          k = k + 1
          If(k.eq.1)Then
            keyTime = j
          End If
          tempStrA = "       "
          tempStrB = "       "
          write(tempStrA,"(I2)") isotopeTallyInt(key,1)
          tempStrB = trim(tempStrA)//elementSymbol(isotopeTallyInt(key,1))
          tempStrA = "       "
          write(tempStrA,"(I2)") isotopeTallyInt(key,2)
          tempStrB = trim(tempStrB)//trim(RemoveSpaces(tempStrA))
          If(isotopeTallyInt(key,3).eq.1)Then
            tempStrB = trim(tempStrB)//"m"
          End If
          fileRow = fileRow(1:(k*16))//tempStrB//"         "
        End If
      End Do
      write(996,*) trim(fileRow)
      Do i=1,isotopeActivityHeight      ! Loop through times
        If(isotopeActivity(i,keyTime,1).eq.-1.0D0)Then
          Exit
        End If
        fileRow = BlankString(fileRow)
        write(tempStrC,"(E14.6)") isotopeActivity(i,keyTime,1)
        fileRow = tempStrC//"  "
        k = 0
        Do j=1,isotopeActivityWidth     ! Loop through isotopes
          key = isotopeActivityKey(1,j)
          If(key.gt.0.and.key.le.70800)Then
            k = k + 1
            If(isotopeActivity(i,j,3).eq.-1.0D0)Then
              tempStrC = " 0.000000E+00 "
              fileRow = fileRow(1:(k*16))//tempStrC//"  "
            Else
              write(tempStrC,"(E14.6)") isotopeActivity(i,j,3)
              fileRow = fileRow(1:(k*16))//tempStrC//"  "
            End If
          End If
        End Do
        write(996,*) trim(fileRow)
      End Do
      close(996)
      If(mpiProcessID.eq.0)Then
! Print gamma data for final tally - open output file
        open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
        write(999,"(A70)") "----------------------------------------------------------------------"
        write(999,"(A33,E20.10)") "Gamma Lines at end of Simulation ",amTime
        write(999,"(A70)") "----------------------------------------------------------------------"
        write(999,"(A9,A9,A9,A4,A20,A4,A20)") "Z        ","A        ","M        ",&
        "    ","Gamma Energy/KeV    ","    ","Gamma count/s      "
      End If
      maxGammaEnergy = 0.0D0
      Do i=1,size(simIsotopeKeys)
        key = simIsotopeKeys(i)
        If(key.gt.0.and.key.le.70800)Then
          z = isotopeTallyInt(key,1)
          a = isotopeTallyInt(key,2)
          m = isotopeTallyInt(key,3)
          If(isotopeTallyActive(key,2).gt.0)Then
            Do j=1,size(gammaLinesKey,1)
              If(gammaLinesKey(j,1).eq.z.and.gammaLinesKey(j,2).eq.a.and.&
                gammaLinesKey(j,3).eq.m)Then
                gammaActivity = 1.0D0*isotopeTallyActive(key,2)*gammaLines(j,2)
                If(mpiProcessID.eq.0)Then
                  gammaEnergy = (gammaLines(j,1)/1.0D3)
                  write(999,"(I8,A1,I8,A1,I8,A4,E20.10,A4,E20.10)")&
                  z," ",a," ",m,"   ",gammaEnergy,"    ",gammaActivity
                  If(gammaEnergy.gt.maxGammaEnergy)Then
                    maxGammaEnergy = gammaEnergy
                  End If
                End If
              End If
            End Do
          End If
        End If
      End Do
! print *,"Max gamma",maxGammaEnergy
! Sum gammas into bin
      gammaLineBins = 0.0D0
      If(gammaResolution.eq.0.0D0)Then
        gammaResolution = ceiling(maxGammaEnergy)
      End If
      Do i=1,size(simIsotopeKeys)
        key = simIsotopeKeys(i)
        If(key.gt.0.and.key.le.70800)Then
          z = isotopeTallyInt(key,1)
          a = isotopeTallyInt(key,2)
          m = isotopeTallyInt(key,3)
          If(isotopeTallyActive(key,2).gt.0)Then
            Do j=1,size(gammaLinesKey,1)
              If(gammaLinesKey(j,1).eq.z.and.gammaLinesKey(j,2).eq.a.and.&
                gammaLinesKey(j,3).eq.m)Then
                gammaActivity = 1.0D0*isotopeTallyActive(key,2)*gammaLines(j,2)
                gammaEnergy = (gammaLines(j,1)/1.0D3)
                gammaBin = ceiling((gammaEnergy/maxGammaEnergy)*gammaResolution)
                gammaLineBins(gammaBin) = gammaLineBins(gammaBin) + gammaActivity
              End If
            End Do
          End If
        End If
      End Do
! Largest bin/max gamma count value
      maxGammaCounts = 0.0D0
      Do i=1,gammaResolution
        If(maxGammaCounts.lt.gammaLineBins(i))Then
          maxGammaCounts = gammaLineBins(i)
        End If
      End Do
      If(mpiProcessID.eq.0)Then
        open(unit=995,file=Trim(outputDirectory)//"/"//"gammaLines.dat")
        write(995,"(A1)") " "
        Do i=1,gammaResolution
          write(995,"(E14.6,I8,E14.6)") &
          ((maxGammaEnergy/gammaResolution)*i),i,gammaLineBins(i)
        End Do
        close(995)
      End If
! order final tally by activity - count active isotopes
      orderedActivityCount = 0
      Do j=1,size(simIsotopeKeys)
        key = simIsotopeKeys(j)
        If(key.gt.0.and.key.le.70800)Then
          If(isotopeTallyActive(key,2).gt.0.0D0)Then
            orderedActivityCount = orderedActivityCount + 1
          End If
        End If
      End Do
      Allocate(orderedActivity(1:orderedActivityCount,1:2))
      orderedActivityCount = 0
      Do j=1,size(simIsotopeKeys)
        key = simIsotopeKeys(j)
        If(key.gt.0.and.key.le.70800)Then
          If(isotopeTallyActive(key,2).gt.0.0D0)Then
            orderedActivityCount = orderedActivityCount + 1
            orderedActivity(orderedActivityCount,1) = key
            orderedActivity(orderedActivityCount,2) = isotopeTallyActive(key,2)
          End If
        End If
      End Do
! Sort orderedActivity array
      sorted = .false.
      Do While(sorted.eqv..false.)
        sorted = .true.
        Do j=1,(size(orderedActivity,1)-1)
          If(orderedActivity(j,2).lt.orderedActivity(j+1,2))Then
            sorted = .false.
! temp store array values
            tempAA = orderedActivity(j,1)
            tempAB = orderedActivity(j,2)
            tempBA = orderedActivity(j+1,1)
            tempBB = orderedActivity(j+1,2)
! swap values
            orderedActivity(j,1) = tempBA
            orderedActivity(j,2) = tempBB
            orderedActivity(j+1,1) = tempAA
            orderedActivity(j+1,2) = tempAB
          End If
        End Do
      End Do
! output
      If(mpiProcessID.eq.0)Then
        write(999,"(A1)") " "
        write(999,"(A70)") "----------------------------------------------------------------------"
        write(999,"(A33)") "Activities Ordered by Most Active"
        write(999,"(A70)") "----------------------------------------------------------------------"
        write(999,"(A1)") " "
      End If
      k = 0
      Do j=1,size(orderedActivity,1)
        key = Int(orderedActivity(j,1))
        If(key.gt.0.and.key.le.70800)Then
          If(mpiProcessID.eq.0)Then
            write(999,"(I8,A1,A2,A1,I8,A1,I8,A1,I8,A4,E20.10)") &
            key," ",elementSymbol(isotopeTallyInt(key,1))," ",&
            isotopeTallyInt(key,1)," ",isotopeTallyInt(key,2)," ",&
            isotopeTallyInt(key,3),"    ",isotopeTallyActive(key,2)
            k = k + 1
            If(k.le.5)Then
              topFiveKeys(k) = key
              print *,k,key
            End If
          End If
        End If
      End Do
! If verbose on
      If(verboseTerminal.eqv..true.)Then
        print "(A4,A21,A5,F8.4)","    ","Most Active Isotopes ","     ",ProgramTime()
        print "(A6,A12,A9,A9,A18)","    ","Symbol    Z ","       A ",&
        "       M ","       Activity/Bq"
        Do j=1,size(orderedActivity,1)
          key = Int(orderedActivity(j,1))
          If(key.gt.0.and.key.le.70800)Then
            Print "(A6,A2,A1,I8,A1,I8,A1,I8,A4,E20.10)", "      ",&
            elementSymbol(isotopeTallyInt(key,1))," ",&
            isotopeTallyInt(key,1)," ",isotopeTallyInt(key,2)," ",&
            isotopeTallyInt(key,3),"    ",isotopeTallyActive(key,2)
          End If
          If(j.eq.5)Then
            Exit
          End If
        End Do
      End If
! -------------------------------------
! PLOTS (if pyplot installed)
! -------------------------------------
! Make pyplot for activities
      open(unit=995,file=trim(tempDirectory)//"/activity.py")
      open(unit=994,file=trim(tempDirectory)//"/atoms.py")
      open(unit=993,file=trim(tempDirectory)//"/gammaLines.py")
! Time sentence
      timeStr = "at                              "
      write(dbleStr,"(F14.6)") amTime
      dbleStr = RemoveSpaces(dbleStr)
      timeStr = Trim(timeStr)//" "//Trim(dbleStr)
! Write python files - activity
      write(995,"(A)") "#!/usr/bin/env python"
      write(995,"(A)") "import numpy as np"
      write(995,"(A)") "import matplotlib"
      write(995,"(A)") "matplotlib.use('Agg')"
      write(995,"(A)") "import matplotlib.pyplot as plt"
! Write python files - atoms
      write(994,"(A)") "#!/usr/bin/env python"
      write(994,"(A)") "import numpy as np"
      write(994,"(A)") "import matplotlib"
      write(994,"(A)") "matplotlib.use('Agg')"
      write(994,"(A)") "import matplotlib.pyplot as plt"
! Write python files - atoms
      write(993,"(A)") "#!/usr/bin/env python"
      write(993,"(A)") "import numpy as np"
      write(993,"(A)") "import matplotlib"
      write(993,"(A)") "matplotlib.use('Agg')"
      write(993,"(A)") "import matplotlib.pyplot as plt"
! Set figure sizes
      write(993,"(A)") "plt.figure(figsize=(1792/144, 1008/144), dpi=144)"
      write(994,"(A)") "plt.figure(figsize=(1792/144, 1008/144), dpi=144)"
      write(995,"(A)") "plt.figure(figsize=(1792/144, 1008/144), dpi=144)"
! ------------
! Activity and Atoms
! ------------
! Prepare and save data
      Do j=1,isotopeActivityWidth     ! Loop through isotopes
        key = isotopeActivityKey(1,j)
        If(key.gt.-1)Then
          plotA = BlankString(plotA)
          plotB = BlankString(plotB)
          plotC = BlankString(plotC)
          plotA = "["
          plotB = "["
          plotC = "["
! Isotope label
          tempStrA = "       "
          tempStrB = "       "
          write(tempStrA,"(I2)") isotopeTallyInt(key,1)
          tempStrB = trim(tempStrA)//elementSymbol(isotopeTallyInt(key,1))
          tempStrA = "       "
          write(tempStrA,"(I2)") isotopeTallyInt(key,2)
          tempStrB = trim(tempStrB)//trim(RemoveSpaces(tempStrA))
          If(isotopeTallyInt(key,3).eq.1)Then
            tempStrB = trim(tempStrB)//"m"
          End If
          tempStrB = RemoveSpaces(tempStrB)
! Store data
          k = 0
          Do i=1,isotopeActivityHeight      ! Loop through times
            If(isotopeActivity(i,j,1).eq.-1.0D0)Then
              Exit
            End If
            write(tempStrC,"(E14.6)") isotopeActivity(i,j,1)
            write(tempStrD,"(E14.6)") isotopeActivity(i,j,2)
            write(tempStrE,"(E14.6)") isotopeActivity(i,j,3)
            If(k.eq.0)Then
              plotA = trim(plotA)//RemoveSpaces(tempStrC)
              plotB = trim(plotB)//RemoveSpaces(tempStrD)
              plotC = trim(plotC)//RemoveSpaces(tempStrE)
            Else
              plotA = trim(plotA)//","//RemoveSpaces(tempStrC)
              plotB = trim(plotB)//","//RemoveSpaces(tempStrD)
              plotC = trim(plotC)//","//RemoveSpaces(tempStrE)
            End If
            k = k + 1
          End Do
          If(k.gt.1)Then
            plotA = trim(plotA)//"]"
            plotB = trim(plotB)//"]"
            plotC = trim(plotC)//"]"
            outputLine = "plt.plot("//trim(plotA)//","//trim(plotB)//&
            ",label='"//trim(tempStrB)//"')"
            write(994,"(A)") outputLine  ! Atoms
            outputLine = "plt.plot("//trim(plotA)//","//trim(plotC)//&
            ",label='"//trim(tempStrB)//"')"
            write(995,"(A)") outputLine  ! Activity
          End If
        End If
      End Do
! Write plot titles
      write(995,"(A)") "plt.title('Target Activity "//Trim(timeStr)//"s')"
      write(994,"(A)") "plt.title('Target Active Nuclei Count "//Trim(timeStr)//"s')"
! Axes Labels
      write(994,"(A)") "plt.ylabel('Nuclei count at time t')"
      write(994,"(A)") "plt.xlabel('Time t after irradiation commenced (s)')"
      write(995,"(A)") "plt.ylabel('Activity at time t (Bq)')"
      write(995,"(A)") "plt.xlabel('Time t after irradiation commenced (s)')"
! End files
      write(994,"(A)") "plt.savefig('"//trim(outputDirectory)//"/atoms',dpi=144)"
      write(995,"(A)") "plt.savefig('"//trim(outputDirectory)//"/activity',dpi=144)"
! Close files
      close(994)
      close(995)
! Run python scripts to make graphs
      Call execute_command_line("python "//trim(tempDirectory)//"/activity.py",&
      exitstat=termExitStat)
      Call execute_command_line("python "//trim(tempDirectory)//"/atoms.py",&
      exitstat=termExitStat)
! Make pyplot for activities
      open(unit=995,file=trim(tempDirectory)//"/activityTop5.py")
      open(unit=994,file=trim(tempDirectory)//"/atomsTop5.py")
! Write python files - activity
      write(995,"(A)") "#!/usr/bin/env python"
      write(995,"(A)") "import numpy as np"
      write(995,"(A)") "import matplotlib"
      write(995,"(A)") "matplotlib.use('Agg')"
      write(995,"(A)") "import matplotlib.pyplot as plt"
! Write python files - atoms
      write(994,"(A)") "#!/usr/bin/env python"
      write(994,"(A)") "import numpy as np"
      write(994,"(A)") "import matplotlib"
      write(994,"(A)") "matplotlib.use('Agg')"
      write(994,"(A)") "import matplotlib.pyplot as plt"
! Set figure sizes
      write(994,"(A)") "plt.figure(figsize=(1792/144, 1008/144), dpi=144)"
      write(995,"(A)") "plt.figure(figsize=(1792/144, 1008/144), dpi=144)"
! Prepare and save data
      Do j=1,isotopeActivityWidth     ! Loop through isotopes
        key = isotopeActivityKey(1,j)
        If(key.gt.-1)Then
          topFive = 0
          Do k=1,5
            If(key.eq.topFiveKeys(k))Then
              topFive = 1
            End If
          End Do
          If(topFive.eq.1)Then
            plotA = BlankString(plotA)
            plotB = BlankString(plotB)
            plotC = BlankString(plotC)
            plotA = "["
            plotB = "["
            plotC = "["
! Isotope label
            tempStrA = "       "
            tempStrB = "       "
            write(tempStrA,"(I2)") isotopeTallyInt(key,1)
            tempStrB = trim(tempStrA)//elementSymbol(isotopeTallyInt(key,1))
            tempStrA = "       "
            write(tempStrA,"(I2)") isotopeTallyInt(key,2)
            tempStrB = trim(tempStrB)//trim(RemoveSpaces(tempStrA))
            If(isotopeTallyInt(key,3).eq.1)Then
              tempStrB = trim(tempStrB)//"m"
            End If
            tempStrB = RemoveSpaces(tempStrB)
! Store data
            k = 0
            Do i=1,isotopeActivityHeight      ! Loop through times
              If(isotopeActivity(i,j,1).eq.-1.0D0)Then
                Exit
              End If
              write(tempStrC,"(E14.6)") isotopeActivity(i,j,1)
              write(tempStrD,"(E14.6)") isotopeActivity(i,j,2)
              write(tempStrE,"(E14.6)") isotopeActivity(i,j,3)
              If(k.eq.0)Then
                plotA = trim(plotA)//RemoveSpaces(tempStrC)
                plotB = trim(plotB)//RemoveSpaces(tempStrD)
                plotC = trim(plotC)//RemoveSpaces(tempStrE)
              Else
                plotA = trim(plotA)//","//RemoveSpaces(tempStrC)
                plotB = trim(plotB)//","//RemoveSpaces(tempStrD)
                plotC = trim(plotC)//","//RemoveSpaces(tempStrE)
              End If
              k = k + 1
            End Do
            If(k.gt.1)Then
              plotA = trim(plotA)//"]"
              plotB = trim(plotB)//"]"
              plotC = trim(plotC)//"]"
              outputLine = "plt.plot("//trim(plotA)//","//trim(plotB)//&
              ",label='"//trim(tempStrB)//"')"
              write(994,"(A)") outputLine  ! Atoms
              outputLine = "plt.plot("//trim(plotA)//","//trim(plotC)//&
              ",label='"//trim(tempStrB)//"')"
              write(995,"(A)") outputLine  ! Activity
            End If
          End If
        End If
      End Do
! Write Titles
      write(994,"(A)") "plt.title('Target Active Nuclei Count  "//Trim(timeStr)//"s - Top 5')"
      write(995,"(A)") "plt.title('Target Activity  "//Trim(timeStr)//"s - Top 5')"
! Axes Labels
      write(994,"(A)") "plt.ylabel('Gamma Activity (Bq)')"
      write(994,"(A)") "plt.xlabel('Gamma Energy (KeV)')"
! Axes Limits/Ranges
      write(994,"(A)") "plt.ylabel('Nuclei count at time t')"
      write(994,"(A)") "plt.xlabel('Time t after irradiation commenced (s)')"
      write(995,"(A)") "plt.ylabel('Activity at time t (Bq)')"
      write(995,"(A)") "plt.xlabel('Time t after irradiation commenced (s)')"
! Legend
      write(994,"(A)") "plt.legend()"
      write(995,"(A)") "plt.legend()"
! ------------
! Gamma Lines
! ------------
! Set axis
      write(dbleStr,"(E14.6)") (1.05D0*maxGammaEnergy)
      write(993,"(A)") "plt.xlim(0.0,"//dbleStr//")"
      write(dbleStr,"(E14.6)") (1.05D0*maxGammaCounts)
      write(993,"(A)") "plt.ylim(0.0,"//dbleStr//")"
! Data
      outputLineA = BlankString(outputLineA)
      outputLineB = BlankString(outputLineB)
      outputLineC = BlankString(outputLineC)
      outputLineA = "x = ["
      outputLineB = "y = ["
      outputLineC = "zero = ["
      Do i=1,Int(gammaResolution)
        If(gammaLineBins(i).gt.0.0D0)Then
          dbleStr = BlankString(dbleStr)
! x val
          gammaEnergy = (i/(1.0D0*gammaResolution))*maxGammaEnergy
          write(dbleStr,"(E14.6)") gammaEnergy
          outputLineA = trim(outputLineA)//trim(dbleStr)
! y val
          write(dbleStr,"(E14.6)") gammaLineBins(i)
          outputLineB = trim(outputLineB)//trim(dbleStr)
! zero
          outputLineC = trim(outputLineC)//"0.0"
! trailing comma
          If(i.lt.gammaResolution)Then
            outputLineA = trim(outputLineA)//","
            outputLineB = trim(outputLineB)//","
            outputLineC = trim(outputLineC)//","
          End If
        End If
      End Do
      outputLineA = trim(outputLineA)//"]"
      outputLineB = trim(outputLineB)//"]"
      outputLineC = trim(outputLineC)//"]"
! write data
      write(993,"(A)") Trim(outputLineA)
      write(993,"(A)") Trim(outputLineB)
      write(993,"(A)") Trim(outputLineC)
! plot chart
      outputLine = "plt.errorbar(x, y, yerr=[y,zero], fmt='o', markersize=0, "//&
      "capsize=0, elinewidth=1)"
      write(993,"(A)") Trim(outputLine)
      write(993,"(A)") "plt.ylabel('Activity/Bq')"
      write(993,"(A)") "plt.xlabel('Gamma Energy/KeV')"
      write(993,"(A)") "plt.title('Predicted Gamma Lines "//Trim(timeStr)//"s')"
! outputLine = "x = arr([1.0,2.0,3.0,4.0,5.0])"
! write(993,"(A)") Trim(outputLine)
! outputLine = "y = arr([1.0,1.0,1.0,2.0,2.0])"
! write(993,"(A)") Trim(outputLine)
! outputLine = "c = arr([[0.0,1.0],[0.0,1.0],[0.0,1.0],[0.0,2.0],[0.0,2.0]]).T"
! write(993,"(A)") Trim(outputLine)
! write(993,"(A)") "plt.errorbar(x, y, yerr=c)"
! outputLine = "xerr = x"
! write(993,"(A)") Trim(outputLine)
! outputLine = "yerr = y"
! write(993,"(A)") Trim(outputLine)
! write(993,"(A)") "plt.errorbar(x, y, yerr=[yerr, 2*yerr])"
! write(993,"(A,A)") "plt.errorbar(x, y, yerr=[0, 2*yerr], ",&
! "xerr=xerr, fmt='o', ecolor='g', capthick=2)"
! write(993,"(A)") "error = 0.1 + 0.2 * x"
! write(993,"(A)") "lower_error = 0.4 * error"
! write(993,"(A)") "upper_error = error"
! write(993,"(A)") "asymmetric_error = [lower_error, upper_error]"
! write(993,"(A)") "fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)"
! write(993,"(A)") "ax0.errorbar(x, y, yerr=[0.1,3.0], fmt='-o')"
! write(993,"(A)") "ax0.set_title('variable, symmetric error')"
! ------------
! Save and close files
! ------------
! End files
      write(993,"(A)") "plt.savefig('"//trim(outputDirectory)//"/gammaLines',dpi=144)"
      write(994,"(A)") "plt.savefig('"//trim(outputDirectory)//"/atomsTop5',dpi=144)"
      write(995,"(A)") "plt.savefig('"//trim(outputDirectory)//"/activityTop5',dpi=144)"
! Close files
      close(993)
      close(994)
      close(995)
! Run python scripts to make graphs
      Call execute_command_line("python "//trim(tempDirectory)//"/gammaLines.py",&
      exitstat=termExitStat)
      Call execute_command_line("python "//trim(tempDirectory)//"/atomsTop5.py",&
      exitstat=termExitStat)
      Call execute_command_line("python "//trim(tempDirectory)//"/activityTop5.py",&
      exitstat=termExitStat)
    End If
! close file
    close(999)
  End Subroutine runInOut

! ------------------------------------------------------------------------!
! SUBROUTINE tallyDecay
!
!
!
! ------------------------------------------------------------------------!

  Subroutine tallyDecay(zP,aP,mP,parentProductionRate,tStart,tEnd,tBeamEnd,tallyColumn)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal) :: parentProductionRate,tStart,tEnd,tBeamEnd
    Integer(kind=StandardInteger) :: keyP, key, tallyColumn
    Integer(kind=StandardInteger) :: zP,aP,mP
    Real(kind=DoubleReal),Dimension(:,:,:),Allocatable::decayChainArray
    Real(kind=DoubleReal),Dimension(:,:),Allocatable::decayDataArray
    Real(kind=DoubleReal),Dimension(:,:),Allocatable::isotopeChange
    Integer(kind=StandardInteger) :: arrayHeight, arrayWidth, chainLength
! Parent/product key
    keyP = makeIsotopeKey (zP, aP, mP)
! make isotope decay chain array
    decayChainArray = IsotopesDecayChains(zP,aP,mP)
! for each chain, run decay and tally
    arrayHeight = size(decayChainArray,1)  !i, max chain length
    arrayWidth = size(decayChainArray,2)    !j, number of decay chains
! Init
    tStart = 0.0D0
! Open file
    If(mpiProcessID.eq.0)Then
      open(unit=995,file=Trim(outputDirectory)//"/"//"decayChains.log",&
      status="old",position="append",action="write")
    End If
! loop each chain
    Do j=1,arrayWidth
      If(mpiProcessID.eq.0)Then
        write(995,"(A7,I4,A4,I4)") "Branch ",j," of ",arrayWidth
      End If
! Allocate array
      If(Allocated(decayDataArray))Then
        Deallocate(decayDataArray)
      End If
! loop through each isotope
      chainLength = 0
      Do i=1,arrayHeight
        If(decayChainArray(i,j,1).ne.-1)Then
          chainLength = chainLength + 1
        Else
          Exit
        End If
      End Do
      If(chainLength.gt.0)Then
! Allocate data array
        Allocate(decayDataArray(1:chainLength,1:6))
! Store data to pass to activity calc function - build decay data array
        Do i=1,chainLength
          key = Int(decayChainArray(i,j,1))
          If(key.gt.0.and.key.le.70800)Then
            decayDataArray(i,1) = key             !Tally key
            decayDataArray(i,2) = isotopeTallyActive(key,tallyColumn)   !No. Atoms
            decayDataArray(i,3) = decayChainArray(i,j,2)      !Half life
            decayDataArray(i,4) = decayChainArray(i,j,3)      !branching factor
            decayDataArray(i,5) = isotopeTallyInt(key,1)    !isotope Z
            decayDataArray(i,6) = isotopeTallyInt(key,2)    !isotope A
          End If
          If(mpiProcessID.eq.0)Then
            If(decayChainArray(i,j,2).lt.0.0D0)Then  !Stable
              write(995,"(I8,A4,I8,I8,E18.8,E18.8,E18.8,E18.8,A9)") &
              i," "//elementSymbol(Int(decayDataArray(i,5)))//" ",&
              Int(decayDataArray(i,5)),Int(decayDataArray(i,6)),&
              decayDataArray(i,1),decayDataArray(i,2),decayDataArray(i,3),&
              decayDataArray(i,4)," *Stable*"
            Else
              write(995,"(I8,A4,I8,I8,E18.8,E18.8,E18.8,E18.8)") &
              i," "//elementSymbol(Int(decayDataArray(i,5)))//" ",&
              Int(decayDataArray(i,5)),Int(decayDataArray(i,6)),&
              decayDataArray(i,1),decayDataArray(i,2),decayDataArray(i,3),&
              decayDataArray(i,4)
            End If
          End If
        End Do
! calculate change in isotope amounts
        isotopeChange=CalcIsotopeAmount(parentProductionRate,decayDataArray,tEnd)
! record decay chain
! apply changes to temp tally
        Do i=1,size(isotopeChange,1)
          key = Int(isotopeChange(i,1))
          If(key.gt.0.and.key.le.70800)Then
            decayTempTally(key) = decayTempTally(key) + isotopeChange(i,2)
          End If
        End Do
      End If
    End Do
    If(mpiProcessID.eq.0)Then
      close(995)
    End If
  End Subroutine tallyDecay

! ------------------------------------------------------------------------!
! SUBROUTINE calcTallyActivity
!
!
!
! ------------------------------------------------------------------------!

  Subroutine calcTallyActivity
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,key
! calculate activities
    totalActivityAtTime = 0.0D0
    Do i=1,size(simIsotopeKeys)
! isotopeTallyActive(i,6) decay constant
! isotopeTallyActive(i,4) number of atoms
! isotopeTallyActive(i,2) activity
      key = simIsotopeKeys(i)
      isotopeTallyActive(key,2) = isotopeTallyActive(key,6) *&
      isotopeTallyActive(key,4)
      totalActivityAtTime = totalActivityAtTime + 1.0D0*isotopeTallyActive(key,2)
    End Do
  End Subroutine calcTallyActivity

! ------------------------------------------------------------------------!
! SUBROUTINE outputTally
!
!
!
! ------------------------------------------------------------------------!

  Subroutine outputTally
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,key
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! save starting isotope tally to output file
      write(999,"(A1)") " "
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
      write(999,"(A15,F8.4)") "Isotope Tally  ",ProgramTime()
      write(999,"(A6,A8,A4,A4,A2,A4,A4,A18,A18,A18,A18,A18,A18,A18)") &
      "Key   ","Element ",&
      "Z   ","A   ","M ",&
      "mk  ","sim ",&
      "Half life         ","Decay Constant    ","Activity          ","Reaction Rate     ",&
      "Start Atoms       ","Beam End Atoms    ","Atoms At Time     "
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
    End If
    Do j=1,size(simIsotopeKeys)
      key = simIsotopeKeys(j)
      If(isotopeTallyActive(key,4).gt.0.0D0)Then
        i = key
        If(mpiProcessID.eq.0)Then
          write(999,"(I6.6,A7,A1,I3.3,A1,I3.3,A1,I1.1,A1,I3.3,A1,I3.3,A1,D17.10,A1,&
          D17.10,A1,D17.10,A1,D17.10,A1,D17.10,A1,D17.10,A1,D17.10)") &
          key,isotopeTallyChar(i)," ",&
          isotopeTallyInt(i,1)," ",isotopeTallyInt(i,2)," ",isotopeTallyInt(i,3)," ",&
          isotopeTallyInt(i,4)," ",isotopeTallyInt(i,5)," ",&
          isotopeTallyActive(i,1)," ",isotopeTallyActive(i,6)," ",isotopeTallyActive(i,2)," ",&
          isotopeTallyActive(i,3)," ",isotopeTallyActive(i,5)," ",isotopeTallyActive(i,7)," ",&
          isotopeTallyActive(i,4)
        End If
      End If
    End Do
    If(mpiProcessID.eq.0)Then
! close output file
      write(999,"(A1)") " "
      close(999)
    End If
  End Subroutine outputTally

! ------------------------------------------------------------------------!
! SUBROUTINE finishProductionLoss
!
!
!
! ------------------------------------------------------------------------!

  Subroutine finishProductionLoss
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,key,z,a,m
    Real(kind=DoubleReal) :: totalActivity
    Real(kind=DoubleReal) :: totalGammaPower
    Real(kind=DoubleReal) :: absorbedDose
    Integer(kind=StandardInteger), Dimension(1:3) :: theTime, theDate
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! save starting isotope tally to output file
      write(999,"(A1)") " "
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
      write(999,"(A15,F8.4)") "Totals         ",ProgramTime()
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
      write(999,"(A1)") " "
    End If
! Total activity
    totalActivity = 0.0D0
    Do i=1,size(simIsotopeKeys)
! isotopeTallyActive(i,6) decay constant
! isotopeTallyActive(i,4) number of atoms
      key = simIsotopeKeys(i)
      If(isotopeTallyActive(key,2).gt.0)Then
        totalActivity = totalActivity + 1.0D0*isotopeTallyActive(key,6)*&
        isotopeTallyActive(key,4)
      End If
    End Do
! Total energy output
    totalGammaPower = 0.0D0
    Do i=1,size(simIsotopeKeys)
      key = simIsotopeKeys(i)
      z = isotopeTallyInt(key,1)
      a = isotopeTallyInt(key,2)
      m = isotopeTallyInt(key,3)
      If(isotopeTallyActive(key,2).gt.0)Then
        Do j=1,size(gammaLinesKey,1)
          If(gammaLinesKey(j,1).eq.z.and.gammaLinesKey(j,2).eq.a.and.&
            gammaLinesKey(j,3).eq.m)Then
            totalGammaPower = totalGammaPower+1.0D0*isotopeTallyActive(key,2)*&
            gammaLines(j,2)*gammaLines(j,1)
          End If
        End Do
      End If
    End Do
! Calculate absorbed dose at 1m, 80kg human, surface area 1m2
    absorbedDose = (totalGammaPower*elementaryCharge)/(4*pi*80)
! Print if verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A4,A19,E20.10,A5,F8.4)","    ","Total Activity/Bq: ",&
      totalActivity,"     ",ProgramTime()
    End If
    If(mpiProcessID.eq.0)Then
! Write to file
      write(999,"(A50,E20.10)") "Total Activity/Bq:                                ",&
      totalActivity
      write(999,"(A50,E20.10)") "Total Gamma Power/eV/s:                           ",&
      totalGammaPower
      write(999,"(A50,E20.10)") "Total Gamma Power/Watts:                          ",&
      (totalGammaPower*elementaryCharge)
      write(999,"(A50,E20.10)") "Absorbed Dose*/Grays/s:                           ",&
      absorbedDose
      write(999,"(A50,E20.10)") "Absorbed Dose*/Grays/hr:                          ",&
      (absorbedDose*3600)
      write(999,"(A50,E20.10)") "Fraction of annual dosage if exposed for 1 hr:    ",&
      ((absorbedDose*3600)/0.001)
      write(999,"(A1)") " "
      write(999,"(A56,A61)") "Absorbed dose, assumes all energy absorbed, 80kg human, ",&
      "1m from point-target, 1m surface area exposed to irradiation."
      write(999,"(A1)") " "
      write(999,"(A12)") "Dose Limits:"
      write(999,"(A35)") "employees 18+ 20 millisieverts/year"
      write(999,"(A33)") "trainees 18+ 6 millisieverts/year"
      write(999,"(A41)") "public and under 18s 1 millisieverts/year"
      write(999,"(A54,E20.10)") "public and under 18s millisieverts averaged per hour: ",&
      (1/(365.25*24))
      write(999,"(A50)") "Dose averaged over area of skin not exceeding 1cm2"
      write(999,"(A55)") "Source: http://www.hse.gov.uk/radiation/ionising/doses/"
! close output file
      write(999,"(A1)") " "
      write(999,"(A1)") " "
      write(999,"(A1)") " "
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
      write(999,"(A30,F8.4)") "Activity calculation complete.",ProgramTime()
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
      write(999,"(A1)") " "
      write(999,"(A1)") " "
      close(999)
    End If
    If(verboseTerminal.eqv..true.)Then
      call idate(theDate)   ! theDate(1)=day, (2)=month, (3)=year
      call itime(theTime)   ! theDate(1)=hour, (2)=minute, (3)=second
      print "(A14,F8.4)","    Run time: ",ProgramTime()
      print "(A14,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)",&
      "    Ended at: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)
    End If
  End Subroutine finishProductionLoss

! ------------------------------------------------------------------------!
! Decay Functions
!
! ------------------------------------------------------------------------!

  Function IsotopesInDecayTree(z,a,m) RESULT (isotopeArray)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,z,a,m,complete
    Real(kind=DoubleReal), Dimension(:,:), &
    Allocatable :: isotopeArray
! Allocate an array to store isotope list
    If(Allocated(decayChainIsotopesArray))Then
      Deallocate(decayChainIsotopesArray)
    End If
    Allocate(decayChainIsotopesArray(1:300,1:3))
! fill with blank data
    Do i=1,300
      Do j=1,3
        decayChainIsotopesArray(i,j) = 0
      End Do
    End Do
! set counter = 1 and store parent isotope
    decayChainIsotopesArray(1,1) = z
    decayChainIsotopesArray(1,2) = a
    decayChainIsotopesArray(1,3) = m
    isotopeCounter = 1
! call function recursively
    complete = IsotopesInDecayTreeR(z,a,m)
! transfer data to new array
    Allocate(isotopeArray(1:isotopeCounter,1:3))
    Do i=1,isotopeCounter
      Do j=1,3
        isotopeArray(i,j) = decayChainIsotopesArray(i,j)
      End Do
    End Do
! Deallocate array
    If(Allocated(decayChainIsotopesArray))Then
      Deallocate(decayChainIsotopesArray)
    End If
  End Function IsotopesInDecayTree
! ---------------------------------------------------------------------
  Recursive Function IsotopesInDecayTreeR(z,a,m) RESULT (lastRow)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,z,a,m,lastRow
    Integer(kind=StandardInteger) :: zP, aP, mP
    Integer(kind=StandardInteger) :: zC, aC, mC
    Logical :: store
! loop through decay chains
    Do i=1,size(decayChar)
      zP = decayInt(i,2)
      aP = decayInt(i,1)
      mP = decayInt(i,3)
      If(z.eq.zP.and.a.eq.aP.and.m.eq.mP)Then
        zC = decayInt(i,5)
        aC = decayInt(i,4)
        mC = decayInt(i,6)
! check not in array already
        store = .true.
        Do j=1,size(decayChainIsotopesArray,1)
          If(decayChainIsotopesArray(j,1).eq.zC.and.&
            decayChainIsotopesArray(j,2).eq.aC.and.&
            decayChainIsotopesArray(j,3).eq.mC)Then
            store = .false.
          End If
        End Do
        If(store.eqv..true.)Then
          isotopeCounter = isotopeCounter + 1
          decayChainIsotopesArray(isotopeCounter,1) = zC
          decayChainIsotopesArray(isotopeCounter,2) = aC
          decayChainIsotopesArray(isotopeCounter,3) = mC
! print *,isotopeCounter,zP,aP,mP,zC,aC,mC
        End If
        lastRow = IsotopesInDecayTreeR(zC,aC,mC)
      End If
    End Do
  End Function IsotopesInDecayTreeR

  Function IsotopesDecayChains(z,a,m) RESULT (decayChains)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,k
    Integer(kind=StandardInteger) :: z, a, m
    Integer(kind=StandardInteger) :: nodes
    Integer(kind=StandardInteger) :: tempArraySize
    Real(kind=DoubleReal),Dimension(:,:,:),Allocatable::decayChains
    Integer(kind=StandardInteger) :: arrayHeight, arrayWidth
    Integer(kind=StandardInteger) :: keyP, keyC, key
    Logical :: arrayEnd
! make temporary square array 50x50
    tempArraySize = 25
    If(Allocated(public3DTempArray))Then
      Deallocate(public3DTempArray)
    End If
    Allocate(public3DTempArray(1:tempArraySize,1:tempArraySize,1:3))
! fill with blank markers (-1)
    Do i=1,tempArraySize
      Do j=1,tempArraySize
        public3DTempArray(i,j,1) = -1
      End Do
    End Do
! run recursive function to fill public3DTempArray
    isotopeCounter = 0
    publicX = 0
    publicY = 1
    nodes = IsotopesDecayChainsR(z,a,m)
! measure public3DTempArray
    arrayHeight = 0
    Do i=1,tempArraySize
      arrayEnd = .true.
      Do j=1,tempArraySize
        If(public3DTempArray(i,j,1).ne.-1)Then
          arrayEnd = .false.
        End If
      End Do
      If(arrayEnd.eqv..true.)Then
        arrayHeight = i - 1
        Exit
      End If
    End Do
    arrayWidth = 0
    Do j=1,tempArraySize
      arrayEnd = .true.
      Do i=1,tempArraySize
        If(public3DTempArray(i,j,1).ne.-1)Then
          arrayEnd = .false.
        End If
      End Do
      If(arrayEnd.eqv..true.)Then
        arrayWidth = j - 1
        Exit
      End If
    End Do
! Allocate decay chains array
    Allocate(decayChains(1:arrayHeight,1:arrayWidth,1:3))
! Transfer into new array of correct size
    Do i=1,arrayHeight
      arrayEnd = .true.
      Do j=1,arrayWidth
        decayChains(i,j,1) = public3DTempArray(i,j,1)
      End Do
    End Do
! Fill in blanks
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
! Fill in half lives
    Do j=1,arrayWidth
      Do i=1,arrayHeight
! Do k=1,size(decayInt,1)
!  If(decayChains(i,j,1).eq.decayInt(k,7))Then
!    decayChains(i,j,2) = decayDouble(k,2)
!  Exit
!  End If
! End Do
        key = Int(decayChains(i,j,1))
        If(key.gt.0.and.key.le.70800)Then
          decayChains(i,j,2) = isotopeTallyActive(key,1)
        End If
      End Do
    End Do
! Fill in branching factors
    Do j=1,arrayWidth
      Do i=1,arrayHeight
        If(i.eq.1)Then
          decayChains(i,j,3) = 0
        Else
          Do k=1,size(decayInt,1)
            keyP = Int(decayChains(i-1,j,1))
            keyC = Int(decayChains(i,j,1))
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
! -----------------------------------------------------------------
  Recursive Function IsotopesDecayChainsR(z,a,m) RESULT (nodes)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: nodes
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: z, a, m
    Integer(kind=StandardInteger) :: zP, aP, mP
    Integer(kind=StandardInteger) :: zC, aC, mC
    Integer(kind=StandardInteger) :: key
    Logical :: unstable
! loop through decay chains
    unstable=.false.
    Do i=1,size(decayChar)
      zP = decayInt(i,2)
      aP = decayInt(i,1)
      mP = decayInt(i,3)
      If(z.eq.zP.and.a.eq.aP.and.m.eq.mP)Then
        unstable=.true.
        zC = decayInt(i,5)
        aC = decayInt(i,4)
        mC = decayInt(i,6)
        isotopeCounter = isotopeCounter + 1
        publicX = publicX + 1
        key = makeIsotopeKey(zP,aP,mP)
        public3DTempArray(publicX,publicY,1) = key
! store half life
! public3DTempArray(publicX,publicY,2) = decayDouble(i,2)
        publicZ = publicX + 1
        nodes = IsotopesDecayChainsR(zC,aC,mC)
        publicX = publicX - 1
      End If
    End Do
    If(unstable.eqv..false.)Then
      key = makeIsotopeKey(z,a,m)
      public3DTempArray(publicZ,publicY,1) = key
      publicY = publicY + 1
    End If
  End Function IsotopesDecayChainsR

  Function BranchProductionLoss&
    (parentProductionRate,decayDataArray,tStart,tEnd,tBeamEnd) RESULT (branchTally)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,decaySteps
    Real(kind=DoubleReal) :: parentProductionRate,tStart,tEnd,tBeamEnd
    Real(kind=DoubleReal) :: timeBPL,timeStepBPL
    Real(kind=DoubleReal) :: lastLoss,sourceAtoms,lossAtoms
    Real(kind=DoubleReal) :: lnTwo
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: branchTally
    Integer(kind=StandardInteger) :: iterations
! decayDataArray(1,1) = Key      !Key
! decayDataArray(1,2) = 100      !Np(tStart)
! decayDataArray(1,3) = 60           !T1/2p
! decayDataArray(1,4) = 0          !
! decayDataArray(2,1) = Key      !Key
! decayDataArray(2,2) = 500      !Na(tStart)
! decayDataArray(2,3) = 120          !T1/2a
! decayDataArray(2,4) = 0.8        !BF P to A
! set variables
    lnTwo = 0.693147180559945
! set iterations
    iterations = 100
! decay steps
    decaySteps = size(decayDataArray,1)
! allocate branch tally array
    Allocate(branchTally(1:decaySteps,1:2))
! fill tally array
    Do i=1,decaySteps
      branchTally(i,1) = decayDataArray(i,1)  !Key
      branchTally(i,2) = decayDataArray(i,2)  !No. atoms
! print *,branchTally(i,1),branchTally(i,2)
    End Do
! set time step
    timeStepBPL = (tEnd-tStart)/iterations
! set time
    timeBPL = tStart
! iterate from start to end time
    Do i=1,iterations
      Do j=1,decaySteps
        sourceAtoms = 0.0D0
        lossAtoms = 0.0D0
! source parent
        If(j.eq.1)Then  !parent
          If(timeBPL.lt.tBeamEnd)Then
            sourceAtoms = timeStepBPL * parentProductionRate
          End If
        End If
! source from previous isotope in decay chain
        If(j.gt.1)Then
          sourceAtoms = sourceAtoms + decayDataArray(j,4)*lastLoss !add bf * parent loss to source
        End If
        lastLoss = 0.0D0
! loss from decay
        If(decayDataArray(1,3).gt.0)Then
          lossAtoms = branchTally(j,2)*timeStepBPL*(lnTwo/decayDataArray(j,3))
        End If
! save last loss amount
        lastLoss = lossAtoms
        If(lastLoss.gt.branchTally(j,2))Then
          lastLoss = branchTally(j,2)
        End If
! tally with source and loss atoms
        branchTally(j,2) = branchTally(j,2) + sourceAtoms - lossAtoms
! check greater than or equal to 0
        If(branchTally(j,2).lt.0)Then
          branchTally(j,2) = 0.0D0
        End If
      End Do
      timeBPL = timeBPL + timeStepBPL
    End Do
! print
    print *,"Prod",parentProductionRate
    print *,"tStart",tStart
    print *,"tEnd",tEnd
    Do i=1,decaySteps
      print *,"S",i,branchTally(i,1),branchTally(i,2),decayDataArray(i,3),decayDataArray(i,4)
    End Do
  End Function BranchProductionLoss

! ------------------------------------------------------------------------!
! SUBROUTINE searchXSa
! Might need to develop, add 3 point interpolation
! ------------------------------------------------------------------------!

  Function searchXSa(energy, i) RESULT (xs)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,startKey,endKey,key,difference
    Real(kind=DoubleReal) :: xs, energy, tempEnergyL, tempEnergyU, factor
    Real(kind=DoubleReal) :: eA,oA,eB,oB,eM,oM
! set start and end points
    startKey = xsKey(i,7)                    !xsKey(i,7)   data row start
    endKey = xsKey(i,7) + xsKey(i,8) - 1     !xsKey(i,8)   data row length
    factor = 0.5
    difference = ceiling(factor*(endKey - startKey))
    print *,xsData(startKey,1),xsData(endKey,1)
    If(energy.lt.xsData(startKey,1))Then
      key = startKey
      eA = 0
      oA = 0
      eB = xsData(key,1)
      oB = xsData(key,2)
      eM = energy
      oM = oA + ((eM-eA)/(eB-eA))*(oB-oA)
      xs = oM
    ElseIf(energy.gt.xsData(endKey,1))Then
      key = endKey
      eA = xsData(key,1)
      oA = xsData(key,2)
      eB = 2*energy
      oB = 0
      eM = energy
      oM = oA + ((eM-eA)/(eB-eA))*(oB-oA)
      xs = oM
    Else
! start point
      key = startKey + difference
      Do j=1,10
! adjust key if too high or too low
        If(key.ge.endKey)Then
          key = endKey-1
        End If
        If(key.lt.startKey)Then
          key = startKey
        End If
        tempEnergyL = xsData(key,1)
        tempEnergyU = xsData(key+1,1)
        If(energy.ge.tempEnergyL.and.energy.le.tempEnergyU)Then
! energy bound found
! Linear interpolate between start and end energy
          xs = 1.0D0 * xsData(key,2) + &
          1.0D0 * ((energy-xsData(key,1))/(xsData(key+1,1)-xsData(key,1))) * &
          (xsData(key+1,2)-xsData(key,2))
          exit !break out of loop
        Else
          If(energy.lt.tempEnergyL)Then
            difference = ceiling(factor * difference)
            key = key - difference
! print *,"Decrease"
          End If
          If(energy.gt.tempEnergyU)Then
            difference = ceiling(factor * difference)
            key = key + difference
! print *,"Decrease"
          End If
        End If
      End Do
    End If
  End Function searchXSa

End Module productionLoss
