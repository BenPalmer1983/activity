Module prep

! Setup Modules
  Use kinds
  Use constants
  Use activityFunctions
  Use calcFunctions
  Use regression
  Use msubs
  Use globals
  Use initialise
  Use input

! force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! declare global variables

! Privacy of functions/subroutines/variables
  Private
  Public :: runPrep          !Subroutine

! Module Subroutines
  contains

! Run all the input subroutines

  Subroutine runPrep()
! Internal subroutine variables
    If(projectileZ.gt.0)Then
      Call fitExyz()
    Else
      Call neutronFit()
    End If
    Call M_synchProcesses()
    If(mpiProcessID.eq.0)Then
      Call materialIsotopes()
      Call calculateNumberDensity()
      Call createMaterial()
      Call addProductIsotopes()
      Call addDecayIsotopes()
      Call makeSimIsotopeKeys()
      Call makeReducedDecayList()
    End If
  End Subroutine runPrep

  Subroutine fitExyz()
! force declaration of all variables
    Implicit None
! Internal subroutine variables
    Integer(kind=StandardInteger) :: i,j,k,ionCount,totalIons,percentPrint
    Integer(kind=StandardInteger) :: intgSub, dataStart, dataLength, dataEnd
    Integer(kind=StandardInteger), Dimension( : ), Allocatable :: ionDataCount
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: specificIonExyzD
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyCoefficients
    Real(kind=DoubleReal), Dimension(1:polyFitOrder+1) :: coefficients
    Real(kind=DoubleReal) :: maxDepth, x, y
    Integer(kind=StandardInteger) :: order
! MPI
    Integer(kind=StandardInteger) :: mpiAssigned
! print if verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A22,F8.4)","    Fitting exyz data ",ProgramTime()
      print "(A22,I8)","    Total data points ",size(exyzKey,1)
    End If
! Init variables
    ionCount = 0
    intgSub = 0
    maxDepth = 0.0D0
    x = 0.0D0
    y = 0.0D0
    If(intgSub.eq.0)Then
! total ions in simulation
      totalIons = size(exyzIonKey,1)
! Set starting values
      percentPrint = 0
! Set coefficient arrays
! set the order of the polynomials to fit
      order = polyFitOrder
! allocate the array
      Allocate(fitCoefficients(0:polyFitOrder))
      Allocate(polyCoefficients(0:polyFitOrder))
      fitCoefficients = 0.0D0  !Fill fitCoefficients with zeros
      fitCoefficientsMpi = 0.0D0
! Loop through ions
      Do i=1,totalIons
! print percentage progress if verbose on
        If(verboseTerminal.eqv..true.)Then
          If((i-1).eq.0.and.percentPrint.eq.0)Then
            print "(A19,I8,A5,F8.4)","    0% exyz fit of ",totalIons," ions",ProgramTime()
            percentPrint = percentPrint + 1
          End If
          If((i-1).ge.(0.25*totalIons).and.percentPrint.eq.1)Then
            print "(A20,I8,A5,F8.4)","    25% exyz fit of ",totalIons," ions",ProgramTime()
            percentPrint = percentPrint + 1
          End If
          If((i-1).ge.(0.50*totalIons).and.percentPrint.eq.2)Then
            print "(A20,I8,A5,F8.4)","    50% exyz fit of ",totalIons," ions",ProgramTime()
            percentPrint = percentPrint + 1
          End If
          If((i-1).ge.(0.75*totalIons).and.percentPrint.eq.3)Then
            print "(A20,I8,A5,F8.4)","    75% exyz fit of ",totalIons," ions",ProgramTime()
            percentPrint = percentPrint + 1
          End If
        End If
! set assigned process id
        mpiAssigned = mod((i-1),mpiProcessCount)
        If(mpiAssigned.eq.mpiProcessID)Then
! Get start end points for data
          dataStart = exyzIonKey(i,1)
          dataLength = exyzIonKey(i,2)
          dataEnd = dataStart + dataLength - 1
! Deallocate arrays
          If(Allocated(specificIonExyzD))Then
            Deallocate(specificIonExyzD)
          End If
          Allocate(specificIonExyzD(1:dataLength,1:2))
! Loop through data points of ion
          k = 0
          Do j=dataStart,dataEnd
            k = k + 1
            specificIonExyzD(k,1) = exyzData(j,1)
            specificIonExyzD(k,2) = exyzData(j,2)
          End Do
! Fit polynomial to the curve
          polyCoefficients = PolyFit(specificIonExyzD,order)
! Add coefficients
          Do j=0,order
            fitCoefficients(j) = fitCoefficients(j) + polyCoefficients(j)
          End Do
        End If
      End Do
! Merge fitCoefficients from all processes
! Use fixed size array
      Do j=0,order
        fitCoefficientsMpi(j+1) = fitCoefficients(j)
      End Do
! Collect, sum and send out array to all processes
! Call MPI_sumData1DDP(fitCoefficientsMpi,10)
      Call M_sumDataDouble1D(fitCoefficientsMpi)
! Reassign to fitCoefficients array
      Do j=0,order
        fitCoefficients(j) = fitCoefficientsMpi(j+1)
      End Do
! Adjust coefficients
      Do j=0,order
        fitCoefficients(j) = fitCoefficients(j) / totalIons
      End Do
    Else
! old integrate sub
      print *,"old"
      Allocate(specificIonExyzD(1:ionDataCount(ionCount),1:2))
    End If
! print if verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A21,I8,A5,F8.4)","    100% exyz fit of ",totalIons," ions",ProgramTime()
    End If
! On root process only
    If(mpiProcessID.eq.0)Then
! Save trajectory to file
      Do i=1,polyFitOrder+1
        coefficients(i) = fitCoefficients(i-1)
      End Do
! Calc max depth
      trajDepth = MaxTrajDepth(coefficients)
      If(verboseTerminal.eqv..true.)Then
        print *,"   Max depth/ang: ",trajDepth
      End If
      Open(unit=961,file=Trim(outputDirectory)//"/"//"ionTraj.dat",&
      status="old",position="append",action="write")
      Do i=0,1000
        x = (i/1000.0D0)*trajDepth
        y = CalcPolynomial(coefficients, x)
        Write(961,"(E18.8,A4,E18.8)") x,"    ",y
      End Do
      Close(961)
    End If
  End Subroutine fitExyz

  Subroutine neutronFit()
! force declaration of all variables
    Implicit None
! Internal subroutine variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: x, y
    Real(kind=DoubleReal), Dimension(1:polyFitOrder+1) :: coefficients
! Set fit coefficients - assumes target too thin to change neutron energy
! Same energy through target
    Allocate(fitCoefficients(0:polyFitOrder))
    fitCoefficients = 0.0D0
    fitCoefficients(0) = 1.0D0*beamEnergy
! On root process only
    If(mpiProcessID.eq.0)Then
! Save trajectory to file
      Do i=1,polyFitOrder+1
        coefficients(i) = fitCoefficients(i-1)
      End Do
! Calc max depth
      trajDepth = MaxTrajDepth(coefficients)
      If(verboseTerminal.eqv..true.)Then
        print *,"   Max depth/ang: ",trajDepth
      End If
      Open(unit=961,file=Trim(outputDirectory)//"/"//"ionTraj.dat",&
      status="old",position="append",action="write")
      Do i=0,1000
        x = (i/1000.0D0)*trajDepth
        y = CalcPolynomial(coefficients, x)
        Write(961,"(E18.8,A4,E18.8)") x,"    ",y
      End Do
      Close(961)
    End If
  End Subroutine neutronFit

  Subroutine materialIsotopes()
! force declaration of all variables
    Implicit None
! Internal subroutine variables
    Integer :: i, j, k
    Real(kind=DoubleReal) :: averageAtomicMass, denominator
! print if verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A30,F8.4)","    Calc isotopic composition ",ProgramTime()
    End If
    k = 0
    denominator = 0.0
! loop through elements in material
    Do i=1,size(elements)
! loop through stable isotopes
      averageAtomicMass = 0.0
      Do j=1,size(isotopesChar)
        If(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.001)Then
          k = k + 1
          averageAtomicMass = averageAtomicMass + isotopesReal(j,2) * isotopesReal(j,1)
        End If
      End Do
      denominator = denominator + materialComposition(i) / averageAtomicMass
! print *,elements(i),averageAtomicMass,denominator
    End Do
! allocate arrays
    Allocate(materialIsotopesChar(1:k))
    Allocate(materialIsotopesInt(1:k,1:2))
    Allocate(materialIsotopesReal(1:k,1:3))
    k = 0
! loop through elements in material
    Do i=1,size(elements)
! loop through stable isotopes
      Do j=1,size(isotopesChar)
        If(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.001)Then
          k = k + 1
          materialIsotopesChar(k) = isotopesChar(j)      !isotope/element symbol
          materialIsotopesInt(k,1) = isotopesInt(j,1)    !isotope Z (atomic number)
          materialIsotopesInt(k,2) = isotopesInt(j,2)    !isotope A (mass number)
          materialIsotopesReal(k,1) = isotopesReal(j,1)    !Relative atomic mass
          materialIsotopesReal(k,2) = isotopesReal(j,2)                !Isotopic Composition
          materialIsotopesReal(k,3) = isotopesReal(j,2) * &
          ((materialComposition(i)/isotopesReal(j,1))/denominator)  !Material Composition by number
! print *,materialIsotopesChar(k),materialIsotopesInt(k,2),materialIsotopesReal(k,3)
        End If
      End Do
    End Do
  End Subroutine materialIsotopes

! ------------------------------------------------------------------------!
! SUBROUTINE calculateNumberDensity
! Calculate number of atoms per cubic metre of material
! ------------------------------------------------------------------------!
  Subroutine calculateNumberDensity()
! atoms per cubic metre
! force declaration of all variables
    Implicit None
! Internal subroutine variables
    Integer :: i,j
    Real(kind=DoubleReal) :: fractionComposition, atomicMass
! loop through elements
    atomicMass = 0.0
    Do i=1,size(elements)
! loop through stable isotopes
      Do j=1,size(isotopesChar)
        If(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.0)Then
          fractionComposition = materialComposition(i) * isotopesReal(j,2)
          atomicMass = atomicMass + fractionComposition * isotopesReal(j,1)
        End If
      End Do
    End Do
! Calculate number density
    numberDensity = ((materialDensity * 1000) / atomicMass) * avogadrosConstant    !density from kgm-3 to gm-3
  End Subroutine calculateNumberDensity

! ------------------------------------------------------------------------
! SUBROUTINE createMaterial
! Make tally arrays for material and mark starting isotopes/amounts
! ------------------------------------------------------------------------

  Subroutine createMaterial()
! force declaration of all variables
    Implicit None
! Internal subroutine variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: z, a, m
    Integer(kind=StandardInteger) :: key, startingIsotopeCount
! print if verbose on
    If(verboseTerminal.eqv..true.)Then
      print "(A28,F8.4)","    Prepare isotope tallies ",ProgramTime()
    End If
    If(mpiProcessID.eq.0)Then
! open output file
      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
      write(999,"(A1)") " "
      write(999,"(A40,F8.4)") "Create material tally                   ",ProgramTime()
    End If
! Make Isotope Tallys
    Allocate(isotopeTallyChar(1:70800))
    Allocate(isotopeTallyInt(1:70800,1:5))
    Allocate(isotopeTallyActive(1:70800,1:7))
! Fill with blank data
    If(mpiProcessID.eq.0)Then
      write(999,"(A40,F8.4)") "Fill tally with blank data              ",ProgramTime()
    End If
    Do i=1,70800
      isotopeTallyChar(i) = "ZZ"
      isotopeTallyInt(i,1) = 0        !Z
      isotopeTallyInt(i,2) = 0        !A
      isotopeTallyInt(i,3) = 0        !M
      isotopeTallyInt(i,4) = 0        !Decay calc marker
      isotopeTallyInt(i,5) = 0        !In sim
      isotopeTallyActive(i,1) = -1.0D0    !half life of isotope (-1 if stable)
      isotopeTallyActive(i,2) = 0.0D0      !activity of isotope in bq
      isotopeTallyActive(i,3) = 0.0D0    !beam creation rate
      isotopeTallyActive(i,4) = 0.0D0    !atom tally
      isotopeTallyActive(i,5) = 0.0D0    !atom tally (starting)
      isotopeTallyActive(i,6) = 0.0D0    !decay constant
      isotopeTallyActive(i,7) = 0.0D0    !atom tally (end beam)
    End Do
! Fill with isotope data
    If(mpiProcessID.eq.0)Then
      write(999,"(A40,F8.4)") "Fill tally with isotope data            ",ProgramTime()
    End If
! Add neutron
    Do j=0,1
      key = makeIsotopeKey(0,1,j)
      isotopeTallyChar(key) = "NN"
      isotopeTallyInt(key,1) = 0
      isotopeTallyInt(key,2) = 1
      isotopeTallyInt(key,3) = j
    End Do
! Add other isotopes
    Do i=1,size(isotopesChar)
      Do j=0,1
! key = 590 * isotopesInt(i,1) + 290 * j + isotopesInt(i,2)
        key = makeIsotopeKey(isotopesInt(i,1),isotopesInt(i,2),j)
        isotopeTallyChar(key) = elementSymbol(isotopesInt(i,1))
        isotopeTallyInt(key,1) = isotopesInt(i,1)
        isotopeTallyInt(key,2) = isotopesInt(i,2)
        isotopeTallyInt(key,3) = j
        z = isotopeTallyInt(key,1)
        a = isotopeTallyInt(key,2)
        m = isotopeTallyInt(key,3)
        Do k=1,size(decayInt,1)
          If(decayInt(k,1).eq.a.and.decayInt(k,2).eq.z.and.decayInt(k,3).eq.m)then  !decayInt A Z wrong way around
            isotopeTallyActive(key,1) = decayDouble(k,2)
            isotopeTallyActive(key,6) = 1.0D0 * (lnTwo / decayDouble(k,2))
            exit
          End If
        End Do
      End Do
    End Do
! Set the composition amounts
    startingIsotopeCount = 0
    Do i=1,size(materialIsotopesChar)
      key = makeIsotopeKey(materialIsotopesInt(i,1),materialIsotopesInt(i,2),0)
! atoms in 1m3
      If(key.gt.0.and.key.le.70800)Then
        isotopeTallyActive(key,4) = 1.0D0 * materialIsotopesReal(i,3)&
        * materialDensity * 1.0D3 * (1.0D0 / materialIsotopesReal(i,1))&
        * avogadrosConstant
! Mark isotope as in simulation
        isotopeTallyInt(key,5) = 1
        startingIsotopeCount = startingIsotopeCount + 1
      End If
    End Do
! save starting isotope tally to output file
    If(mpiProcessID.eq.0)Then
      write(999,"(A1)") " "
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
      write(999,"(A39)") "Starting Isotope Tally - Input Material"
      write(999,"(A8,A4,A4,A2,A4,A4,A18,A18,A18,A18,A21)") &
      "Element ",&
      "Z   ","A   ","M ",&
      "mk  ","sim ",&
      "Half life         ","Decay Constant    ","Activity          ","Reaction Rate     ",&
      "Atoms/mg            "
      write(999,"(A60,A80)") "------------------------------------------------------------",&
      "--------------------------------------------------------------------------------"
    End If
    Do i=1,size(isotopeTallyChar)
      If(isotopeTallyInt(i,5).eq.1)Then
        If(mpiProcessID.eq.0)Then
          write(999,&
          "(A7,A1,I3.3,A1,I3.3,A1,I1.1,A1,I3.3,A1,I3.3,A1,d17.10,A1,d17.10,A1,d17.10,A1,d17.10,A1,d17.10)") &
          isotopeTallyChar(i)," ",&
          isotopeTallyInt(i,1)," ",&
          isotopeTallyInt(i,2)," ",isotopeTallyInt(i,3)," ",&
          isotopeTallyInt(i,4)," ",isotopeTallyInt(i,5)," ",&
          isotopeTallyActive(i,1)," ",isotopeTallyActive(i,6)," ",&
          isotopeTallyActive(i,2)," ",&
          isotopeTallyActive(i,3)," ",isotopeTallyActive(i,4)
        End If
      End If
    End Do
    If(mpiProcessID.eq.0)Then
      write(999,"(A25,I8)") "Total starting isotopes: ",startingIsotopeCount
    End If
! initialise tally
    totalSimulationAtoms = 0.0D0
! close output file
    If(mpiProcessID.eq.0)Then
      write(999,"(A1)") " "
      close(999)
    End If
  End Subroutine createMaterial

! ------------------------------------------------------------------------
! SUBROUTINE addProductIsotopes
! Make tally arrays for material and mark starting isotopes/amounts
! ------------------------------------------------------------------------

  Subroutine addProductIsotopes()
! force declaration of all variables
    Implicit None
! Internal subroutine variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: keyP
!    If(mpiProcessID.eq.0)Then
! open output file
!      open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! write to output file
!      write(999,"(A1)") " "
!      write(999,"(A40,F8.4)") "Add Product Isotopes to Tally           ",ProgramTime()
!    End If
! output tally to file
! Add product isotopes to tally
    Do i=1,size(xsKey,1)
! Get key and data for cross section
      keyP = makeIsotopeKey(xsKey(i,4),xsKey(i,5),xsKey(i,6))
      If(keyP.gt.0.and.keyP.le.70800)Then
        isotopeTallyInt(keyP,5) = 1
      End If
    End Do
!    If(mpiProcessID.eq.0)Then
! output tally to file
!      write(999,"(A1)") " "
! close output file
!      close(999)
!    End If
  End Subroutine addProductIsotopes

! ------------------------------------------------------------------------
! SUBROUTINE addDecayIsotopes
! Make tally arrays for material and mark starting isotopes/amounts
! ------------------------------------------------------------------------

  Subroutine addDecayIsotopes()
! force declaration of all variables
    Implicit None
! Internal subroutine variables
    Integer(kind=StandardInteger) :: i,j
    Integer(kind=StandardInteger) :: z,a,m
    Integer(kind=StandardInteger) :: key
    Real(kind=DoubleReal), Dimension(:,:), Allocatable :: isotopeArray
    Do i=1,size(isotopeTallyChar)
      If(isotopeTallyInt(i,5).eq.1.and.isotopeTallyActive(i,1).gt.(-1.0D0))Then
        isotopeArray = IsotopesInDecayTree(isotopeTallyInt(i,1),&
        isotopeTallyInt(i,2),isotopeTallyInt(i,3))
        Do j=1,size(isotopeArray,1)
          z = Int(isotopeArray(j,1))
          a = Int(isotopeArray(j,2))
          m = Int(isotopeArray(j,3))
          key = makeIsotopeKey(z,a,m)
          If(key.gt.0.and.key.le.70800)Then
            isotopeTallyInt(key,5) = 1
          End If
        End Do
      End If
    End Do
  End Subroutine addDecayIsotopes

! ------------------------------------------------------------------------
! SUBROUTINE makeSimIsotopeKeys
! Make tally arrays for material and mark starting isotopes/amounts
! ------------------------------------------------------------------------

  Subroutine makeSimIsotopeKeys()
! force declaration of all variables
    Implicit None
! Internal subroutine variables
    Integer(kind=StandardInteger) :: i,j
    Integer(kind=StandardInteger) :: key, simIsotopesCount
! Count sim isotopes
    simIsotopesCount = 0
    Do i=1,size(isotopeTallyChar,1)
      If(isotopeTallyInt(i,5).eq.1)Then
! print *,i,isotopeTallyInt(i,1),isotopeTallyInt(i,2),isotopeTallyInt(i,3)
        simIsotopesCount = simIsotopesCount + 1
      End If
    End Do
! Allocate array
    Allocate(simIsotopeKeys(1:simIsotopesCount))
    j = 0
! store sim isotope keys
    Do i=1,size(isotopeTallyChar,1)
      If(isotopeTallyInt(i,5).eq.1)Then
        key = makeIsotopeKey(isotopeTallyInt(i,1),isotopeTallyInt(i,2),isotopeTallyInt(i,3))
        If(key.gt.0.and.key.le.70800)Then
          j = j + 1
          simIsotopeKeys(j) = key
        End If
      End If
    End Do
  End Subroutine makeSimIsotopeKeys

! ------------------------------------------------------------------------!
! SUBROUTINE makeReducedDecayList
! Reduce the size of the decay list
! ------------------------------------------------------------------------!
  Subroutine makeReducedDecayList()
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,j,k,decayReducedCount,key
    Integer, Dimension( : , : ), Allocatable :: decayIntReduced
    Double Precision, Dimension( : , : ), Allocatable :: decayDoubleReduced
! If(mpiProcessID.eq.0)Then
! open output file
!  open(unit=999,file=trim(fileOutputData),status="old",position="append",action="write")
! End If
! count reduced array size
    decayReducedCount = 0
! loop through isotopes in tally
    Do i=1,size(decayInt,1)
      key=makeIsotopeKey(decayInt(i,2),decayInt(i,1),decayInt(i,3))
      If(isotopeTallyInt(key,5).eq.1.and.isotopeTallyActive(key,1).ge.0)Then
        decayReducedCount = decayReducedCount + 1
      End If
    End Do
! write to output file
! If(mpiProcessID.eq.0)Then
!  write(999,"(A1)") " "
!  write(999,"(A30,I4,A4,I4,A10,F8.4)") "Decay data array reduced from ",size(decayInt,1),&
!  " to ",decayReducedCount,"          ",ProgramTime()
! End If
! Allocate arrays
    Allocate(decayIntReduced(1:decayReducedCount,1:6))
    Allocate(decayDoubleReduced(1:decayReducedCount,1:2))
! loop through isotopes in tally
    j = 0
    Do i=1,size(decayInt,1)
      key=makeIsotopeKey(decayInt(i,2),decayInt(i,1),decayInt(i,3))
      If(key.gt.0.and.key.le.70800)Then
        If(isotopeTallyInt(key,5).eq.1.and.isotopeTallyActive(key,1).ge.0)Then
          j = j + 1
          Do k=1,6
            decayIntReduced(j,k) = decayInt(i,k)
          End Do
          Do k=1,2
            decayDoubleReduced(j,k) = decayDouble(i,k)
          End Do
        End If
      End If
    End Do
! write to file
! If(mpiProcessID.eq.0)Then
!  Do i=1,size(decayIntReduced,1)
!    write(999,"(I8,A2,I8,I8,I8,A2,I8,I8,I8,E20.10,E20.10)") &
!    i,"  ",decayIntReduced(i,1),decayIntReduced(i,2),decayIntReduced(i,3),&
!    "  ",decayIntReduced(i,4),decayIntReduced(i,5),decayIntReduced(i,6),&
!    decayDoubleReduced(i,1),decayDoubleReduced(i,2)
!  End Do
! End If
! Loop and transfer
! Do i=1,decayReducedCount
! key=makeIsotopeKey(decayInt(i,2),decayInt(i,1),decayInt(i,3))
! Do k=1,6
! decayInt(i,k) = decayIntReduced(i,k)
! End Do
! Do k=1,2
! decayDouble(i,k) = decayDoubleReduced(i,k)
! End Do
! End Do
! close output file
! If(mpiProcessID.eq.0)Then
! write(999,"(A1)") " "
! close(999)
! End If
  End Subroutine makeReducedDecayList

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
    If(Allocated(decayChainIsotopesArr))Then
      Deallocate(decayChainIsotopesArr)
    End If
    Allocate(decayChainIsotopesArr(1:300,1:3))
! fill with blank data
    Do i=1,300
      Do j=1,3
        decayChainIsotopesArr(i,j) = 0
      End Do
    End Do
! set counter = 1 and store parent isotope
    decayChainIsotopesArr(1,1) = z
    decayChainIsotopesArr(1,2) = a
    decayChainIsotopesArr(1,3) = m
    isotopeCounterA = 1
! call function recursively
    complete = IsotopesInDecayTreeR(z,a,m)
! transfer data to new array
    Allocate(isotopeArray(1:isotopeCounterA,1:3))
    Do i=1,isotopeCounterA
      Do j=1,3
        isotopeArray(i,j) = decayChainIsotopesArr(i,j)
      End Do
    End Do
! Deallocate array
    If(Allocated(decayChainIsotopesArr))Then
      Deallocate(decayChainIsotopesArr)
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
! Initialise
    lastRow = 0
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
        Do j=1,size(decayChainIsotopesArr,1)
          If(decayChainIsotopesArr(j,1).eq.zC.and.&
            decayChainIsotopesArr(j,2).eq.aC.and.&
            decayChainIsotopesArr(j,3).eq.mC)Then
            store = .false.
          End If
        End Do
        If(store.eqv..true.)Then
          isotopeCounterA = isotopeCounterA + 1
          decayChainIsotopesArr(isotopeCounterA,1) = zC
          decayChainIsotopesArr(isotopeCounterA,2) = aC
          decayChainIsotopesArr(isotopeCounterA,3) = mC
        End If
        lastRow = IsotopesInDecayTreeR(zC,aC,mC)
      End If
    End Do
  End Function IsotopesInDecayTreeR

End Module prep
