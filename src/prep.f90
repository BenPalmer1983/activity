Module prep

! Setup Modules
  Use kinds
  Use constants
  Use stringfunctions		!string functions
  Use maths
  Use initialise				
  Use input


!force declaration of all variables
  Implicit None
  
!declare global variables  
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: fitCoefficients
  Character(len=2), Dimension( : ), Allocatable :: materialIsotopesChar
  Integer, Dimension( : , : ), Allocatable :: materialIsotopesInt
  Real, Dimension( : , : ), Allocatable :: materialIsotopesReal
  Character(len=2), Dimension( : ), Allocatable :: isotopeTallyChar
  Integer, Dimension( : , : ), Allocatable :: isotopeTallyInt
  Real(kind=DoubleReal) :: totalSimulationAtoms
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeTallyActive
  Integer(kind=StandardInteger), Dimension( : , :), Allocatable :: decayChainIsotopesArr
  Integer(kind=StandardInteger) :: isotopeCounterA
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: simIsotopeKeys
  
  
!Privacy of functions/subroutines/variables
  Private
  Public :: runPrep					!Subroutine
  Public :: fitCoefficients			!Variable
  Public :: materialIsotopesChar	!Variable
  Public :: materialIsotopesInt		!Variable
  Public :: materialIsotopesReal	!Variable
  Public :: isotopeTallyChar		!Variable
  Public :: isotopeTallyInt			!Variable
  Public :: isotopeTallyActive      !Variable
  Public :: totalSimulationAtoms	!Variable
  Public :: decayChainIsotopesArr	!Variable	
  Public :: isotopeCounterA	        !Variable	
  Public :: simIsotopeKeys	        !Variable	
  
!Module Subroutines  
contains 

!Run all the input subroutines

  Subroutine runPrep()
	
	!Internal subroutine variables
	Integer :: i, j, k
	
	Call fitExyz()
	Call materialIsotopes()
	Call calculateNumberDensity()
	Call createMaterial()
	Call addProductIsotopes()
	Call addDecayIsotopes()
	Call makeSimIsotopeKeys()
	Call makeReducedDecayList()

  End Subroutine runPrep

  
  Subroutine fitExyz()
    
!force declaration of all variables
	Implicit None

!Internal subroutine variables
	Integer :: i, j, k, ionCount, dataPointCount, totalIons
	Integer, Dimension( : ), Allocatable :: ionDataCount
	Real, Dimension( : , : ), Allocatable :: specificIonExyz	
	double precision, Dimension( : , : ), Allocatable :: specificIonExyzD
	double precision, Dimension( : ), Allocatable :: polyCoefficients
	integer :: order
!number of ions
	ionCount = 0
	do i=1,size(exyzKey)
	  if(i.eq.1)then
	    ionCount = ionCount + 1  
	  else
	    if(exyzKey(i-1).ne.exyzKey(i))then
		  ionCount = ionCount + 1
		endif
	  endif
	enddo
	totalIons = ionCount
!Allocate arrays
	Allocate(ionDataCount(1:totalIons))
!save ion counts
	ionCount = 0
	do i=1,size(exyzKey)
	  if(i.eq.1)then
	    ionCount = ionCount + 1 !set ion
	    dataPointCount = 0 !set counter
	  else
	    if(exyzKey(i-1).ne.exyzKey(i))then
          ionDataCount(ionCount) = dataPointCount !save data point count
		  dataPointCount = 0 !reset counter
		  ionCount = ionCount + 1 !set ion
		endif
	  endif
	  dataPointCount = dataPointCount + 1	!increment counter
	  if(i.eq.size(exyzKey)) then
        ionDataCount(ionCount) = dataPointCount !save (final) data point count
	  endif	
	enddo
!set the order of the polynomials to fit
	order = polyFitOrder
!allocate the array
	Allocate(fitCoefficients(0:order))
	Allocate(polyCoefficients(0:order))
!Fill fitCoefficients with zeros
	do j=0,order
	  fitCoefficients(j) = 0
	enddo
	ionCount = 0
	do i=1,size(exyzKey)
	  if(i.eq.1)then
	    ionCount = ionCount + 1 !set ion
	    dataPointCount = 0 !set counter
		Allocate(specificIonExyz(1:ionDataCount(ionCount),1:2))
		Allocate(specificIonExyzD(1:ionDataCount(ionCount),1:2))		
	  else
	    if(exyzKey(i-1).ne.exyzKey(i))then
!Process data before next ion
		  polyCoefficients = PolyFit(specificIonExyzD,order)
		  do j=0,order
		    fitCoefficients(j) = fitCoefficients(j) + (polyCoefficients(j) / (1.0 * totalIons))
		  enddo
!Next ion
		  dataPointCount = 0 !reset counter
		  ionCount = ionCount + 1 !set ion
!deallocate and allocate array
		  Deallocate(specificIonExyz)
		  Allocate(specificIonExyz(1:ionDataCount(ionCount),1:2))
		  Deallocate(specificIonExyzD)
		  Allocate(specificIonExyzD(1:ionDataCount(ionCount),1:2))
		endif
	  endif	  
!increment counter
	  dataPointCount = dataPointCount + 1	
!Store data
	  specificIonExyz(dataPointCount,1) = exyzData(i,1)
	  specificIonExyz(dataPointCount,2) = exyzData(i,2)
	  specificIonExyzD(dataPointCount,1) = 1.0D0*exyzData(i,1)
	  specificIonExyzD(dataPointCount,2) = 1.0D0*exyzData(i,2)
	  if(i.eq.size(exyzKey)) then
	    !Process data for last ion
		!Process specificIonExyz array
		!Call polyFit(specificIonExyz,order,polyCoefficients)
		polyCoefficients = PolyFit(specificIonExyzD,order)
		do j=0,order
		  fitCoefficients(j) = fitCoefficients(j) + (polyCoefficients(j) / (1.0 * totalIons))
		enddo  
	  endif	
	enddo
	  
  End Subroutine fitExyz
  
  
  
  
  Subroutine materialIsotopes()

!force declaration of all variables
	Implicit None

!Internal subroutine variables
    Integer :: i, j, k
	Real :: averageAtomicMass, denominator

	k = 0
	denominator = 0.0
	!loop through elements in material
	do i=1,size(elements)
	  !loop through stable isotopes
	  averageAtomicMass = 0.0
	  do j=1,size(isotopesChar)
	    if(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.001)then
	      k = k + 1
		  averageAtomicMass = averageAtomicMass + isotopesReal(j,2) * isotopesReal(j,1)
	    endif
      enddo	
	  denominator = denominator + materialComposition(i) / averageAtomicMass
	  !print *,elements(i),averageAtomicMass,denominator 
	enddo	
	
	
	!allocate arrays
	Allocate(materialIsotopesChar(1:k))
	Allocate(materialIsotopesInt(1:k,1:2))
	Allocate(materialIsotopesReal(1:k,1:3))	
    
	k = 0
	!loop through elements in material
	do i=1,size(elements)
	  !loop through stable isotopes
	  do j=1,size(isotopesChar)
	    if(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.001)then
	      k = k + 1
		  materialIsotopesChar(k) = isotopesChar(j)			!isotope/element symbol
		  materialIsotopesInt(k,1) = isotopesInt(j,1)		!isotope Z (atomic number)
		  materialIsotopesInt(k,2) = isotopesInt(j,2)		!isotope A (mass number)
		  materialIsotopesReal(k,1) = isotopesReal(j,1)		!Relative atomic mass
		  materialIsotopesReal(k,2) = isotopesReal(j,2)								!Isotopic Composition
		  materialIsotopesReal(k,3) = isotopesReal(j,2) * &
		  ((materialComposition(i)/isotopesReal(j,1))/denominator)	!Material Composition by number
		  !print *,materialIsotopesChar(k),materialIsotopesInt(k,2),materialIsotopesReal(k,3)
	    endif
      enddo	
	enddo	
  
  End Subroutine materialIsotopes
  
    
 

!------------------------------------------------------------------------!
! SUBROUTINE calculateNumberDensity 
! Calculate number of atoms per cubic metre of material
!------------------------------------------------------------------------! 
  Subroutine calculateNumberDensity()      
!atoms per cubic metre	 
!force declaration of all variables
	Implicit None
!Internal subroutine variables
    Integer :: i, j, k
	Real :: fractionComposition, atomicMass
!loop through elements
	atomicMass = 0.0
    do i=1,size(elements)
!loop through stable isotopes
	  do j=1,size(isotopesChar)
	    if(elements(i).eq.isotopesChar(j).and.isotopesReal(j,2).gt.0.0)then
          fractionComposition = materialComposition(i) * isotopesReal(j,2)
		  atomicMass = atomicMass + fractionComposition * isotopesReal(j,1)
		endif
	  enddo
	enddo
!Calculate number density  
	numberDensity = ((materialDensity * 1000) / atomicMass) * avogadrosConstant		!density from kgm-3 to gm-3
  End Subroutine calculateNumberDensity
  
  
!------------------------------------------------------------------------
! SUBROUTINE createMaterial                                                    
! Make tally arrays for material and mark starting isotopes/amounts                                    
!------------------------------------------------------------------------
  
  Subroutine createMaterial()      
!force declaration of all variables
	Implicit None
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: z, a, m
	Integer(kind=StandardInteger) :: key, startingIsotopeCount
	Real :: tempReal
	Double Precision :: tempDouble
!open output file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	write(999,"(A1)") " "
	write(999,"(A40,F8.4)") "Create material tally                   ",ProgramTime()
!Make Isotope Tallys
    Allocate(isotopeTallyChar(1:70800))
    Allocate(isotopeTallyInt(1:70800,1:5))
	Allocate(isotopeTallyActive(1:70800,1:7))
!Fill with blank data	
	write(999,"(A40,F8.4)") "Fill tally with blank data              ",ProgramTime()
	do i=1,70800
	  isotopeTallyChar(i) = "ZZ"
	  isotopeTallyInt(i,1) = 0				!Z
	  isotopeTallyInt(i,2) = 0				!A
	  isotopeTallyInt(i,3) = 0				!M
	  isotopeTallyInt(i,4) = 0				!Decay calc marker
	  isotopeTallyInt(i,5) = 0				!In sim 
	  isotopeTallyActive(i,1) = -1.0D0		!half life of isotope (-1 if stable)
	  isotopeTallyActive(i,2) = 0.0D0    	!activity of isotope in bq
	  isotopeTallyActive(i,3) = 0.0D0		!beam creation rate
	  isotopeTallyActive(i,4) = 0.0D0		!atom tally
	  isotopeTallyActive(i,5) = 0.0D0		!atom tally (starting)
	  isotopeTallyActive(i,6) = 0.0D0		!decay constant
	  isotopeTallyActive(i,7) = 0.0D0		!atom tally (end beam)
	enddo	
!Fill with isotope data		
    write(999,"(A40,F8.4)") "Fill tally with isotope data            ",ProgramTime()
!Add neutron
    Do j=0,1
      key = makeIsotopeKey(0,1,j)
	  isotopeTallyChar(key) = "NN"
	  isotopeTallyInt(key,1) = 0
	  isotopeTallyInt(key,2) = 1
	  isotopeTallyInt(key,3) = j	
	End Do
!Add other isotopes
    do i=1,size(isotopesChar)
	  do j=0,1
	    !key = 590 * isotopesInt(i,1) + 290 * j + isotopesInt(i,2)
		key = makeIsotopeKey(isotopesInt(i,1),isotopesInt(i,2),j)
	    isotopeTallyChar(key) = elementSymbol(isotopesInt(i,1))
	    isotopeTallyInt(key,1) = isotopesInt(i,1)
	    isotopeTallyInt(key,2) = isotopesInt(i,2)
	    isotopeTallyInt(key,3) = j	
		z = isotopeTallyInt(key,1)
		a = isotopeTallyInt(key,2)
		m = isotopeTallyInt(key,3)
		do k=1,size(decayInt,1)		
		  if(decayInt(k,1).eq.a.and.decayInt(k,2).eq.z.and.decayInt(k,3).eq.m)then  !decayInt A Z wrong way around
		    isotopeTallyActive(key,1) = decayDouble(k,2)	
            isotopeTallyActive(key,6) = 1.0D0 * (lnTwo / decayDouble(k,2))
			exit
		  endif
		enddo			
	  enddo
	enddo
!Set the composition amounts
	startingIsotopeCount = 0
    do i=1,size(materialIsotopesChar)
	  !key = 590 * materialIsotopesInt(i,1) + 290 * 0 + materialIsotopesInt(i,2)
	  key = makeIsotopeKey(materialIsotopesInt(i,1),materialIsotopesInt(i,2),0)
!atoms in 1m3	  
	  isotopeTallyActive(key,4) = 1.0D0 * materialIsotopesReal(i,3)&
        * materialDensity * 1.0D3 * (1.0D0 / materialIsotopesReal(i,1))&
        * avogadrosConstant
!Mark isotope as in simulation
	  isotopeTallyInt(key,5) = 1
	  startingIsotopeCount = startingIsotopeCount + 1
	enddo	
!save starting isotope tally to output file
	write(999,"(A1)") " "	
	write(999,"(A140)") "-----------------------------------------------------------------&
	---------------------------------------------------------------------------"
	write(999,"(A39)") "Starting Isotope Tally - Input Material"
	write(999,"(A8,A4,A4,A2,A4,A4,A18,A18,A18,A18,A21)") &
	"Element ",&
	"Z   ","A   ","M ",& 
	"mk  ","sim ",& 
	"Half life         ","Decay Constant    ","Activity          ","Reaction Rate     ",&
	"Atoms/mg            "
	write(999,"(A140)") "-----------------------------------------------------------------&
	---------------------------------------------------------------------------"
	Do i=1,size(isotopeTallyChar)
      !if(isotopeTallyChar(i).ne."ZZ".and.isotopeTallyInt(i,5).eq.1)then
      !if(isotopeTallyChar(i).ne."ZZ")then
	  If(isotopeTallyInt(i,5).eq.1)Then
	    write(999,"(A7,A1,&
		I3.3,A1,I3.3,A1,I1.1,A1,&
		I3.3,A1,I3.3,A1,&
		d17.10,A1,d17.10,A1,d17.10,A1,d17.10,A1,d17.10)") &
		isotopeTallyChar(i)," ",&
		isotopeTallyInt(i,1)," ",&
		isotopeTallyInt(i,2)," ",isotopeTallyInt(i,3)," ",&
		isotopeTallyInt(i,4)," ",isotopeTallyInt(i,5)," ",&
		isotopeTallyActive(i,1)," ",isotopeTallyActive(i,6)," ",&
		isotopeTallyActive(i,2)," ",&
		isotopeTallyActive(i,3)," ",isotopeTallyActive(i,4)		
	  End If
    End Do
	write(999,"(A25,I8)") "Total starting isotopes: ",startingIsotopeCount	
!initialise tally 	
	totalSimulationAtoms = 0.0D0
!close output file
	write(999,"(A1)") " "
    close(999)
  End Subroutine createMaterial
  

!------------------------------------------------------------------------
! SUBROUTINE addProductIsotopes                                                    
! Make tally arrays for material and mark starting isotopes/amounts                                    
!------------------------------------------------------------------------
  
  Subroutine addProductIsotopes()      
!force declaration of all variables
	Implicit None
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: keyP
!open output file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	write(999,"(A1)") " "
	write(999,"(A40,F8.4)") "Add Product Isotopes to Tally           ",ProgramTime()
!output tally to file	
    !write(999,"(A1)") " "
    !write(999,"(A70)") "----------------------------------------------------------------------"
    !write(999,"(A28)") "Material Isotopes"
	!write(999,"(A6,F8.4)") "Time: ",ProgramTime()
    !write(999,"(A70)") "----------------------------------------------------------------------"
	!Do i=1,size(isotopeTallyChar)
    !  If(isotopeTallyInt(i,5).eq.1)Then
	!	write(999,"(I8)") i
	!  End If
	!End Do
	!write(999,"(A1)") " "
!Add product isotopes to tally
    Do i=1,size(xsKey,1)		
!Get key and data for cross section	
      !print *,xsKey(i,1),xsKey(i,2),xsKey(i,3),"   ",xsKey(i,4),xsKey(i,5),xsKey(i,6)
	  keyP = makeIsotopeKey(xsKey(i,4),xsKey(i,5),xsKey(i,6))
	  isotopeTallyInt(keyP,5) = 1
    End Do	
!output tally to file	
    !write(999,"(A1)") " "
    !write(999,"(A70)") "----------------------------------------------------------------------"
    !write(999,"(A28)") "Material Isotopes + Products"
	!write(999,"(A6,F8.4)") "Time: ",ProgramTime()
    !write(999,"(A70)") "----------------------------------------------------------------------"
	!Do i=1,size(isotopeTallyChar)
    !  If(isotopeTallyInt(i,5).eq.1)Then
	!	write(999,"(I8)") i
	!print *,isotopeTallyInt(i,1),isotopeTallyInt(i,2),isotopeTallyInt(i,3)
	!  End If
	!End Do
    write(999,"(A1)") " "
!close output file
    close(999) 
  End Subroutine addProductIsotopes
  
  
  
  
!------------------------------------------------------------------------
! SUBROUTINE addDecayIsotopes                                                    
! Make tally arrays for material and mark starting isotopes/amounts                                    
!------------------------------------------------------------------------
  
  Subroutine addDecayIsotopes()      
!force declaration of all variables
	Implicit None
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: z,a,m
    Integer(kind=StandardInteger) :: key
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: isotopeArray
!open output file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	write(999,"(A1)") " "
	write(999,"(A40,F8.4)") "Add Product Isotopes to Tally           ",ProgramTime()
	Do i=1,size(isotopeTallyChar)
      If(isotopeTallyInt(i,5).eq.1.and.isotopeTallyActive(i,1).gt.(-1.0D0))Then  
		isotopeArray = IsotopesInDecayTree(isotopeTallyInt(i,1),&
		isotopeTallyInt(i,2),isotopeTallyInt(i,3))
		Do j=1,size(isotopeArray,1)
		  z = isotopeArray(j,1)
		  a = isotopeArray(j,2)
		  m = isotopeArray(j,3)
		  key = makeIsotopeKey(z,a,m)
		  isotopeTallyInt(key,5) = 1
		End Do
	  End If
    End Do		
!output tally to file	
    !write(999,"(A1)") " "
    !write(999,"(A70)") "----------------------------------------------------------------------"
    !write(999,"(A59)") "Material Isotopes + Products/Decay Parents + Decay Children"
	!write(999,"(A6,F8.4)") "Time: ",ProgramTime()
    !write(999,"(A70)") "----------------------------------------------------------------------"
	!Do i=1,size(isotopeTallyChar)
    !  If(isotopeTallyInt(i,5).eq.1)Then
	!	write(999,"(I8)") i
	!  End If
	!End Do
    write(999,"(A1)") " "
!close output file
    close(999)   
  
  End Subroutine addDecayIsotopes  
  
  
  
!------------------------------------------------------------------------
! SUBROUTINE makeSimIsotopeKeys                                                    
! Make tally arrays for material and mark starting isotopes/amounts                                    
!------------------------------------------------------------------------
  
  Subroutine makeSimIsotopeKeys()      
!force declaration of all variables
	Implicit None
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: key, simIsotopesCount
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: isotopeArray 
!Count sim isotopes
    simIsotopesCount = 0
    Do i=1,size(isotopeTallyChar,1)
      If(isotopeTallyInt(i,5).eq.1)Then 
        !print *,i,isotopeTallyInt(i,1),isotopeTallyInt(i,2),isotopeTallyInt(i,3)  
	    simIsotopesCount = simIsotopesCount + 1
	  End If
	End Do
!Allocate array
    Allocate(simIsotopeKeys(1:simIsotopesCount))
    j = 0
!store sim isotope keys
    Do i=1,size(isotopeTallyChar,1)
      If(isotopeTallyInt(i,5).eq.1)Then  
	    j = j + 1
        key = makeIsotopeKey(isotopeTallyInt(i,1),isotopeTallyInt(i,2),isotopeTallyInt(i,3)) 
        simIsotopeKeys(j) = key		
	  End If
    End Do	  
  
  
  End Subroutine makeSimIsotopeKeys   
  
  
!------------------------------------------------------------------------!
! SUBROUTINE makeReducedDecayList                                         
! Reduce the size of the decay list       
!------------------------------------------------------------------------!    
  Subroutine makeReducedDecayList()  
!force declaration of all variables
	Implicit None	
!declare variables
 	Integer(kind=StandardInteger) :: i,j,k,decayReducedCount,key
	Integer, Dimension( : , : ), Allocatable :: decayIntReduced
    Double Precision, Dimension( : , : ), Allocatable :: decayDoubleReduced
!open output file
    open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!count reduced array size
	decayReducedCount = 0
!loop through isotopes in tally
	do i=1,size(decayInt,1)
	  key=makeIsotopeKey(decayInt(i,2),decayInt(i,1),decayInt(i,3))	
      if(isotopeTallyInt(key,5).eq.1.and.isotopeTallyActive(key,1).ge.0)then		
		decayReducedCount = decayReducedCount + 1
	  endif
	enddo
!write to output file
	write(999,"(A1)") " "
	write(999,"(A30,I4,A4,I4,A10,F8.4)") "Decay data array reduced from ",size(decayInt,1),&
	  " to ",decayReducedCount,"          ",ProgramTime()	
!Allocate arrays
  Allocate(decayIntReduced(1:decayReducedCount,1:6))
  Allocate(decayDoubleReduced(1:decayReducedCount,1:2))
!loop through isotopes in tally
    j = 0
	do i=1,size(decayInt,1)
	  key=makeIsotopeKey(decayInt(i,2),decayInt(i,1),decayInt(i,3))	
      if(isotopeTallyInt(key,5).eq.1.and.isotopeTallyActive(key,1).ge.0)then		
		j = j + 1
		Do k=1,6
		  decayIntReduced(j,k) = decayInt(i,k) 
		End Do
		Do k=1,2
		  decayDoubleReduced(j,k) = decayDouble(i,k) 
		End Do
	  endif
	enddo
!write to file
    Do i=1,size(decayIntReduced,1)
	  write(999,"(I8,A2,I8,I8,I8,A2,I8,I8,I8,E20.10,E20.10)") &
	    i,"  ",decayIntReduced(i,1),decayIntReduced(i,2),decayIntReduced(i,3),&
	    "  ",decayIntReduced(i,4),decayIntReduced(i,5),decayIntReduced(i,6),&
	    decayDoubleReduced(i,1),decayDoubleReduced(i,2)
	End Do
!Transfer to decay arrays
    !Deallocate(decayInt)
    !Deallocate(decayDouble)
!Allocate arrays
    !Allocate(decayInt(1:decayReducedCount,1:6))
    !Allocate(decayDouble(1:decayReducedCount,1:2))
!Loop and transfer
	Do i=1,decayReducedCount
	  !key=makeIsotopeKey(decayInt(i,2),decayInt(i,1),decayInt(i,3))		
	  Do k=1,6
		!decayInt(i,k) = decayIntReduced(i,k) 
	  End Do
	  Do k=1,2
		!decayDouble(i,k) = decayDoubleReduced(i,k) 
	  End Do
	End Do
	
!close output file
	write(999,"(A1)") " "
    close(999)
  End Subroutine makeReducedDecayList
  
  
  
  
  
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
    If(Allocated(decayChainIsotopesArr))Then
	  Deallocate(decayChainIsotopesArr)
	End If
	Allocate(decayChainIsotopesArr(1:300,1:3))
!fill with blank data
    do i=1,300
	  do j=1,3
	    decayChainIsotopesArr(i,j) = 0
	  enddo	
	enddo
!set counter = 1 and store parent isotope
    decayChainIsotopesArr(1,1) = z
    decayChainIsotopesArr(1,2) = a
    decayChainIsotopesArr(1,3) = m
    isotopeCounterA = 1
!call function recursively
    complete = IsotopesInDecayTreeR(z,a,m)	
!transfer data to new array	
	Allocate(isotopeArray(1:isotopeCounterA,1:3))
    do i=1,isotopeCounterA
	  do j=1,3
	    isotopeArray(i,j) = decayChainIsotopesArr(i,j)
	  enddo
	enddo
!Deallocate array
	If(Allocated(decayChainIsotopesArr))Then
	  Deallocate(decayChainIsotopesArr)
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
		do j=1,size(decayChainIsotopesArr,1)
		  if(decayChainIsotopesArr(j,1).eq.zC.and.&
		  decayChainIsotopesArr(j,2).eq.aC.and.&
		  decayChainIsotopesArr(j,3).eq.mC)then
		    store = .false.
		  endif
		enddo
		if(store.eqv..true.)then
		  isotopeCounterA = isotopeCounterA + 1
		  decayChainIsotopesArr(isotopeCounterA,1) = zC
		  decayChainIsotopesArr(isotopeCounterA,2) = aC
		  decayChainIsotopesArr(isotopeCounterA,3) = mC
		endif
		lastRow = IsotopesInDecayTreeR(zC,aC,mC)
	  endif	  
	enddo  
  End Function IsotopesInDecayTreeR
  

End Module prep