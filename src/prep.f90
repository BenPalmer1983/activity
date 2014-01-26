Module prep

! Setup Modules
  Use kinds
  Use constants
  Use stringfunctions		!string functions
  Use maths
  Use input


!force declaration of all variables
  Implicit None
  
!declare global variables  
  double precision, Dimension( : ), Allocatable :: fitCoefficients
  Character(len=2), Dimension( : ), Allocatable :: materialIsotopesChar
  Integer, Dimension( : , : ), Allocatable :: materialIsotopesInt
  Real, Dimension( : , : ), Allocatable :: materialIsotopesReal
  Character(len=2), Dimension( : ), Allocatable :: isotopeTallyChar
  Integer, Dimension( : , : ), Allocatable :: isotopeTallyInt
  Double Precision, Dimension( : , : ), Allocatable :: isotopeTallyActivity
  Integer(kind=ib), Dimension( : , : ), Allocatable :: isotopeTallyAtoms
  Integer(kind=ib) :: totalSimulationAtoms
  Double Precision, Dimension( : , : ), Allocatable :: isotopeTallyActive
  
  
!Privacy of functions/subroutines/variables
  Private
  Public :: runPrep					!Subroutine
  Public :: fitCoefficients			!Variable
  Public :: materialIsotopesChar	!Variable
  Public :: materialIsotopesInt		!Variable
  Public :: materialIsotopesReal	!Variable
  Public :: isotopeTallyChar		!Variable
  Public :: isotopeTallyInt			!Variable
  Public :: isotopeTallyAtoms		!Variable	
  Public :: isotopeTallyActivity	!Variable	
  Public :: totalSimulationAtoms	!Variable	
  Public :: isotopeTallyActive      !Variable
  
!Module Subroutines  
contains 

!Run all the input subroutines

  Subroutine runPrep()
	
	!Internal subroutine variables
	Integer :: i, j, k
	
	Call fitExyz()
	Call materialIsotopes()
	Call createMaterial()
	Call calculateNumberDensity()
	
	
	

  End Subroutine runPrep

  
  Subroutine fitExyz()
    
!force declaration of all variables
	Implicit None

!Internal subroutine variables
	Integer :: i, j, k, ionCount, dataPointCount, totalIons
	Integer, Dimension( : ), Allocatable :: ionDataCount
	Real, Dimension( : , : ), Allocatable :: specificIonExyz	
	double precision, Dimension( : ), Allocatable :: polyCoefficients
	integer :: order
    
	!Allocate array size (could be reduced)
	!Allocate(tempExyzArrayLarge(1:size(exyzKey),1:2))
	
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
	
	!Allocate array
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
	
	!Fill fitCoefficients with zeros
	do j=0,order
	  fitCoefficients(j) = 0
	enddo
	
	ionCount = 0
	do i=1,size(exyzKey)
	  if(i.eq.1)then
	    ionCount = ionCount + 1 !set ion
	    dataPointCount = 0 !set counter
		!deallocate and allocate array
		!Deallocate(specificIonExyz)
		Allocate(specificIonExyz(1:ionDataCount(ionCount),1:2))
	  else
	    if(exyzKey(i-1).ne.exyzKey(i))then
		  !Process data before next ion
		  !Process specificIonExyz array
		  !if(ionCount.eq.1)then
		  !Call polyFitVerb(specificIonExyz,order,polyCoefficients)
		  !else
		  !  Call polyFit(specificIonExyz,order,polyCoefficients)
		  !endif
		  Call polyFit(specificIonExyz,order,polyCoefficients)
		  do j=0,order
		    fitCoefficients(j) = fitCoefficients(j) + (polyCoefficients(j) / (1.0 * totalIons))
		  enddo
		  !print
		  !if(ionCount.eq.1)then
		    !do j=1,(size(specificIonExyz)/2)
			!  print *,specificIonExyz(j,1),specificIonExyz(j,2)
			!enddo
			!do j=0,order
			!  print *,polyCoefficients(j)
			!enddo
		  !endif
		  !Next ion
		  dataPointCount = 0 !reset counter
		  ionCount = ionCount + 1 !set ion
		  !deallocate and allocate array
		  Deallocate(specificIonExyz)
		  Allocate(specificIonExyz(1:ionDataCount(ionCount),1:2))
		endif
	  endif
	  !
	  !increment counter
	  dataPointCount = dataPointCount + 1	
	  !
	  !Store data
	  specificIonExyz(dataPointCount,1) = exyzData(i,1)
	  specificIonExyz(dataPointCount,2) = exyzData(i,2)
	  if(i.eq.size(exyzKey)) then
	    !Process data for last ion
		!Process specificIonExyz array
		Call polyFit(specificIonExyz,order,polyCoefficients)
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
	!Real, Dimension( : , : ), Allocatable :: materialIsotopesReal
	
	!Character(len=2), Dimension( : ), Allocatable :: materialIsotopesChar
    !Integer, Dimension( : , : ), Allocatable :: materialIsotopesInt
    !Real, Dimension( : , : ), Allocatable :: materialIsotopesReal

	!load all decay chain steps
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
  
  
  Subroutine createMaterial()
      
!force declaration of all variables
	Implicit None

!Internal subroutine variables
    Integer :: i, j, k
	Integer :: key
	Real :: tempReal
	Double Precision :: tempDouble
!Temporary isotope arrays
    !Character(len=2), Dimension( : ), Allocatable :: isotopeTallyCharTemp
    !Integer, Dimension( : , : ), Allocatable :: isotopeTallyIntTemp
    !Double Precision, Dimension( : ), Allocatable :: isotopeTallyDoubleTemp
	
!Create an isotope tally array
!Allocate arrays
!Z 0-120 A 0-295 M0-1 
!Key = 590Z + 290M + A
    Allocate(isotopeTallyChar(1:70800))
    Allocate(isotopeTallyInt(1:70800,1:4))
    Allocate(isotopeTallyAtoms(1:70800,1:2))	
	Allocate(isotopeTallyActivity(1:70800,1:2))
	Allocate(isotopeTallyActive(1:70800,1:3))
	
	
!Fill with blank data	
	do i=1,70800
	  isotopeTallyChar(i) = "ZZ"
	  isotopeTallyInt(i,1) = -1		!Z
	  isotopeTallyInt(i,2) = -1		!A
	  isotopeTallyInt(i,3) = -1		!M
	  isotopeTallyInt(i,4) = -1		!Decay calc marker
	  isotopeTallyActivity(i,1) = -1.0		!half life of isotope (-1 if stable)
	  isotopeTallyActivity(i,2) = 0.0		!activity of isotope in bq
	  isotopeTallyActive(i,1) = -1.0		!half life of isotope (-1 if stable)
	  isotopeTallyActive(i,2) = 0.0	    !activity of isotope in bq
	  isotopeTallyActive(i,3) = 0.0		!beam creation rate
	enddo
	
!Fill with isotope data		
    do i=1,size(isotopesChar)
	  do j=0,1
	    key = 590 * isotopesInt(i,1) + 290 * j + isotopesInt(i,2)
	    isotopeTallyChar(key) = isotopesChar(i)
	    isotopeTallyInt(key,1) = isotopesInt(i,1)
	    isotopeTallyInt(key,2) = isotopesInt(i,2)
	    isotopeTallyInt(key,3) = j		
	  enddo
	enddo
	
!Set the composition amounts
    do i=1,size(materialIsotopesChar)
	  key = 590 * materialIsotopesInt(i,1) + 290 * 0 + materialIsotopesInt(i,2)
!atoms in a milligram of the material (max 6.02E17 atoms)
	  tempDouble = 1.0 * (materialIsotopesReal(i,3) / materialIsotopesReal(i,1))&
	  * (avogadrosConstant/1E6)	
	  isotopeTallyAtoms(key,1) = tempDouble	
	  isotopeTallyAtoms(key,2) = 0
	enddo
	
!load half lives of unstable isotopes
    !  decayChar(n) = Element e.g. FE, CR  
    !  decayInt(n,1) = parent Z     
    !  decayInt(n,2) = parent A    
    !  decayInt(n,3) = parent meta      
    !  decayInt(n,4) = child Z     
    !  decayInt(n,5) = child A    
    !  decayInt(n,6) = child meta    
    !  decayDouble(n,1) = branching factor   
    !  decayDouble(n,2) = half life   	
	do i=1,size(decayChar)
	  key = 590 * decayInt(i,2) + 290 * decayInt(i,3) + decayInt(i,1)
	  if(decayDouble(i,2).gt.0)then
	    isotopeTallyActive(key,1) = decayDouble(i,2)   !store half life
	    !print *,decayInt(i,1),decayInt(i,2),decayInt(i,3),decayDouble(i,2),isotopeTallyActive(i,1)
	  endif
	enddo
	
!tally all the atoms in the 	
	totalSimulationAtoms = 0
	do i=1,size(isotopeTallyChar)
	  totalSimulationAtoms = totalSimulationAtoms + isotopeTallyAtoms(i,1)
	enddo

!print out table	
    !do i=1,size(isotopeTallyChar)
    !  if(isotopeTallyChar(i).ne."ZZ")then
	!    print *,i,isotopeTallyChar(i),isotopeTallyInt(i,1),isotopeTallyInt(i,2), &
    !		isotopeTallyInt(i,3),isotopeTallyAtoms(i,1),isotopeTallyAtoms(i,2),&
	!	isotopeTallyActive(i,1),isotopeTallyActive(i,2)
	!  endif
    !enddo
  
  End Subroutine createMaterial
  
  
  
  Subroutine calculateNumberDensity()
      
!atoms per cubic metre	  
	  
!force declaration of all variables
	Implicit None

!Internal subroutine variables
    Integer :: i, j, k
	Real :: fractionComposition, atomicMass
	
  
!Calculate number density  
	
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
	
	numberDensity = ((materialDensity * 1000) / atomicMass) * avogadrosConstant		!density from kgm-3 to gm-3
	
	
	!print *,avogadrosConstant,atomicMass,materialDensity,numberDensity
		  
  
  End Subroutine calculateNumberDensity
  
  
  
  
  
  
  

End Module prep