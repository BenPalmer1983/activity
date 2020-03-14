Module globals

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Declare all global variables

! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------

! Setup Modules
  Use kinds

! Force declaration of all variables
  Implicit None

! Declare global variables

! ------------------------
! MPI variables
! ------------------------
  Integer(kind=StandardInteger) :: mpiProcessCount, mpiProcessID
! ------------------------
! Initialise variables
! ------------------------
  Character(len=64) :: compileLine
  Real(kind=DoubleReal) :: programStartTime
  Character(len=255) :: currentWorkingDirectory
  Character(len=255) :: outputFile
  Character(len=255) :: activityHistoryFile
  Character(len=255) :: inOutFile
  Character(len=255) :: isotopeActivityFile
  Character(len=255) :: outputDirectory
  Character(len=255) :: tempDirectory
  Character(len=255) :: fileOutputData, fileActivityHistory
  Character(len=255) :: fileIsotopeActivity, fileIsotopeActivityG
  Character(len=255) :: fileInOut

! ------------------------
! Input variables
! ------------------------
  Character(len=2), Dimension( : ), Allocatable :: isotopesChar
  Integer, Dimension( : , : ), Allocatable :: isotopesInt
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopesReal
  Character(len=2), Dimension( : ), Allocatable :: decayChar
  Integer, Dimension( : , : ), Allocatable :: decayInt
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDouble
  Character(len=2), Dimension( : ), Allocatable :: elements
  Character(len=2), Dimension( : ), Allocatable :: elementSymbol
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: materialComposition
  Integer, Dimension( : ), Allocatable :: exyzKey
  Integer, Dimension( : , : ), Allocatable :: exyzIonKey
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: exyzData
  Integer :: polyFitOrder
  Integer, Dimension( : , : ), Allocatable :: xsKey
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xsData
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xsMCData
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: gammaLines
  Integer, Dimension( : , : ), Allocatable :: gammaLinesKey
  Integer :: integrationGranularity
  Real(kind=DoubleReal) :: beamEnergy
  Real(kind=DoubleReal) :: beamFlux
  Real(kind=DoubleReal) :: beamDuration
  Real(kind=DoubleReal) :: beamArea
  Real(kind=DoubleReal) :: amTime
  Real(kind=DoubleReal) :: timeStep
  Real(kind=DoubleReal) :: targetThickness
  Real(kind=DoubleReal) :: numberDensity
  Real(kind=DoubleReal) :: materialDensity
  Real(kind=DoubleReal) :: vpi
  Real(kind=DoubleReal) :: activityTimeStep
  Integer :: projectileZ
  Integer :: projectileA
  Character(len=1) :: individualIsotopeActivity
  Character(len=255) :: isotopeFile
  Character(len=255) :: activityFile
  Character(len=255) :: decayModesFile
  Character(len=255) :: gammaEnergiesFile
  Character(len=255) :: xsFiles
  Character(len=255) :: trajFile
  Logical :: verboseTerminal
  Real(kind=DoubleReal) :: targetDPA

! ------------------------
! Prep variables
! ------------------------
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: fitCoefficients
  Real(kind=DoubleReal), Dimension(1:10) :: fitCoefficientsMpi
  Real(kind=DoubleReal) :: trajDepth
  Character(len=2), Dimension( : ), Allocatable :: materialIsotopesChar
  Integer, Dimension( : , : ), Allocatable :: materialIsotopesInt
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: materialIsotopesReal
  Character(len=2), Dimension( : ), Allocatable :: isotopeTallyChar
  Integer, Dimension( : , : ), Allocatable :: isotopeTallyInt
  Real(kind=DoubleReal) :: totalSimulationAtoms
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeTallyActive
  Integer(kind=StandardInteger), Dimension( : , :), Allocatable :: decayChainIsotopesArr
  Integer(kind=StandardInteger) :: isotopeCounterA
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: simIsotopeKeys

! ------------------------
! Production Loss
! ------------------------
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
  Real(kind=DoubleReal) :: maxBeamDepth
  Real(kind=DoubleReal) :: simulationScalingFactor
  Real(kind=DoubleReal) :: affectedMaterialDepth
  Real(kind=DoubleReal) :: affectedMaterialVolume
  Real(kind=DoubleReal) :: dpa
  Real(kind=DoubleReal) :: totalActivityAtTime
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: gammaLineSpectra
  Real(kind=DoubleReal) :: totalActivity, totalActivityPowerOutput, affectedAtoms
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: decayTempTally
  Real(kind=DoubleReal), Dimension(1:10000) :: gammaLineBins
  Integer(kind=StandardInteger) :: gammaResolution  ! Number of bins

! Privacy of variables
  Public :: initGlobals
! ------------------------
! MPI variables
! ------------------------
  Public :: mpiProcessCount, mpiProcessID
! ------------------------
! Initialise variables
! ------------------------
  Public :: compileLine
  Public :: programStartTime
  Public :: outputFile
  Public :: activityHistoryFile
  Public :: inOutFile
  Public :: isotopeActivityFile
  Public :: currentWorkingDirectory
  Public :: outputDirectory
  Public :: tempDirectory
  Public :: fileOutputData, fileActivityHistory
  Public :: fileIsotopeActivity, fileIsotopeActivityG
  Public :: fileInOut

! ------------------------
! Input variables
! ------------------------
  Public :: isotopesChar      !Variable
  Public :: isotopesInt        !Variable
  Public :: isotopesReal      !Variable
  Public :: decayChar        !Variable
  Public :: decayInt        !Variable
  Public :: decayDouble        !Variable
  Public :: elements        !Variable
  Public :: elementSymbol      !Variable
  Public :: materialComposition      !Variable
  Public :: exyzKey            !Variable
  Public :: exyzIonKey
  Public :: exyzData        !Variable
  Public :: polyFitOrder      !Variable
  Public :: xsKey            !Variable
  Public :: xsData          !Variable
  Public :: xsMCData        !Variable
  Public :: gammaLines        !Variable
  Public :: gammaLinesKey        !Variable
  Public :: integrationGranularity  !Variable
  Public :: beamEnergy        !Variable
  Public :: beamFlux          !Variable
  Public :: beamDuration        !Variable
  Public :: beamArea           !Variable
  Public :: amTime          !Variable
  Public :: timeStep          !Variable
  Public :: targetThickness      !Variable
  Public :: numberDensity          !Variable
  Public :: materialDensity      !Variable
  Public :: vpi                  !Variable
  Public :: projectileZ        !Variable
  Public :: projectileA        !Variable
  Public :: activityTimeStep        !Variable
  Public :: individualIsotopeActivity !Variable
  Public :: verboseTerminal         !Variable
  Public :: targetDPA               !Variable

! ------------------------
! Prep variables
! ------------------------
  Public :: fitCoefficients, fitCoefficientsMpi
  Public :: trajDepth
  Public :: materialIsotopesChar
  Public :: materialIsotopesInt
  Public :: materialIsotopesReal
  Public :: isotopeTallyChar
  Public :: isotopeTallyInt
  Public :: isotopeTallyActive
  Public :: totalSimulationAtoms
  Public :: decayChainIsotopesArr
  Public :: isotopeCounterA
  Public :: simIsotopeKeys

! ------------------------
! Production Loss
! ------------------------
  Public :: reactionRate            !Variable
  Public :: maxBeamDepth            !Variable
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
  Public :: gammaLineBins
  Public :: gammaResolution

  Contains

! Init global variables
  Subroutine initGlobals()
    Implicit None
! Initialise Subroutine Variable
    compileLine = "15:08:41  11/01/2016"
! ------------------------
! Prep variables
! ------------------------
! isotopeTallyChar = "ZZ"
! isotopeTallyInt = 0
! isotopeTallyActive = 0.0D0
! ------------------------
! Production Loss
! ------------------------
    gammaResolution = 0
  End Subroutine initGlobals

End Module globals
