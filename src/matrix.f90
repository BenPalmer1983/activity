Module matrix
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use arrayFunctions
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
! --- Functions
  Public :: InvertMatrix
  Public :: TransposeMatrix
  Public :: IdentityMatrix
  Public :: DiagMatrix
  Public :: Trace
  Public :: MatAdd
  Public :: MatMult
  Public :: ScalarMult  
  Public :: PivotMatrix
  Public :: PositiveYMatrix
  Public :: LnYMatrix
! --- Subroutines
  Public :: LUDecomp
! Interfaces  
  Interface Trace
    Module Procedure Trace_R, Trace_I
  End Interface Trace
  Interface PivotMatrix
    Module Procedure PivotMatrix_1D, PivotMatrix_2D
  End Interface PivotMatrix
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------

  Function InvertMatrix(xMatrix) RESULT (xMatrixInverse)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,1)) :: xMatrixInverse
! Private: Declare variables
    Integer(kind=StandardInteger) :: row,col,rowb
    Integer(kind=StandardInteger) :: matrixSize
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:2*size(xMatrix,1)) :: xMatrixWorking
    Real(kind=DoubleReal), Dimension(1:2*size(xMatrix,1)) :: xMatrixRow
! matrix(row,column)  
! Initialise variables
    row = 0
    rowb = 0
    col = 0
    matrixSize = size(xMatrix,1)
    xMatrixWorking = 0.0D0
    xMatrixInverse = 0.0D0
    xMatrixRow = 0.0D0
! if a square matrix
    If(size(xMatrix,1).eq.size(xMatrix,2))Then
! Fill working array
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
        End Do
      End Do
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.eq.col)Then
            xMatrixWorking(row,col+matrixSize) = 1.0D0
          End If
        End Do
      End Do
! make lower triangle of zeros
      Do row=1,matrixSize-1
        Do rowb=row+1,matrixSize
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the lower triangle
        Do rowb=row+1,matrixSize
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! re-force zeros in the lower triangle
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.gt.col)Then
            xMatrixWorking(row,col) = 0.0D0
          End If
        End Do
      End Do
! make upper triangle of zeros
      Do row=matrixSize,2,-1
        Do rowb=row-1,1,-1
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the upper triangle
        Do rowb=row-1,1,-1
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! Divide rhs by diagonal on lhs and store in inverse
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixInverse(row,col) = 1.0D0*&
          xMatrixWorking(row,col+matrixSize)/xMatrixWorking(row,row)
        End Do
      End Do
    End If 
  End Function InvertMatrix

  Function TransposeMatrix(xMatrix) RESULT (xMatrixTranspose)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,2),1:size(xMatrix,1)) :: xMatrixTranspose
    Integer(kind=StandardInteger) :: row,col
! Transpose
    Do row=1,size(xMatrix,1)
      Do col=1,size(xMatrix,2)
        xMatrixTranspose(col,row) = xMatrix(row,col)
      End Do
    End Do
  End Function TransposeMatrix

  Function IdentityMatrix(iMatrix) RESULT (oMatrix)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: row,col
    Real(kind=DoubleReal), Dimension(:,:) :: iMatrix
    Real(kind=DoubleReal), Dimension(1:size(iMatrix,1),1:size(iMatrix,2)) :: oMatrix
! Transpose
    Do row=1,size(iMatrix,1)
      Do col=1,size(iMatrix,2)
        If(row.eq.col)Then
          oMatrix(row,col) = 1.0D0
        Else
          oMatrix(row,col) = 0.0D0
        End If
      End Do
    End Do
  End Function IdentityMatrix

  Function DiagMatrix(iMatrix) RESULT (oMatrix)
! Takes diagonal of a matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: row,col
    Real(kind=DoubleReal), Dimension(:,:) :: iMatrix
    Real(kind=DoubleReal), Dimension(1:size(iMatrix,1),1:size(iMatrix,2)) :: oMatrix
! Transpose
    If(size(iMatrix,1).ne.size(iMatrix,2))Then
      oMatrix = iMatrix
    Else
      Do row=1,size(iMatrix,1)
        Do col=1,size(iMatrix,2)
          If(row.eq.col)Then
            oMatrix(row,col) = 1.0D0*iMatrix(row,col)
          Else
            oMatrix(row,col) = 0.0D0
          End If
        End Do
      End Do
    End If
  End Function DiagMatrix

  Function Trace_R(iMatrix) RESULT (output)
! Trace of matrix (sum of diagonal)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: row,col
    Real(kind=DoubleReal), Dimension(:,:) :: iMatrix
    Real(kind=DoubleReal) :: output
! Trace
    output = 0.0D0
    If(size(iMatrix,1).eq.size(iMatrix,2))Then
      Do row=1,size(iMatrix,1)
        col = row
        output = output + iMatrix(row,col)
      End Do
    End If
  End Function Trace_R

  Function Trace_I(iMatrix) RESULT (output)
! Trace of matrix (sum of diagonal)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: row,col
    Integer(kind=StandardInteger), Dimension(:,:) :: iMatrix
    Integer(kind=StandardInteger) :: output
! Trace
    output = 0
    If(size(iMatrix,1).eq.size(iMatrix,2))Then
      Do row=1,size(iMatrix,1)
        col = row
        output = output + iMatrix(row,col)
      End Do
    End If
  End Function Trace_I

  Function MatAdd(xMatrix,yMatrix) RESULT (oMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix, yMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,2)) :: oMatrix
! Initialise variables
    oMatrix = 0.0D0
! Add matrices
    If(size(xMatrix,1).eq.size(yMatrix,1).and.size(xMatrix,2).eq.size(yMatrix,2))Then
      Do i=1,size(xMatrix,1)
        Do j=1,size(xMatrix,2)
          oMatrix(i,j) = xMatrix(i,j) + yMatrix(i,j)
        End Do
      End Do
    End If
  End Function MatAdd

  Function MatMult(xMatrix,yMatrix) RESULT (oMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j,k,xRows,xCols,yRows,yCols
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix, yMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(yMatrix,2)) :: oMatrix
! matrix(row,column)
! Initialise variables
    oMatrix = 0.0D0
    xRows = size(xMatrix,1)  ! Output rows
    xCols = size(xMatrix,2)
    yRows = size(yMatrix,1)
    yCols = size(yMatrix,2)  ! Output cols
! Multiply matrices
    If(xCols.eq.yRows)Then
      Do i=1,xRows
        Do j=1,yCols
          Do k=1,xCols
            oMatrix(i,j) = oMatrix(i,j) + xMatrix(i,k)*yMatrix(k,j)
          End Do
        End Do
      End Do
    End If
  End Function MatMult

  Function ScalarMult(scalar,xMatrix) RESULT (oMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal) :: scalar
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,2)) :: oMatrix
! Mult 2D
    Do i=1,size(xMatrix,1)
      Do j=1,size(xMatrix,2)
        oMatrix(i,j) = scalar*xMatrix(i,j)
      End Do
    End Do
  End Function ScalarMult  

! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------  
  
  Subroutine LUDecomp(xMatrix, lMatrix, uMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: row, rowB, col, matrixSize
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,2)) :: lMatrix, uMatrix
    Real(kind=DoubleReal) :: factor
! matrix(row,column)
! Initialise variables
    lMatrix = 0.0D0
    uMatrix = 0.0D0
    matrixSize = size(xMatrix,1)
    factor = 0.0D0
! Decompose - if square matrix
    If(size(xMatrix,1).eq.size(xMatrix,2))Then
      Do row=1,matrixSize-1                  ! Loop through rows in matrix
        Do rowB=row+1,matrixSize             ! Loop through rows rowB below row
          If(xMatrix(rowb,row).ne.0.0D0)Then
            factor = (1.0D0*xMatrix(rowb,row))/(1.0D0*xMatrix(row,row))
            Do col=1,matrixSize
              xMatrix(rowb,col) = xMatrix(rowb,col)-factor*xMatrix(row,col)
            End Do
          Else
            factor = 0.0D0
          End If
          lMatrix(rowB,row) = factor
        End Do
      End Do
! Identity line in lower
      Do row=1,matrixSize
        lMatrix(row,row) = 1.0D0
      End Do
! transfer upper triangle to upper
      Do row=1,matrixSize
        Do col=row,matrixSize
          uMatrix(row,col) = xMatrix(row,col)
        End Do
      End Do
    End If
  End Subroutine LUDecomp
    
  Subroutine PivotMatrix_1D(xMatrix, pivotMap, operationIn)
! Force declaration of all variables
    Implicit None
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:) :: xMatrix
    Integer(kind=StandardInteger), Dimension(:) :: pivotMap
    Character(Len=1), Optional :: operationIn
! Private: Declare variables    
    Integer(kind=StandardInteger) :: row
    Logical :: sortFlag
    Character(Len=1) :: operation
    Logical, Dimension(1:size(xMatrix,1)) :: swapFlags
! Optional     
    operation = "M"  ! M make, A apply, R reverse
    If(Present(operationIn))Then
      operation = operationIn
    End If
! M: Make pivot map and pivot the input matrix    
    If(operation.eq."M")Then
! Sort matrix
      sortFlag = .true.
      Do while(sortFlag)
        sortFlag = .false.
        Do row=1,(size(xMatrix,1)-1)
          If(xMatrix(row).lt.xMatrix(row+1))Then
            sortFlag = .true.
            Call swapRows(xMatrix,row,row+1)
            Call swapRows(pivotMap,row,row+1)
          End If
        End Do
      End Do
    End If  
! A: Apply pivot map to the input matrix (Reverse is just the same)
    If(operation.eq."A".or.operation.eq."R")Then
      swapFlags = .true.
      Do row=1,size(xMatrix,1)
        If(swapFlags(row))Then
          swapFlags(row) = .false.
          swapFlags(pivotMap(row)) = .false.
          Call swapRows(xMatrix,row,pivotMap(row))
        End If
      End Do      
    End If    
  End Subroutine PivotMatrix_1D
    
  Subroutine PivotMatrix_2D(xMatrix, pivotMap, operationIn)
! Force declaration of all variables
    Implicit None
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Integer(kind=StandardInteger), Dimension(:) :: pivotMap
    Character(Len=1), Optional :: operationIn
! Private: Declare variables    
    Integer(kind=StandardInteger), Dimension(1:size(xMatrix,1)) :: sortMatrix
    Integer(kind=StandardInteger) :: row, col
    Logical :: sortFlag
    Character(Len=1) :: operation
    Logical, Dimension(1:size(xMatrix,1)) :: swapFlags
! Optional     
    operation = "M"  ! M make, A apply, R reverse
    If(Present(operationIn))Then
      operation = operationIn
    End If
! M: Make pivot map and pivot the input matrix    
    If(operation.eq."M")Then
! Init  
      sortMatrix = 0
! Loop
      Do row=1,size(xMatrix,1)   
        Do col=1,size(xMatrix,2)       
          If(xMatrix(row,col).ne.0.0D0)Then
            sortMatrix(row) = sortMatrix(row) + 2**(size(xMatrix,2)-col+1)
          End If  
        End Do
        pivotMap(row) = row
      End Do
! Sort matrix
      sortFlag = .true.
      Do while(sortFlag)
        sortFlag = .false.
        Do row=1,(size(xMatrix,1)-1)
          If(sortMatrix(row).lt.sortMatrix(row+1))Then
            sortFlag = .true.
            Call swapRows(xMatrix,row,row+1)
            Call swapRows(sortMatrix,row,row+1)
            Call swapRows(pivotMap,row,row+1)
          End If
        End Do
      End Do
    End If  
! A: Apply pivot map to the input matrix (Reverse is just the same)
    If(operation.eq."A".or.operation.eq."R")Then
      swapFlags = .true.
      Do row=1,size(xMatrix,1)
        If(swapFlags(row))Then
          swapFlags(row) = .false.
          swapFlags(pivotMap(row)) = .false.
          Call swapRows(xMatrix,row,pivotMap(row))
        End If
      End Do      
    End If    
  End Subroutine PivotMatrix_2D
  
  
  Subroutine PositiveYMatrix(X,Y)
! Make Y side positive
! Used when taking ln of Y matrix to fit exp(P(x)) spline
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: X
    Real(kind=DoubleReal), Dimension(:) :: Y
! Private: Declare variables
    Integer(kind=StandardInteger) :: row, col
! Check X and Y size
    If(Size(X,1).eq.Size(Y,1))Then
      Do row=1,Size(Y,1)
        If(Y(row).lt.0.0D0)Then
          Y(row) = -1.0D0*Y(row)
          Do col=1,size(X,1)
            X(row,col) = -1.0D0*X(row,col)
          End Do
        End If
      End Do
    End If  
  End Subroutine PositiveYMatrix  
  
  Subroutine LnYMatrix(Y)
! Take LN of Y matrix
! Used to fit exp(P(x)) spline
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:) :: Y
! Private: Declare variables
    Integer(kind=StandardInteger) :: row
! Check X and Y size
    Do row=1,Size(Y,1)
      If(Y(row).eq.0.0D0)Then
        Y(row) = 1.0D-99 ! stop 0.0D0 breaking it
      End If
      Y(row) = log(Y(row))
    End Do 
  End Subroutine LnYMatrix  
    
End Module matrix