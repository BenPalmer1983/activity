Module linearAlgebra
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: SolveLinearSet
  Public :: SolveLinearMatrix
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------

  Function SolveLinearSet(aMatrixIn, yMatrixIn) RESULT (xMatrix)
! Force declaration of all variables
    Implicit None
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: aMatrixIn
    Real(kind=DoubleReal), Dimension(:) :: yMatrixIn
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:size(yMatrixIn,1)) :: xMatrix
! Private: Declare variables
    Real(kind=DoubleReal), Dimension(1:size(aMatrixIn,1),1:size(aMatrixIn,2)) :: aMatrix
    Real(kind=DoubleReal), Dimension(1:size(yMatrixIn,1)) :: yMatrix
    Integer(kind=StandardInteger) :: row, col, matrixSize
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: sumTerm
    Real(kind=DoubleReal), Dimension(1:size(aMatrixIn,1),1:size(aMatrixIn,2)) :: lMatrix, uMatrix
    Real(kind=DoubleReal), Dimension(1:size(yMatrixIn,1)) ::  bMatrix
    Integer(kind=StandardInteger), Dimension(1:size(yMatrixIn,1)) :: pivot
! Transfer data    
    aMatrix = aMatrixIn
    yMatrix = yMatrixIn
! If matrix sizes are correct
    If(size(aMatrix,1).eq.size(aMatrix,2).and.size(aMatrix,1).eq.size(yMatrix,1))Then
! Init matrices
      bMatrix = 0.0D0
! size of square matrix
      matrixSize = size(aMatrix,1)
! Pivot Matrix
      Call PivotMatrix(aMatrix, pivot, 'M')
      Call PivotMatrix(yMatrix, pivot, 'A')    
! LU decomposition
      Call LUDecomp(aMatrix, lMatrix, uMatrix)
! Ax=(LU)x=L(Ux)=y   Ux=b   Lb=y
! Solve Lb=y first
      Do i=1,matrixSize        ! Loop through rows from top to bottom
        row = i
        sumTerm = 0.0D0
        Do j=1,row             ! Loop through columns left to right
          col = j
          If(col.lt.row)Then
            sumTerm = sumTerm + lMatrix(row,col) * bMatrix(col)
          ElseIf(col.eq.row)Then
            bMatrix(row) = (yMatrix(row)-sumTerm)/lMatrix(row,col)
            Exit  ! break
          End If
        End Do
      End Do
! Solve Ux=b
      Do i=1,matrixSize        ! Loop through rows from bottom
        row = matrixSize-i+1
        sumTerm = 0.0D0
        Do j=1,matrixSize             ! Loop through columns right to left
          col = matrixSize-j+1
          If(col.gt.row)Then
            sumTerm = sumTerm + uMatrix(row,col) * xMatrix(col)
          ElseIf(col.eq.row)Then
            xMatrix(row) = (bMatrix(row)-sumTerm)/uMatrix(row,col)
            Exit  ! break
          End If
        End Do
      End Do 
    End If
  End Function SolveLinearSet  
  
  Function SolveLinearMatrix(aMatrixIn, yMatrixIn) RESULT (xMatrix)
! Force declaration of all variables
    Implicit None
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: aMatrixIn
    Real(kind=DoubleReal), Dimension(:) :: yMatrixIn
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:size(aMatrixIn,1)) :: xMatrix
! Private: Declare variables
    Real(kind=DoubleReal), Dimension(1:size(aMatrixIn,1),1:size(aMatrixIn,2)) :: aMatrix
    Real(kind=DoubleReal), Dimension(1:size(yMatrixIn,1)) :: yMatrix
    Integer(kind=StandardInteger), Dimension(1:size(aMatrixIn,1)) :: pivot
! Transfer data    
    aMatrix = aMatrixIn
    yMatrix = yMatrixIn
! If matrix sizes are correct
    If(size(aMatrix,1).eq.size(aMatrix,2).and.size(aMatrix,1).eq.size(yMatrix,1))Then
! Pivot
      Call PivotMatrix(aMatrix, pivot, 'M')
      Call PivotMatrix(yMatrix, pivot, 'A')    
      aMatrix = InvertMatrix(aMatrix)    
      xMatrix = matmul(aMatrix,yMatrix)
    End If  
  End Function SolveLinearMatrix
End Module linearAlgebra