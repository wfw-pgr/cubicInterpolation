program main
  use genMatrixMod
  implicit none
  double precision :: BMatrix(64,64), DMatrix(64,64)
  
  call generate__BMatrix( BMatrix )
  call generate__DMatrix( DMatrix )
  call save__Matrices( BMatrix, DMatrix )
  
  return
end program main
