! ====================================================== !
! === main.f90 :: cubic interpolation (1D) test      === !
! ====================================================== !
program main
  use cubInterpMod
  implicit none
  integer                     :: i, j
  integer         , parameter :: lun  = 50
  integer         , parameter :: nRef = 21
  integer         , parameter :: nItp = 201
  double precision, parameter :: xMin = 0.0
  double precision, parameter :: xMax = 3.14159265358979d0 / 4.d0
  double precision            :: dx
  double precision            :: xg(nRef), yg(nRef), xi(nItp), yi(nItp), ya(nItp)

  ! ------------------------------------------------------ !
  ! --- [1] xGrid / yGrid making                       --- !
  ! ------------------------------------------------------ !
  dx = ( xMax-xMin ) / dble( nRef-1 )
  do i=1, nRef
     xg(i) = dx * dble(i-1) + xMin
     yg(i) = testFunc( xg(i) )
  enddo
  ! ------------------------------------------------------ !
  ! --- [2] x_at_interpolation point making            --- !
  ! ------------------------------------------------------ !
  dx = ( xMax-xMin ) / dble( nItp-1 )
  do i=1, nItp
     xi(i) = dx * dble(i-1) + xMin
     yi(i) = 0.d0
     ya(i) = testFunc( xi(i) )
  enddo
  ! ------------------------------------------------------ !
  ! --- [3] cubic interpolation                        --- !
  ! ------------------------------------------------------ !
  call cubicInterpolation_1D( xg, yg, xi, yi, nRef, nItp )

  ! ------------------------------------------------------ !
  ! --- [4] output result                              --- !
  ! ------------------------------------------------------ !
  !  -- [4-1] grid data output                         --  !
  open( lun, file=trim("dat/xRef.dat"), form="formatted" )
  do i=1, nRef
     write(lun,*) xg(i), yg(i)
  enddo
  close(lun)
  !  -- [4-2] interpolated data output                 --  !
  open( lun, file=trim("dat/xItp.dat"), form="formatted" )
  do i=1, nItp
     write(lun,*) xi(i), yi(i), ya(i)
  enddo
  close(lun)
  
  return
contains

  
  ! ====================================================== !
  ! === Function to be interpolated                    === !
  ! ====================================================== !
  function testFunc( xval )
    implicit none
    double precision, intent(in) :: xval
    double precision             :: testFunc
    
    testFunc = sin( xval )
    return
  end function testFunc
  
  
end program main
