! ====================================================== !
! === main.f90 :: cubic interpolation (3D) test      === !
! ====================================================== !
program main
  use cubInterpMod
  implicit none
  integer                       :: i, j, k, m
  integer         , parameter   :: lun  = 50
  integer         , parameter   :: LI   = 21
  integer         , parameter   :: LJ   = 21
  integer         , parameter   :: LK   = 21
  integer         , parameter   :: nItp = 201
  double precision, parameter   :: xMin = - 3.14159265358979d0 / 4.d0
  double precision, parameter   :: xMax = + 3.14159265358979d0 / 4.d0 * 7.d0
  double precision, parameter   :: yMin = - 3.14159265358979d0 / 4.d0
  double precision, parameter   :: yMax = + 3.14159265358979d0 / 4.d0 * 7.d0
  double precision, parameter   :: zMin = - 3.14159265358979d0 / 4.d0
  double precision, parameter   :: zMax = + 3.14159265358979d0 / 4.d0 * 7.d0
  double precision              :: dx, dy, dz, val
  double precision, allocatable :: xRef(:,:,:,:), xItp(:,:)
  integer         , parameter   :: x_=1, y_=2, z_=3, v_=4

  ! ------------------------------------------------------ !
  ! --- [1] xGrid / yGrid making                       --- !
  ! ------------------------------------------------------ !
  allocate( xRef(4,LI,LJ,LK), xItp(4,nItp) )
  dx = ( xMax-xMin ) / dble( LI-1 )
  dy = ( yMax-yMin ) / dble( LJ-1 )
  dz = ( zMax-zMin ) / dble( LK-1 )
  do k=1, LK
     do j=1, LJ
        do i=1, LI
           xRef(x_,i,j,k) = dx*dble(i-1) + xMin
           xRef(y_,i,j,k) = dy*dble(j-1) + yMin
           xRef(z_,i,j,k) = dz*dble(k-1) + zMin
           xRef(v_,i,j,k) = exp( 0.1d0 * ( xRef(x_,i,j,k)**2 + xRef(y_,i,j,k)**2 + xRef(z_,i,j,k)**2 ) )
        enddo
     enddo
  enddo
  ! ------------------------------------------------------ !
  ! --- [2] x_at_interpolation point making            --- !
  ! ------------------------------------------------------ !
  dx = ( xMax-xMin ) / dble( nItp-1 )
  dy = ( yMax-yMin ) / dble( nItp-1 )
  dz = ( zMax-zMin ) / dble( nItp-1 )
  do m=1, nItp
     xItp(x_,m) = dx*dble(m-1) + xMin
     xItp(y_,m) = dy*dble(m-1) + yMin
     xItp(z_,m) = dz*dble(m-1) + zMin
     xItp(v_,m) = 0.d0
  enddo
  ! ------------------------------------------------------ !
  ! --- [3] cubic interpolation                        --- !
  ! ------------------------------------------------------ !
  ! call cubicInterpolation_3D( xRef, xItp, LI, LJ, LK, nItp )

  ! ------------------------------------------------------ !
  ! --- [4] output result                              --- !
  ! ------------------------------------------------------ !
  !  -- [4-1] grid data output                         --  !
  open( lun, file=trim("xRef.dat"), form="formatted" )
  do k=1, LK
     do j=1, LJ
        do i=1, LI
           write(lun,*) xRef(x_,i,j,k), xRef(y_,i,j,k), xRef(z_,i,j,k), xRef(v_,i,j,k)
        enddo
     enddo
  enddo
  close(lun)
  !  -- [4-2] interpolated data output                 --  !
  open( lun, file=trim("xItp.dat"), form="formatted" )
  do m=1, nItp
     val = exp( 0.1d0 * ( xItp(x_,m)**2 + xItp(y_,m)**2 + xItp(z_,m)**2 ) )
     write(lun,*) xItp(x_,m), xItp(y_,m), xItp(z_,m), xItp(v_,m), val
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
