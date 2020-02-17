module cubInterpMod
contains

  ! ====================================================== !
  ! === 1D cubic interpolation ( Lagrange )            === !
  ! ====================================================== !
  subroutine cubicInterpolation_1D( xg, yg, xi, yi, nRef, nItp )
    implicit none
    integer         , intent(in)  :: nRef, nItp
    double precision, intent(in)  :: xg(nRef), yg(nRef)
    double precision, intent(in)  :: xi(nItp)
    double precision, intent(out) :: yi(nItp)
    integer         , parameter   :: lun = 50
    integer                       :: i, j, ki
    integer         , allocatable :: xi_index(:)
    double precision              :: bMatrix(4,4), bInvMat(4,4)
    double precision              :: avec(4), pvec(4)
    double precision              :: dx, dxInv, xMin, xi_norm, fi

    ! ------------------------------------------------------ !
    ! --- [1] Prepare Bmatrix Coefficients               --- !
    ! ------------------------------------------------------ !
    open( lun, file=trim("dat/bInvMat.dat"), form="formatted" )
    do i=1, 4
       read(lun,*) bMatrix(i,1:4)
    enddo
    close(lun)

    ! ------------------------------------------------------ !
    ! --- [2] interpolation pt. ==> grid number          --- !
    ! ------------------------------------------------------ !
    allocate( xi_index(nItp) )
    xMin  = xg(1)
    dx    = xg(2) - xg(1)
    dxInv = 1.d0 / dx
    do ki=1, nItp
       xi_index(ki) = ( ceiling( ( xi(ki)-xMin ) * dxInv ) - 1 ) + 1
    enddo
    
    ! ------------------------------------------------------ !
    ! --- [3] cubic interpolation loop                   --- !
    ! ------------------------------------------------------ !
    do ki=1, nItp
       !  -- [3-1] coefficient calculation                  --  !
       avec(1:4) = 0.d0
       pvec(1:4) = yg( (xi_index(ki)-1):( xi_index(ki)+2 ) )
       do i=1, 4
          do j=1, 4
             avec(i) = avec(i) + bMatrix(i,j)*pvec(j)
          enddo
       enddo
       !  -- [3-2] polynomial interpolation                 --  !
       xi_norm = ( xi(ki) - xg( xi_index(ki) ) )*dxInv
       fi      = 0.d0
       yi(ki)  =   avec(1)*xi_norm**3 + avec(2)*xi_norm**2 &
            &    + avec(3)*xi_norm**1 + avec(4)
    enddo
    
    return
  end subroutine cubicInterpolation_1D
  
end module cubInterpMod
