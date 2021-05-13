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
    double precision              :: bInvMat(4,4)
    double precision              :: avec(4), pvec(4)
    double precision              :: dx, dxInv, xMin, xi_norm, fi

    ! ------------------------------------------------------ !
    ! --- [1] Prepare Bmatrix Coefficients (analytic )   --- !
    ! ------------------------------------------------------ !
    bInvMat(1,1) = - 0.5d0
    bInvMat(1,2) = + 1.5d0
    bInvMat(1,3) = - 1.5d0
    bInvMat(1,4) = + 0.5d0

    bInvMat(2,1) = + 1.0d0
    bInvMat(2,2) = - 2.5d0
    bInvMat(2,3) = + 2.0d0
    bInvMat(2,4) = - 0.5d0

    bInvMat(3,1) = - 0.5d0
    bInvMat(3,2) = + 0.0d0
    bInvMat(3,3) = + 0.5d0
    bInvMat(3,4) = + 0.0d0

    bInvMat(4,1) = + 0.0d0
    bInvMat(4,2) = + 1.0d0
    bInvMat(4,3) = + 0.0d0
    bInvMat(4,4) = + 0.0d0

    
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
             avec(i) = avec(i) + bInvMat(i,j)*pvec(j)
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
