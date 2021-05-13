module genMatrixMod
contains

  ! ====================================================== !
  ! === generate Interpolation (B) Matrix              === !
  ! ====================================================== !
  subroutine generateBMatrix( BMatrix )
    implicit none
    integer                       :: i, j, k, idx, ex, ey, ez
    integer                       :: exyz(64,3)
    double precision              :: xn, yn, zn
    double precision              :: corners(8,3)
    double precision, intent(out) :: BMatrix(64,64)

    ! ------------------------------------------------------ !
    ! --- [1] index preparation                          --- !
    ! ------------------------------------------------------ !
    !  -- [1-1] normalized corner points [ 0.0, 1.0 ]    --  !
    do k=1, 2
       do j=1, 2
          do i=1, 2
             idx            = 4*(k-1) + 2*(j-1) + (i-1) + 1
             corners(idx,1) = i
             corners(idx,2) = j
             corners(idx,3) = k
          enddo
       enddo
    enddo
    !  -- [1-2] exponential number x^e                   --  !
    do k=1, 4
       do j=1, 4
          do i=1, 4
             idx         = 16*(k-1) + 4*(j-1) + (i-1) + 1
             exyz(idx,1) = i
             exyz(idx,2) = j
             exyz(idx,3) = k
          enddo
       enddo
    enddo
    !  -- [1-3] Bmatrix Calculation                      --  !
    do i=1, 64
       ex = exyz(i,1)
       ey = exyz(i,2)
       ez = exyz(i,3)
       do k=1, 8
          xn               = corners(k,1)
          yn               = corners(k,2)
          zn               = corners(k,3)
          BMatrix(0*8+k,i) = xn**ex             * yn**ey             * zn**ez
          BMatrix(1*8+k,i) = ex*xn**(abs(ex-1)) * yn**ey             * zn**ez
          BMatrix(2*8+k,i) = xn**ex             * ey*yn**(abs(ey-1)) * zn**ez
          BMatrix(3*8+k,i) = xn**ex             * yn**ey             * ez*zn**(abs(ez-1))
          BMatrix(4*8+k,i) = ex*xn**(abs(ex-1)) * ey*yn**(abs(ey-1)) * zn**ez
          BMatrix(5*8+k,i) = ex*xn**(abs(ex-1)) * yn**ey             * ez*zn**(abs(ez-1))
          BMatrix(6*8+k,i) = xn**ex             * ey*yn**(abs(ey-1)) * ez*zn**(abs(ez-1))
          BMatrix(7*8+k,i) = ex*xn**(abs(ex-1)) * ey*yn**(abs(ey-1)) * ez*zn**(abs(ez-1))
       enddo
    enddo
    return
  end subroutine generateBMatrix


  ! ====================================================== !
  ! === save__Matrices                                 === !
  ! ====================================================== !
  subroutine save__Matrices
    implicit none
    integer         , parameter  :: lun = 50
    double precision, intent(in) :: Bmatrix(64,64)
    integer                      :: i, j
    
    open(lun,file='BMatrix.dat',form='formatted')
    do j=1, 64
       write(lun,*) ( Bmatrix(i,j), i=1, 64 )
    enddo
    close(lun)
    return
  end subroutine save__Matrices

  
end module genMatrixMod
