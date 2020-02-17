module genMatrixMod
contains

  ! ====================================================== !
  ! === generate Interpolation (B) Matrix              === !
  ! ====================================================== !
  subroutine generate__BMatrix( BMatrix )
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
             corners(idx,1) = i-1
             corners(idx,2) = j-1
             corners(idx,3) = k-1
          enddo
       enddo
    enddo
    !  -- [1-2] exponential number x^e                   --  !
    do k=1, 4
       do j=1, 4
          do i=1, 4
             idx         = 16*(k-1) + 4*(j-1) + (i-1) + 1
             exyz(idx,1) = i-1
             exyz(idx,2) = j-1
             exyz(idx,3) = k-1
          enddo
       enddo
    enddo
    ! ------------------------------------------------------ !
    ! --- [2] BMatrix Calculation                        --- !
    ! ------------------------------------------------------ !
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
  end subroutine generate__BMatrix


  ! ====================================================== !
  ! ===  generate Differential Matrix ( D )            === !
  ! ====================================================== !
  subroutine generate__DMatrix( DMatrix )
    implicit none
    integer                       :: i, j, m
    integer                       :: cuboid(8)
    double precision, intent(out) :: DMatrix(64,64)

    ! ------------------------------------------------------ !
    ! --- [1] preparation cboid / DMatrix                --- !
    ! ------------------------------------------------------ !
    data cuboid/22,23,26,27,38,39,42,43/
    do j=1, 64
       do i=1, 64
          DMatrix(i,j) = 0.d0
       enddo
    enddo
    
    ! ------------------------------------------------------ !
    ! --- [1] preparation 'cuboid' Coordinates           --- !
    ! ------------------------------------------------------ !
    do m=1, 8
       DMatrix( m+0 , cuboid(m)    ) = 1.d0
    enddo
    do m=1, 8
       DMatrix( m+8 , cuboid(m)-1  ) = -0.5d0
       DMatrix( m+8 , cuboid(m)+1  ) = +0.5d0
    enddo
    do m=1, 8
       DMatrix( m+16, cuboid(m)-4  ) = -0.5d0
       DMatrix( m+16, cuboid(m)+4  ) = +0.5d0
    enddo
    do m=1, 8
       DMatrix( m+24, cuboid(m)-16 ) = -0.5d0
       DMatrix( m+24, cuboid(m)+16 ) = +0.5d0
    enddo
    do m=1, 8
       DMatrix( m+32, cuboid(m)+5  ) = +0.25d0
       DMatrix( m+32, cuboid(m)-3  ) = -0.25d0
       DMatrix( m+32, cuboid(m)+3  ) = -0.25d0
       DMatrix( m+32, cuboid(m)-5  ) = +0.25d0
    enddo
    do m=1, 8
       DMatrix( m+40, cuboid(m)+17 ) = +0.25d0
       DMatrix( m+40, cuboid(m)-15 ) = -0.25d0
       DMatrix( m+40, cuboid(m)+15 ) = -0.25d0
       DMatrix( m+40, cuboid(m)-17 ) = +0.25d0
    enddo
    do m=1, 8
       DMatrix( m+48, cuboid(m)+20 ) = +0.25d0
       DMatrix( m+48, cuboid(m)-12 ) = -0.25d0
       DMatrix( m+48, cuboid(m)+12 ) = -0.25d0
       DMatrix( m+48, cuboid(m)-20 ) = +0.25d0
    enddo
    do m=1, 8
       DMatrix( m+56, cuboid(m)+21 ) = +0.125d0
       DMatrix( m+56, cuboid(m)+13 ) = -0.125d0
       DMatrix( m+56, cuboid(m)+19 ) = -0.125d0
       DMatrix( m+56, cuboid(m)+11 ) = +0.125d0
       DMatrix( m+56, cuboid(m)-11 ) = -0.125d0
       DMatrix( m+56, cuboid(m)-19 ) = +0.125d0
       DMatrix( m+56, cuboid(m)-13 ) = +0.125d0
       DMatrix( m+56, cuboid(m)-21 ) = -0.125d0
    enddo
    return
  end subroutine generate__DMatrix

  
  ! ====================================================== !
  ! === save__Matrices                                 === !
  ! ====================================================== !
  subroutine save__Matrices( BMatrix, DMatrix )
    implicit none
    integer         , parameter  :: lun = 50
    double precision, intent(in) :: BMatrix(64,64), DMatrix(64,64)
    integer                      :: i, j
    
    ! ------------------------------------------------------ !
    ! --- [1] save BMatrix                               --- !
    ! ------------------------------------------------------ !
    open(lun,file='BMatrix.dat',form='formatted')
    do i=1, 64
       write(lun,*) ( BMatrix(i,j), j=1, 64 )
    enddo
    close(lun)
    ! ------------------------------------------------------ !
    ! --- [2] save DMatrix                               --- !
    ! ------------------------------------------------------ !
    open(lun,file='DMatrix.dat',form='formatted')
    do i=1, 64
       write(lun,*) ( DMatrix(i,j), j=1, 64 )
    enddo
    close(lun)

    return
  end subroutine save__Matrices

  
end module genMatrixMod
