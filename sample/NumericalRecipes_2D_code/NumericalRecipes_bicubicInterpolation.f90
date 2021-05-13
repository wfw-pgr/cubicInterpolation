module cubInterpMod
contains

  subroutine bcucof(y,y1,y2,y12,d1,d2,c)
    REAL d1,d2,c(4,4),y(4),y1(4),y12(4),y2(4)
    ! Given arrays y,y1,y2, and y12,
    ! each of length 4, containing the function, gradients, and cross derivative at the four grid points of a rectangular grid cell (numbered counterclockwise from the lower left), and given d1 and d2, the length of the grid cell in the 1- and 2- directions, this routine returns the table c(1:4,1:4) that is used by routine bcuint for bicubic interpolation.
    INTEGER i,j,k,l
    REAL d1d2,xx,cl(16),wt(16,16),x(16)
    SAVE wt

    DATA wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4 &
         & ,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4 &
         & ,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2&
         & ,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2 &
         & ,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2&
         & ,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2 &
         & ,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1&
         & ,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/ 
    
    d1d2=d1*d2
    d1d2=d1*d2
    do 11 i=1,4 Pack a temporary vector x.
       x(i)=y(i)
       x(i+4)=y1(i)*d1 x(i+8)=y2(i)*d2 x(i+12)=y12(i)*d1d2
    enddo
    ! Matrix multiply by the stored table.
    do 13 i=1,16
       xx=0.
       do 12 k=1,16
          xx=xx+wt(i,k)*x(k) 
       enddo
       cl(i)=xx
    enddo
    ! Unpack the result into the output table.
    l=0
    do 15 i=1,4
       do 14 j=1,4 l=l+1
          c(i,j)=cl(l) enddo 14
       enddo
    enddo
    return
  end subroutine bcucof

  
  SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy, ansy1,ansy2)
    REAL ansy,ansy1,ansy2,x1,x1l,x1u,x2,x2l,x2u,y(4),y1(4), y12(4),y2(4)
    
    ! USES bcucof
    ! Bicubic interpolation within a grid square. Input quantities are y,y1,y2,y12 (as described in bcucof); x1l and x1u, the lower and upper coordinates of the grid square in the 1- direction; x2l and x2u likewise for the 2-direction; and x1,x2, the coordinates of the desired point for the interpolation. The interpolated function value is returned as ansy, and the interpolated gradient values as ansy1 and ansy2. This routine calls bcucof.
    INTEGER i
    REAL t,u,c(4,4)
    call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c) !Get the c’s.
    if(x1u.eq.x1l.or.x2u.eq.x2l) pause 'bad input in bcuint'
    ! Equation (3.6.4).
    t=(x1-x1l)/(x1u-x1l)
    u=(x2-x2l)/(x2u-x2l)
    ansy=0.
    ansy2=0.
    ansy1=0.
    ! Equation (3.6.6).
    do 11 i=4,1,-1
       ansy  = t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
       ansy2 = t*ansy2+(3.*c(i,4)*u+2.*c(i,3))*u+c(i,2)
       ansy1 = u*ansy1+(3.*c(4,i)*t+2.*c(3,i))*t+c(2,i)
    enddo
    ansy1=ansy1/(x1u-x1l)
    ansy2=ansy2/(x2u-x2l)
    return
  END SUBROUTINE bcuint
  
