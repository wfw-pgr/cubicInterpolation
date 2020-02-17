! ====================================================== !
! === program borrowed from Shiki's Note             === !
! ====================================================== !
program main
  implicit none
  integer::i,N,M
  double precision::x,f,df,df2,s
  double precision,allocatable::xdata(:),fdata(:)
  double precision::pi=dacos(-1d0)
  double precision,external::cubic

  N=30
  M=100
  allocate(xdata(0:N),fdata(0:N))
  xdata=0d0
  fdata=0d0

  do i=0,N
     xdata(i)=dble(i)*0.1d0*pi
     fdata(i)=sin(xdata(i))
     write(10,*)xdata(i),fdata(i)
  enddo

  ! Cubic-spline interpolation given position as point
  do i=0,M
     x=dble(i)*0.03d0*pi-1d0
     write(11,'(2e20.7e2)')x,cubic(x,N,xdata,fdata)
  enddo

  stop
end program main

double precision function cubic(x,N,x0,f0)
  implicit none
  integer,intent(in)::N
  double precision,intent(in)::x,f0(0:N),x0(0:N)

  integer::i,i0,i1,i2,i3
  double precision::tx  
  double precision::a,b,c,d,p,q,r,s,t,u

  tx = x-x0(0)
  i1 = 0
  do i=1,N-2
     tx = x-x0(i)
     if(tx.gt.0d0)then
        i1 = i
     else
        exit
     endif
  enddo
  if(i1.eq.0)i1=1
  i0=i1-1
  i2=i1+1
  i3=i1+2

  a = x-x0(i0)
  b = x-x0(i1)
  c = x-x0(i2)
  d = x-x0(i3)  

  p = x0(i1)-x0(i0)
  q = x0(i2)-x0(i1)
  r = x0(i3)-x0(i2)

  s = x0(i2)-x0(i0)
  t = x0(i3)-x0(i1)

  u = x0(i3)-x0(i0)

  cubic = a*c*( d*f0(i1)/(p*q) + b*f0(i3)/(u*r) )/t &
       - b*d*( c*f0(i0)/(p*u) + a*f0(i2)/(q*r) )/s

  return
end function cubic


