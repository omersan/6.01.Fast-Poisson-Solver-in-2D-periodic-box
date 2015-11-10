!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Fast Poisson Solver
!     Fast Fourier Transform (FFT) Solver for solving 2D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Periodic b.c.
!
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012)
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Nov. 11, 2015
!-----------------------------------------------------------------------------!
program poisson_solver
implicit none
integer::nx,ny
real*8 ::pi,dx,dy,x,y,emax
real*8,allocatable::u(:,:),f(:,:),v(:,:),e(:,:)
integer::i,j

pi = 4.0d0*datan(1.0d0)

nx = 64  !should be power of 2
ny = nx

dx = 1.0d0/dfloat(nx)
dy = 1.0d0/dfloat(ny)

allocate(u(0:nx,0:ny))
allocate(f(0:nx,0:ny))

!exact solution
do j=0,ny
do i=0,nx
x = dfloat(i)*dx
y = dfloat(j)*dy
u(i,j) = dcos(2.0d0*pi*x) + dcos(2.0d0*pi*y)
f(i,j) =-4.0d0*pi*pi*u(i,j)
end do
end do

!compute numerical solution
allocate(v(0:nx,0:ny))
call poissonS(nx,ny,dx,dy,f,v) !spectral solution
!call poisson2(nx,ny,dx,dy,f,v) !second-order solution


!compute error
allocate(e(0:nx,0:ny))
do j=0,ny
do i=0,nx
e(i,j) = dabs(u(i,j) - v(i,j))
end do
end do
emax = maxval(e)
print*,nx,emax


open(100,file='field.plt')
write(100,*)'variables ="x","y","u","v"'
write(100,*)'zone f=point i=',nx+1,',j=',ny+1
do j=0,ny
do i=0,nx
write(100,*) dfloat(i)*dx,dfloat(j)*dy,u(i,j),v(i,j)
end do
end do
close(100)

end


!---------------------------------------------------------------------------!
!Spectral accurate Poisson solver
!Periodic, equidistant grid
!---------------------------------------------------------------------------!
subroutine poissonS(nx,ny,dx,dy,f,u)
implicit none
integer::nx,ny
real*8 ::dx,dy
real*8 ::f(0:nx,0:ny),u(0:nx,0:ny)
real*8 ::pi,Lx,Ly,den
real*8 ::kx(0:nx-1),ky(0:ny-1) 
real*8 ::data1d(2*nx*ny) 
integer::i,j,k,isign,ndim,nn(2)

!2d data
ndim =2
nn(1)=nx
nn(2)=ny

!1.Find the f coefficient in Fourier space
!assign 1d data array
k=1
do j=0,ny-1  
do i=0,nx-1   
	data1d(k)   =  f(i,j)
	data1d(k+1) =  0.0d0    
k = k + 2
end do
end do
!normalize
do k=1,2*nx*ny
data1d(k)=data1d(k)/dfloat(nx*ny)
end do
!inverse fourier transform
isign= -1
call fourn(data1d,nn,ndim,isign)

!2.Solve for u coeeficient in Fourier space
!coefficients
Lx = dfloat(nx)*dx
Ly = dfloat(ny)*dy

!wave numbers (scaled)
pi = 4.0d0*datan(1.0d0)
do i=0,nx/2-1
kx(i)      = (2.0d0*pi/Lx)*dfloat(i)
kx(i+nx/2) = (2.0d0*pi/Lx)*dfloat(i-nx/2)
end do
do j=0,ny/2-1
ky(j)      = (2.0d0*pi/Ly)*dfloat(j)
ky(j+ny/2) = (2.0d0*pi/Ly)*dfloat(j-ny/2)
end do
kx(0) = 1.0d-6 !to eleminate zero division
ky(0) = 1.0d-6 !to eleminate zero division
data1d(1) = 0.0d0
data1d(2) = 0.0d0

!Fourier coefficients for u
k=1
do j=0,ny-1
do i=0,nx-1   
    den = -(kx(i)*kx(i))-(ky(j)*ky(j))
	data1d(k)   =  data1d(k)/den
	data1d(k+1) =  data1d(k+1)/den
k = k + 2
end do
end do

!3. Find u values on physical space
!forward fourier transform
isign= 1
call fourn(data1d,nn,ndim,isign)

!assign 2d array
k=1
do j=0,ny-1
do i=0,nx-1
u(i,j)=data1d(k)
k=k+2
end do
end do

! periodicity
do j=0,ny-1
u(nx,j)=u(0,j)
end do
do i=0,nx
u(i,ny)=u(i,0)
end do

return
end

!---------------------------------------------------------------------------!
!2nd-order accurate Poisson solver
!Periodic, equidistant grid
!---------------------------------------------------------------------------!
subroutine poisson2(nx,ny,dx,dy,f,u)
implicit none
integer::nx,ny
real*8 ::dx,dy
real*8 ::f(0:nx,0:ny),u(0:nx,0:ny)
real*8 ::pi,a,b,c,den
real*8 ::kx(0:nx-1),ky(0:ny-1) 
real*8 ::data1d(2*nx*ny) 
integer::i,j,k,isign,ndim,nn(2)

!2d data
ndim =2
nn(1)=nx
nn(2)=ny

!1.Find the f coefficient in Fourier space
!assign 1d data array
k=1
do j=0,ny-1  
do i=0,nx-1   
	data1d(k)   =  f(i,j)
	data1d(k+1) =  0.0d0    
k = k + 2
end do
end do
!normalize
do k=1,2*nx*ny
data1d(k)=data1d(k)/dfloat(nx*ny)
end do
!inverse fourier transform
isign= -1
call fourn(data1d,nn,ndim,isign)

!2.Solve for u coeeficient in Fourier space
!coefficients for 2nd-order discretization
a =-2.0d0/(dx*dx) - 2.0d0/(dy*dy)  
b = 1.0d0/(dx*dx)
c = 1.0d0/(dy*dy)  

!wave numbers (scaled)
pi = 4.0d0*datan(1.0d0)
do i=0,nx/2-1
kx(i)      = (2.0d0*pi/dfloat(nx))*dfloat(i)
kx(i+nx/2) = (2.0d0*pi/dfloat(nx))*dfloat(i-nx/2)
end do
do j=0,ny/2-1
ky(j)      = (2.0d0*pi/dfloat(ny))*dfloat(j)
ky(j+ny/2) = (2.0d0*pi/dfloat(ny))*dfloat(j-ny/2)
end do
kx(0) = 1.0d-6 !to eleminate zero division
ky(0) = 1.0d-6 !to eleminate zero division
data1d(1) = 0.0d0
data1d(2) = 0.0d0

!Fourier coefficients for u
k=1
do j=0,ny-1
do i=0,nx-1   
    den = a + 2.0d0*b*dcos(kx(i)) + 2.0d0*c*dcos(ky(j))
	data1d(k)   =  data1d(k)/den
	data1d(k+1) =  data1d(k+1)/den
k = k + 2
end do
end do

!3. Find u values on physical space
!forward fourier transform
isign= 1
call fourn(data1d,nn,ndim,isign)

!assign 2d array
k=1
do j=0,ny-1
do i=0,nx-1
u(i,j)=data1d(k)
k=k+2
end do
end do

! periodicity
do j=0,ny-1
u(nx,j)=u(0,j)
end do
do i=0,nx
u(i,ny)=u(i,0)
end do

return
end




!---------------------------------------------------------------------------!
!N-dimensional FFT routine 
!Power of 2
!by Numerical Recipes
!For two dimensional problems ndim = 2
!nn(1)=nx
!nn(2)=ny
!---------------------------------------------------------------------------!
subroutine fourn(data,nn,ndim,isign)
implicit none
integer::ndim,isign
integer::nn(ndim)
real*8:: wr,wi,wpr,wpi,wtemp,theta
integer::i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3
integer::k1,k2,n,nprev,nrem,ntot
real*8:: tempr,tempi
real*8:: data(*)

ntot=1
do idim=1,ndim
ntot=ntot*nn(idim)
end do
nprev=1
do idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do i1=i2,i2+ip1-2,2
              do i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
			  end do
			end do
          endif
          ibit=ip2/2
		  1 continue
          if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          go to 1
          endif
          i2rev=i2rev+ibit
		end do
        ifp1=ip1
		2 continue
        if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do i3=1,ifp1,ip1
            do i1=i3,i3+ip1-2,2
              do i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
			  end do
			end do
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
		  end do
          ifp1=ifp2
        go to 2
        endif
nprev=n*nprev
end do
return
end




