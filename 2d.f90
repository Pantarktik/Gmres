	program poisson2d
	implicit none
	integer, parameter:: nr=10,nteta=10
	integer,parameter:: nrow=nr*nteta
	integer,parameter:: nnzero=5*nrow-2*nr-2*nteta
	real*8,parameter :: g=6.67d-7,pi=3.1415926d0	
	
	real*8 x(nrow),r(nr),dr,rhs(nrow),a(nrow,nrow),r_in,r_out,rpop(nr+1)
	real*8 rho(nrow),out_bound,var(nnzero),v(nrow),mus,y(nrow),x0(nrow)
	real*8 teta(nteta),dteta,tetapop(nteta+1),teta_bound(nr),teta_niz,teta_verh
	integer ncols(nnzero),nrows(nnzero),point(nrow+1),nvarinrow
	
	integer i,j,k,ivar,irow
	
	!!!!!!!!!!!!SETKA
	
	r_in=0.
	r_out=100.
	
	dr=(r_out-r_in)/nr
	r(1)=r_in+dr/2.
	do i=2,nr
	r(i)=r(i-1)+dr
	enddo
	
	do i=1,nr
	rpop(i)=r(i)-dr/2.
	enddo
	rpop(nr+1)=r_out
		

	
	!!!!!!!!!!!TETA
	
	teta_niz=-pi/2.
	teta_verh=pi/2.
	
	tetapop(1)=teta_niz
	tetapop(nteta+1)=teta_verh
	
	dteta=(teta_verh-teta_niz)/(nteta)

	do j=2,nteta
	tetapop(j)=tetapop(j-1)+dteta
	enddo


	do j=1,nteta
	teta(j)=tetapop(j)+dteta/2.
	enddo
	


	!open(1,file='teta2d')
	!do j=1,nteta!!
	!write(1,*) teta(j),tetapop(j)
	!enddo
	!close(1)






	!!!!!!!!RHO
	
	
	rho=1.
	
	!!!!!!!!Boundary
	
	open(1,file='1dx.dat')
	do i=1,nr
	read(1,*)mus,teta_bound(i)
	enddo
	
	out_bound=-4./3.*pi*r_out**3*rho(1)*g/r_out	
	
	!!!!MATRIX
	
	a=0.
	ivar=1
	var=0.
	point(1)=1
	
!	open(1,file='ij')
	do irow=1,nrow
	
		i=irow/nteta+1
		j=mod(irow,nteta)
		if(mod(irow,nteta).eq.0)then
			i=i-1
			j=nteta
		endif
	
	nvarinrow=5
	if(i.eq.1.or.i.eq.nr)nvarinrow=nvarinrow-1
	if(j.eq.1.or.j.eq.nteta)nvarinrow=nvarinrow-1	
	point(irow+1)=point(irow)+nvarinrow



	if(i.gt.1)then
	var(ivar)=rpop(i)**2/dr
	ncols(ivar)=irow-nteta
	nrows(ivar)=irow
	ivar=ivar+1
	endif



	if(j.gt.1)then
	var(ivar)=dr*sin(tetapop(j))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))
	ncols(ivar)=irow-1
	nrows(ivar)=irow
	ivar=ivar+1
	endif


		
!	a(i,i)=-(rpop(i)**2+rpop(i+1)**2)/dr   !!!AA
		var(ivar)=-(rpop(i)**2+rpop(i+1)**2)/dr &
			 -dr*(sin(tetapop(j+1))+sin(tetapop(j)))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))
		if(j.eq.1)then
		var(ivar)=var(ivar)+dr*sin(tetapop(j))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))
		elseif(j.eq.nteta)then
		var(ivar)=var(ivar)+dr*sin(tetapop(j+1))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))
		endif		
		ncols(ivar)=irow
		nrows(ivar)=irow
		ivar=ivar+1



	if(j.lt.nteta)then
	var(ivar)=dr*sin(tetapop(j+1))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))
	ncols(ivar)=irow+1
	nrows(ivar)=irow
	ivar=ivar+1
	endif
	
	if(i.lt.nr)then
	var(ivar)=rpop(i+1)**2/dr!*(cos(tetapop(j+1))-cos(tetapop(j)))    !!!BB
	ncols(ivar)=irow+nteta
	nrows(ivar)=irow
	ivar=ivar+1
	endif


	
	!!!!RHS
	
	v(irow)=2*pi/3.*(rpop(i+1)**3-rpop(i)**3)!*abs(cos(tetapop(j))-cos(tetapop(j+1)))
	rhs(irow)=-2*g*rho(irow)*v(irow)     !4*pi/3.*(rpop(i+1)**3-rpop(i)**3)     !!DD
	if(i.eq.nr)then
	rhs(irow)=rhs(irow)-out_bound*rpop(nr+1)**2/dr	
	endif
	if(j.eq.1)then
	rhs(irow)=rhs(irow)!-teta_bound(i)*dr*sin(tetapop(j))/dteta
	elseif(j.eq.nteta)then
	rhs(irow)=rhs(irow)!-teta_bound(i)*dr*sin(tetapop(j+1))/dteta
	endif


	enddo
	

	
 	open(1,file='2dmtx.dat')
 	do i=1,nnzero
 	write(1,*) nrows(i),ncols(i),var(i)
 	enddo
 	close(1)
 	
 	open(1,file='2drhs')
 	do i=1,nrow
 	write(1,*) rhs(i)
 	enddo
 	close(1)
	
	!stop


	
!	write(*,*) point
!	stop
	

	
	do i=1,nrow
	x(i)=teta_bound((i-1)/nteta+1)
	x0(i)=x(i)
	enddo

	call amux(nrow,x,y,var,ncols,point)

	open(1,file='2dx0.dat')
	do i=1,nrow
	write(1,*) x(i),y(i),rhs(i)
	if(mod(i,nteta).eq.0)then
	write(1,*)
	endif
	enddo
		
 	call mgmres_st (nrow, nnzero, nrows, ncols, var, x, rhs, nrow, nrow-1, 0.,1.e-5 )
	

	call amux(nrow,x,y,var,ncols,point)	
	!write(*,*) out_bound
	
	open(1,file='2dx.dat')
	do i=1,nrow
	write(1,*) x(i),y(i),rhs(i)
	if(mod(i,nteta).eq.0)then
	write(1,*)
	endif
	enddo
	!write(1,*) r(nr)+dr,out_bound
	close(1)
	
	open(1,file='2dr.dat')
	do i=1,nrow

	if(mod(i,nteta).eq.1)then
	write(1,*) r(i/nteta+1),x(i),-2./3.*pi*r(i/nteta+1)**2*rho(i)*g-3./5.*pi*r_out**3*rho(1)*g/r_out
	endif
	enddo
	write(1,*) r(nr)+dr,out_bound,-2./3.*pi*(r(nr)+dr)**2*rho(nr)*g-3./5.*pi*r_out**3*rho(1)*g/r_out
	close(1)	
	end
	
	subroutine amux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(*)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(n)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
    do k = ia(i), ia(i+1)-1
      t = t + a(k) * x(ja(k))
    end do

    y(i) = t

  end do

  return
end	
