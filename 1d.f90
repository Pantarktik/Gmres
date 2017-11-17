	program poisson1d
	implicit none
	integer, parameter:: nr=10
	integer,parameter:: nnzero=3*nr-2
	real*8,parameter :: g=6.67d-7,pi=3.1415926d0	
	
	real*8 x(nr),r(nr),dr(nr),rhs(nr),a(nr,nr),r_in,r_out,rpop(nr+1)
	real*8 rho(nr),out_bound,ddr,var(nnzero)
	integer ncols(nnzero),nrows(nnzero)
	
	integer i,j,k,ivar
	
	!!!!!!!!!!!!SETKA
	
	r_in=0.
	r_out=100.
	
	dr=(r_out-r_in)/nr
	r(1)=r_in+dr(1)/2.
	do i=1,nr
	r(i)=r(i-1)+dr(i-1)/2.+dr(i)/2.
	enddo
	
	do i=1,nr
	rpop(i)=r(i)-dr(i)/2.
	enddo
	rpop(nr+1)=r_out
	
	
	ddr=dr(1)
	
	rho=1.
	out_bound=-4./3.*pi*r_out**3*rho(1)*g/r_out	
	
	!!!!MATRIX
	
	a=0.
	ivar=1
	
	
	do i=1,nr
	if(i.gt.1)then
	a(i,i-1)=rpop(i)**2/ddr       !!!!!!CC
	var(ivar)=a(i,i-1)
	ncols(ivar)=i-1
	nrows(ivar)=i
	ivar=ivar+1
	endif	
	a(i,i)=-(rpop(i)**2+rpop(i+1)**2)/ddr   !!!AA
	var(ivar)=a(i,i)
	ncols(ivar)=i
	nrows(ivar)=i
	ivar=ivar+1
	if(i.lt.nr)then
	a(i,i+1)=rpop(i+1)**2/ddr    !!!BB
	var(ivar)=a(i,i+1)
	ncols(ivar)=i+1
	nrows(ivar)=i
	ivar=ivar+1
	endif
	enddo
	
	!!!!RHS
	
	do i=1,nr
	rhs(i)=-g*rho(i)*4*pi/3.*(rpop(i+1)**3-rpop(i)**3)     !!DD
	enddo
	
	rhs(nr)=rhs(nr)-out_bound*rpop(nr+1)**2/ddr
	
	open(1,file='1dmtx.dat')
	do i=1,nnzero
	write(1,*) nrows(i),ncols(i),var(i)
	enddo
	close(1)
	
	open(1,file='1drhs')
	do i=1,nr
	write(1,*) rhs(i)
	enddo
	close(1)
		
 	call mgmres_st (nr, nnzero, nrows, ncols, var, x, rhs, nr, 10, 0.,1.e-5 )
	
	write(*,*) out_bound
	
	open(1,file='1dx.dat')
	do i=1,nr
	write(1,*) r(i),x(i)
	enddo
	write(1,*) r(nr)+ddr,out_bound
	close(1)
	end
