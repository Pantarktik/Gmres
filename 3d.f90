	program poisson3d
	implicit none
	integer, parameter:: nr=30,nteta=30,nphi=20
	integer,parameter:: nrow=nr*nteta*nphi
	integer,parameter :: np=nrow+1	
	integer,parameter:: nnzero=7*nrow-2*nr*nphi-2*nteta*nphi!-2*nr*nteta
	real*8,parameter :: g=6.67d-7,pi=3.1415926d0	
	
	real*8 x(nrow),r(nr),dr,rhs(nrow),r_in,r_out,rpop(nr+1),phi(nphi),dphi
	real*8 rho(nrow),out_bound,var(nnzero),v(nrow),mus,y(nrow),x0(nrow)
	real*8 teta(nteta),dteta,tetapop(nteta+1),teta_bound(nr),teta_niz,teta_verh
	integer ncols(nnzero),nrows(nnzero),point(nrow+1),nvarinrow,iwk(nnzero)
	
	integer i,j,k,ivar,irow,ii,jj,kk,ierr,nstr_p,nstr_r,nstr_v,nstr,m
	
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


	!!!!!!!!!PHI
	
	dphi=2.*pi/nphi
	phi(1)=dphi/2.
	do k=2,nphi
	phi(i)=phi(i-1)+dphi
	enddo



	!!!!!!!!RHO
	
	
	rho=1.
	
	!!!!!!!!Boundary
	

	
	out_bound=-4./3.*pi*(r_out+dr/2.)**2*rho(1)*g	
	
	!!!!MATRIX
	

	ivar=1
	var=0.
	point(1)=1
	
!	open(1,file='ij')
	do irow=1,nrow
	
		i=irow/(nteta*nphi)+1
		if(mod(irow,(nteta*nphi)).eq.0)i=i-1
		
		j=mod(irow,(nteta*nphi))/nphi+1
		if(mod(mod(irow,(nteta*nphi)),nphi).eq.0)j=j-1
		if(j.eq.0)j=nteta
		
		k=mod(mod(irow,(nteta*nphi)),nphi)
		if(mod(mod(irow,(nteta*nphi)),nphi).eq.0)k=nphi
	
		nvarinrow=7
		if(i.eq.1.or.i.eq.nr)then
			nvarinrow=nvarinrow-1
		end if
		if(j.eq.1.or.j.eq.nteta)then
			nvarinrow=nvarinrow-1
		end if
!		if(k(irow).eq.1.or.k(irow).eq.nphi)then
!			nvarinrow(irow)=nvarinrow(irow)-1
!		end if

	point(irow+1)=point(irow)+nvarinrow
	
	if(irow.eq.nrow.and.point(nrow+1).ne.nnzero+1)then
	write(*,*)'error 1, point(nrow+1)=',point(nrow+1),nnzero+1
	stop
	end if



	if(i.gt.1)then
	var(ivar)=rpop(i)**2/dr*dphi
	ncols(ivar)=irow-nteta*nphi
	nrows(ivar)=irow
	ivar=ivar+1
	endif



	if(j.gt.1)then
	var(ivar)=dr*sin(tetapop(j))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))*dphi
	ncols(ivar)=irow-nphi
	nrows(ivar)=irow
	ivar=ivar+1
	endif

	if(k.eq.nphi)then
	var(ivar)=dr/dphi*(log( abs(tan(tetapop(j+1)/2.)) ) - log( abs(tan(tetapop(j)/2.)) ) )/(cos(tetapop(j+1))-cos(tetapop(j)))!DLYA K+1
	ncols(ivar)=irow-k+1
	nrows(ivar)=irow
	ivar=ivar+1	
	endif
	
	if(k.gt.1)then
	var(ivar)=dr/dphi*(log( abs(tan(tetapop(j+1)/2.)) ) - log( abs(tan(tetapop(j)/2.)) ) )/(cos(tetapop(j+1))-cos(tetapop(j)))
	ncols(ivar)=irow-1
	nrows(ivar)=irow
	ivar=ivar+1	
	end if	
		
!	a(i,i)=-(rpop(i)**2+rpop(i+1)**2)/dr   !!!AA
		var(ivar)=-(rpop(i)**2+rpop(i+1)**2)/dr*dphi &
			 -dr*(sin(tetapop(j+1))+sin(tetapop(j)))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))*dphi &
			 -2*dr/dphi*(log( abs(tan(tetapop(j+1)/2.)) ) - log( abs(tan(tetapop(j)/2.)) ) )/(cos(tetapop(j+1))-cos(tetapop(j)))
		!!!!! GRANUSLOVIA DLYA TETA	 
		if(j.eq.1)then
		var(ivar)=var(ivar)+dr*sin(tetapop(j))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))*dphi
		elseif(j.eq.nteta)then
		var(ivar)=var(ivar)+dr*sin(tetapop(j+1))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))*dphi
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		ncols(ivar)=irow
		nrows(ivar)=irow
		ivar=ivar+1

	if(k.lt.nphi)then
	var(ivar)=dr/dphi*(log( abs(tan(tetapop(j+1)/2.)) ) - log( abs(tan(tetapop(j)/2.)) ) )/(cos(tetapop(j+1))-cos(tetapop(j)))
	ncols(ivar)=irow+1
	nrows(ivar)=irow
	ivar=ivar+1	
	end if

	if(k.eq.1)then
	var(ivar)=dr/dphi*(log( abs(tan(tetapop(j+1)/2.)) ) - log( abs(tan(tetapop(j)/2.)) ) )/(cos(tetapop(j+1))-cos(tetapop(j)))!DLYA K-1
	ncols(ivar)=irow-k+nphi
	nrows(ivar)=irow
	ivar=ivar+1	
	endif
	
	if(j.lt.nteta)then
	var(ivar)=dr*sin(tetapop(j+1))/dteta/(cos(tetapop(j+1))-cos(tetapop(j)))*dphi
	ncols(ivar)=irow+nphi
	nrows(ivar)=irow
	ivar=ivar+1
	endif
	
	if(i.lt.nr)then
	var(ivar)=rpop(i+1)**2/dr*dphi!*(cos(tetapop(j+1))-cos(tetapop(j)))    !!!BB
	ncols(ivar)=irow+nteta*nphi
	nrows(ivar)=irow
	ivar=ivar+1
	endif


	
	!!!!RHS
	
	v(irow)=dphi/3.*(rpop(i+1)**3-rpop(i)**3)!*abs(cos(tetapop(j))-cos(tetapop(j+1)))
	rhs(irow)=-4*pi*g*rho(irow)*v(irow)     !4*pi/3.*(rpop(i+1)**3-rpop(i)**3)     !!DD
	if(i.eq.nr)then
	rhs(irow)=rhs(irow)-out_bound*rpop(nr+1)**2/dr*dphi	
	endif
	if(j.eq.1)then
	rhs(irow)=rhs(irow)!-teta_bound(i)*dr*sin(tetapop(j))/dteta
	elseif(j.eq.nteta)then
	rhs(irow)=rhs(irow)!-teta_bound(i)*dr*sin(tetapop(j+1))/dteta
	endif


	enddo
	

	
 	open(1,file='3dmtx.dat')
 	do i=1,nnzero
 	write(1,*) nrows(i),ncols(i),var(i)
 	enddo
 	close(1)
 	
 	open(1,file='3drhs')
 	do i=1,nrow
 	write(1,*) rhs(i)
 	enddo
 	close(1)
	
	!stop


	
!	write(*,*) point
!	stop
	
!	open(1,file='1dx.dat')
!	do i=1,nr
!	read(1,*)mus,x0(i)
!	enddo
	x0=0.
	x=x0

	call amux(nrow,x,y,var,ncols,point)

	open(1,file='3dx0.dat')
	do i=1,nrow
	write(1,*) x(i),y(i),rhs(i)
	if(mod(i,nteta).eq.0)then
	write(1,*)
	endif
	enddo


	nstr_p=np/10
	
	nstr_r=nnzero/10
	
	nstr_v=nnzero/5			
	
	open(1,file='phi_matrix.rua')
	write(1,'(A72,A8)') 'TEST MATRIX','example'
	write(1,'(5i14)') nstr,nstr_p,nstr_r,nstr_v,0
	write(1,'(A3,11x,4i14)') 'RUA',nrow,nrow,nnzero,0
	write(1,'(2A16,2A20)') '(10I8)','(10I8)','(1P5E16.8)','(1P5E16.8)'        
	

	
	do ii=1,nstr_p
	write(1,'(10i8)') (point(m),m=10*(ii-1)+1,10*ii)
	enddo
	if(mod(np,10).gt.0)then
	write(1,'(10i8)') (point(m),m=10*nstr_p+1,np)
	end if

!	write(1,*) 'ENDMASSIVE1'

	do ii=1,nstr_r
	write(1,'(10i8)') (ncols(m),m=10*(ii-1)+1,10*ii)
	enddo	
	if(mod(nnzero,10).gt.0)then
	write(1,'(10i8)') (ncols(m),m=10*nstr_r+1,nnzero)	
	end if
	
!	write(1,*) 'ENDMASSIVE2'

	do ii=1,nstr_v
	write(1,'(1P5E16.8)') (var(m),m=5*(ii-1)+1,5*ii)
	enddo	
	if(mod(nnzero,5).gt.0)then
	write(1,'(1P5E16.8)') (var(m),m=5*nstr_v+1,nnzero)
	endif
	
	close(1)

	call transp(nrow,nrow,var,ncols,point,iwk,ierr)

!	write(*,*) 'ttt2'
			
	open(1,file='phi_matrix_t.rua')
	write(1,'(A72,A8)') 'TEST MATRIX','example'
	write(1,'(5i14)') nstr,nstr_p,nstr_r,nstr_v,0
	write(1,'(A3,11x,4i14)') 'RUA',nrow,nrow,nnzero,0
	write(1,'(2A16,2A20)') '(10I8)','(10I8)','(1P5E16.8)','(1P5E16.8)'        

	
	do ii=1,nstr_p
	write(1,'(10i8)') (point(m),m=10*(ii-1)+1,10*ii)
	enddo
	if(mod(np,10).gt.0)then
	write(1,'(10i8)') (point(m),m=10*nstr_p+1,np)
	end if

!	write(1,*) 'ENDMASSIVE1'

	do ii=1,nstr_r
	write(1,'(10i8)') (ncols(m),m=10*(ii-1)+1,10*ii)
	enddo	
	if(mod(nnzero,10).gt.0)then
	write(1,'(10i8)') (ncols(m),m=10*nstr_r+1,nnzero)	
	end if
	
!	write(1,*) 'ENDMASSIVE2'

	do ii=1,nstr_v
	write(1,'(1P5E16.8)') (var(m),m=5*(ii-1)+1,5*ii)
	enddo	
	if(mod(nnzero,5).gt.0)then
	write(1,'(1P5E16.8)') (var(m),m=5*nstr_v+1,nnzero)
	endif
	
	close(1)		
	
	write(*,*) 'ttt',nrow,nnzero
	call transp(nrow,nrow,var,ncols,point,iwk,ierr)	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		
! 	call mgmres_st (nrow, nnzero, nrows, ncols, var, x, rhs, 10, nrow-1, 0.,1.e-5 )
	
	open(1,file='x1c.dat')
	do i=1,nrow
	read(1,*) x(i)
	enddo
	close(1)

	call amux(nrow,x,y,var,ncols,point)	
	!write(*,*) out_bound
	
	open(1,file='3dx.dat')
	do i=1,nrow
	write(1,*) x(i),y(i),rhs(i)
	if(mod(i,nteta*nphi).eq.0)then
	write(1,*)
	endif
	enddo
	!write(1,*) r(nr)+dr,out_bound
	close(1)
	
	open(1,file='3dr.dat')
	do i=1,nrow

	if(mod(i,nteta*nphi).eq.1)then
	write(1,*) r(i/(nteta*nphi)+1),x(i),-2./3.*pi*r(i/(nteta*nphi)+1)**2*rho(i)*g!-2./3.*pi*(r(nr)+dr)**2*rho(1)*g
	endif
	enddo
	write(1,*) r(nr)+dr,out_bound,-2./3.*pi*(r(nr)+dr)**2*rho(nr)*g!-2./3.*pi*(r(nr)+dr)**2*rho(1)*g
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

subroutine transp ( nrow, ncol, a, ja, ia, iwk, ierr )

!*****************************************************************************80
!
!! TRANSP carries out in-place transposition routine.
!
!  Discussion:
!
!    This routine transposes a matrix stored in compressed sparse row
!    format.  The transposition is done in place in that the arrays 
!    A, JA, and IA of the transpose are overwritten onto the original arrays.
!
!    If you do not need the transposition to be done in place
!    it is preferrable to use the conversion routine csrcsc
!    (see conversion routines in formats).
!
!    The entries of the output matrix are not sorted (the column
!    indices in each are not in increasing order).  Use CSRCSC
!    if you want them sorted.
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
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Workspace, integer ( kind = 4 ) IWK(*), of the same length as JA.
!
! on return:
!
!
! ncol      = actual row dimension of the transpose of the input matrix.
!         Note that this may be <= the input value for ncol, in
!         case some of the last columns of the input matrix are zero
!         columns. In the case where the actual number of rows found
!         in transp(A) exceeds the input value of ncol, transp will
!         return without completing the transposition. see ierr.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NCOL+1), the transposed
!    matrix in CSR Compressed Sparse Row format.
!
! ierr      = integer ( kind = 4 ). error message. If the number of rows for the
!         transposed matrix exceeds the input value of ncol,
!         then ierr is  set to that number and transp quits.
!         Otherwise ierr is set to 0 (normal return).
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) inext
  integer ( kind = 4 ) init
  integer ( kind = 4 ) iwk(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnz
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
	

	
  ierr = 0
  nnz = ia(nrow+1) - 1
!
!  Determine the column dimension.
!
  jcol = 0


  do k = 1, nnz
    jcol = max ( jcol, ja(k) )
  end do

  if ( ncol < jcol ) then
     ierr = jcol
     return
  end if
!	write(*,*) 'ttt3',jcol,ncol

!
!  Convert to coordinate format.  Use IWK for row indices.
!

!  ncol = jcol
!	write(*,*) 'ttt4'
  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      iwk(k) = i
    end do
  end do
!
!  Find pointer array for transpose.
!
  ia(1:ncol+1) = 0

  do k = 1, nnz
    i = ja(k)
    ia(i+1) = ia(i+1) + 1
  end do
  ia(1) = 1

  do i = 1, ncol
    ia(i+1) = ia(i) + ia(i+1)
  end do
!
!  Loop for a cycle in chasing process.
!
  init = 1
  k = 0

 5    continue

  t = a(init)
  i = ja(init)
  j = iwk(init)
  iwk(init) = -1

 6 continue

   k = k + 1
!
!  Current row number is I.  Determine where to go.
!
  l = ia(i)
!
!  Save the chased element.
!
  t1 = a(l)
  inext = ja(l)
!
!  Then occupy its location.
!
  a(l) = t
  ja(l) = j
!
!  Update pointer information for next element to be put in row i.
!
  ia(i) = l + 1
!
!  Determine next element to be chased.
!
  if ( iwk(l) < 0 ) then
    go to 65
  end if

  t = t1
  i = inext
  j = iwk(l)
  iwk(l) = -1

  if ( k < nnz ) then
    go to 6
  end if

  do i = ncol, 1, -1
    ia(i+1) = ia(i)
  end do

  ia(1) = 1

  return

 65   continue

  init = init + 1

  if ( nnz < init ) then

    do i = ncol, 1, -1
      ia(i+1) = ia(i)
    end do

    ia(1) = 1

    return
  end if

  if ( iwk(init) < 0 ) then
    go to 65
  end if
!
!  Restart chasing.
!
  go to 5
end
