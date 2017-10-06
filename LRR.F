C_______________________________________________________________
c	this program is logit.3 LRR(the Logit Regression Recurrence)
c	Cn=1/(XnT*Wn-1*Xn)
c	wn+1=y1n+1*(1-y1n+1)
c	dn+1=(wn+1+Xn+1*Cn*Xn+1)
c	Cn+1=Cn-Cn*Xn+1T*Xn+1*Cn/dn+1
c	bn+1=bn+Cn*Xn+1*dn+1*wn+1*(yn+1-y1n+1)
c	yn+1 -obtation
c	y1n+1 -estimator
c	y1n+1=xn+1T*bn+1
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine LRR(cn,bn,x,y,yf,m,n,ll,bs0)
C       (M1=M)
	PARAMETER(M1=12)
	real cn(m,m),x(n,m),bn(m),yf(n),y(n)
	REAL x1(m1),x2(m1),cc(m1,m1),b1(m1),xn(m1)
	bs=0.0
	do 100 k=1,n
		do j=1,m
			xn(j)=x(k,j)
		enddo
		yn=y(k)
c==
		y1=0.0
		do i=1,m
			y1=y1+xn(i)*bn(i)
		enddo
c===
	do i=1,m
		x1(i)=0.0
		do j=1,m
			x1(i)=x1(i)+cn(i,j)*xn(j)
		enddo
	enddo
c===
	dn=0.0
	do i=1,m
		dn=dn+xn(i)*x1(i)
	enddo
	if(ll.eq.0)then
		y1=1.0/(1.0+exp(-y1))
		wn=y1*(1.0-y1)
	else
		wn=1.0
	endif
	dn1=wn/(dn*wn+1.0)
	dn2=1.0/(dn*wn+1.0)
c 	write(*,'(2e12.3)')dn1,dn2
c===	
	do i=1,m
		x2(i)=0.0
		do j=1,m
			x2(i)=x2(i)+cn(j,i)*xn(j)
		enddo
	enddo

c===	
	do i=1,m
		do j=1,m
			cc(i,j)=x1(i)*dn1*x2(j)
		enddo
	enddo
	do i=1,m
		b1(i)=x1(i)*dn2*(yn-y1)
	enddo
c===	
	do i=1,m
		bn(i)=bn(i)+b1(i)
	enddo
c===	
	do i=1,m
		do j=1,m
			cn(i,j)=cn(i,j)-cc(i,j)
		enddo
	enddo
c===	
	y1=0.0
	do i=1,m
		y1=y1+xn(i)*bn(i)
	enddo
	y1=1.0/(1.0+exp(-y1))
	dy=yn-y1
	yf(k)=y1
	bs=bs+dy*dy
100	continue
	bs=bs/n
	bs0=bs
	write(*,'(f8.6)')bs0
	return
	end

c	this is the main of logit3.for
	parameter(ln=184,lm=12,lm0=11)
	real x(ln,lm),b(lm),y(ln),xw(lm,lm),yf(ln),x0(lm),
     1y2(ln)
	character*12 fa,fa1,f1,f2,f3,f4
c===========read data==============
	kk=0
	bs0=1.0
	eta=0.3
	aph=0.7
	open(2,file='mmj.dat')
	read(2,'(3i5)')m1,N,NT
	read(2,'(a)')f1
	READ(2,'(A)')F2
	READ(2,'(A)')F3
	READ(2,'(A)')F4
	close(2)
	N=N-NT
	m=m1+1
	m0=m1
	write(f1(2:3),'(i2.2)')m
	write(*,'(2i5)')n,m1
 100	open(2,file=F3)
	READ(2,'(A)')F1
	do j=1,n
		read(2,f1,end=111)(x(j,i),i=1,m1),Y(J)
	enddo
 111	close(2)

	do i=1,n
		x(i,m)=1.0
	enddo

	open(2,file='b00.dat')
	read(2,'(a)')fa1
	read(2,'(a)')
	do i=1,m-1
		read(2,fa1,end=112)ij,x0(i)
	enddo
 112	close(2)

	m1=m-1
	do i=1,n
		do j=1,m1
			x(i,j)=x(i,j)*10.0
			if(x0(j).gt.10.0)x(i,j)=x(i,j)*10.0
		enddo
		x(i,m)=1.0
	enddo

	open(2,file='b2.dat')
	read(2,'(a)')fa
	read(2,fa)(b(i),i=1,m)
	close(2)
	do i=1,n
		yf(i)=0.0
		do j=1,m
			yf(i)=yf(i)+x(i,j)*b(j)
		enddo
	enddo
	do i=1,n
		yf(i)=1.0/(1.0+exp(-yf(i)))
		if(yf(i).lt.0.0)yf(i)=0.0
	enddo
	do i=1,n
		y2(i)=y(i)-yf(i)
	enddo

	bs=0.0
	do i=1,n
		bs=bs+y2(i)*y2(i)
	enddo
	bs=bs/n
	write(*,'(a,f12.8)')' bs=',bs

	open(2,file='xw.dat')
	read(2,'(a)')fa
	read(2,fa)((xw(i,j),i=1,m),j=1,m)
	close(2)
	bs0=bs
	open(2,file='xw.dat')
	write(2,'(a)')fa
	call lrr(xw,b,x,y,yf,m,n,0,bs)
	write(2,fa)((xw(i,j),i=1,m),j=1,m)
	close(2)
	kk=kk+1
	do i=1,n
		yf(i)=0.0
		do j=1,m
			yf(i)=yf(i)+x(i,j)*b(j)
		enddo
	enddo
	do i=1,n
		yf(i)=1.0/(1.0+exp(-yf(i)))
	enddo
	do i=1,n
		y2(i)=y(i)-yf(i)
	enddo

	bs=0.0
	do i=1,n
		bs=bs+y2(i)*y2(i)
	enddo
	bs=bs/n
	write(*,'(a,f12.8)')' bs=',bs
	if(bs.gt.bs0)goto 1
	bs0=bs
	open(2,file='b2.dat')
	write(2,'(a)')fa
	write(2,fa)(b(i),i=1,m)
	close(2)
	goto 100
 1	write(*,'(a,i5,a,f8.6)')' k=',kk,' BS=',bs
	open(2,file='plj3.dat')
	write(2,'(a)')'(2f8.4)'
	do i=1,n
		write(2,'(2f8.4)')yf(i),y(i)
	enddo
	close(2)
	END

