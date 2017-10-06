	program main
	use msimsl
	parameter(N=5)
	real ryx1,ryx2
	dimension yt(5),x1t(5),x2t(5),dety(4),detx1(4),detx2(4),
     &     yx1(4),yx2(4)   
	DATA  (yt(i), i=1,5) / 90.,97.63,104.88,118.84,156.34/
	DATA  (x1t(i), i=1,5) / 64.19,70.67,74.57,86.76,109.95/
	DATA  (x2t(i), i=1,5) / 25.81,26.96,30.31,32.08,46.39/
C   Calculating medium varibles
	do i=1,N-1
  	dety(i)=yt(i+1)-yt(i)
  	detx1(i)=x1t(i+1)-x1t(i)
  	detx2(i)=x2t(i+1)-x2t(i)
	end do
C   Calculating corelation function
	do i=1,N-1
  	yx1(i)=1/(1+abs(detx1(i)/x1t(i)-dety(i)/yt(i)))
  	yx2(i)=1/(1+abs(detx2(i)/x2t(i)-dety(i)/yt(i)))
	end do
C   Calculating corelation
	ryx1=0.0
	ryx2=0.0
	do i= 1,N-1
  	ryx1=ryx1+yx1(i)
  	ryx2=ryx2+yx2(i)
	end do
	ryx1=ryx1/(N-1)
	ryx2=ryx2/(N-1)
C   Write calculated result
	WRITE(6,*) ryx1,ryx2
	pause	
	end
