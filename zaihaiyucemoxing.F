c  The data is the rainfall of JINHUA from 1953-1986,and follow is the main program.
	program main
	use msimsl
	parameter(N=34,valve=700.0,M=7,Long=8)
	REAL x0(N),p(M),pp(M),p1(M),forcastp(Long), y(6,1),ainv(2,2),
     &medium(2,6),b(M-1,2),a(2,M-1),u(2,1)
	data (x0(i),i=1,34)/336,1171,792,444,311,663,553,501,554,710,
     &509,424,536,437,720,645,723,607,573,384,
     &917,533,535,627,773,333,360,402,344,439,
     &631,660,411,392/
c   Take the catastrophic sequence   
	k=1
	do i=1,N
	if(x0(i)>=700.0)	then
	p(k)=i
	k=k+1
	end if
	end do
c   Take the first accumulation and make the matrix B and Y.
	p1(1)=p(1)
	do i=2,M
		p1(i)=p1(i-1)+p(i)
	end do
	do i=2,M
		y(i-1,1)=p(i)
	end do
	do i=1,M-1
		b(i,1)=-(p1(i)+p1(i+1))/2
		b(i,2)=1
	end do
c   Calculating matrix u and forcasting 
CALL MXTXF (6, 2, b, 6, 2, a, 2)
CALL LINRG (2, A, 2, AINV, 2)
	medium=MATMUL(AINV, TRANSPOSE(b))
	u=matmul(medium,y)
	forcastp(1)=p1(1)
	do i=1,Long
		forcastp(i+1)=(p1(1)-u(2,1)/u(1,1))*exp(-(1,1)*i)+u(2,1)/u(1,1)
	end do
	pp(1)=P1(1)
	do i=2,M
		pp(i)=forcastp(i)-forcastp(i-1)
	end do
	write(*,*)forcastp
	write(*,*)pp
	write(*,*)'The next catastrophic year is:'
	write(*,*) 1952+(forcastp(Long)-forcastp(Long-1))
	pause
	end
