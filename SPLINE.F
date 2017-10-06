! SPLINE.FOR -- Natural Cubic Spline
! Call SplineCalculate with the two arrays of control points (and
! the number of elements in each array).  To interpolate the spline,
! Call SplineEvaluate for each 'X' for which you want a 'Y' along the 
! spline's curve.
!

	subroutine SplineCalculate( inx, ina, numx )
!ms$if .not. defined(LINKDIRECT)
!ms$attributes dllexport :: SplineCalculate
!ms$endif

	implicit none
	real*4 ina(*), inx(*)
	integer numx

	integer MAXPOINTS
	parameter( MAXPOINTS = 1000 )
	integer nxs
	real*4 x(0:MAXPOINTS)
	real*4 a(0:MAXPOINTS), b(0:MAXPOINTS)
	real*4 c(0:MAXPOINTS), d(0:MAXPOINTS)
	common /SplineData/ a, b, c, d, nxs, x

	integer i
	real*4 alpha(0:MAXPOINTS), l(0:MAXPOINTS)
	real*4 mu(0:MAXPOINTS), z(0:MAXPOINTS)
	real*4 h(0:MAXPOINTS)

	! ===== Step 1
	nxs = numx-1
	do i = 0, nxs
		x(i) = inx(i+1)
		a(i) = ina(i+1)
	end do
	do i = 0, nxs-1
		h(i) = x(i+1) - x(i)
	end do
	x(nxs+1) = 40000

	! ===== Step 2
	do i = 1, nxs-1
		alpha(i) = 3.0 * (a(i+1) * h(i-1)-a(i) *
     &	(x(i+1)-x(i-1)) + a(i-1) * h(i)) /
     &	(h(i-1) * h(i))
	end do

	! ===== Step 3
	l(0) = 1.0
	mu(0) = 0.0
	z(0) = 0.0

	! ===== Step 4
	do i = 1, nxs-1
		l(i) = 2.0 * (x(i+1)-x(i-1))-h(i-1) * mu(i-1)
		mu(i) = h(i) / l(i)
		z(i) = (alpha(i)-h(i-1) * z(i-1)) / l(i)
	end do

	! ===== Step 5
	l(nxs) = 1.0
	z(nxs) = 0.0
	c(nxs) = 0.0

	! ===== Step 6
	do i = nxs-1, 0, -1
		c(i) = z(i) - mu(i)*c(i+1)
		b(i) = (a(i+1)-a(i))/h(i) - h(i)*
     &	(c(i+1)+2.0 * c(i)) / 3.0
		d(i) = (c(i+1)-c(i)) / (3.0 * h(i))
	end do

	end

	subroutine SplineEvaluate( evalx, evaly )
!ms$if .not. defined(LINKDIRECT)
!ms$attributes dllexport :: SplineEvaluate
!ms$endif

	real*4 evalx, evaly

	integer MAXPOINTS
	parameter( MAXPOINTS = 1000 )
	integer nxs
	real*4 x(0:MAXPOINTS)
	real*4 a(0:MAXPOINTS), b(0:MAXPOINTS)
	real*4 c(0:MAXPOINTS), d(0:MAXPOINTS)
	common /SplineData/ a, b, c, d, nxs, x

	real*4 term
	integer i

	i = 0
	do while( .not. ((x(i) .le. evalx) .and. (x(i+1) .gt. evalx)) )
		i = i + 1
	end do

	term = evalx - x(i)
	evaly = a(i) + b(i)*term + c(i)*term*term +
     &		  d(i)*term*term*term

	end
