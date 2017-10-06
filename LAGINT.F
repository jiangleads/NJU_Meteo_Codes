!根据给定点的函数值,应用拉格朗日插值公式,计算指定插值点处的
!函数值。
!在指定点的前后各取四个结点,用七次拉格朗日插值公式。
!X,Y为输入参数,各有N个点
!T为指定插入点
!Z为输出参数,是插值点T处的近似值。
	SUBROUTINE  LAGINT(X,Y,N ,T,Z)
	DIMENSION X(N),Y(N)
	DOUBLE PRECISION X,Y,T,S
	Z=0.0
	IF (N.LE.0) RETURN
	IF(N.EQ.1) THEN
		Z=Y(1)
	ENDIF
	IF(N.EQ.2) THEN
		Z=(Y(1)*(T-X(2))-Y(2)*(T-X(1)))/(X(1)-X(2))
		RETURN
	ENDIF
	I=1
10	IF(X(1) .LT. T) THEN !//人工加了then，可能会出错
	I=I+1
	IF(I .LE. N) GOTO 10
	ENDIF
	K=I-4
	IF(K .LT. 1) K=1
	M=I+3
	IF( M .GT. N) M=N
	DO 30 I=K,M
		S=1.0
		DO 20 J=K,M
			IF(J .NE. I) THEN
				S=S*(T-X(J))/(X(I)-X(J))
			ENDIF
20		CONTINUE
		Z=Z+S*Y(I)
30	CONTINUE
	RETURN
	END
