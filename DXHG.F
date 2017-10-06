	SUBROUTINE DXHG(M,N,X,Y,A,Q,S,R,V,U,B)
	REAL(KIND=8),DIMENSION(M,N)::X
	REAL(KIND=8),DIMENSION(N)::Y
	REAL(KIND=8),DIMENSION(M+1)::A
	REAL(KIND=8),DIMENSION(M+1,M+1)::B
	REAL(KIND=8),DIMENSION(M)::V
	REAL(KIND=8) Q,S,R,U,YY,DYY,P,PP
	MM=M+1
	B(1,1)=N
	DO J=2,MM
		B(1,J)=0
		DO I=1,N
			B(1,J)=B(1,J)+X(J-1,I)
		END DO
		B(J,1)=B(1,J)
	END DO
	DO I=2,MM
		DO J=I,MM
			B(I,J)=0
			DO K=1,N
				B(I,J)=B(I,J)+X(I-1,K)*X(J-1,K)
			END DO
			B(J,I)=B(I,J)
		END DO
	END DO
	A(1)=0
	DO I=1,N
		A(1)=A(1)+Y(I)
	END DO
	DO I=2,MM
		A(I)=0
		DO J=1,N
			A(I)=A(I)+X(I-1,J)*Y(J)
		END DO
	END DO
	CALL CHOLESKY(B,MM,1,A,L)
	YY=0
	DO I=1,N
		YY=YY+Y(I)
	END DO
	YY=YY/N
	Q=0
	DYY=0
	U=0
	DO I=1,N
		P=A(1)
		DO J=1,M
			P=P+A(J+1)*X(J,I)
		END DO
		Q=Q+(Y(I)-P)*(Y(I)-P)
		DYY=DYY+(Y(I)-YY)*(Y(I)-YY)
		U=U+(YY-P)*(YY-P)
	END DO
	S=SQRT(Q/N)
	R=SQRT(1-Q/DYY)
	DO J=1,M
		P=0
		DO I=1,N
			PP=A(1)
			DO K=1,M
				IF(K/=J)PP=PP+A(K+1)*X(K,I)
			END DO
			P=P+(Y(I)-PP)*(Y(I)-PP)
		END DO
		V(J)=SQRT(1-Q/P)
	END DO
	END

	SUBROUTINE CHOLESKY(C,N,M,D,L)
	REAL(KIND=8),DIMENSION(N,N)::C
	REAL(KIND=8),DIMENSION(N,M)::D
	L=1
	IF(ABS(C(1,1))<1.0E-10)THEN
		L=0
		WRITE(*,'(" FAIL")')
		RETURN
	END IF
	C(1,1)=SQRT(C(1,1))
	DO J=2,N
		C(1,J)=C(1,J)/C(1,1)
	END DO
	DO I=2,N
		DO J=2,I
			C(I,I)=C(I,I)-C(J-1,I)*C(J-1,I)
		END DO
		IF(ABS(C(I,I))<1.0E-10)THEN
			L=0
			WRITE(*,'(" FAIL")')
			RETURN
		END IF
		C(I,I)=SQRT(C(I,I))
		IF(I/=N)THEN
			DO J=I+1,N
				DO K=2,I
					C(I,J)=C(I,J)-C(K-1,I)*C(K-1,J)
				END DO
				C(I,J)=C(I,J)/C(I,I)
			END DO
		END IF
	END DO
	DO J=1,M
		D(1,J)=D(1,J)/C(1,1)
		DO I=2,N
			DO K=2,I
				D(I,J)=D(I,J)-C(K-1,I)*D(K-1,J)
			END DO
			D(I,J)=D(I,J)/C(I,I)
		END DO
	END DO
	DO J=1,M
		D(N,J)=D(N,J)/C(N,N)
		DO K=N,2,-1
			DO I=K,N
				D(K-1,J)=D(K-1,J)-C(K-1,I)*D(I,J)
			END DO
			D(K-1,J)=D(K-1,J)/C(K-1,K-1)
		END DO
	END DO
	END
