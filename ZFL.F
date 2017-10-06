	SUBROUTINE ZFL(X,N,P,V,B,G)
	INTEGER::P
	REAL(8),DIMENSION(P,N)::X
	REAL(8),DIMENSION(P)::XV
	REAL(8),DIMENSION(P)::B,G
	REAL(8),DIMENSION(P,P)::V
	REAL(8),DIMENSION(P,P)::S
!   求每一点的均值
	XV=0
	DO J=1,P
		DO I=1,N
			XV(J)=XV(J)+X(I,J)
		END DO
		XV(J)=XV(J)/N
	END DO
!	求协方差阵,该协方差阵为对称方阵
	S=0
	DO I=1,P
		DO J=1,P
			DO K=1,N
				S(I,J)=S(I,J)+(X(I,K)-XV(I))*(X(J,K)-XV(J))
			END DO
			S(I,J)=S(I,J)/N
		END DO
	END DO
	CALL JCB(S,P,1.0E-8,V,B,L)  !	调用雅可比法求解矩阵的特征值和特征向量
	CALL MSORT(B,P,V)  !   将特征值按大小排列
	GP=0
	DO I=1,P
		GP=GP+B(I)
	END DO
	DO I=1,P
		G(I)=B(I)/GP
	END DO
	END

	SUBROUTINE JCB(A,N,EPS,V,B,L)
!   A:调用时存放实对称矩阵
!   B:返回时存放矩阵的全部特征值
!   V:存放特征向量,其中第i列为与第i个特征值相对应的特征向量
!	EPS:存放精度要求
	REAL(8),DIMENSION(N,N)::A,V
	REAL(8),DIMENSION(N)::B
	REAL(8)::FM,CN,SN,OMEGA,X,Y
	INTEGER::P,Q
	L=1
	V=0
	DO I=1,N
		V(I,I)=1
	END DO
10	FM=0
	DO I=1,N
		DO J=1,N
			IF(I/=J.AND.ABS(A(I,J))>FM)THEN
				FM=ABS(A(I,J))
				P=I
				Q=J
			END IF
		END DO
	END DO
	IF(FM<EPS)THEN
		L=1
		RETURN
	END IF
	IF(L>100)THEN
		L=0
		RETURN
	END IF
	L=L+1
	X=-A(P,Q)
	Y=(A(Q,Q)-A(P,P))/2
	OMEGA=X/SQRT(X*X+Y*Y)
	IF(Y<0)OMEGA=-OMEGA
	SN=1+SQRT(1-OMEGA*OMEGA)
	SN=OMEGA/SQRT(2*SN)
	CN=SQRT(1-SN*SN)
	FM=A(P,P)
	A(P,P)=FM*CN*CN+A(Q,Q)*SN*SN+A(P,Q)*OMEGA
	A(Q,Q)=FM*SN*SN+A(Q,Q)*CN*CN-A(P,Q)*OMEGA
	A(P,Q)=0
	A(Q,P)=0
	DO J=1,N
		IF(J/=P.AND.J/=Q)THEN
			FM=A(P,J)
			A(P,J)=FM*CN+A(Q,J)*SN
			A(Q,J)=-FM*SN+A(Q,J)*CN
		END IF
	END DO
	DO I=1,N
		IF(I/=P.AND.I/=Q)THEN
			FM=A(I,P)
			A(I,P)=FM*CN+A(I,Q)*SN
			A(I,Q)=-FM*SN+A(I,Q)*CN
		END IF
	END DO
	DO I=1,N
		FM=V(I,P)
		V(I,P)=FM*CN+V(I,Q)*SN
		V(I,Q)=-FM*SN+V(I,Q)*CN
	END DO
	DO I=1,N
		B(I)=A(I,I)
	END DO
	GOTO 10
	END

	SUBROUTINE MSORT(B,N,V)
!   将特征值按大小排列
	REAL(8),DIMENSION(N)::B
	REAL(8),DIMENSION(N,N)::V
	REAL(8),DIMENSION(N)::V1
	REAL(8)::B1
	M=N
20	IF(M>0)THEN
		J=M-1
		M=0
		DO I=1,J
			IF(B(I).LT.B(I+1))THEN
				B1=B(I)
				B(I)=B(I+1)
				B(I+1)=B1
				M=I
				DO K=1,N
					V1(K)=V(K,I)
					V(K,I)=V(K,I+1)
					V(K,I+1)=V1(K)
				END DO
			ENDIF
		ENDDO
		GOTO 20
	ENDIF
	END


	PROGRAM MAIN
	INTEGER,PARAMETER::N=50
	INTEGER,PARAMETER::P=4
	REAL(8),DIMENSION(N,P)::X
	REAL(8),DIMENSION(P)::B,G
	REAL(8),DIMENSION(P,P)::V
	OPEN(10,FILE='DXHG.DAT')
	DO I=1,N
		READ(10,*)M,A,(X(I,J),J=1,P)
	END DO
	CALL ZFL(X,N,P,V,B,G)
	OPEN(12,FILE='ZFL.DAT')
	DO I=1,P
		WRITE(12,10)I,B(I),I,(V(J,I),J=1,P)
	END DO
 10	FORMAT(1X,'第',I1,'特征值=',D11.5,2X,'V(',I1,')=(',3D12.5)
	DO I=1,P
		WRITE(12,'(1X,"第",I1,"主分量解释的方差百分比=",F7.3,"%")')I,G(I)*100
	END DO
	END
