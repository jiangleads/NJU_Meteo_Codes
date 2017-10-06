	SUBROUTINE EOF(X,P,N,LW,XF)
	INTEGER::P
	INTEGER::N
	INTEGER::LW
	REAL(8),DIMENSION(P,N)::X,XF
	REAL(8),DIMENSION(P,P)::A,V,V1
	REAL(8),DIMENSION(P,N)::T
	REAL(8),DIMENSION(P)::B,GM,GA
	REAL(8),DIMENSION(P,LW)::VF
	REAL(8),DIMENSION(LW,N)::TF
!	求X乘以X的转置,即A=XXˊ
	DO I=1,P
	  DO J=1,P
	    A(I,J)=0
	    DO K=1,N
	      A(I,J)=A(I,J)+X(I,K)*X(J,K)
	    END DO
	  END DO
	END DO
!	用Jacobi法求A的特征值和特征向量
!   返回时B存放矩阵的全部特征值,V存放特征向量为列组成的矩阵
	CALL JCB(A,P,1.0E-6,V,B,L)
	DO I=1,P
	  GA(I)=0
	  DO J=1,I
	    GA(I)=GA(I)+B(J)
	  END DO
	END DO
	DO I=1,P
	  GM(I)=GA(I)/GA(P)
	END DO
	DO I=1,P
	  DO J=1,P
	    V1(I,J)=V(J,I)
	  END DO
	END DO
	T=MATMUL(V1,X)	
	WRITE(12,'(" 特征值")')
	WRITE(12,'(<P>I10)')(I,I=1,P)
	WRITE(12,'(3X,<P>D10.4)')B
	WRITE(12,'(" 解释的方差(%)")')
	WRITE(12,'(<P>I7)')(I,I=1,P)
	WRITE(12,'(3X,<P>F7.2)')GM*100
	WRITE(12,'(" 特征向量为列组成的矩阵,即空间函数V")')
	WRITE(12,'(<P>F7.4)')((V(I,J),J=1,P),I=1,P)
	WRITE(12,'("  时间函数T")')
	WRITE(12,'(<P>F10.4)')((T(I,J),J=1,P),I=1,N)
	DO I=1,P
	  DO J=1,LW
	    VF(I,J)=V(I,J)
	  END DO
	END DO
	DO I=1,LW
	  DO J=1,N
	    TF(I,J)=T(I,J)
	  END DO
	END DO
	XF=MATMUL(VF,TF)
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
	  GOTO 15
	END IF
	IF(L>100)THEN
	  L=0
	  WRITE(12,'("L=",I3,2X,"没有达到精度要求")')L
	  GOTO 15
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
!   将特征值按大小排列
15	M=N
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
	        V1=V(K,I)
	        V(K,I)=V(K,I+1)
	        V(K,I+1)=V1
	      END DO
	    ENDIF
	  ENDDO
	  GOTO 20
	  ENDIF
	END

 	PROGRAM MAIN
 	INTEGER,PARAMETER::P=20 	!资料的空间点数
 	INTEGER,PARAMETER::N=123	!资料的时间长度
 	REAL(8),DIMENSION(P,N)::X,XF,ERROR
 	REAL(8),DIMENSION(P)::XV  !X的平均值
 	INTEGER::LP=P;LW=3
 	OPEN(12,FILE='F58Z850.DAT')
 	DO IT=1,123
 	READ(12,'(<LP>F8.1)')(X(I,IT),I=1,P)
 	END DO
 	CLOSE(12)
 	XV=0
 	DO I=1,P
 	  DO J=1,N
 	    XV(I)=XV(I)+X(I,J)
 	  END DO
 	  XV(I)=XV(I)/N
 	END DO
 	DO I=1,P
 	  DO J=1,N
 	    X(I,J)=X(I,J)-XV(I)
 	  END DO
 	END DO
 	OPEN(12,FILE='EOF.DAT')
 	CALL EOF(X,P,N,LW,XF)
 	ERROR=X-XF
 	DO I=1,P
 	  DO J=1,N
 	    XF(I,J)=XF(I,J)+XV(I)
 	  END DO
 	END DO
 	WRITE(12,'("  恢复场")')
 	WRITE(12,'(<LP>F8.1)')((XF(I,J),J=1,P),I=1,N)
 	WRITE(12,'("  误差场")')
 	WRITE(12,'(<LP>F8.1)')((ERROR(I,J),J=1,P),I=1,N)
 	CLOSE(12)
 	END
