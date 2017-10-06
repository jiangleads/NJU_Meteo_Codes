	SUBROUTINE ARMA(X,N,P,Q,M,FAI,THITA,EPS)
	INTEGER::TAO	                   !落后时间
	INTEGER::P	                   !自回归阶数
	INTEGER::Q	                   !滑动平均阶数
	INTEGER::M	                   !M=P+Q
	REAL(8),DIMENSION(N)::X		   !输入序列
	REAL(8),DIMENSION(0:P)::FAI	   !自回归系数
	REAL(8),DIMENSION(P,P)::A      !工作数组
	REAL(8),DIMENSION(P)::B        !工作数组
	REAL(8),DIMENSION(0:M)::S	   !协方差,S(0)即为方差
	REAL(8),DIMENSION(0:Q)::SC	   !自回归后的协方差
	REAL(8),DIMENSION(Q)::THITA    !滑动平均系数
	REAL(8),DIMENSION(Q)::THIT     !迭代中用的滑动系数,中间变量
	REAL(8)::A1,A2,A3  	           !A1,A2,A3:中间变量
	REAL(8)::S2A                   !S2A:自回归后的序列a(t)的方差
	REAL(8)::EPS,EP1,EP2           !EPS:迭代的精度
	S=0
	DO TAO=0,M
		DO I=1,N-TAO
			S(TAO)=S(TAO)+X(I)*X(I+TAO)
		END DO
		S(TAO)=S(TAO)/(N-TAO)
	END DO
	DO I=1,P
		DO J=1,P
			A(I,J)=S(ABS(Q+I-J))
		END DO
		B(I)=S(Q+I)
	END DO
	CALL GASJDN(A,B,P)
	FAI(1:P)=B(1:P)
	FAI(0)=-1
	A1=0
	DO I=0,P
		A1=A1+FAI(I)*FAI(I)
	END DO
	DO K=0,Q
		A2=0
		DO I=1,P
			A3=0
			DO J=0,P-I
				A3=A3+FAI(J)*FAI(J+I)
			END DO
			A2=A2+A3*(S(K+I)+S(ABS(K-I)))
		END DO
		SC(K)=A1*S(K)+A2
	END DO
	S2B=0
	THIT=0
	NN=0
	DO 				  !迭代
		NN=NN+1
		WRITE(*,'(" NN=",I3)')NN
		A1=1
		DO I=1,Q
			A1=A1+THIT(I)*THIT(I)
		END DO
		S2A=SC(0)/A1
		DO K=1,Q
			THITA(K)=-SC(K)/S2A
			DO I=1,Q-K
				THITA(K)=THITA(K)+THIT(I)*THIT(K+I)
			END DO
		END DO
		EP1=ABS(S2A-S2B)
		EP2=MAXVAL(ABS(THIT-THITA))
		WRITE(*,*)S2A,EP1,EP2
		IF(EP1<EPS.AND.EP2<EPS)EXIT
		THIT=THITA
		S2B=S2A
	END DO
	END


! 全选主元高斯——约当法(Gauss-Jordan)求解n阶线性代数方程组
	SUBROUTINE GASJDN(A,B,N)
	REAL(8),DIMENSION(N,N)::A
	REAL(8),DIMENSION(N)::B
	REAL(8),DIMENSION(N)::JA
	REAL(8)::DMAX,DD
	LL=1
	DO K=1,N
		DMAX=0
		DO I=K,N
			DO J=K,N
				IF(ABS(A(I,J))>DMAX)THEN
					DMAX=ABS(A(I,J))
					JA(K)=J
					IA=I
				END IF
			END DO
		END DO
		IF(DMAX+1==1)THEN
			WRITE(*,'(" 主元为0,求解失败 ")')
			LL=0
			RETURN
		END IF
		DO J=K,N
			DD=A(K,J)
			A(K,J)=A(IA,J)
			A(IA,J)=DD
		END DO
		DD=B(K)
		B(K)=B(IA)
		B(IA)=DD
		DO I=1,N
			DD=A(I,K)
			A(I,K)=A(I,JA(K))
			A(I,JA(K))=DD
		END DO
		DO J=K+1,N
			A(K,J)=A(K,J)/A(K,K)
		END DO
		B(K)=B(K)/A(K,K)
		DO J=K+1,N
			DO I=1,N
				IF(I/=K)A(I,J)=A(I,J)-A(I,K)*A(K,J)
			END DO
		END DO
		DO I=1,N
			IF(I/=K)THEN
				B(I)=B(I)-A(I,K)*B(K)
			ENDIF
		END DO
	END DO
	DO K=N,1,-1
		DD=B(K)
		B(K)=B(JA(K))
		B(JA(K))=DD
	END DO
	END
	
 	!PROGRAM ARMAPQ
 	!INTEGER,PARAMETER::N=30  
 	!INTEGER,PARAMETER::P=2,Q=1,M=P+Q !P自回归阶数,Q滑动平均阶数
 	!REAL(8),DIMENSION(N)::X
 	!REAL(8),DIMENSION(0:P)::FAI
 	!REAL(8),DIMENSION(Q)::THITA      !滑动系数
 	!REAL(8)::XV                      !X的平均值
 	!REAL(8)::EPS
 	!EPS=1.E-4
 	!OPEN(10,FILE='BEIJING.DAT')
 	!READ(10,*)X
 	!CLOSE(10)
 	!XV=0
 	!DO I=1,N
 	!  XV=XV+X(I)
 	!END DO
 	!XV=XV/N
 	!X=X-XV
 	!CALL ARMA(X,N,P,Q,M,FAI,THITA,EPS)
 	!OPEN(12,FILE='ARMA.DAT')
 	!WRITE(12,'(2X,"XV=",F8.4)')XV
 	!DO I=1,P
 	!  WRITE(12,'("   FAI(",I2,")=",F8.4)')I,FAI(I)
 	!END DO
 	!DO I=1,Q
 	!  WRITE(12,'(" THITA(",I2,")=",F8.4)')I,THITA(I)
 	!END DO
 	!CLOSE(12)
 	!END
