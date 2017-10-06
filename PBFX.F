!	判别分析程序
	! M个因子,序列总长度为N,分为G类,每类的个数为NG(G),N=NG(1)+NG(2)+...+NG(G)
	! 每类的均值为XV(M,G),总的均值为XVV(M),
	SUBROUTINE PBFX(X,N,M,G,KG)
	! M个因子,序列总长度为N,分为G类,每类的个数为NG(G),N=NG(1)+NG(2)+...+NG(G)
	! 每类的均值为XV(M,G),总的均值为XVV(M),
	IMPLICIT NONE
	INTEGER::G,I,J,K,N,M,L,MG
	REAL(8),DIMENSION(N,M)::X
	REAL(8),DIMENSION(M)::XVV
	INTEGER,DIMENSION(G)::NG
	INTEGER,DIMENSION(N)::KG
	REAL(8),DIMENSION(M,G)::XV
	REAL(8),DIMENSION(M,M)::T,S,B
	REAL(8),DIMENSION(M)::ST,DT,XX1
	REAL(8),DIMENSION(M,M)::ST12,ST_12,SS_12,SS,V,S12,S_12,VD,D,VS
	REAL(8),DIMENSION(:),ALLOCATABLE::AA,X2,X22
	REAL(8),DIMENSION(:,:,:),ALLOCATABLE::XG
	REAL(8),DIMENSION(:,:,:),ALLOCATABLE::DS
	INTEGER,DIMENSION(:,:),ALLOCATABLE::KGJ
	INTEGER::L0,NP,NP1,N0,NN
	REAL(8)::AA1,XX2,XX3,DMIN,ALF,H
	NG=0
	DO I=1,N
		NG(KG(I))=NG(KG(I))+1
	END DO
	MG=MAXVAL(NG)	!求数组的最大值
	ALLOCATE(XG(M,G,MG))
	XG=0
	NG=0
	DO I=1,N
		NG(KG(I))=NG(KG(I))+1
		DO J=1,M
			XG(J,KG(I),NG(KG(I)))=X(I,J)
		END DO
	END DO
!   求每个因子总的均值
	XVV=0
	DO K=1,M
		DO I=1,N
			XVV(K)=XVV(K)+X(I,K)
		END DO
		XVV(K)=XVV(K)/N
	END DO
	WRITE(12,'(4X,"总的均值=",<M>F8.2)')XVV
!   求每个因子各组的均值
	XV=0
	DO K=1,M
		DO I=1,G
			DO J=1,NG(I)
				XV(K,I)=XV(K,I)+XG(K,I,J)
			END DO
			XV(K,I)=XV(K,I)/NG(I)
		END DO
	END DO
	WRITE(12,*)
	DO I=1,G
		WRITE(12,'("第",I2,"组的均值=",<M>F8.2)')I,(XV(J,I),J=1,M)
	END DO
	!   求总的离差平方和阵T和组内间离差平方和阵S
	T=0
	DO K=1,M
		DO L=1,M
			DO I=1,G
				DO J=1,NG(I)
					T(K,L)=T(K,L)+(XG(K,I,J)-XVV(K))*(XG(L,I,J)-XVV(L))
					S(K,L)=S(K,L)+(XG(K,I,J)-XV(K,I))*(XG(L,I,J)-XV(L,I))
				END DO
			END DO
		END DO
	END DO
	WRITE(12,*)
	WRITE(12,'("总的离差阵T")')
	WRITE(12,'(<M>F9.2)')T
	WRITE(12,*)
	WRITE(12,'("组内的离差阵S")')
	WRITE(12,'(<M>F9.2)')S
!   求组间离差平方和阵B,可以采用两种方法：1. B=T-S；2. 直接计算法。两种方法计算结果完全一致。
!	第一种方法:B=T-S
	B=T-S
!	第二种方法:直接计算法
!    B=0
!	DO K=1,M
!	  DO L=1,M
!	    DO I=1,G
!		  B(K,L)=B(K,L)+NG(I)*(XV(K,I)-XVV(K))*(XV(L,I)-XVV(L))
!		END DO
!	  END DO
!	END DO
	WRITE(12,*)
	WRITE(12,'("组间离差阵B")')
	WRITE(12,'(<M>F9.2)')B
!   求特征方程(B-λS)C=0的特征值和特征向量
!	第一步,求S的特征值即特征向量(用Jacobi法),S可以分解为S=(VS)(ST)(VS)ˊ
!   返回时ST存放矩阵的全部特征值,VS存放特征向量为列组成的矩阵
	CALL JCB(S,M,1.0E-8,VS,ST,L0)
	ST12=0
	ST_12=0
	DO I=1,M
		ST12(I,I)=SQRT(ST(I))  !ST12为ST的开方
		ST_12(I,I)=1/ST12(I,I)  !ST_12为ST的(-1/2)次方
	END DO
!  求S12=(VS)(ST12)			      ! S12为S的(1/2)次方
!  求S_12=(VS)(ST_12)			  ! S_12为S的(-1/2)次方
	S12=MATMUL(VS,ST12)		  !求矩阵的乘积
	S_12=MATMUL(VS,ST_12)
!	求 D=(S_12)ˊB(S_12)
!	先S_12求的转置SS_12
	DO I=1,M
		DO J=1,M
			SS_12(I,J)=S_12(J,I)
		END DO
	END DO
	SS=MATMUL(SS_12,B)
	D=MATMUL(SS,S_12)
!	求D的特征值及特征向量
	CALL JCB(D,M,1.0E-8,VD,DT,L0)
!	求D的特征值即特征向量(用Jacobi法),D可以分解为D=(VD)(DD)(VD)ˊ,
!   而DD的对角线值放在DT中
!   返回时DT存放矩阵的全部特征值,VD存放特征向量为列组成的矩阵
	V=MATMUL(S_12,VD)      !V和DT为方程(B-λS)C=0的特征向量和特征值
	NP=MIN(G-1,M)
	WRITE(12,*)
	WRITE(12,'(5X,"特征值")')
	WRITE(12,'(<NP>E12.4)')(DT(I),I=1,NP)
	WRITE(12,*)
	WRITE(12,'(5X,"特征向量")')
	WRITE(12,'(<NP>E12.4)')((V(I,J),J=1,NP),I=1,M)
	ALLOCATE(AA(NP))
	ALLOCATE(X2(NP))
	ALLOCATE(X22(NP))
	NP1=0
	DO I=1,NP
		AA1=1
		DO J=I,NP
			AA(J)=AA1/(1+DT(J))
			AA1=AA(J)
		END DO
		X2(I)=-(N-1-0.5*(M+G))*LOG(AA(I))
		N0=(M+1-I)*(G-I)	!检验的自由度
		ALF=0.05   !检验的置信度
		H=-0.005	 !检验时的积分步长
		NN=(1.-ALF)/ABS(H)
		CALL EULER1(H,NN,N0,X22(I))  !哑元依次为积分步长、积分步数、自由度及χ2值
		IF(X2(I)>X22(I))THEN
			NP1=I
		ELSE
			GOTO 30
		END IF
	END DO
30	CONTINUE
	WRITE(12,'("NP1=",I2)')NP1
	WRITE(12,'(5X," χ2检验 ")')
	WRITE(12,'("X2=",<NP1>F8.4)')(X2(I),I=1,NP1)
	WRITE(12,'("X22="<NP1>F8.4)')(X22(I),I=1,NP1)
!   经判别,有NP1个判别函数显著
!	计算每一样品与各组判别函数重心的距离
	ALLOCATE(DS(N,G,NP1))
	DO L=1,NP1
		DO I=1,N
			DO J=1,G
				DO K=1,M
					XX1(K)=X(I,K)-XV(K,J)
				END DO
				XX2=0
				XX3=0
				DO K=1,M
					XX2=XX2+XX1(K)*V(K,L)
					XX3=XX3+V(K,L)*XX1(K)
				END DO
				DS(I,J,L)=XX2*XX3
			END DO
		END DO
	END DO
	ALLOCATE(KGJ(N,NP1))
	DO L=1,NP1
		DO I=1,N
			DMIN=1.0E+30
			DO J=1,G
				IF(DS(I,J,L)<DMIN)THEN
					DMIN=DS(I,J,L)
					KGJ(I,L)=J
				ENDIF
			END DO
		END DO
	END DO
	WRITE(12,'(4X,"样品与判别函数重心的距离")')
	DO L=1,NP1
		DO J=1,G
			WRITE(12,'(" 组别=",I2,<N>F6.3)')J,(DS(I,J,L),I=1,N)
		END DO
		WRITE(12,'("最后分组",<N>I6)')(KGJ(I,L),I=1,N)
	END DO
	DEALLOCATE(AA)
	DEALLOCATE(X2)
	DEALLOCATE(X22)
	DEALLOCATE(DS)
	DEALLOCATE(KGJ)
	END

	SUBROUTINE EULER1(H,N,N0,Y)  !此子程序用于显著性检验
	REAL(8)::Y,D,X,T,H,B
	INTEGER::N,N0
	REAL(8),DIMENSION(N)::Z
	B=2./N0
	T=1
	X=T+H*(J-1)
	IF(N0==1)THEN
		Y=0.000005
	ELSE
		Y=0
	ENDIF
	Z(1)=Y
	DO J=2,N
		X=T+(J-2)*H
		CALL F(X,Y,D,N0)
		Y=Z(J-1)+H*D
		X=T+(J-1)*H
		CALL F(X,Y,D,N0)
		D=Z(J-1)+H*D
		Y=(Y+D)/2.
		Z(J)=Y
	END DO
	DO J=1,N
		IF(N0>2)THEN
			Z(J)=Z(J)**B
		ENDIF
	END DO
	Y=Z(N)
	END

	SUBROUTINE F(X,Y,D,N0)
	REAL(8)::Y,D,X,A,GAMMA,B
	INTEGER::N0
	A=N0/2.
	B=2./N0
	IF(N0<=2)THEN
		D=-2**A*GAMMA(A)*EXP(0.5*Y)*Y**(1-A)
	ELSE
		D=-A*2**A*GAMMA(A)*EXP(0.5*Y**B)
	ENDIF
	X=X
	END

	FUNCTION GAMMA(X)      ! 求Γ函数值
	REAL(8),DIMENSION(0:10)::A
	REAL(8)::S,T,Y,GAMMA,X
	DATA A/6.77106D-5,-3.442342D-4,1.5397681D-3,-2.446748D-3,1.09736958D-2,
     &-2.109075D-4,7.42379071D-2,8.15782188D-2,4.118402518D-1
     &,4.22784337D-1,1.0/
	IF(X<=0)THEN
		PRINT*,'X<0;error!'
		RETURN
	ENDIF
	Y=X
	IF(Y<=1.0)THEN
	T=1/(Y*(Y+1))
	Y=Y+2
	ELSEIF(Y<=2)THEN
		T=1/Y
		Y=Y+1
	ELSEIF(Y<=3)THEN
		T=1
	ELSEIF(Y>3)THEN
		T=1
		DO WHILE(Y>3)
			Y=Y-1
			T=T*Y
		ENDDO
	ELSE
		GAMMA=0
		RETURN
	ENDIF
	S=A(0)
	DO I=1,10
		S=S*(Y-2)+A(I)
		GAMMA=T*S
	ENDDO
	END
	PROGRAM MAIN
 	INTEGER,PARAMETER::N=17
 	INTEGER,PARAMETER::M=4
 	INTEGER,PARAMETER::G=3
 	INTEGER,DIMENSION(N)::KG
 	REAL(8),DIMENSION(N,M)::X
 	OPEN(10,FILE='PB.DAT')
 	DO I=1,N
 	  READ(10,*)I0,(X(I,J),J=1,M),KG(I)
 	END DO
 	CLOSE(10)
 	OPEN(12,FILE='PBFX.DAT')
 	CALL PBFX(X,N,M,G,KG)
 	CLOSE(12)
 	END 

