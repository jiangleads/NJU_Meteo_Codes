	SUBROUTINE JULEI(X,N,M)
	REAL(4),DIMENSION(N,M)::X
	REAL(4),DIMENSION(N,N)::D2 !D2开始存放样品间距离,然后为类与类间距离
	INTEGER,DIMENSION(N)::NP
	REAL(4),DIMENSION(N-1)::DMIN  !存放每次聚类的最小距离
	INTEGER::P,Q
	OPEN(12,FILE='JULEI.DAT')
	WRITE(6,'(<N>F7.1)')((X(I,J),I=1,N),J=1,M)
	DO
		PRINT*," 输入资料的处理方法(提供了3种):   "
		PRINT*," ISB=1,数据不作变换			  "
		PRINT*," ISB=2,将数据中心化,即减去其均值"
		PRINT*," ISB=3,对数据标准化"
		READ(*,'(I1)')ISB
		IF(ISB<=3.AND.ISB>=1)EXIT
		PRINT*,"输入资料的处理方法错误,请重新输入ISB"
	END DO
	CALL ZLCL(X,N,M,ISB)  ! 调用资料处理子程序
	DO
		PRINT*," 输入计算样品间距离的方法(提供了5种):"
		PRINT*,"JL=1,取绝对值距离               "
		PRINT*,"JL=2,取欧式距离				 "
		PRINT*,"JL=3,取切比雪夫距离			 "
		PRINT*,"JL=4,取马氏距离				 "
		PRINT*,"JL=5,取兰氏距离                 "
		READ(*,'(I1)')JL
		IF(JL<=5.AND.JL>=1)EXIT
		PRINT*,"输入的JL错误,请重新输入JL"
	END DO
	CALL JULI(X,N,M,D2,JL)   !调用求距离子程序
	DO
		PRINT*, " 输入计算类间距离的方法(提供了6种)  "
		PRINT*,"LJJ=1,最短距离法					 "
		PRINT*,"LJJ=2,最长距离法					 "
		PRINT*,"LJJ=3,中间距离法					 "
		PRINT*,"LJJ=4,重心法(要求JL=2)			 "
		PRINT*,"LJJ=5,类平均法						 "
		PRINT*,"LJJ=6,离差平方和法(要求JL=2)		 "
		READ(*,'(I1)')LJJ
		IF(LJJ<=6.AND.LJJ>=1)EXIT
		PRINT*,"输入的LJJ错误,请重新输入LJJ"
	END DO
	WRITE(12,'(5X,"ISB=",I2,"; JL=",I2,"; LJJ=",I2)')ISB,JL,LJJ
	WRITE(12,'("样本序号 ",<N>I3)')(I,I=1,N)
	NP=1
	DO II=1,N-1
!	  WRITE(12,'("距离阵",2X,<N>E11.3)')((D2(I,J),I=1,N),J=1,N)
			DMIN(II)=1.0E30
			DO I=2,N					!选最小距离
				DO J=1,I-1
					IF(D2(I,J)>=0)THEN
						IF(D2(I,J)<DMIN(II))THEN
							DMIN(II)=D2(I,J)
							P=I
							Q=J
						END IF
					END IF
				END DO
			END DO
			WRITE(6,'(3I4,E12.4)')II,P,Q,DMIN(II)
			NR=NP(P)+NP(Q)
			DO I=1,N
				IF(NP(I)==0.OR.I==P.OR.I==Q) GOTO 100
				CALL XISHU(AP,AQ,BT,GM,N,NP,P,Q,I,NR,LJJ)  !调用求类间距离的公式系数
				DD=AP*D2(P,I)+AQ*D2(Q,I)+BT*D2(P,Q)+GM*ABS(D2(P,I)-D2(Q,I))
				D2(P,I)=DD
				D2(I,P)=DD
			END DO
100			D2(P,P)=0
			DO I=1,N
				D2(I,Q)=-999
				D2(Q,I)=-999
			END DO
			NP(P)=NR
			NP(Q)=0
			WRITE(12,'("第",I3," 步 ",<N>I3)')II,NP
		END DO
		CLOSE(12)
	END

	SUBROUTINE ZLCL(X,N,M,ISB)
	REAL(4),DIMENSION(N,M)::X
	REAL(4),DIMENSION(M)::XV
	REAL(4),DIMENSION(M)::S
	IF(ISB==1)RETURN                !ISB=1,资料不作处理
	IF(ISB==2)THEN					!ISB=2,中心化处理
		DO J=1,M
			XV(J)=0
			DO I=1,N
				XV(J)=XV(J)+X(I,J)
			END DO
			XV(J)=XV(J)/N
		END DO
		DO I=1,N
			DO J=1,M
				X(I,J)=X(I,J)-XV(J)
			END DO
		END DO
	ELSEIF(ISB==3)THEN				 !ISB=3,标准化处理
		DO J=1,M
			S(J)=0
			DO I=1,N
				S(J)=S(J)+X(I,J)*X(I,J)
			END DO
			X(I,J)=X(I,J)/S(J)
		END DO
	END IF
	END

	SUBROUTINE XISHU(AP,AQ,BT,GM,N,NP,P,Q,I,NR,LJJ)
	REAL(4)::AP,AQ,BT,GM
	INTEGER,DIMENSION(N)::NP
	INTEGER::NR,LJJ,P,Q,I
	SELECT CASE(LJJ)
	CASE (1)	   !最短距离法
		AP=0.5
		AQ=0.5
		BT=0
		GM=-0.5
	CASE (2)	   !最长距离法
		AP=0.5
		AQ=0.5
		BT=0
		GM=0.5
	CASE (3)      !中间距离法
		AP=0.5
		AQ=0.5
		BT=-0.25
		GM=0
	CASE (4)	   !重心法
		AB=1.0/NR
		AP=NP(P)*AB
		AQ=NP(Q)*AB
		BT=-AP*AQ
		GM=0
	CASE (5)	   !类平均法
		AB=1.0/NR
		AP=NP(P)*AB
		AQ=NP(Q)*AB
		BT=0
		GM=0
	CASE (6)	   !离差平方和法
		AB=1.0/(NP(I)+NR)
		AP=(NP(I)+NP(P))*AB
		AQ=(NP(I)+NP(Q))*AB
		BT=-NP(I)*AB
		GM=0
	END SELECT
	END

	SUBROUTINE JULI(X,N,M,D,JL)   !注意:这里采用的是距离的平方
	REAL(4),DIMENSION(N,M)::X
	REAL(4),DIMENSION(M,N)::X1
	REAL(4),DIMENSION(N,M)::X2
	REAL(4),DIMENSION(N,N)::D,S,S1
	SELECT CASE(JL)
	CASE (1)  !绝对值距离
		DO I=1,N
			DO J=1,I
				D(I,J)=0
				IF(I/=J)THEN
					DO K=1,M
						D(I,J)=D(I,J)+ABS(X(I,K)-X(J,K))
					END DO
					D(I,J)=D(I,J)*D(I,J)
					D(J,I)=D(I,J)
				END IF
			END DO
		END DO
	CASE(2)	  !欧式距离
		DO I=1,N
			DO J=1,I
				D(I,J)=0
				IF(I/=J)THEN
					DO K=1,M
						D(I,J)=D(I,J)+(X(I,K)-X(J,K))**2
					END DO
					D(J,I)=D(I,J)
				END IF
			END DO
		END DO
	CASE(3)   !切比雪夫距离
		DO I=1,N
			DO J=1,I
				IF(I==J)THEN
					D(I,J)=0
				ELSE
					DMIN=1.0E30
					DO K=1,M
						F=ABS(X(I,K)-X(J,K))
						IF(F<DMIN)THEN
							DMIN=F
						END IF
						D(I,J)=DMIN
					END DO
				END IF
			END DO
		END DO
	CASE(4)	  !Mahalanobis距离
		DO I=1,N         !求协方差阵
			DO J=1,I
				S(I,J)=0
				DO K=1,M
					S(I,J)=S(I,J)+(X(I,K)-X(J,K))**2
				END DO
				S(I,J)=S(I,J)/N
				S(J,I)=S(I,J)
			END DO
		END DO
		S1=S
		CALL NIHS(S1,N)		 !求协方差阵的逆S1
		DO J1=1,N
			DO J2=1,N
				DO K=1,M
					DO I=1,N
						X1(K,I)=X(J1,K)-X(J2,K)
						X2(I,K)=X(J1,K)-X(J2,K)
					END DO
				END DO
			END DO
		END DO
		D=MATMUL(X1,S1)
		D=MATMUL(D,X2)  !注意:这里的D实际是D的平方
	CASE(5)	  !求兰氏距离
		DO I=1,N
			DO J=1,I
				D(I,J)=0
				IF(I/=J)THEN
					DO K=1,M
						D(I,J)=D(I,J)+ABS(X(I,K)-X(J,K))/(X(I,K)+X(J,K))
					END DO
					D(I,J)=D(I,J)*D(I,J)
					D(J,I)=D(I,J)
				END IF
			END DO
		END DO
	END SELECT
	END

	SUBROUTINE NIHS(A,N)  !求矩阵的逆子程序
	REAL(4),DIMENSION(N,N)::A
	INTEGER,DIMENSION(N)::IS,JS
	REAL(4)::B,AMAX
	DO K=1,N
		AMAX=0
		DO I=K,N
			DO J=K,N
				IF(ABS(A(I,J))>AMAX)THEN
					AMAX=ABS(A(I,J))
					IS(K)=I
					JS(K)=J
				END IF
			END DO
		END DO
		IF(AMAX+1.0==1.0)THEN
			WRITE(*,'("矩阵奇异,无逆阵")')
			STOP
		ENDIF
		DO J=1,N
			B=A(K,J)
			A(K,J)=A(IS(K),J)
			A(IS(K),J)=B
		END DO
		DO I=1,N
			B=A(I,K)
			A(I,K)=A(I,JS(K))
			A(I,JS(K))=B
		END DO
		A(K,K)=1/A(K,K)
		DO J=1,N
			IF(J/=K)THEN
				A(K,J)=A(K,J)*A(K,K)
			END IF
		END DO
		DO I=1,N
			IF(I/=K)THEN
				DO J=1,N
					IF(J/=K)THEN
						A(I,J)=A(I,J)-A(I,K)*A(K,J)
					END IF
				END DO
			END IF
		END DO
		DO I=1,N
			IF(I/=K)THEN
				A(I,K)=-A(I,K)*A(K,K)
			END IF
		END DO
	END DO
	DO K=N,1,-1
		DO J=1,N
			B=A(K,J)
			A(K,J)=A(JS(K),J)
			A(JS(K),J)=B
		END DO
		DO I=1,N
			B=A(I,K)
			A(I,K)=A(I,IS(K))
			A(I,IS(K))=B
		END DO
	END DO
	END
 
 	PROGRAM MAIN
 	INTEGER,PARAMETER::N=20,M=3	!N:样品数;M:指标数
 	REAL(4),DIMENSION(N,M)::X
 	OPEN(10,FILE='BEIJING.DAT')
 	READ(10,*)(X(I,1),I=1,N)
 	READ(10,*)(X(I,2),I=1,N)
 	READ(10,*)(X(I,3),I=1,N)
 	CLOSE(10)
 	CALL JULEI(X,N,M)	
 	END
