	SUBROUTINE DIVER(U,V,P,DL,DM,F0,L,M,N,DIV)
!   U、V分别为风速的东西和南北分量；P为气压值；DX、DY分别为X、Y方向的格距(单位为弧度)；
!   F0为J=1处的纬度值(弧度)；L、M分别为X、Y方向的格点数；N为垂直方向的层数；DIV为散度。
	REAL(8),DIMENSION(L,M,N)::U,V
	REAL(8),DIMENSION(L,M,N)::DIV
	REAL(8),DIMENSION(L,M)::DV,DVV,ESP
	REAL(8),DIMENSION(N)::P
	REAL(8)::DL,DM,F0
	REAL(8),PARAMETER::AA=6.371E6	 !地球半径
!   计算未订正的散度
	L1=L-1
	M1=M-1
	DO K=1,N
		DO J=2,M1
			DO I=2,L1
				DIV(I,J,K)=(U(I+1,J,K)-U(I-1,J,K))/(AA*COS(DM*(J-1)+F0)*2*DL)
     &			+(V(I,J+1,K)-V(I,J-1,K))/(AA*2*DM)
     &			-V(I,J,K)*TAN(DM*(J-1)+F0)/AA
			END DO
		END DO
	END DO
	CALL BOUND(DIV,L,M,N)
!   进行散度订正
	DO J=1,M
		DO I=1,L
			DV(I,J)=0.0
			DVV(I,J)=0.0
			DO K=1,N-1
				DV(I,J)=DV(I,J)+(DIV(I,J,K)+DIV(I,J,K+1))*(P(K+1)-P(K))/2.
		DVV(I,J)=DVV(I,J)+(ABS(DIV(I,J,K))+ABS(DIV(I,J,K+1)))*(P(K+1)-P(K))/2.
			END DO
		END DO
	END DO
	DO J=1,M
		DO I=1,L
			ESP(I,J)=-DV(I,J)/DVV(I,J)
			DO K=1,N
				DIV(I,J,K)=DIV(I,J,K)+ESP(I,J)*ABS(DIV(I,J,K))
			END DO
		END DO
	END DO
	RETURN
	END

	SUBROUTINE BOUND(A,L,M,N)
	REAL(8),DIMENSION(L,M,N)::A
	L1=L-1
	M1=M-1
	DO K=1,N
		DO I=2,L1
			A(I,1,K)=2*A(I,2,K)-A(I,3,K)
			A(I,M,K)=2*A(I,M-1,K)-A(I,M-2,K)
		END DO
		DO J=2,M1
			A(1,J,K)=2*A(2,J,K)-A(3,J,K)
			A(L,J,K)=2*A(L-1,J,K)-A(L-2,J,K)
		END DO
		A(1,1,K)=A(1,2,K)+A(2,1,K)-(A(1,3,K)+A(3,1,K))*0.5
		A(L,1,K)=A(L,2,K)+A(L-1,1,K)-(A(L,3,K)+A(L-2,1,K))*0.5
		A(1,M,K)=A(1,M-1,K)+A(2,M,K)-(A(1,M-2,K)+A(3,M,K))*0.5
		A(L,M,K)=A(L,M-1,K)+A(L-1,M,K)-(A(L,M-2,K)+A(L-2,M,K))*0.5
	END DO
	END

	SUBROUTINE VORTICITY(U,V,DL,DM,F0,L,M,N,VOR)
	REAL(8),DIMENSION(L,M,N)::U,V
	REAL(8),DIMENSION(L,M,N)::VOR
!   U、V分别为风速的东西和南北分量；DX、DY分别为X、Y方向的格距(单位为弧度)；
!   F0为J=1处的纬度值(弧度)；L、M分别为X、Y方向的格点数；N为垂直方向的层数；
!   VOR为涡度。
	REAL(8),PARAMETER::AA=6.371E6	  !地球半径
	REAL(8)::DL,DM,F0
	L1=L-1
	M1=M-1
	DO K=1,N
		DO J=2,M1
			DO I=2,L1
				VOR(I,J,K)=(V(I+1,J,K)-V(I-1,J,K))/(AA*COS(DM*(J-1)+F0)*2*DL)
     &			-(U(I,J+1,K)-U(I,J-1,K))/(AA*2*DM)
     &			+U(I,J,K)*TAN(DM*(J-1)+F0)/AA
			END DO
		END DO
	END DO
	CALL BOUND(VOR,L,M,N)
	RETURN
	END
