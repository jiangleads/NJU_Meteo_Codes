	SUBROUTINE CDLRR(B0,P0,W,V,X,Y1,Y,DY,B1,P1
     ,PK,EK,GA,EKIN,XT,GAT,P,PXT,N,M)
C       Y1 = X*B0
C       PK = P0	 
C       EK = X*PK*XT+I
C       GA = PK*XT*EKIN
C       B1 = B0+GA*(Y-Y1)
C       P1 = PK-GA*EK*GAT
C______________________________________________________
	DIMENSION B0(M),P0(M,M),W(M,M),V(N,N)
	DIMENSION X(N,M),Y1(N),Y(N),DY(N)
	DIMENSION B1(M),P1(M,M)
	DIMENSION PK(M,M),EK(N,N),GA(M,N),EKIN(N,N)
	DIMENSION XT(M,N),GAT(N,M),P(M,M),PXT(M,N)
C______________________________________________________
	CALL VLT(PK,P0,M,M)
C ______________________________________________________
	CALL TRN(XT,X,N,M)
	CALL MLT(PXT,PK,XT,M,M,N)
	CALL MLT(EK,X,PXT,N,M,N)
	CALL ADD(EK,V,N,N)
C______________________________________________________
	CALL VLT(EKIN,EK,N,N)
	CALL INVR(EKIN,N)

	CALL MLT(GA, XT,EKIN,M,N,N)
	CALL MLT(GA,PK,GA,M,M,N)
C	CALL MLT(EK,XT,PK,N,M,N)
C______________________________________________________
	CALL VLT(DY,Y,N,1)
	CALL SBT(DY,Y1,N,1)
	CALL MLT(B1,GA,DY,M,N,1)
	CALL ADD(B1,B0,M,1)
C______________________________________________________
	CALL TRN(GAT,GA,M,N)
	CALL MLT(PXT,GA,EK,M,N,N)
	CALL MLT(P,PXT,GAT,M,N,M)
	CALL VLT(P1,PK,M,M)
	CALL SBT(P1,P,M,M)
C______________________________________________________
	RETURN
	END
C______________________________________________________

C     THE PROGRAM **DLRR.FOR** (1999.8.20)
c      -------------------------------------------------------
c     这里的M=M1即因子数加上预报量的总个数
	PARAMETER (N=1,M=5)
	CHARACTER*32 FLNM,BFMT,VFMT,FYFMT,XYWENJIAN
	CHARACTER*2 CM(20)
	DIMENSION B0(M),P0(M,M),W(M,M),V(N,N)
	DIMENSION X(N,M),Y1(N),Y(N),DY(N)
	DIMENSION B1(M),P1(M,M)
	DIMENSION PK(M,M),EK(N,N),GA(M,N),EKIN(N,N)
	DIMENSION XT(M,N),GAT(N,M),P(M,M),PXT(M,N)
	DIMENSION XY(50),BB(M,2000)
	DIMENSION PAF(M,M,M),AF(M ,M)
	INTEGER  NDAY,NUM
cQ-----------------------------
	REAL   B5(M),P5(M,M),B6(M),P6(M,M)
     *,B7(M),P7(M,M),B8(M),P8(M,M),
     *B9(M),P9(M,M)
cQ-----------------------------
	DATA CM/' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10'
     &,'11','12','13','14','15','16','17','18','19','20'/
C
 501	FORMAT(A)
	WRITE(*,*)'请输入存放旧的B、P、W和V的文件名'
	READ(*,501)FLNM
	OPEN (UNIT=1,FILE=FLNM,STATUS='OLD')
C
	WRITE(*,*)'请输入存放X和Y的文件名'
	READ(*,501)FLNM
	OPEN (UNIT=2,FILE=FLNM,STATUS='OLD')
C
	WRITE(*,*)'请输入存放新计算的B和P的文件名'
	READ(*,501)FLNM
	OPEN (UNIT=3,FILE=FLNM)
	WRITE(*,*)'请输入存放预报量估计值的文件名'
	READ(*,501)FLNM
 105	OPEN (UNIT=4,FILE=FLNM)
CC
C       ***********读入   B0,P0,W,V ***************
C
	READ(1,666)(B0(I),I=1,M),((P0(I,J),J=1,M),I=1,M)
     &,((W(I,J),J=1,M),I=1,M)
CC
	READ(1,669)((V(I,J),J=1,N),I=1,N)
 666	FORMAT(5E12.5)
 669	FORMAT(E12.5)
CC       ********************************************
	DO I=1,M
		DO J=1,M
			W(I,J)=0.0
		ENDDO
	ENDDO
	DO I=1,N
		DO J=1,N
			V(I,J)=0.0
			IF(I.EQ.J)V(I,J) =1.0
		ENDDO
	ENDDO

	BFMT(1:4)='(1X,'
	BFMT(5:6)=CM(M)
	BFMT(7:13)='F12.5/)'
	FYFMT(1:4)='(1X,'
	FYFMT(5:6)=CM(M+N)
	FYFMT(7:13)='F12.5/)'
	VFMT(1:7)='(F12.5)'
CC
	NUM=0
	YMSE=0.0
	AMSE=0.0
	NX=M
	M1=NX
C       *********读入因子X和预报量Y *******

 10	NUM=NUM+1
	DO I=1,M
		BB(I,NUM)=B0(I)
	ENDDO
	READ(2,15,END=888)(XY(J),J=1,NX),NDAY
C---------------------------------------------------------
 15	FORMAT(15X,4F11.5,F6.1,3X,I9)
	X(1,1)=1.0
	J=2
	DO 1010 I=1,NX
		X(1,J)=XY(I)
		J=J+1
 1010	CONTINUE
	Y(1)=XY(NX)
CC        ******************************************
	CALL MLT(Y1,X,B0,N,M1,1)
	DY(1)=Y(1)-Y1(1)
	CALL CDLRR(B0,P0,W,V,X,Y1,Y,DY,B5,P5
     *,PK,EK,GA,EKIN,XT,GAT,P,PXT,N,M1)
CC
 642	FORMAT(1X,'F B0 :')
	WRITE(*,642)
	WRITE(*,FYFMT)(X(1,J),J=1,M1),Y(1)
	WRITE(*,BFMT)(B0(I),I=1,M1)
CC
 640	FORMAT(1X,'YO',F12.5,5X,'YF=',F12.5,5X,'DY=',F12.5,3X,I9,2X,I3)
 644	WRITE(*,640)Y(1),Y1(1),DY(1),NDAY,NUM
CC ****************************************************
 641	FORMAT(1X,'B5 P5')
	WRITE(*,641)
	WRITE(*,BFMT)(B5(I),I=1,M1)
	WRITE(*,BFMT)((P5(I,J),J=1,M1),I=1,M1)
	WRITE(3,641)
	WRITE(3,BFMT)(B5(I),I=1,M1)
	WRITE(3,BFMT)((P5(I,J),J=1,M1),I=1,M1)
	YMSE=YMSE+(DY(1)*DY(1))
	AMSE=AMSE+ABS(DY(1))
	CALL VLT(B0,B5,M1,1)
	CALL VLT(P0,P5,M1,M1)
cc---------------------------------------------------
	GOTO 10
 888	CLOSE(2)
ccQ----------------------------------------------------
	WRITE(*,*)'请输入存放XY的文件名'
	READ(*,501)FLNM
	OPEN (UNIT=2,FILE=FLNM,STATUS='OLD')
 5677	DO II=1,M1
		DO JJ=1,M1
			DO I1=1,M1
				P6(I1,JJ)=0.0
				IF(I1.EQ.JJ)P6(I1,JJ)=1.0
			ENDDO
			B6(JJ)=1.0/M1
		ENDDO

		DO KK=M1+1,NUM
			DO JJ=1,M1
				X(1,JJ)=BB(II,KK-JJ)
			ENDDO
			Y(1)=BB(II,KK)
ccq----------------------------------------
			CALL MLT(Y1,X,B6,N,M1,1)
			DY(1)=Y(1)-Y1(1)
ccq----------------------------------------
			CALL CDLRR(B6,P6,W,V,X,Y1,Y,DY,B7,P7
     &		,PK,EK,GA,EKIN,XT,GAT,P,PXT,N,M1)
			DO JJ=1,M1
				AF(II,JJ)=B7(JJ)
			ENDDO
ccQ-----------------------
		DO IO=1,M1
			B6(IO)=B7(IO)
		ENDDO
		DO IO=1,M1
			DO IP=1,M1
				P6(IO,IP)=P7(IO,IP)
			ENDDO
		ENDDO
ccQ------------------------
	ENDDO
	ENDDO
C-----------------------------------------------
C	计算预报值
 224	DO I=1,M
		B8(I)=0.0
		DO J1=1,M
			B8(I) = B8(I)+BB(I,NUM-M+J1)*AF(I,J1)
		ENDDO
	ENDDO


	READ(2,15,END=8889)(XY(J),J=1,NX),NDAY
	X(1,1)=1.0
	J=2
	DO 1013 I=1,NX
		X(1,J)=XY(I)
		J=J+1
 1013	CONTINUE
	Y(1)=XY(NX)
CC  ******************************************
	CALL MLT(Y1,X,B8,N,M1,1)
	DY(1)=Y(1)-Y1(1)
	WRITE(4,640)Y(1),Y1(1),DY(1),NDAY
	CALL CDLRR(B5,P5,W,V,X,Y1,Y,DY,B9,P9
     &,PK,EK,GA,EKIN,XT,GAT,P,PXT,N,M1)

	NUM=NUM+1
	DO I=1,M
		BB(I,NUM)=B9(I)
	ENDDO
	DO I=1,M
		B5(I)=B9(I)
	ENDDO
	DO I=1,M
		DO J=1,M
			P5(I,J)=P9(I,J)
		ENDDO
	ENDDO
	GOTO 5677
 8889	WRITE(*,*)'请输入存放最新的B、P、W和V的文件名'
	READ(*,501)FLNM
	OPEN (UNIT=13,FILE=FLNM,STATUS='NEW')
	WRITE(13,BFMT)(B5(I),I=1,M1)
	WRITE(13,BFMT)((P5(I,J),J=1,M1),I=1,M1)
	WRITE(13,BFMT) ((W(I,J),J=1,M1),I=1,M1)
	WRITE(13,669)((V(I,J),J=1,N),I=1,N)
	CLOSE(1)
	CLOSE(2)
	CLOSE(3)
	CLOSE(4)
	CLOSE(13)
	END
