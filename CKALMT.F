	SUBROUTINE CKALMT(B0,P0,W,V,X,Y1,Y,DY,B1,P1
     ,PK,EK,GA,EKIN,XT,GAT,P,PXT,N,M)
C       PK  = P0+W	 P0=C0,PK=R,X=F,P1=C1,A=GA
C       EK = X*PK*XT+V
C       GA  = PK*XT*EKIN
C       B1 = B0+GA*(Y-Y1)
C       P1 = PK-GA*EK*GAT
C       or P1 = PK-GA*X*PK
C________________________________________________
	DIMENSION B0(M),P0(M,M),W(M,M),V(N,N)
	DIMENSION X(N,M),Y1(N),Y(N),DY(N)
	DIMENSION B1(M),P1(M,M)
	DIMENSION PK(M,M),EK(N,N),GA(M,N),EKIN(N,N)
	DIMENSION XT(M,N),GAT(N,M),P(M,M),PXT(M,N)
C___________________________________________________
	CALL VLT(PK,P0,M,M)
	CALL ADD(PK,W,M,M)
C_________________________________________________
	CALL TRN(XT,X,N,M)
	CALL MLT(PXT,PK,XT,M,M,N)
	CALL MLT(EK,X,PXT,N,M,N)
	CALL ADD(EK,V,N,N)
C__________________________________________________
	CALL VLT(EKIN,EK,N,N)
	CALL INVR(EKIN,N)
	CALL MLT(GA,PXT,EKIN,M,N,N)
C__________________________________________________
	CALL VLT(DY,Y,N,1)
	CALL SBT(DY,Y1,N,1)
	CALL MLT(B1,GA,DY,M,N,1)
	CALL ADD(B1,B0,M,1)
C__________________________________________________
	CALL TRN(GAT,GA,M,N)
	CALL MLT(PXT,GA,EK,M,N,N)
	CALL MLT(P,PXT,GAT,M,N,M)
	CALL VLT(P1,PK,M,M)
	CALL SBT(P1,P,M,M)
C__________________________________________________	
	RETURN
	END
C__________________________________________________

C       THE PROGRAM **KALMA.FOR** (1999.8.20)
	PARAMETER (N=1,M=5)
	CHARACTER*32 FLNM,BFMT,VFMT,FYFMT
	CHARACTER*2 CM(20)
	DIMENSION B0(M),P0(M,M),W(M,M),V(N,N)
	DIMENSION X(N,M),Y1(N),Y(N),DY(N)
	DIMENSION B1(M),P1(M,M)
	DIMENSION PK(M,M),EK(N,N),GA(M,N),EKIN(N,N)
	DIMENSION XT(M,N),GAT(N,M),P(M,M),PXT(M,N)
	DIMENSION XY(50)
	INTEGER*4  NDAY
C       -------------------------------------------
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
	OPEN (UNIT=3,FILE=FLNM,STATUS='NEW')
	WRITE(*,*)'请输入存放预报量估计值的文件名'
	READ(*,501)FLNM
 105	OPEN (UNIT=4,FILE=FLNM,STATUS='NEW')
CC
C       ***********读入   B0,P0,W,V ***************
C
	READ(1,666)(B0(I),I=1,M),((P0(I,J),J=1,M),I=1,M)
     &  ,((W(I,J),J=1,M),I=1,M)
CC
	READ(1,669)((V(I,J),J=1,N),I=1,N)
 666	FORMAT(5E12.5)
 669	FORMAT(E12.5)
CC       ********************************************
	BFMT(1:4)='(1X,'
	BFMT(5:6)=CM(M)
	BFMT(7:13)='F12.5/)'
	FYFMT(1:4)='(1X,'
	FYFMT(5:6)=CM(M+N)
	FYFMT(7:13)='F12.5/)'
	VFMT(1:7)='(F12.5)'
	WRITE(*,*)'请输入因子和预报量的总数'
	READ(*,'(I2)')NX
CC
	NUM=0
	YMSE=0.0
	AMSE=0.0
	M1=NX
C       *********读入因子X和预报量Y *******
 10	NUM=NUM+1
	READ(2,15,END=888)(XY(J),J=1,NX),NDAY
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
	CALL CKALMT(B0,P0,W,V,X,Y1,Y,DY,B1,P1
     &,PK,EK,GA,EKIN,XT,GAT,P,PXT,N,M1)
CC
 642	FORMAT(1X,'F B0 :')
	WRITE(*,642)
	WRITE(*,FYFMT)((X(I,J),I=1,N),J=1,M1),(Y(I),I=1,N)
	WRITE(*,BFMT)(B0(I),I=1,M1)
CC
 640	FORMAT(/1X,'YO',F12.5,5X,'YF=',F12.5,5X,'DY=',F12.5,3X,I9,2X,I4)
 644	WRITE(*,640)(Y(I),Y1(I),DY(I),I=1,N),NDAY,NUM
	WRITE(4,640)(Y(I),Y1(I),DY(I),I=1,N),NDAY,NUM
CC      **********************************************************
 641	FORMAT(1X,'B1 P1 W和V:')
	WRITE(*,641)
	WRITE(*,BFMT)(B1(I),I=1,M1)
	WRITE(*,BFMT)((P1(I,J),J=1,M1),I=1,M1)
	WRITE(3,641)
	WRITE(3,BFMT)(B1(I),I=1,M1)
	WRITE(3,BFMT)((P1(I,J),J=1,M1),I=1,M1)
C	WRITE(3,BFMT) ((W(I,J),J=1,M1),I=1,M1)
C	WRITE(3,669)((V(I,J),J=1,N),I=1,N)
CC
	YMSE=YMSE+(DY(1)*DY(1))
	AMSE=AMSE+ABS(DY(1))
CC
	CALL VLT(B0,B1,M1,1)
	CALL VLT(P0,P1,M1,M1)
C	PAUSE
	GO TO 10
 888	YMSE=SQRT(YMSE/NUM)
	AMSE=AMSE/NUM
	WRITE(*,770)YMSE
	WRITE(4,770)YMSE
	WRITE(*,776)AMSE
	WRITE(4,776)AMSE
 770	FORMAT(1X,7H均方差:,F6.2)
 776	FORMAT(1X,9H绝对误差:,F6.2)
	WRITE(*,*)'请输入存放最新的B、P、W和V的文件名'
	READ(*,501)FLNM
	OPEN (UNIT=13,FILE=FLNM,STATUS='NEW')
	WRITE(13,BFMT)(B1(I),I=1,M1)
	WRITE(13,BFMT)((P1(I,J),J=1,M1),I=1,M1)
	WRITE(13,BFMT) ((W(I,J),J=1,M1),I=1,M1)
	WRITE(13,669)((V(I,J),J=1,N),I=1,N)
	CLOSE(1)
	CLOSE(2)
	CLOSE(3)
	CLOSE(4)
	CLOSE(13)
	STOP
	END
