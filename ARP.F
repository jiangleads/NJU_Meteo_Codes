      SUBROUTINE ARP(X,N,M,R,FAI)
	INTEGER::TAO	                !落后时间
	REAL(4),DIMENSION(N)::X
	REAL(4),DIMENSION(M,M)::FAI
	REAL(4),DIMENSION(M)::R
	REAL(4),DIMENSION(M)::S		!协方差
	REAL(4)::S2,A1,A2	            !S2:方差, A1,A2:中间变量
	S=0
	DO TAO=1,M
	  DO I=1,N-TAO
	    S(TAO)=S(TAO)+X(I)*X(I+TAO)
	  END DO
	  S(TAO)=S(TAO)/(N-TAO)
	END DO
	S2=0
	DO I=1,N
	  S2=S2+X(I)*X(I)
	END DO
	S2=S2/N
	DO TAO=1,M
	  R(TAO)=0
	  DO I=1,N-TAO
	    R(TAO)=R(TAO)+X(I)*X(I+TAO)/S2
	  END DO
	  R(TAO)=R(TAO)/(N-TAO)
	END DO
	FAI(1,1)=R(1)
	FAI(2,2)=(R(2)-R(1)*R(1))/(1-R(1)*R(1))
	FAI(1,2)=FAI(1,1)-FAI(2,2)*FAI(1,1)
	DO J=3,M
	  A1=0              
	  A2=0				
	  DO K=1,J-1
	    A1=A1+FAI(K,J-1)*R(J-K)
	    A2=A2+FAI(K,J-1)*R(K)
	  END DO
	  FAI(J,J)=(R(J)-A1)/(1-A2)
	  DO K=1,J-1
	    FAI(K,J)=FAI(K,J-1)-FAI(J,J)*FAI(J-K,J-1)
	  END DO
	END DO
	END

 	PROGRAM MAIN
 	INTEGER,PARAMETER::N=264
 	INTEGER,PARAMETER::M=12
 	REAL(4),DIMENSION(N)::X
 	REAL(4),DIMENSION(M,M)::FAI
 	REAL(4),DIMENSION(M)::R
 	REAL(4)::XV     !X的平均值
 	OPEN(10,FILE='AA2.DAT')
 	DO I=1,N
 	READ(10,'(F8.2)')X(I)
 	END DO
 	CLOSE(10)
 	XV=0
 	DO I=1,N
 	  XV=XV+X(I)
 	END DO
 	XV=XV/N
 	X=X-XV
 	CALL ARP(X,N,M,R,FAI)
 	OPEN(12,FILE='ARP.DAT')
 	WRITE(12,'(2X,"XV=",F8.4)')XV
 	DO I=1,M
	WRITE(12,'("R(",I2,")=",F8.4,"  FAI(",I2,")=",F8.4)')I,R(I),I,FAI(I,M)
 	END DO
 	CLOSE(12)
 	END
