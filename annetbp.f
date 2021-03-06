	SUBROUTINE annetbp(LR,A,B,BT1,BT2)
	PARAMETER (NUM=3005,NFA=31,NF=256,NUMT=546,NUMTF=60)
	DIMENSION X(NUM,NFA),W1(40,20),W2(20),YY(NUM),ess(20),
     *HMAX(NFA),HMIN(NFA),POSITION(NFA),EY(20),
     *X0(NUM,NFA),CT(NUMT,NFA),YYT(NUMT)
	DIMENSION A(NFA,2749),B(NFA,NF),BT1(NFA,NUMT-NUMTF),
     *BT2(NFA,NUMTF)
	INTEGER S,SIGN,SIGN1,POSITION,LR
	F(z)=1.0/(1+EXP(-z))
	F1(Z)=Z
C	 计算网络连接权矢量
	DO J=1,20
		DO I=1,40
			CALL RANDOM_NUMBER(W1(I,J))
		ENDDO
		CALL RANDOM_NUMBER(W2(J))
	ENDDO

	DO 209 II=1,NUMT-NUMTF
		DO 109 JJ=1,NFA
			CT(II,JJ)=BT1(JJ,II)
			YYT(II)=0
109		CONTINUE
209	CONTINUE
	do 219 II=NUMT-NUMTF,NUMT
		DO 119 JJ=1,NFA
			CT(II,JJ)=BT2(JJ,II-NUMT-NUMTF)
			YYT(II)=1
119		CONTINUE
219	CONTINUE

	print*,'进入BP神经网络训练'
	DO 11 J=1,NFA
		DO 12 I=1,NF
			X(11*I-10,J)=A(J,I)
			X(11*I-9,J)=A(J,I)
			X(11*I-8,J)=A(J,I)
			X(11*I-7,J)=A(J,I)
			X(11*I-6,J)=A(J,I)
			X(11*I-5,J)=A(J,I)
			X(11*I-4,J)=A(J,I)
			X(11*I-3,J)=A(J,I)
			X(11*I-2,J)=A(J,I)
			X(11*I-1,J)=A(J,I)
			X(11*I,J)=B(J,I)
			YY(11*I-10)=0
			YY(11*I-9)=0
			YY(11*I-8)=0
			YY(11*I-7)=0
			YY(11*I-6)=0
			YY(11*I-5)=0
			YY(11*I-4)=0
			YY(11*I-3)=0
			YY(11*I-2)=0
			YY(11*I-1)=0
			YY(11*I)=1
12		CONTINUE
		do 15 i=2817,NUM
			x(i,j)=a(j,i)
			yy(i)=0.0
15		continue
11	CONTINUE

	DO 20 J=1,NFA
		POSITION(J)=J
20	CONTINUE
	DO 55 J=1,NFA
		CALL QQ(X,NUM,J,HMAX,HMIN)

		DO 65 I=1,NUM
			X(I,J)=(X(I,J)-HMIN(J))/(HMAX(J)-HMIN(J))
65		CONTINUE
55	CONTINUE


	do 23 i=1,NUM
		do 24 j=1,NFA
			x0(i,j)=x(i,j)
24		continue
23	continue

	DO 50 II=6,15
		NYZ=II
		S=1
		DO 103 I=1,3
			CALL	INTRU(X,W1,W2,YY,S,SIGN,POSITION,nyz,NUM,NFA,LR)
			IF(SIGN.EQ.1)THEN
				S=S+1
				OPEN(UNIT=101,FILE='POSITION.DAT')
				WRITE(UNIT=101,FMT=FORM2)(POSITION(J),J=1,S)
				CLOSE(101)
			ENDIF
			IF(SIGN.EQ.0)THEN
				GOTO	1110
			ENDIF
103		CONTINUE
		DO 110 I=1,NFA
			CALL DISCARD(X,W1,W2,YY,S,SIGN1,POSITION,NYZ,NUM,NFA,LR)
			IF(SIGN1.EQ.1)THEN
				S=S
			ENDIF
			IF(SIGN1.EQ.0)THEN
				S=S-1
				OPEN(UNIT=101,FILE='POSITION.DAT')
				WRITE(UNIT=101,FMT=FORM2)(POSITION(J),J=1,S)
				CLOSE(101)
			ENDIF
			CALL	INTRU(X,W1,W2,YY,S,SIGN,POSITION,NYZ,NUM,NFA,LR)
			IF(SIGN.EQ.1)THEN
				S=S+1
				OPEN(UNIT=101,FILE='POSITION.DAT')
				WRITE(UNIT=101,FMT=FORM2)(POSITION(J),J=1,S)
				CLOSE(101)
			ENDIF
			IF(SIGN.EQ.0.AND.SIGN1.EQ.1)THEN
				CALL STUDY(X,CT,w1,w2,YY,S,NYZ,NUM,NFA,NUMT,NUMTF,LR)
				GOTO	1110
			ENDIF
110		CONTINUE

1110		CALL TXYZ(X,W1,W2,YY,S,NYZ,ESS,NUM,LR)
		WRITE(*,'(A,2I5)')' S=,NYZ= ',S,NYZ

		EY(II)=ESS(nyz)

50	CONTINUE
	CALL MINWC(EY,EYM,6,15,NYZ0)
	print*,'nyz0=',nyz0

	DO 22 J=1,NFA
		POSITION(J)=J
22	CONTINUE

c	print*, x0(1,16),x0(1,8),x0(1,7),x0(1,6)
	S=3
	DO 121 I=1,28
		CALL DISCARD(X0,W1,W2,YY,S,SIGN1,POSITION,NYZ0,NUM,NFA,LR)
		IF(SIGN1.EQ.1)THEN
			S=S
		ENDIF
		IF(SIGN1.EQ.0)THEN
			S=S-1
			OPEN(UNIT=101,FILE='POSITION.DAT')
			WRITE(UNIT=101,FMT=FORM2)(POSITION(J),J=1,S)
			CLOSE(101)
		ENDIF
		CALL	INTRU(X0,W1,W2,YY,S,SIGN,POSITION,NYZ0,NUM,NFA,LR)
		IF(SIGN.EQ.1)THEN
			S=S+1
			OPEN(UNIT=101,FILE='POSITION.DAT')
			WRITE(UNIT=101,FMT=FORM2)(POSITION(J),J=1,S)
			CLOSE(101)
		ENDIF
		IF(SIGN.EQ.0.AND.SIGN1.EQ.1)THEN
c		print*, x0(1,1),x0(1,16),x0(1,12),'begin1'
				CALL STUDY(X,CT,w1,w2,YY,S,NYZ0,NUM,NFA,NUMT,NUMTF,LR)
				GOTO	1000
			ENDIF
 121	CONTINUE
c		print*, x0(1,1),x0(1,2),x0(1,3),x0(1,4)
c	pause
	CALL STUDY(X,CT,w1,w2,YY,S,NYZ0,NUM,NFA,NUMT,NUMTF,LR)

	OPEN(UNIT=50,FILE='WW1.DAT')
	WRITE(UNIT=50,FMT='(40F10.6)')((W1(I,J),I=1,S+1),J=1,nyz0+1)
	CLOSE(50)
	OPEN(UNIT=60,FILE='WW2.DAT')
	WRITE(UNIT=60,FMT='(20F10.6)')(W2(J),J=1,S+1)
	CLOSE(60)
 115	CONTINUE
1000	OPEN(UNIT=50,FILE='WW1.DAT')
	WRITE(UNIT=50,FMT='(40F10.6)')((W1(I,J),I=1,S+1),J=1,nyz0+1)
	CLOSE(50)
	OPEN(UNIT=60,FILE='WW2.DAT')
	WRITE(UNIT=60,FMT='(20F10.6)')(W2(J),J=1,S+1)
	CLOSE(60)
	CLOSE(101)
	write(*,111)((w1(i,j),i=1,S+1),j=1,NYZ0+1)
	write(*,111)(w2(i),J=1,NYZ0+1)
111	format(20F10.5)

	RETURN
	END

	!在本程序中还应用了如下子程序：
	!SUBROUTINE INTRU() 引入因子子程序
	!SUBROUTINE DISCARD() 剔除因子子程序
	!SUBROUTINE STUDY () 学习子程序
	!SUBROUTINE TXYZ() 隐层结点数选择子程序
	!SUBROUTINE MINWC() 求极小值子程序
	!SUBROUTINE QQ() 求极大小值子程序
	!SUBROUTINE TEST() 网络测试子程序
!	PROGRAM Artificial Neural Network LP
!	!NUM=3005为总样本数,2749个晴天样本,256个雾天样本
!	!NFA=31为因子数(这是芜湖机场雾预报的例子)
!	PARAMETER (NUM=3005,NFA=31,NF=256,NUMT=546,NUMTF=60)
!	DIMENSION ,A(NFA,2749),B(NFA,NF),BT1(NFA,NUMT-NUMTF),
!     *BT2(NFA,NUMTF)
!	INTEGER LR
!	CHARACTER FORM1*9,FORM2*6
!	F(z)=1.0/(1+EXP(-z))
!	F1(Z)=Z
!
!	LR=1 !线性(LR=0为非线性)
!	NC=NUM-NF
!	FORM1='(31F10.1)'
!	FORM2='(31I4)'
!	WRITE(FORM1(2:3),'(I2.2)')NFA
!	WRITE(FORM2(2:3),'(I2.2)')NFA
!C 	读网络训练样本
!	!clear Days.
!	OPEN(UNIT=10,FILE='c24dat.dat',STATUS='OLD')
!	READ(UNIT=10,FMT=FORM1)((a(I,J),i=1,NFA),j=1,NC)
!	CLOSE(10)
!	!Foggy !days
!	OPEN(UNIT=20,FILE='f24dat.dat',STATUS='OLD')
!	READ(UNIT=20,FMT=FORM1)((b(I,J),i=1,NFA),j=1,NF)
!	CLOSE(10)
!C	读检验样本
!	OPEN(UNIT=108,FILE='cj24DAT.DAT',STATUS='OLD')
!	READ(UNIT=108,FMT=FORM1)((BT1(II,JJ),II=1,NFA),JJ=1,NUMT-NUMTF)
!	CLOSE(108)
!
!	OPEN(UNIT=208,FILE='fj24DAT.DAT',STATUS='OLD')
!	READ(UNIT=208,FMT=FORM1)((BT2(II,JJ),II=1,NFA),JJ=1,NUMTF)
!	CLOSE(208)
!
!	call annetbp(LR,A,B,BT1,BT2)
!	stop
!	end

! ! 误差反向传播学习算法就是一种有效的多层网络学习算法,通常把它称为Back Propagation,简写为BP算法,也叫做多层感知器算法。
	!也是一种建立信息与信息之间联系的数学模型的方法,和单层线性感知器网络不同的是它是多层非线性网络自学习,因此,它是一种用途更加广泛的建模方法。
! ! 
! 
!	PROGRAM Artificial Neural Network BP
!	!NUM=3005为总样本数,2749个晴天样本,256个雾天样本
!	!NFA=31为因子数(这是芜湖机场雾预报的例子)
!	PARAMETER (NUM=3005,NFA=31,NF=256,NUMT=546,NUMTF=60)
!	DIMENSION ,A(NFA,2749),B(NFA,NF),BT1(NFA,NUMT-NUMTF),
!	* 	       BT2(NFA,NUMTF)
!	INTEGER LR
!	CHARACTER FORM1*9,FORM2*6
!	F(z)=1.0/(1+EXP(-z))
!	F1(Z)=Z
!
!	LR=0 !非线性(LR=1为线性)
!	NC=NUM-NF
!	FORM1='(31F10.1)'
!	FORM2='(31I4)'
!	WRITE(FORM1(2:3),'(I2.2)')NFA
!	WRITE(FORM2(2:3),'(I2.2)')NFA
!C 	读网络训练样本
!	!clear Days.
!	OPEN(UNIT=10,FILE='c24dat.dat',STATUS='OLD')
!	READ(UNIT=10,FMT=FORM1)((a(I,J),i=1,NFA),j=1,NC)
!	CLOSE(10)
!	!Foggy days
!	OPEN(UNIT=20,FILE='f24dat.dat',STATUS='OLD')
!	READ(UNIT=20,FMT=FORM1)((b(I,J),i=1,NFA),j=1,NF)
!	CLOSE(10)
!C	读检验样本
!	OPEN(UNIT=108,FILE='cj24DAT.DAT',STATUS='OLD')
!	READ(UNIT=108,FMT=FORM1)((BT1(II,JJ),II=1,NFA),JJ=1,NUMT-NUMTF)
!	CLOSE(108)
!
!	OPEN(UNIT=208,FILE='fj24DAT.DAT',STATUS='OLD')
!	READ(UNIT=208,FMT=FORM1)((BT2(II,JJ),II=1,NFA),JJ=1,NUMTF)
!	CLOSE(208)
!
!	call annetbp(LR,A,B,BT1,BT2)
!	stop
!	end

