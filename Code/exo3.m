clear all
close all
clc

Ke = 3.6/1000*60/(2*pi);
Ks=10;
Kg=0.105;
Te = 0.01;
Km=10;
Tm=0.3;
Kc=3.5/100;


A=[0 Ks/(9*Kg);0 (-1/Tm)]
B=[0; (Km*Kg)/Tm]
C=[1 0]
D=[0]

sys=ss(A,B,C,D)
Vp=eig(A)

Co=ctrb(sys);
rang_co=rank(Co);

Obs=obsv(sys);
rang_obs=rank(Obs);

K=acker(A,B,[-2.4+5.5*i -2.4-5.5*i])
% K1=1.3890;
% K2=0.5986;

N=1/(C*inv(-A+B*K)*B)


% N=0.9723

Abf=A-B*K;

sys=ss(Abf,B,C,D);
Vp1=eig(Abf);

%**********Observateur identité *************

P=[-2.4+5.5*i -2.4-5.5*i]
G=acker(A',C',P)'
G2=acker(A',C',2*P)'
G4=acker(A',C',4*P)'
G8=acker(A',C',8*P)'

F=A-G*C
F2=A-G2*C
F4=A-G4*C
F8=A-G8*C

%***** Observateur identité minimal*****

A11=A(1,1)
A12=A(1,2)
A21=A(2,1)
A22=A(2,2) 

 
G1=acker(A22',A12',4*(-2.4))'
F1=A22-G1*A12;
%G1=-(F1-A22)/A12
G1_tild=(F1*G1)-(G1*A11)+A21
H1_tild=3.5
M=eye(1)
%N=G

%****** observateur fonctionel******

  % observateur fonctionnel
 Pmin = 2*(-2.4) 
 T1 = acker (A22',-A12',Pmin )
 Ffc = A22 +T1* A12 
 Gfc = A21 + T1* A11 - Ffc *T1
 Nfc = K(1) - K(2)*T1
 M = K(2) 
 T = [T1 1]
 Hfc = T*B
 Afc = Ffc 
 Bfc = [ Gfc Hfc ]
 Cfc = M
 Dfc = [ Nfc 0]















