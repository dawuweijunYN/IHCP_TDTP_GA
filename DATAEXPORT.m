%function GA
clc
clear all
format long eng

DATACONDUC(1)=14.1;         % K1
DATACONDUC(2)=.0166;      % K2
DATACONDUC(3)=0;           % K3
DATACONDUC(4)=448.;         % CP1
DATACONDUC(5)=.291;       % CP2
DATACONDUC(6)=0;           % CP3
DATACONDUC(7)=6288.25;        % Ro
DATACONDUC(8)=23.;         % T Initial
DATACONDUC(9)=1000/2+1;         % NT
DATACONDUC(10)=11;          % NX
DATACONDUC(11)=21;          % Ny
DATACONDUC(12)=2*2;        % DT
DATACONDUC(13)=120*10^3;    %q
DATACONDUC(14)=.3;           % th
DATACONDUC(15)=.1;           % Lx
DATACONDUC(16)=.2;           % Ly
DATACONDUC(17)=7;           % Error
NT=DATACONDUC(9);
NX=DATACONDUC(10);
NY=DATACONDUC(11);
DT=DATACONDUC(12);
TIME=DATACONDUC(12)*(NT-1.)
TIMEPLUS=TIME*(DATACONDUC(1)/(DATACONDUC(7)*DATACONDUC(4)))/(DATACONDUC(15)^2)

% DATACONDUC(1)=14.14731903419716;         % K1
% DATACONDUC(2)=0.016458479016857;      % K2
% DATACONDUC(3)=0;           % K3
% DATACONDUC(4)=447.0254445212983;         % CP1
% DATACONDUC(5)=0.309626632851004;       % CP2
% DATACONDUC(6)=0;           % CP3


T=DIRECTSOLUTION(DATACONDUC);
zigma=.01*0;
TE=REFRENCESOLUTION(T,DATACONDUC,zigma);
Tmax=max(max(max(TE)))
ij=0;
for m=1:NT
    TIME=DATACONDUC(12)*(m-1.);
    TIMEPLUS=TIME*(DATACONDUC(1)/(DATACONDUC(7)*DATACONDUC(4)))/(DATACONDUC(15)^2);
    AA(m,1)=TIMEPLUS;
    AA(m,2)=TE(m,1);
    AA(m,3)=TE(m,2);
    AA(m,4)=TE(m,3);
    AA(m,5)=TE(m,4);
    AA(m,6)=TE(m,5);
    AA(m,7)=TE(m,6);
end
Ia=1;
A2=DATACONDUC(Ia)+.05*DATACONDUC(Ia)
A1=DATACONDUC(Ia)-.05*DATACONDUC(Ia)
DATACONDUC1=DATACONDUC;
qlk=DATACONDUC(13)*(.08)/DATACONDUC(1)
for jj=1:-5
    if (jj==1);zigma=.01*0;end
    if (jj==2);zigma=.01;end
    if (jj==3);zigma=.01;end
    if (jj==4);zigma=.01;end
    if (jj==5);zigma=.01;end
    T=DIRECTSOLUTION(DATACONDUC);
    TE=REFRENCESOLUTION(T,DATACONDUC,zigma);
    ij=0;
    for A=A1:(A2-A1)/50:A2
        DATACONDUC1(Ia)=A;
        T=DIRECTSOLUTION(DATACONDUC1);
        q=ARMS(T,TE,DATACONDUC1)/qlk^2;
        FITTNES=1/(.001+sqrt(q));
        ij=ij+1;
        TT1(ij,1)=A;
        TT1(ij,jj+1)=(FITTNES);
    end
    q(1:51,1)= TT1(1:51,1)
end