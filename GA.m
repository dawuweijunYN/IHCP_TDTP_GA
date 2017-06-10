%function GA
clc
clear all
format short eng
%========== GENETIC  PARAMETER =====================
NS=40;
NG=1000;
NP=input('NP Dimensional?')
zigma=input('sigma?')
q=input('Constant[0]  Linear[1]  Parabolic[2]?')
if (NP == 6);NN=40;end;
if (NP == 4);NN=30;end;
if (NP == 2);NN=20;end;
%NN=1
if (NP == 6);Rc=.99;end;
if (NP == 4);Rc=.99;end;
if (NP == 2);Rc=.98;end;
PR=.05;
PM=.12;
PC=.9999;

if (NP==2)
    AU(1)=100;AL(1)=0;
    AU(2)=1000;AL(2)=0;
    AU(1)=15;AL(1)=13.5;
    AU(2)=455;AL(2)=440;
 end
if (NP==4)
    AU(1)=100/2;AL(1)=0;
    AU(2)=1/4;AL(2)=-1/4;
    AU(3)=1000;AL(3)=250;
    AU(4)=10/4;AL(4)=-10/4;

    AU(1)=15;AL(1)=13.5;
    AU(2)=.03;AL(2)=.015;
    AU(3)=455;AL(3)=440;
    AU(4)=.63;AL(4)=0.45;

end
if (NP==6)
    AU(1)=100/4;AL(1)=5;
    AU(2)=1/10;AL(2)=-1/20;        %.0166
    AU(3)=0.01/100;AL(3)=-0.01/100;  %-0.000021
    AU(4)=550;AL(4)=350;
    AU(5)=10/20;AL(5)=-10/50;      %.291
    AU(6)=0.1/100;AL(6)=-0.1/100;    %-0.0005597

    AU(1)=17.;AL(1)=12.5;
    AU(2)=.01;AL(2)=-.01;        %.02
    AU(3)=+0.0001;AL(3)=-0.0001;  %-0.000021
    AU(4)=470;AL(4)=430;
    AU(5)=.1;AL(5)=-.1;      %.54
    AU(6)=0.01;AL(6)=-0.01;    %-0.0005597

end
NSPR=round(NS*PR);
%=====================================
LP1=1;LP2=1;
if (q == 0)
    DATACONDUC(1)=14.1;         % K1
    DATACONDUC(2)=0;      % K2
    DATACONDUC(3)=0;           % K3
    DATACONDUC(4)=448.;         % CP1
    DATACONDUC(5)=0;       % CP2
    DATACONDUC(6)=0;           % CP3

end
if (q == 1)
    DATACONDUC(1)=14.1;         % K1
    DATACONDUC(2)=.0166;      % K2
    DATACONDUC(3)=0;           % K3
    DATACONDUC(4)=448.;         % CP1
    DATACONDUC(5)=.291;       % CP2
    DATACONDUC(6)=0;           % CP3

end
if (q == 2)
    DATACONDUC(1)=14.1;         % K1
    DATACONDUC(2)=0.0200;      % K2
    DATACONDUC(3)=-0.000021;           % K3
    DATACONDUC(4)=448.;         % CP1
    DATACONDUC(5)=0.5455;       % CP2
    DATACONDUC(6)=-0.0005597;           % CP3

end

% DATACONDUC(1)=14.1;         % K1
% DATACONDUC(2)=.0166*LP1;      % K2
% DATACONDUC(3)=-0.000021*LP2;           % K3
% DATACONDUC(4)=448.;         % CP1
% DATACONDUC(5)=.291*LP1;       % CP2
% DATACONDUC(6)=-0.0005597*LP2;           % CP3
DATACONDUC(7)=6288.25;        % Ro
DATACONDUC(8)=23.;         % T Initial
DATACONDUC(9)=1000/2;         % NT
DATACONDUC(10)=11;          % NX
DATACONDUC(11)=21;          % Ny
DATACONDUC(12)=2*2;        % DT
DATACONDUC(13)=120*10^3;    %q
DATACONDUC(14)=.3;           % th
DATACONDUC(15)=.1;           % Lx
DATACONDUC(16)=.2;           % Ly
DATACONDUC(17)=6;           % Error
DATACONDUC(18)=DATACONDUC(13)*(.08)/DATACONDUC(1);
NT=DATACONDUC(9);
NX=DATACONDUC(10);
NY=DATACONDUC(11);
TIME=DATACONDUC(12)*(NT-1.)
TIMEPLUS=TIME*(DATACONDUC(1)/(DATACONDUC(7)*DATACONDUC(4)))/(DATACONDUC(15)^2)
T=DIRECTSOLUTION(DATACONDUC);

TE=REFRENCESOLUTION(T,DATACONDUC,zigma);
Tmax=max(max(max(TE)))
ij=0;
for m=1:NT
    AA(m,1)=TE(m,1);
    AA(m,2)=TE(m,2);
    AA(m,3)=TE(m,3);
    AA(m,4)=TE(m,4);
    AA(m,5)=TE(m,5);
    AA(m,6)=TE(m,6);
end
for I=1:NX
    for J=1:NY
        ij=ij+1;
        TT1(ij,1)=(I-1)*(DATACONDUC(15)/(NX-1));
        TT1(ij,2)=(J-1)*(DATACONDUC(16)/(NY-1));
        TT1(ij,3)=T(round(NT/2),I,J);
    end
end


%GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
A=INTITALPOPULATION(TE,DATACONDUC,NS,NP,NN,AU,AL);
for K=1:NG
    A=REALTIVEFITTNESS(A,NS,NP);
    A=TOURNOMENTSELECTION(NS,NP,A);
    ANEW=CROOSOVER(A,NS,NP,PC,AU,AL);
    ANEW=MOTTATION(ANEW,AU,AL,PM,NS,NP);
    ANEW=FITTNESSCACULATION(TE,DATACONDUC,ANEW,NS,NSPR,NP,AU,AL);
    ANEW=FINALSELECTION(A,ANEW,NS,NP,NSPR);
    %==========================================================================
    q=0;q(1:NP)=A(1,1:NP);MAX=A(1,NP+3);AVE=A(2,NP+3); q=[K,q,AVE,MAX,(1/MAX-.001)^2];q=q'
    A=0;A=ANEW;
    if (K ==1)
        DATAPRINT=q';
    else
        DATAPRINT=[DATAPRINT;q'];
    end
    %======================================================================
    %    if (K > 50); NS=100;end
    if (K > 200)
        format long eng
    end
    %     if (PR < .1)
    %         PR=PR*1.005;
    %         NSPR=round(NS*PR);
    %         if (NSPR < 1);NSPR=1;end;
    %     end
    %=======
    if ( K > 1 );PM=.13;end;
    if ( K > 25);PM=.12;end;
    if ( K > 100);PM=.11;end;
    if ( K > 150);PM=.10;end;
    if ( K > 200);PM=.1;end;
    %     if ( K > 100);PM=.1;end;
    %==========
    if (K > 2 )
        AU(1:NP)=(1-Rc)*A(1,1:NP)+Rc*AU(1:NP);
        AL(1:NP)=(1-Rc)*A(1,1:NP)+Rc*AL(1:NP);
        qq=[AU;AL]
    end
    save
    %======================================================================
end
%GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

q