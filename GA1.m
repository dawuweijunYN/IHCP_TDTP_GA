%function GA
clc
clear all
format short eng
%========== GENETIC  PARAMETER =====================
NS=100;
NG=500;
NP=4;
NN=1;

RC=.9;
PR=.05;
PM=.12;
PC=.9999;
AU(1)=25;AL(1)=7;
AU(2)=.02;AL(2)=.01;
%AU(3)=-0.00003;AL(3)=-0.00001;
AU(3)=550;AL(3)=350;
AU(4)=.4;AL(4)=.2;
%AU(6)=-0.0006;AL(6)=-0.0005;
NSPR=round(NS*PR);
if (NSPR < 1); NSPR=1; end;
%=====================================
DATACONDUC(1)=14.1;         % K1
DATACONDUC(2)=.0166;      % K2
DATACONDUC(3)=-0.000021*0;           % K3
DATACONDUC(4)=448.;         % CP1
DATACONDUC(5)=.291;       % CP2
DATACONDUC(6)=-0.0005597*0;           % CP3
DATACONDUC(7)=6288.25;        % Ro
DATACONDUC(8)=300.;         % T Initial
DATACONDUC(9)=800;         % NT
DATACONDUC(10)=11;          % NX
DATACONDUC(11)=21;          % Ny
DATACONDUC(12)=2;        % DT
DATACONDUC(13)=120*10^3;    %q
DATACONDUC(14)=.3;           % th
DATACONDUC(15)=.1;           % Lx
DATACONDUC(16)=.2;           % Ly
DATACONDUC(17)=4;           % Error
NT=DATACONDUC(9);
NX=DATACONDUC(10);
NY=DATACONDUC(11);

net1=0;
T=DIRECTSOLUTION(DATACONDUC);
ij=0;
for I=1:NX
    for J=1:NY
        ij=ij+1;
        TT1(ij,1)=(I-1)*(DATACONDUC(15)/(NX-1));
        TT1(ij,2)=(J-1)*(DATACONDUC(16)/(NY-1));
        TT1(ij,3)=T(round(NT/2),I,J);
    end
end


zigma=.01
TE=REFRENCESOLUTION(T,DATACONDUC,zigma);
Tmax=max(max(max(TE)))
for m=1:NT
    AA(m,1)=TE(m,1);
    AA(m,2)=TE(m,2);
    AA(m,3)=TE(m,3);
    AA(m,4)=TE(m,4);
    AA(m,5)=TE(m,5);
    AA(m,6)=TE(m,6);
end



%GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
ijk=0;
NNN=1;
A=INTITALPOPULATION(net1,NNN,TE,DATACONDUC,NS,NP,NN,AU,AL);

for K=1:NG
    A=REALTIVEFITTNESS(A,NS,NP);
    A=TOURNOMENTSELECTION(NS,NP,A);
    ANEW=CROOSOVER(A,NS,NP,PC,AU,AL);
    ANEW=MOTTATION(ANEW,AU,AL,PM,NS,NP);
    ANEW=FITTNESSCACULATION(net1,NNN,TE,DATACONDUC,ANEW,NS,NSPR,NP,AU,AL);
    ANEW=FINALSELECTION(A,ANEW,NS,NP,NSPR);

    q=0; q(1:NP)=A(1,1:NP);MAX=A(1,NP+3);AVE=A(2,NP+3); q=[K,q,AVE,MAX,(1/MAX-.001)];q=q'
    A=0;A=ANEW;
    %============================================================
    if (K==1)
        ANN=A
    else
        ANN=[ANN;A];
    end

    if  (K == 22)
        NNN=2;
        B=sortrows(ANN,-(NP+1));
        for i=1:22*NS-1
            if (B(i,1)==B(i+1,1))
                if (B(i,2)==B(i+1,2))
                    if (B(i,3)==B(i+1,3))
                        B(i+1,NP+1)=0;
                    end
                end
            end
        end
        B=sortrows(B,-(NP+1));
        b=B(1:10*NS,1:NP+1);
        B=0;B=b;
       
        UL=minmax(B')       
        for i=1:10*NS
            for j=1:NP+1
                B(i,j)=(B(i,j)-UL(j,1))/(UL(j,2)-UL(j,1));
            end
        end
        B;
        X(1:20*NS,1:4)=B(1:20*NS,1:4);
        Y(1:20*NS,1)=B(1:20*NS,NP+1);
        X=X';
        Y=Y';
        net1 = newff(minmax(X),[20  20 1],{'tansig' 'tansig' 'tansig'},'trainlm', 'learngdm','mse');
        %net1 = newff([10 20],[10  10  1],{'tansig' 'tansig' 'tansig'},'trainlm', 'learngdm','mse')
        %net.trainParam.lr = 0.85;
        net1.trainParam.show =100;
        net1.trainParam.epochs =1000;
        net1.trainParam.goal = 0.;
        net1=train(net1,X,Y);
        an = sim(net1,X)';
        Y=Y';
        AA=[Y,an];

        save net1
        AUN(1:4)=0;ALN(1:4)=1;
        B=INTITALPOPULATION(net1,NNN,TE,DATACONDUC,NS,NP,NN,AUN,ALN);
        for KK=1:NG
            B=REALTIVEFITTNESS(B,NS,NP);
            B=TOURNOMENTSELECTION(NS,NP,B);
            BNEW=CROOSOVER(B,NS,NP,PC,AUN,ALN);
            BNEW=MOTTATION(BNEW,AUN,ALN,PM,NS,NP);
            BNEW=FITTNESSCACULATION(net1,NNN,TE,DATACONDUC,BNEW,NS,NSPR,NP,AUN,ALN);
            BNEW=FINALSELECTION(B,BNEW,NS,NP,NSPR);
            q=0;q(1:NP)=B(1,1:NP);
            for J=1:NP
            q(J)=q(J)*(UL(J,2)-UL(J,1))+UL(J,1);
            end
            MAX=B(1,NP+3);AVE=B(2,NP+3); q=[KK,q,AVE,MAX];q=q
            B=0;B=BNEW;
        end
    NNN=1;
    q
    GENDATA(1:NP)=ANEW(I,1:NP);
    MAX=FITTNESFUNCTION(net1,NNN,TE,DATACONDUC,GENDATA,NP);

    end
    
    %============================================================
    ijk=ijk+1;
    DATAPRINT(ijk,1)=ijk;
    if (ijk==1)
        DATAPRINT(ijk,2)=NS*(NN+1);
    else
        DATAPRINT(ijk,2)=DATAPRINT(ijk-1,2)+NS;
    end
    DATAPRINT(ijk,3)=AVE;
    DATAPRINT(ijk,4)=MAX;
    DATAPRINT(ijk,4+1:5+NP)=q(1+1:NP+2);
    DATAPRINT(ijk,4+1:5+NP)=q(1+1:NP+2);
    %    save
    %============================================================
end
%GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
format long eng
q