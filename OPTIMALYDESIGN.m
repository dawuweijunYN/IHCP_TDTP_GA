clc
clear all
%format long
format short eng
DATACONDUC(1)=14.1;         % K1
DATACONDUC(2)=.0200;      % K2
DATACONDUC(3)=-0.000021;           % K3
DATACONDUC(4)=448.;         % CP1
DATACONDUC(5)=.5455;       % CP2
DATACONDUC(6)=-0.0005597;           % CP3
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
%TIMEPLUS=TIME*(DATACONDUC(1)/(DATACONDUC(7)*DATACONDUC(4)))/(DATACONDUC(15)^2)
T=DIRECTSOLUTION(DATACONDUC);
%TIMEPLUS=TIME*(DATACONDUC(1)/(DATACONDUC(7)*DATACONDUC(4)))/(DATACONDUC(15)^2)
epsilon=10^(-6)
zigma=.01*0;

%TT(1:4,1:4)=TE(900,1:4,1:4)
iii=0;

for v=.3:.1:.3

    DATACONDUC(14)=v;
    T=DIRECTSOLUTION(DATACONDUC);
    TE=REFRENCESOLUTION(T,DATACONDUC,zigma);
    Tmax=max(max(max(TE)));
    %Tmax=800;
    X=0;XN=0;XM=0;
    for L=1:6
        %         if (L==1);LL=1; end;
        %         if (L==2);LL=2; end;
        %         if (L==3);LL=4; end;
        %         if (L==4);LL=5; end;
        LL=L;
        DATACONDUC1=DATACONDUC;
        DATACONDUC1(LL)=DATACONDUC(LL)*epsilon+DATACONDUC(LL);
        T=DIRECTSOLUTION(DATACONDUC1);
        TE2=REFRENCESOLUTION(T,DATACONDUC1,zigma);

        for K=2:NT
            ij=0;
            for i=1:6
                ij=ij+1;
                x=(TE2(K,i)-TE(K,i))/(epsilon*DATACONDUC1(LL));
                x=x*DATACONDUC1(LL)/(Tmax-DATACONDUC(8));
                XM(ij,1)=x;
            end

            if (K==2)
                XN=XM;
            else
                XN=[XN;XM];
            end
            XM=0;
        end
        if (L==1)
            X=XN;
        else
            X=[X,XN];
        end
        XN=0;

    end
    F=X'*X;
    iii=iii+1;
    D(iii,2)=det(F);
    D(iii,1)=v

    for ii=1:1
        jjj=0;
        for jk=2:NT
            jjj=jjj+1;
            BB(jjj,ii)=X((jk-2)*2+3,ii);
        end
    end
end
%=======================================

for K=1:NT-1
    for m=1:6
        for n=1:1
            AX(n,1)=X((K-1)*6+n,m);
        end

        F=AX'*AX
        FNP(K,m)=det(F);
    end
end