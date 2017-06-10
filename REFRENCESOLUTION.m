function TE=REFRENCESOLUTION(T,DATACONDUC,zigma)


NT=DATACONDUC(9);
NX=DATACONDUC(10);
NY=DATACONDUC(11);
Lx=DATACONDUC(15);
Ly=DATACONDUC(16);

DX=Lx/(NX-1.);
DY=Ly/(NY-1.);

for I=1:NX
    for J=1:NY
        X(I)=(I-1)*DX;
        Y(J)=(J-1)*DY;
    end
end
XE(1)=0;XE(2)=0;XE(3)=.5;XE(4)=.5;XE(5)=.5;XE(6)=.5;
YE(1)=.5+(1/2)*(1/3);YE(2)=.5+(1/2)*(2/3);YE(3)=.5-1/4;YE(4)=.5;YE(5)=.5+1/4;YE(6)=1;
XE=XE*Lx;
YE=YE*Ly;
for K=1:NT
    for II=1:6
        XX=XE(II);
        YY=YE(II);
        TE(K,II)=-.1;
        for I=1:NX-1
            QA=1;
            for J=1:NY-1
                QA=1;
                if ((XX >= X(I))&&(XX <= X(I+1)))
                    if ((YY >= Y(J))&&(YY <= Y(J+1)))
                        %Q=0;
                        DX=X(I+1)-X(I);
                        DY=Y(J+1)-Y(J);
                        TX1=(T(K,I+1,J)-T(K,I,J))/(DX)*(XX-X(I))+T(K,I,J);
                        TX2=(T(K,I+1,J+1)-T(K,I,J+1))/(DX)*(XX-X(I))+T(K,I,J+1);
                        TXY=(TX2-TX1)/(DY)*(YY-Y(J))+TX1;
                        TE(K,II)=TXY;
                        QA=-1;
                    end
                end
                if (QA==-1);break;end
            end
            if (QA==-1);break;end
        end
    end
end
Tmax=max(max((TE)));
%Tmax=
TE2=random('norm',TE,zigma*Tmax);
a=zigma;
if  (a==0)
    TE=TE2;
else
    for L=1:6
        for II=1:6
            for K=2:NT-1
                TE3(K,II)=(TE2(K-1,II)+2*TE2(K,II)+TE2(K+1,II))/4;
            end
            TE3(1,II)=2*TE3(2,II)-TE3(3,II);
            TE3(NT,II)=2*TE3(NT-1,II)-TE3(NT-2,II);
        end
        TE2=TE3;
    end
    TE=TE3;
end
%                             Q(1,1)=T(K,I,J);Q(1,2)=T(K,I+1,J);Q(1,3)=TX1;
%                             Q(2,1)=T(K,I,J+1);Q(2,2)=T(K,I+1,J+1);Q(2,3)=TX2;
%                             Q(3,1)=TXY;Q(3,2)=II;Q(3,3)=JJ;
%                             Q(4,1)=I;Q(4,2)=J;Q(4,3)=K;
%                             Q(5,1)=X(I);Q(5,2)=X(I+1);Q(5,3)=XX;
%                             Q(6,1)=Y(J);Q(6,2)=Y(J+1);Q(6,3)=YY;
%                             Q;
