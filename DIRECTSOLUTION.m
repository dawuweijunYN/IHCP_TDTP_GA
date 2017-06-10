function T=DIRECTSOLUTION(DATACONDUC)

K1=DATACONDUC(1);
K2=DATACONDUC(2);
K3=DATACONDUC(3);
CP1=DATACONDUC(4);
CP2=DATACONDUC(5);
CP3=DATACONDUC(6);
RO=DATACONDUC(7);
TINTIAL=DATACONDUC(8);
NT=DATACONDUC(9);
NX=DATACONDUC(10);
NY=DATACONDUC(11);
DT=DATACONDUC(12);
Q=DATACONDUC(13);
TH=DATACONDUC(14);
Lx=DATACONDUC(15);
Ly=DATACONDUC(16);
QERROR=DATACONDUC(17);

TIME=DT*(NT-1.);
TH=TH*TIME;

UFF=.9;
TETA=.5;
DX=Lx/(NX-1.);
DY=Ly/(NY-1.);
FAMOUR=0.;

for I=1:NX
    for J=1:NY
        X(I)=(I-1)*DX;
        Y(J)=(J-1)*DY;
    end
end

%intial
T=0.;
T(1,1:NX,1:NY)=TINTIAL;
% SOlVING PROCEDURE
if (FAMOUR < 20000)

    for K=2:NT


        TT(1:NX,1:NY)=T(K-1,1:NX,1:NY);
        T(K,1:NX,1:NY)=T(K-1,1:NX,1:NY);
        ITERATION=0;UF=UFF;


        %=========      SPACAE SOlVING	  ================
        AERROR=-1;
        while AERROR == -1
            if (ITERATION/5 == floor(ITERATION/5)); UF=UF*.8;end;
            ITERATION=ITERATION+1;
            FAMOUR=FAMOUR+1;
            TTP=TT;
            if (FAMOUR > 20000);break;end;
            for J=2:NY-1
                for I=2:NX-1

                    XT=TTP(I,J);AK=K1+K2*XT+K3*XT^2;CP=CP1+CP2*XT+CP3*XT^2;
                    %                 XT=TTP(I+1,J);AKE=K1+K2*XT+K3*XT^2;%CPE=CP1+CP2*XT+CP3*XT^2;
                    %                 XT=TTP(I-1,J);AKW=K1+K2*XT+K3*XT^2;%CPW=CP1+CP2*XT+CP3*XT^2;
                    %                 XT=TTP(I,J+1);AKN=K1+K2*XT+K3*XT^2;%CPN=CP1+CP2*XT+CP3*XT^2;
                    %                 XT=TTP(I,J-1);AKS=K1+K2*XT+K3*XT^2;%CPS=CP1+CP2*XT+CP3*XT^2;

                    XT=T(K-1,I,J);AKn=K1+K2*XT+K3*XT^2;CPn=CP1+CP2*XT+CP3*XT^2;
                    %                 XT=T(K-1,I+1,J);AKnE=K1+K2*XT+K3*XT^2;%CPnE=CP1+CP2*XT+CP3*XT^2;
                    %                 XT=T(K-1,I-1,J);AKnW=K1+K2*XT+K3*XT^2;%CPnW=CP1+CP2*XT+CP3*XT^2;
                    %                 XT=T(K-1,I,J+1);AKnN=K1+K2*XT+K3*XT^2;%CPnN=CP1+CP2*XT+CP3*XT^2;
                    %                 XT=T(K-1,I,J-1);AKnS=K1+K2*XT+K3*XT^2;%CPnS=CP1+CP2*XT+CP3*XT^2;

                    ALFAn=AKn/(CPn*RO);
                    ALFA=AK/(CP*RO);

                    %LANDAn=1./(RO*CPn)*((AKnE-AKnW)/(2.*DX)*(T(K-1,I+1,J)-T(K-1,I-1,J))/(2*DX)+(AKnN-AKnS)/(2.*DY)*(T(K-1,I,J+1)-T(K-1,I,J-1))/(2*DY));
                    %LANDA=1./(RO*CP)*((AKE-AKW)/(2.*DX)*(T(K,I+1,J)-T(K,I-1,J))/(2*DX)+(AKN-AKS)/(2.*DY)*(T(K,I,J+1)-T(K,I,J-1))/(2*DY));
                    LANDAn=0.;LANDA=0;
                    An=ALFAn/DX^2.*(T(K-1,I+1,J)-2.*T(K-1,I,J)+T(K-1,I-1,J))+ALFAn/DY^2.*(T(K-1,I,J+1)-2.*T(K-1,I,J)+T(K-1,I,J-1));
                    A=ALFA/DX^2.*(TTP(I+1,J)+TTP(I-1,J))+ALFA/DY^2.*(TTP(I,J+1)+TTP(I,J-1));

                    B=(1-TETA)*(LANDAn+An)+TETA*(LANDA+A)+T(K-1,I,J)/DT;
                    C=1./DT-TETA*ALFA*(-2/DX^2-2/DY^2);

                    TT(I,J)=TTP(I,J)+UF*(-TTP(I,J)+B/C);
                end
            end

            for J=1:NY
                TT(NX,J)=TT(NX-1,J);
                XT=TT(1,J);AK=K1+K2*XT+K3*XT^2;
                if ((K-1)*DT >= TH)
                    bQ=0;
                else
                    if ((J-1)*DY >= .1)
                        bQ=Q;
                    else
                        bQ=0;
                    end
                end
                TT(1,J)=(18*TT(2,J)-9*TT(3,J)+2*TT(4,J)-(-bQ)*6*DX/AK)/11;
            end

            for I=1:NX
                TT(I,NY)=TT(I,NY-1);
                XT=TT(I,1);AK=K1+K2*XT+K3*XT^2;
                if ((K-1)*DT >= TH)
                    bQ=0;
                else
                    bQ=0;
                end
                TT(I,1)=(18*TT(I,2)-9*TT(I,3)+2*TT(I,4)-(-bQ)*6*DX/AK)/11;
            end
            %=========     ERROR EVAlUATION  ================
            AERROR=1;
            for J=1:NY
                for I=1:NX
                    ERROR=0.;
                    if (abs(TT(I,J))> 10^(-16))
                        ERROR=abs(TTP(I,J)-TT(I,J))/abs(TT(I,J));
                    end
                    %                if (ITERATION < 5); AERROR = -1; end;
                    if (ERROR > 10^(-QERROR));AERROR=-1;end;
                    if (AERROR == -1) ; break; end;
                end
                if (AERROR == -1) ; break; end;
            end

            %=================================================

        end
        %=========     SUBSTITUTING     ================
        T(K,1:NX,1:NY)=TT(1:NX,1:NY);


    end
else
    T(1:NT,1:NX,1:NY)=TINTIAL;
    FAMOUR
    DATACONDUC
end
