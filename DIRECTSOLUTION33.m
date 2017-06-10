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

TETA=.5;
DX=Lx/(NX-1.);
DY=Ly/(NY-1.);

%intial
T=0.;
T(1,1:NX,1:NY)=TINTIAL;
% SOlVING PROCEDURE
for K=2:NT
    K;
    TN(1:NX,1:NY)=T(K-1,1:NX,1:NY);
    ij=0;
    %==========  i=2:NX-1   j=2:NY-1  ======   Governing Equation  ==========================
    for I=2:NX-1
        for J=2:NY-1
            ij=ij+1;
            XT=TN(I,J);
            KK=K1+K2*XT+K3*XT^2;
            CP=CP1+CP2*XT+CP3*XT^2;
            ALFA=KK/(CP*RO);

            LANDA=(1-TETA)*((TN(I+1,J)-2*TN(I,J)+TN(I-1,J))/DX^2+(TN(I,J+1)-2*TN(I,J)+TN(I,J-1))/DY^2);
            LANDA=LANDA+RO*CP/(KK*DT)*TN(I,J);

            A(ij,(I-1)*NY+J)=RO*CP/(KK*DT)+2*TETA/(DX^2)+2*TETA/(DY^2);
            A(ij,((I+1)-1)*NY+J)=-TETA/(DX^2);
            A(ij,((I-1)-1)*NY+J)=-TETA/(DX^2);
            A(ij,(I-1)*NY+(J+1))=-TETA/(DY^2);
            A(ij,(I-1)*NY+(J-1))=-TETA/(DY^2);
            C(ij,1)=LANDA;
        end
    end

    %==========  i=2:NX-1   j=1     ======   X-direction Bottom  ==========================
    J=1;
    for I=2:NX-1
        ij=ij+1;
        XT=TN(I,J);
        KK=K1+K2*XT+K3*XT^2;
        bQ=0;
        A(ij,(I-1)*NY+J)=-11/(6*DX);
        A(ij,(I-1)*NY+(J+1))=18/(6*DX);
        A(ij,(I-1)*NY+(J+2))=-9/(6*DX);
        A(ij,(I-1)*NY+(J+2))=2/(6*DX);
        C(ij,1)=-bQ/KK;
    end
    %==========  i=2:NX   j=NY    ======   X-direction Top     ==========================
    J=NY;
    for I=2:NX-1
        ij=ij+1;
        XT=TN(I,J);
        KK=K1+K2*XT+K3*XT^2;
        bQ=0;
        A(ij,(I-1)*NY+J)=11/(6*DX);
        A(ij,(I-1)*NY+(J-1))=-18/(6*DX);
        A(ij,(I-1)*NY+(J-2))=9/(6*DX);
        A(ij,(I-1)*NY+(J-3))=-2/(6*DX);
        C(ij,1)=-bQ/KK;
    end
    %==========  i=1      j=2:NY  ======   Y-direction Left    ==========================
    I=1;
    for J=1:NY
        ij=ij+1;
        XT=TN(I,J);
        KK=K1+K2*XT+K3*XT^2;
        if ((K-1)*DT >= TH)
            bQ=0;
        else
            if ((J-1)*DY >= .1)
                bQ=Q;
            else
                bQ=0;
            end
        end
        A(ij,(I-1)*NY+J)=-11/(6*DX);
        A(ij,((I+1)-1)*NY+J)=18/(6*DX);
        A(ij,((I+2)-1)*NY+J)=-9/(6*DX);
        A(ij,((I+3)-1)*NY+J)=2/(6*DX);
        C(ij,1)=-bQ/KK;
    end
    %==========  i=NX     j=2:NY  ======   Y-direction Right   ==========================
    I=NX;
    for J=1:NY
        ij=ij+1;
        XT=TN(I,J);
        KK=K1+K2*XT+K3*XT^2;
        bQ=0;
        A(ij,(I-1)*NY+J)=11/(6*DX);
        A(ij,((I-1)-1)*NY+J)=-18/(6*DX);
        A(ij,((I-2)-1)*NY+J)=9/(6*DX);
        A(ij,((I-3)-1)*NY+J)=-2/(6*DX);
        C(ij,1)=-bQ/KK;
    end
    %======================================================================

    % =====   Substituting  ========
    X=A^(-1)*C;
    n=max(size(X));
    for ij=1:n
        I=floor(ij/NY-10^(-10))+1;
        J=ij-(I-1)*NY;
        TT(I,J)=X(ij) ;
    end
    T(K,1:NX,1:NY)=TT(1:NX,1:NY);
end
