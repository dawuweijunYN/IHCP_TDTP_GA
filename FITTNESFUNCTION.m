function FITTNES=FITTNESFUNCTION(TE,DATACONDUC,GENDATA,NP)
qlk=DATACONDUC(18);
if (NP==2)
    DATACONDUC(1)=GENDATA(1);         % K1
    DATACONDUC(4)=GENDATA(2);         % CP1
end
if (NP==4)
    DATACONDUC(1)=GENDATA(1);         % K1
    DATACONDUC(2)=GENDATA(2);      % K2
    DATACONDUC(4)=GENDATA(3);         % CP1
    DATACONDUC(5)=GENDATA(4);       % CP2

end
if (NP==6)
    DATACONDUC(1)=GENDATA(1);         % K1
    DATACONDUC(2)=GENDATA(2);      % K2
    DATACONDUC(3)=GENDATA(3);           % K3
    DATACONDUC(4)=GENDATA(4);         % CP1
    DATACONDUC(5)=GENDATA(5);       % CP2
    DATACONDUC(6)=GENDATA(6);           % CP3
end
% ======================================================================
T=DIRECTSOLUTION(DATACONDUC);
B1=ARMS(T,TE,DATACONDUC)/qlk^2;
FITTNES=1/(.001+sqrt(B1));
%==========================================================================
if ((FITTNES > -1000) && (FITTNES < 10000))
    c=1;
else
    FITTNES=0;
    %    DATACONDUC'
end