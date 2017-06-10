function BARGHRAFTGA
format short eng
load
AK=K+1
%GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
R=random('unif',100,400);
rand('twister', sum(R*clock));
for K=AK:NG
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
     if ( K > 1 );PM=.15;end;
     if ( K > 25);PM=.14;end;
     if ( K > 100);PM=.13;end;
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