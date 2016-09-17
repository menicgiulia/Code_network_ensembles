 function [W]=Wupdate(B,class,vecM)
 Nd=length(B);
 B2=zeros(Nd,1);
 W=zeros(Nd,1);

 parfor d=1:Nd
                 	
    B2(d)=sum(vecM.*(class==d));%puoi fare a meno del diag diag
    if (B2(d)*B(d))>0.0
       W(d)=B(d)/(B2(d));
       W(d)=max(W(d),10^(-15));
       W(d)=min(W(d),10^15);
    else
         W(d)=0;
    end
 end