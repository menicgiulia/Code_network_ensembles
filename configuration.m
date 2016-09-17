%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
% pp unweighted undirected adjacency matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
% S_c canonical entropy
% P link probability matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S_c,P]=configuration(pp)
precision=10^(-5);
loops=10000;

n=size(pp,1);
connectivity=sum(pp,2);

avg_conn=sum(connectivity)/n;
% suggested starting conditions
z=connectivity/(sqrt(n*avg_conn));
oldz=zeros(n,1);

for kk=1:loops
        
            U=ones(n,1)*z';
            D=ones(n,n) + z*z';  
            z=connectivity ./ sum(( U./D - diag(diag(U./D))),2);
            z=max(z,10^(-15));
        
        if max(abs((z>0).*(1-z./(oldz+(oldz==0)))))<precision
            break
        end
        oldz=z;
	end

    
%Compute link probability
P=(z*z')./(ones(n,n)+z*z');
P=P-diag(diag(P));

%Compute Shannon entropy
P1=P.*log(P+(P==0));
P2=(ones(n,n)-P).*log(ones(n,n)-P +((ones(n,n)-P)==0));
Smatrix=-triu(P1+P2,1);
S_c=sum(sum(Smatrix));
display(kk)

return
