%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
% ppvec unweighted undirected adjacency matrix (vector=squareform(matrix))
% n number of nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
% S_c canonical entropy 
% P link probability matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S_c,P]=configurationp(ppvec,n)

precision=10^(-5);
loops=100000;

k=sum(squareform(ppvec), 2);
avg_conn=sum(k)/n;

    
clear ppvec

% compute exp(Lagrangian multipliers)
% suggested starting conditions
zold=k/(sqrt(n*avg_conn));




         for kk=1:loops
             
             
             % update z
             z=zupdateconf(k, zold);
  

         if (max(abs((z>0).*(1-z./(zold+(zold==0)))))<precision)
         break
         end

         zold=z;
         end

      
%Compute link probability
P=(z*z')./(1+z*z');
P=P-diag(diag(P));
P=squareform(P);

%Compute Shannon entropy per node
P1=P.*log(P+(P==0));
P2=(1-P).*log(1-P +((1-P)==0));

Svec=-(P1+P2);

S_c=sum(Svec);
display(kk)


return