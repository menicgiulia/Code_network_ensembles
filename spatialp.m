%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
% ppvec unweighted undirected adjacency matrix (vector=squareform(matrix))
% ddvec matrix of distances between each node i and node j of the network (vector=squareform(matrix))
% n number of nodes
% Nd number of bins Nd that we consider in the distance classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
% S_d canonical entropy 
% P link probability matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S_d,P]=spatialp(ppvec,ddvec,n,Nd)

precision=10^(-5);
loops=10000;


dmax=max(ddvec);
dmin=min(ddvec);

%linear binning
l=linspace(dmin,dmax,Nd+1);
l(Nd+1)=dmax+2*eps;



class=zeros(1,n*(n-1)/2);


	for ind = 2:length(l)
		if ind==length(l)
            class(ddvec>=l(ind-1) & ddvec<=l(ind))=ind-1;
		break
		end
   		class(ddvec>=l(ind-1) & ddvec<l(ind))=ind-1;
	end



k=sum(squareform(ppvec), 2);
avg_conn=sum(k)/n;


B=zeros(Nd,1);
	for d=1:Nd
        B(d)=sum(ppvec.*(class==d));
    end
    
clear ppvec

% compute exp(Lagrangian multipliers)
% suggested starting conditions
zold=k/(sqrt(n*avg_conn));
Wold=B/(n*avg_conn);




         for kk=1:loops
             
             bigWvec=zeros(n*(n-1)/2,1);

	         for d=1:Nd,
                bigWvec(class==d)=Wold(d);
             end
             
             bigW=squareform(bigWvec);
             clear bigWvec
             
             % update z
             z=zupdate(k, zold, bigW);
               
             % update W
             M= (z*z') ./ ( 1 + ((z*z') .* bigW)) ;
             M=M-diag(diag(M));
             vecM=squareform(M);
             clear M
              
             W=Wupdate(B,class,vecM);
          
             
             

         if (max(abs((W>0).*(1-W./(Wold+(Wold==0)))))<precision) && (max(abs((z>0).*(1-z./(zold+(zold==0)))))<precision)
         break
         end

         Wold=W; zold=z;
         end


bigWvec=zeros(n*(n-1)/2,1);
for d=1:Nd,
  bigWvec(class==d)=W(d);
end
             
bigW=squareform(bigWvec);
clear bigWvec

    
    
%Compute link probability
P=(z*z'.*bigW)./(1+z*z'.*bigW);
P=P-diag(diag(P));
P=squareform(P);

%Compute Shannon entropy
P1=P.*log(P+(P==0));
P2=(1-P).*log(1-P +((1-P)==0));
Svec=-(P1+P2);
S_d=sum(Svec);
display(kk)


return

