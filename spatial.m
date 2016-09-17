%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
% pp undirected adjacency matrix
% dd matrix of distances between each node i and node j of the network
% Nd number of bins Nd that we consider in the distance classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
% S_d canonical entropy
% P link probability matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S_d, P]=spatial(pp,dd,Nd)

precision=10^(-5);
loops=100000;



n=size(pp,1);
dd=dd-diag(diag(dd));

dmax=max(squareform(dd));
dmin=min(squareform(dd));

%linear binning
l=linspace(dmin,dmax,Nd+1);
l(Nd+1)=dmax+2*eps;


class=zeros(n,n);


	for ind = 2:length(l)
		if ind==length(l)
            class(dd>=l(ind-1) & dd<=l(ind))=ind-1;
		break
		end
   		class(dd>=l(ind-1) & dd<l(ind))=ind-1;
    end
class=class-diag(diag(class));



connectivity=sum(pp,2);
avg_conn=sum(connectivity)/n;


B=zeros(Nd,1);
	for d=1:Nd
        B(d)=sum(sum(pp.*(class==d)))/2;
    end

% compute exp(Lagrangian multipliers)
% suggested starting conditions
z=connectivity/(sqrt(n*avg_conn));
W=B/(n*avg_conn);

oldW=zeros(Nd,1); 
oldz=zeros(n,1);
         
         for kk=1:loops
         bigW=zeros(n,n);

	         for d=1:Nd,
	         bigW=bigW+W(d)*(class==d);
	         end

         U=ones(n,1)*z' .* bigW;
         D=ones(n,n) + ( z*z'.* bigW);
         z=connectivity ./ sum(U./D,2);
         z=max(z,10^(-15));
     
         B2=zeros(Nd,1);

	         for d=1:Nd
                 	M=(class==d) .* (z*z') ./ ( ones(n,n) + ((z*z') .* bigW)) ;
                 	B2(d)=(sum(sum(M)) )/2;
                 	if (B2(d)*B(d))>0.0
                      	W(d)=B(d)/(B2(d));
                     	W(d)=max(W(d),10^(-15));
                      	W(d)=min(W(d),10^15);
                 	else
                      	W(d)=0;
                 	end
             end
          
             
             

         if (max(abs((W>0).*(1-W./(oldW+(oldW==0)))))<precision) && (max(abs((z>0).*(1-z./(oldz+(oldz==0)))))<precision)
         break
         end

         oldW=W; oldz=z;
         end


bigW=zeros(n,n);
for d=1:Nd,
bigW=bigW+W(d)*(class==d);
end

    
    
%Compute link probability
P=(z*z'.*bigW)./(ones(n,n)+z*z'.*bigW);
P=P-diag(diag(P));


%Compute Shannon entropy
P1=P.*log(P+(P==0));
P2=(ones(n,n)-P).*log(ones(n,n)-P +((ones(n,n)-P)==0));
Smatrix=-triu(P1+P2,1);
S_d=sum(sum(Smatrix));
display(kk)


return