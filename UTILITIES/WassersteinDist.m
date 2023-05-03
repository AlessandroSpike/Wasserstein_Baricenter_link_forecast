function WD = WassersteinDist(P1,P2,xi,xit,C,maxiter)

% inizializzazione
b = ones(size(xi,1),1);
iter = 1;
while iter<maxiter  
    a = P1 ./ (xi*b);
    b = P2 ./ (xit*a);    
   iter = iter+1;
end

Pi = diag(a)*xi*diag(b);
WD = sum(sum(Pi.*C)); 


