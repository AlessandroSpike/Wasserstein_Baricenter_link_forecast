function [P,w]=MainOT(X,y,N,K,KT,S,L,lambda)
%% initialization
b = ones(N,S,L);
Phi = zeros(N,S,L);
lcs = logsumexp(lambda);
lambda = log(lambda)+lcs;
%% SINKHORN LOOP
for l = 2:L        
   Phi_lambda = zeros(N,S);
   for s = 1:S            
       % auxiliary
       Phi(:,s,l) =  KT*(X(:,s)./(K*b(:,s,l-1))); 
       Phi_lambda(:,s) = Phi(:,s,l).^lambda(s);
    end
    % baricenter
    P = prod(Phi_lambda,2);
    if sum(isnan(P))>0
        error('Computation blew up, epsilon too small');
    end
    % update
    b(:,:,l) = repmat(P,1,S)./Phi(:,:,l);   
end
%% GRADIENT OF THE LOSS
g = (P-y).*P;
%% REVERSE LOOP
w = zeros(S,1);
r = zeros(N,S);
for l = L:-1:2
    for s = 1:S
        % auxiliary
        w(s)=w(s) + trace(log(Phi(:,s,l))'*g);
        aus1 = (lambda(s)*g-r(:,s))./Phi(:,s,l);
        aus2 = X(:,s)./((K*b(:,s,l-1)).^2);
        r(:,s)=-KT*(K*aus1.*aus2).*b(:,s,l-1);        
    end
    g = sum(r,2);
end


