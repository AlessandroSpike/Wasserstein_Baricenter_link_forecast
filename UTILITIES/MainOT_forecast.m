function P=MainOT_forecast(X,N,K,KT,S,L,lambda)
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
    % update
    b(:,:,l) = repmat(P,1,S)./Phi(:,:,l);   
end