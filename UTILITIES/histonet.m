function [count,centri,pos] = histonet(X,Cap_bin,Nodi)
[NodiSort,pos] = sort(Nodi);
centri{1}=Cap_bin(1:end-1)+ diff(Cap_bin)/2;
centri{2}=NodiSort;
centri{3}=NodiSort;

count = zeros(length(centri{1}),length(centri{2}),length(centri{3}));
for i = 1:length(NodiSort)
    for j = 1:length(NodiSort)
        chiNodi = intersect(find(X(:,2)==NodiSort(i)),find(X(:,3)==NodiSort(j)));
        if isempty(chiNodi)==0
            Cap = X(chiNodi,1);
            for h = 2:length(Cap_bin)
                chiCap = length(intersect(find(Cap>Cap_bin(h-1)),find(Cap<=Cap_bin(h))));
                if isempty(chiCap)==0
                    count(h-1,i,j)=chiCap;
                end
            end
        end
    end
end

