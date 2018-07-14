
function pj = mixweight(assignment,K,data)

for j = 1:K 
    pj(:,j) = length(find(assignment==j))/length(data); 
end


