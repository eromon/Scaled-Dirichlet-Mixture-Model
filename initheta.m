function [alpha,alphasum,beta] = initheta(mean_data ,var_data,K,data,D);
for j = 1:K
    %Added an edit on the negtive mean square
    alphasum(j,:)= (((-(mean_data(j,:).^2)) + mean_data(j,:))/var_data(j,:))-1; 
    alpha(j,:) = (mean_data(j,:) .* alphasum(j,:))./sum(mean_data); 
end
beta = ones(j,size(data,2));

end

