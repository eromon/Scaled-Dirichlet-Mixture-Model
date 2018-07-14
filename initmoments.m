function [mean_data ,var_data] = initmoments(assignment,data,K);
for j = 1:K
    mean_data(j,:) = mean(data(find(assignment==j),:));
    var_data(j,:) = var(data(find(assignment==j),:));
end
end

    


