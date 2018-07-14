function [h_prior] = priorMML(K,D,alpha,beta);

%To calculate prior of Mixing Weight
for i =K-1
    if i == 0
        h=log(realmin);
    else
    h(i)=log(i);
    end
end
%Prior of Mixing Weight
HP = sum(h);

%To calculate Prior of Alpha parameter
for j = 1:K
    A(j) = log(norm(alpha(j,:)));
end
%AA is for the norm of alpha j vector
AA = D*sum(A);

%BB is for each alpha parameter
BB = sum(sum(log(alpha),2));



% CC is for the beta parameter
 %JJ = K*D*log(1);
for j = 1:K
    C(j) = log(norm(beta(j,:)));
end
CC = D*sum(C);

DD = sum(sum(log(beta),2));

% (6*K*D) is for the Exponent
h_prior = HP-(6*K*D)-AA+BB-(6*K*D)-CC+DD;


    