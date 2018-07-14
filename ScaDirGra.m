%Gradient is a D by 1 matrix
%Our components are J
%Code is written by Eromonsele Oboh on the 1st of February 2016

%Code is wriiten in normal space and no log space

function G = ScaDirGra(post,alpha ,beta,data)
    %aa = (psi(1,sum(alpha))-psi(1,alpha))+log(beta);
    %bb = repmat(aa,N,1)+ log(data+eps) -
    %repmat(log(beta(1,:)*(data+eps)')',1,D);
    %%a1(sum((log(beta*((data+eps)'))),2)*sum(post))   a2%%((log(beta*((data+eps)'))*post)); a3((log(beta*sum(((data+eps)'),2)))*sum(post))
    %G11 = sum(post).* sum(bb);
    %G11 = sum(post).*((psi(1,sum(alpha))-psi(1,alpha))+((log(beta)+log(sum(data+eps)))-log(beta*sum(data+eps)')));

    %G11 = ((sum(post).*psi(0,sum(alpha)))-((sum(post).*psi(0,alpha))+(sum(post).*log(beta))))+(((log(data+eps)')*post)')-(sum(log((beta)*((data+eps)')))*sum(post));   

    G11 = (sum(post).*(psi(0,sum(alpha))-psi(0,alpha)))+((sum(post).*log(beta))+(((log(data+eps)')*post)'))-((log(beta*(data+eps)'))*post); 
    G12 = ((sum(post).*alpha)./beta)-(((sum(alpha).*(data+eps)')*post)'./(sum(data+eps)*beta'));
    G = [G11,G12]';
end

%(sum(log((beta)*((data+eps)')))*sum(post))

%(log((beta)*((data+eps)'))*post)

