%Data is a N by D matrix
% alpha and beta are 1 by D vectors
%posterior is a 1 by N vector
%Hessian is a 2D by 2D matrix
%Gradient is a D by 1 matrix
%Our components are J
%Code is written by Eromonsele Oboh on the 1st of February 2016


function [H,Hiv] = ScaDirHes(post,alpha ,beta,data,D)

H11 = sum(post).*(psi(1,(sum(alpha)))- diag(psi(1,alpha)));    

                        %aa = repmat((sum(data+eps)./(beta*(sum(data+eps)'))),D,1);
                        % aa = repmat(((data+eps)'*post)'./(sum(post).*sum(beta*(data+eps)')),D,1);
%aa = sum(post).*repmat(sum(((data+eps)'),2)./sum(beta*(data+eps)'),1,D);
aa = repmat((((data+eps)'*post)./((beta)*sum(data+eps)'))',D,1);
%                         % bb = diag(1./beta);
%                         % bb = sum(post)*diag(1./beta); 
bb = diag(sum(post)./(beta));
H12 = bb-aa;
H12 = H12*0;
 %H12 = (aa+bb -(2*(diag(diag(aa)))));
                        %H12 = sum(post).*((diag(1./beta))-repmat(sum(data+eps)./(beta*sum(data+eps)'),D,1));

H21 = H12';
H21 = H21*0;
                        %cc = ((sum(data+eps)'*sum(data+eps)*(sum(alpha)))./((beta*sum(data+eps)').^2));
                        % cc = ((((data+eps)'*post)*sum(data+eps))*(sum(alpha)))./(sum(post).*(sum(beta*((data+eps)')).^2));
                        %cc = (sum(post).*((((data+eps)')*(data+eps))*(sum(alpha))))./(sum(beta*((data+eps)')).^2);
%  cc =((((data+eps)')*(data+eps))*sum(alpha))./((((beta)*sum(data+eps)').^2));
 cc =((((data+eps)')*(data+eps))*sum(alpha))./sum((((beta)*(data+eps)').^2),2);
%cc =repmat((((((data+eps)').^2)*post)./((beta*sum(data+eps)').^2)),1,D);
                        % dd = (-1*diag(alpha./(beta.^2)));
                        % dd = (-1*diag(sum(post).*(alpha./(beta.^2))));

dd = (diag((sum(post).*(alpha))./(beta.^2)));
H22 = cc-dd;
                  %H22 = (cc+dd -(2*(diag(diag(cc)))));
                  %H22 = sum(post).*(diag((alpha./beta.^2))-((sum(data+eps)'*sum(data+eps)*(sum(alpha)))./(beta*sum(data+eps)'.^2)));

H = [H11,H12;H21,H22];
%H = diag(diag(H));

Hiv = pinv(H);

%A is invertible formula
% Hiv11  = pinv(H11)+(pinv(H11)*H12)*pinv(H22-(H21*pinv(H11)*H12))*H21*pinv(H11);
% 
% Hiv12 = -(pinv(H11)*H12*pinv(H22-(H21*pinv(H11)*H12)));
% 
% Hiv21 = -(pinv(H22-(H21*pinv(H11)*H12))*H21*pinv(H11));
% 
% Hiv22 = pinv(H22-(H21*pinv(H11)*H12));
% 
% M= [Hiv11,Hiv12;Hiv21,Hiv22];
% Hiv = diag(diag(M));

%D is invertible formula
% Hiv11 = pinv(H11-(H12*pinv(H22))*H21);
% Hiv12 = -(pinv(H11)*H12*pinv(H22-(H21*pinv(H11))*H12));
% Hiv21 = -(pinv(H22)*H21*(H11-(H12*pinv(H22))*H21));
% Hiv22 = pinv(H22-(H21*pinv(H11))*H12);






%H11 = sum(post1).*(psi(2,(sum(alphaJ1)))- diag(psi(2,alphaJ1)));
%H12 = sum(post1).*((diag(1./betaJ1))+(-(repmat((sum(data)),3,1)./sum(betaJ1.*sum(data)))));
%H21 = H12';
%H22 = sum(post1).*(((sum(data)'*sum(data)*psi(2,(sum(alphaJ1))))./sum(betaJ1.*sum(data)))-diag((alphaJ1./betaJ1.^2)));

%H = [H11,H12;H21,H22];
%Hiv = pinv(H);







%D = size(data,2);
%DD = 2*D;
%kk = [zeros(306,3),(1-data)];
%rr = [zeros(3,306);(1-data)'];
%tt = [zeros(1,3),(-alphaJ1+sum(alphaJ1)./(betaJ1.^2))];
%zz = [zeros(1,3),betaJ1];
%yy = [zeros(3,1);betaJ1'];
%H = zeros(6,6);

%for i = 1:(DD);   
 %   for j =1:(DD);
        
  %      if (i<=D)&&(j<=D);
   %         H(i,j)= sum((post1*psi(2,(sum(alphaJ1))))) ;
    %    elseif (i<=D)&&(j<=DD);
            %H(i,j)= sum(post1'.*kk)./betaJ1;
            
     %       H(i,j)= sum(post1'.*kk(:,j))./zz(:,j); 
      %  elseif (i<=DD)&&(j<=D);
            %H(i,j)= sum(post1.*rr)./yy; 
       %     H(i,j)= sum(post1.*rr(i,:))./yy(i,:);  
        %else
        % H(i,j)= sum(post1*tt(:,j));
        %end
    %end   
     
%end

