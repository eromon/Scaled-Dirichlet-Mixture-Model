
%Initialization
%preprocessing
clear;
%1. load data dimension (N by D)
        %load('haberdata.mat');
        %data = haberdata./repmat(sum(haberdata,2),1,3);

%data1mix is a 1 mixture dirichlet distributed data with parameters 8 2 8 

load('synthetic4.mat')
%data3mix is a 3 mixture dirichlet, para 16-65-30;65-15-30;30-34-35

D = size(data,2);
N = size(data,1);
K = 4;
maxIter =5;
%MML Criterion Starts From Here
%for K = 1:K
  rng('shuffle'); 
% implement Kmeans
assignment = kmeans(data+eps,K);

% 2a. Initializing the moments
[mean_data ,var_data] = initmoments(assignment,data,K);


% 2b. Calculating Mixing Weight
pj = mixweight(assignment,K,data);

%2c. Initializing theta
[alpha,alphasum,beta] = initheta(mean_data ,var_data,K,data,D);

post=zeros(N,K);
pdf=zeros(N,K);

alphaprev = zeros(K,D);
betaprev = zeros(K,D);
% for (iter =1:2000)
%figure;
LLK=[];
%while 1
 for iter = 1:100
    fprintf(1,'%d',iter);

%3. Expectation  Step
%for each cluster
    for j = 1 : K
    %Evaluate pdf for all data points belonging to cluster J
        pdf(:,j) = computepdf(alpha(j,:) ,beta(j,:),data);
    end

%Multiply each pdf value by its prior probability according to cluster
            %pdf_weight = pdf + repmat(log(pj),306,1);
%pdf_weight = pdf.*repmat(pj,N,1); %bsxfun(@times, pdf, pj);
        pdf_weight = repmat(log(pj),N,1)+pdf;

%To calculate the posterior, we didive the weighted pdf by its sum
                    %post = pdf_weight - repmat((sum(pdf_weight,2)),1,3);
                    %post = bsxfun(@rdivide, pdf_weight, sum(pdf_weight, 2));
%post = pdf_weight ./repmat((sum(pdf_weight,2)),1,K);
        post =exp(pdf_weight);
        post(post==inf)=realmax;
        post(post==0)=realmin;

        post = post ./ repmat((sum(post,2)),1,K);

%4. Maximization Step
%Updating parameters
        alphaprev = alpha;
        betaprev= beta;
    for j = 1:K 
    % Calculate the prior probability for cluster 'j'.
        pj(:,j) = sum(post(:,j))./size(data,1);  %mean(post(:, j), 1);

        GG=[];
%         cost = 0;
%         counter(j) = 0;
%         while 1
%              counter(j)=counter(j)+1;
         for counter = 1:maxIter  
             thetaold(j,:)= [alpha(j,:),beta(j,:)];
             G(:,:,j) = ScaDirGra(post(:,j),alpha(j,:) ,beta(j,:),data);
             %[H(:,:,j),Hiv(:,:,j)] = ScaDirHes(post(:,j),alpha(j,:) ,beta(j,:),data,D);
        
             %cost = (Hiv(:,:,j)*G(:,:,j))'; 
             thetanew(j,:) = [alpha(j,:),beta(j,:)] - [.000001*G(:,:,j)'];
             thetanew(thetanew<=0)= .00001;
             thetanew(thetanew>=realmax)=10000;  
%              absError(j) = norm(thetanew-thetaold)/norm(thetanew);
       
             alpha(j,:) = thetanew(j,1:D);
             beta(j,:) = ones(1,size(data,2));
%this code is to normalise beta to be less than %
             %beta(j,:) = thetanew(j,D+1:2*D);
%                   if (counter(j)>maxIter);
%                        haveWeFoundSolution = false;
%                        break
%                   end
%                   if (absError(j)<=1e-10);
%                        haveWeFoundSolution = true;
%                        break
%                   end       


            Gnorm = norm(G(:,:,j));
            GG=[GG,Gnorm];
        end
%         plot(GG,'b.')
%         if (haveWeFoundSolution)
%             fprintf('Converging\n');
%         else
%             fprintf('Not Converging\n');
%         end

 %Just Added 7/4/16 This is to discard pj and go to Expectation step

   end 


    for j = 1:K
    newlog(j,:) = sum(repmat(log(pj(:,j)),N,1)+(computepdf(alpha(j,:) ,beta(j,:),data)));
    oldlog(j,:) = sum(repmat(log(pj(:,j)),N,1)+(computepdf(alphaprev(j,:) ,betaprev(j,:),data)));
    end
    LLK=[LLK,sum(newlog)];

    if ((sum(newlog)- sum(oldlog)<= 1e-5)&&((norm(alpha)-norm(alphaprev))<=1e-5));
        
        break
    end

end
LLK;
% figure;
 plot(LLK,'r.')

%Post_final assigns the data objects with its class labels
 %[value(:,K), Post_final(:,K)]= newp(post,N);
 
[value, Post_final]= newp(post,N);

%Calculating Determinant of Fisher
% nj = num_clus(Post_final(:,K),K);
nj= num_clus(Post_final,K);

%Fdet(K) = FisherHes(nj,alpha ,beta,data,D,K,N,pj);
Fdet = FisherHes(nj,alpha ,beta,data,D,K,N,pj);

%Calculating prior
 %h_prior(K) = priorMML(K,D,alpha,beta);
h_prior = priorMML(K,D,alpha,beta);

%Calculate log pdf
 %pdf_final(K) = sum(sum(pdf));

% for j = 1:K
%    pdf_f(j)= sum(computepdf(alpha(j,:) ,beta(j,:),data));
% end
%  pdf_final = sum(pdf_f);
%pdf_final = log(sum(sum(exp(pdf))));
%pdf_final = sum(sum(pdf)); %made this change 5/26/16
%pdf_final = sum(sum(log(post)+pdf_weight));
pdf_final = sum(sum(post.*pdf_weight));


%Calculate MML constant
 %NP_final(K) = ((K*(2*D-1))/2)*(1-log(12));
NP_final = ((K*(2*D-1))/2)*(1-log(12));

%Calculate MML
 %MML(K) = -h_prior(K)-pdf_final(K)+ (.5*log(Fdet(K)))+NP_final(K);
MML(K) = -h_prior-pdf_final+(.5*log(Fdet))+NP_final;
%MML(K) = -pdf_final+NNN;

%    end  
% %  figure
%   plot(MML)
%plot(www)     
    
    
    

    