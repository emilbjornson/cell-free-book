function SE = functionDownlinkSE_maxmin_dist(bk,Ck,preLogFactor,L,K,D,rhomax)
%Compute downlink SE according to Corollary 6.3 with max-min power allocation
%in Algorithm 7.5 
%
%INPUT:
%bk             = Matrix with dimension L x K where (1:La,k) is the non-zero
%                 portion of \tilde{b}_k in (7.25) (La is the number of
%                 APs serving UE k)
%Ck             = Matrix with dimension L x L x K x K where (1:La,1:La,k,i) is the 
%                 non-zero portion of \tilde{C}_{ki} in (7.26) (La is the
%                 number of APs serving UE i)
%preLogFactor   = Pre-log factor in the SE expression of Corollary 6.3
%L              = Number of APs
%K              = Number of UEs in the network
%rhomax         = Maximum allowable AP transmit power
%
%OUTPUT:
%SE             = SEs achieved with max-min power allocation algorithm
%
%
%This Matlab function was developed to generate simulation results to:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO", 
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.

%Implement Algorithm 7.5

%Prepare array to store SEs
SE = zeros(K,1);

%Solution accuracy of bisection search
delta = 0.1;
tLow = 0; %Line 3 of Algorithm 7.5
%Compute the upper bound on t according to Line 4 of Algorithm 7.5
tUpp = inf; 
for k = 1:K
    servingAPs = find(D(:,k)==1);
    La = length(servingAPs);
    tUpp = min(tUpp, La*rhomax*norm(bk(1:La,k))^2);
    
end

%Initialize the optimization variables
tilrho = zeros(L,K);

%Prepare array to store the number of APs serving a specficic UE
La = zeros(K,1);
%Prepare cell to store the AP indices serving a specficic UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specficic UE
NoServ = cell(K,1);

%Construc the above array and cells
for k = 1:K
    servingAPs = find(D(:,k)==1);
    NoservingAPs = find(D(:,k)==0);
    
    Serv{k} = servingAPs;
    NoServ{k} = NoservingAPs;
    
    La(k) = length(servingAPs);
end

%Compute the concatenated matrix whose block-diagonals are square-roots of
%the non-zero portions of of the matrices \tilde{C}_{ki} in (7.30), and 
%the concatenated vectors whose portions are the non-zero elements of
%\tilde{b}_k in (7.30)
Ck2 = zeros(sum(La),sum(La),K);
bk2 = zeros(sum(La),K);
for k = 1:K
    bk2(sum(La(1:k-1))+1:sum(La(1:k)),k) = bk(1:La(k),k);
    for i = 1:K
        
        Ck2(sum(La(1:i-1))+1:sum(La(1:i)),sum(La(1:i-1))+1:sum(La(1:i)),k) = sqrtm(Ck(1:La(i),1:La(i),k,i));
    end
end


while tUpp-tLow > delta %Condition in Line 6
    
    tMid = (tLow + tUpp)/2; %Line 7
    S = sqrt(1+1/tMid); %The constant term appearing in the constraints of (7.30)
    
    %Solve the problem in (7.30) with CVX
    cvx_begin quiet
    variable rho(L,K)
    variable rho2(sum(La),1) %non-zero portions of the main optimization variable rho
    minimize norm(rho2)
    subject to
    
    rho2 >= zeros(sum(La),1);
    for k=1:K
        norm([Ck2(:,:,k)*rho2; 1]) <= S*bk2(:,k)'*rho2;
        
        rho(Serv{k},k) == rho2(sum(La(1:k-1))+1:sum(La(1:k)),1);
        rho(NoServ{k},k) == zeros(length(NoServ{k}),1);
    end
    
    for l = 1:L
        norm(rho(l,:)) <= sqrt(rhomax);
    end
    
    cvx_end
    if cvx_status(1)=='S' %Line 9
        tLow = tMid; %Line 10
        tilrho = rho; %Line 11
    else %Line 12
        tUpp = tMid; %Line 13
    end
end

%Go through all UEs
for k = 1:K
    
    %Compute the numerator and denominator of (7.23)
    numm = abs(bk(1:La(k),k)'*tilrho(Serv{k},k))^2;
    denomm = 1-numm;
    
    for i = 1:K
        
        denomm = denomm+tilrho(Serv{i},i)'*Ck(1:La(i),1:La(i),k,i)*tilrho(Serv{i},i);
       
    end
    
    %Compute SE using Corollary 6.3 and SINRs in (7.23)
    SE(k) = preLogFactor*log2(1+numm/denomm);
    
end