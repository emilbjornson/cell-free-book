function SE = functionDownlinkSE_sumSE_dist(bk,Ck,preLogFactor,L,K,D,rhomax,tau_p)
%Compute downlink SE according to Corollary 6.3 with sum SE maximizing power allocation
%in Algorithm 7.6
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
%tau_p          = Length of pilot sequences
%
%OUTPUT:
%SE             = SEs achieved with sum SE maximizing power allocation algorithm
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

%Implement Algorithm 7.6

%Initialize the current objective function
objec_new = inf;

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
%the non-zero portions of of the matrices \tilde{C}_{ki} in (7.33), and 
%the concatenated vectors whose portions are the non-zero elements of
%\tilde{b}_k in (7.33)
Ck2 = zeros(sum(La),sum(La),K);
bk2 = zeros(sum(La),K);
for k = 1:K
    bk2(sum(La(1:k-1))+1:sum(La(1:k)),k) = bk(1:La(k),k);
    for i = 1:K
        
        Ck2(sum(La(1:i-1))+1:sum(La(1:i)),sum(La(1:i-1))+1:sum(La(1:i)),k) = Ck(1:La(i),1:La(i),k,i);
    end
end

%Initialize the power allocation coefficients randomly (\tilde{rho})
rho = sqrt(rhomax/tau_p)*rand(sum(La),1);

%Prepare arrays to store e_k, u_k, and d_k in Algoirthm 7.6
eee = zeros(K,1);
uuu = zeros(K,1);
ddd = zeros(K,1);
%Intialize the difference between current and previous objective values 
diff = 100;

%Initizalize the iteration counter to zero
iterr = 0;
%Go through the algorithm steps if the objective function is improved
%more than 0.1 or not improved at all
while (diff>0.1) || (diff<0)
    %Increase iteration counter by one
    iterr = iterr+1;
    %Update the previous objective value by the current objective value
    objec_old = objec_new;
    
    %Go through all UEs
    for k = 1:K
        %Compute the numerator and denominator in Line 4 of Algorithm 7.6
        numm = bk2(:,k)'*rho;
        denomm = 1+rho'*Ck2(:,:,k)*rho;
        
        %Update u_k according to Line 4
        uuu(k,1) = numm/denomm;
        %Update e_k and d_k as in Line 5
        eee(k,1) = 1-abs(numm)^2/denomm;
        ddd(k,1) = 1/eee(k,1);
        
    end
    
    %Solve the convex problem in (7.33) with CVX
    cvx_begin quiet
    variable rho3(L,K)
    variable rho2(sum(La),1) %non-zero portions of the main optimization variable rho3
    variable sss(K,1)
    minimize sum(sss)
    subject to
    
    for k=1:K
            quad_form(rho2,ddd(k,1)*uuu(k,1)^2*Ck2(:,:,k))...
            -2*ddd(k,1)*uuu(k,1)*bk2(:,k)'*rho2<=sss(k,1);
        
        rho3(Serv{k},k) == rho2(sum(La(1:k-1))+1:sum(La(1:k)),1);
        rho3(NoServ{k},k) == zeros(length(NoServ{k}),1);
    end
    for l = 1:L
        norm(rho3(l,:)) <= sqrt(rhomax);
    end
    
    rho2 >= zeros(sum(La),1);
    
    cvx_end
    
    %Update the power allocation coefficients 
    %obtained by CVX
    rho = rho2;
    %Update the current objective value
    objec_new = sum(ddd.*eee-log(ddd));
    %Obtain the difference between current and previous objective values
    diff = objec_old - objec_new;
    
end
%Compute SEs
SE = preLogFactor*log2(ddd);