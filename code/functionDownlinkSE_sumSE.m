function SE = functionDownlinkSE_sumSE(bk,ck,wk,preLogFactor,L,K,rhomax,tau_p)
%Compute downlink SE according to Theorem 6.1 with sum SE maximizing 
%power allocation in Algorithm 7.4 
%
%INPUT:
%bk             = Vector of length K with elements \tilde{b}_k in (7.13)
%ck             = Matrix with dimension K x K where (i,k) is the 
%                 \tilde{c}_{ki} in (7.14)-(7.15)
%wk             = Matrix with dimension L x K where (l,k) is the
%                 expected value of the norm square of the portion
%                 of the normalized centralized transmit precoder of UE k
%                 corresponding to AP l in (7.17)
%preLogFactor   = Pre-log factor in the SE expression of Theorem 6.1
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

%Implement Algorithm 7.4

%Initialize the current objective function
objec_new = inf;
%Initialize the power allocation coefficients randomly
eta = rhomax/tau_p*rand(K,1);

%Prepare arrays to store e_k, u_k, and d_k in Algoirthm 7.4
eee = zeros(K,1);
uuu = zeros(K,1);
ddd = zeros(K,1);
%Intialize the difference between current and previous objective values 
diff = 100;

%Initizalize the iteration counter to zero
iterr = 0;

%Go through the algorithm steps if the objective function is improved
%more than 0.2 or not improved at all
while (diff>0.2) || (diff<0)
    %Increase iteration counter by one
    iterr = iterr+1;
    %Update the previous objective value by the current objective value
    objec_old = objec_new;
    
    %Go through all UEs
    for k = 1:K
        %Compute the numerator and denominator in Line 4 of Algorithm 7.4
        numm = sqrt(eta(k,1)*bk(k,1));
        denomm = ck(:,k)'*eta+1+eta(k,1)*bk(k,1);
        
        %Update u_k according to Line 4
        uuu(k,1) = numm/denomm;
        %Update e_k and d_k as in Line 5
        eee(k,1) = 1-abs(numm)^2/denomm;
        ddd(k,1) = 1/eee(k,1);
        
    end
    
    %Solve the convex problem in (7.21) with CVX
    cvx_begin quiet
    variable rho(K,1) %square root of the power allocation coefficients
    variable sss(K,1)
    minimize sum(sss)
    subject to
    
    
    for k = 1:K
        quad_form(rho(k,1),ddd(k,1)*uuu(k,1)^2*bk(k,1))...
            +quad_form(rho,ddd(k,1)*uuu(k,1)^2*diag(ck(:,k)))...
            -2*ddd(k,1)*uuu(k,1)*sqrt(bk(k,1))*rho(k,1)<=sss(k,1);
    end
    for l = 1:L
        sum_square(sqrt(wk(l,:).').*rho)<=rhomax;
        
    end
    
    
    cvx_end
    
    %Compute the power allocation coefficients by using their square roots
    %obtained by CVX
    eta = rho.^2;
    %Update the current objective value
    objec_new = sum(ddd.*eee-log(ddd));
    %Obtain the difference between current and previous objective values
    diff = objec_old - objec_new;
    
end
%Compute SEs
SE = preLogFactor*log2(ddd);