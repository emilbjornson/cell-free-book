function SE = functionUplinkSE_sumSE(bk,ck,sigma2,preLogFactor,K,pmax)
%Compute uplink SE according to Theorem 5.2 with sum SE maximizing 
%power control in Algorithm 7.2 
%
%INPUT:
%bk             = Vector of length K with elements b_k in (7.2)
%ck             = Matrix with dimension K x K where (i,k) is the c_{ki} in
%                 (7.3)-(7.4)
%sigma2         = Vector of length K with elements \sigma_k^2 in (7.5)
%preLogFactor   = Pre-log factor in the SE expression of Theorem 5.2
%K              = Number of UEs in the network
%pmax           = Maximum allowable UE transmit power
%
%OUTPUT:
%SE             = SEs achieved with sum SE maximizing power control algorithm
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

%Implement Algorithm 7.2

%Initialize the current objective function
objec_new = inf;
%Initialize the power control coefficients randomly
eta = pmax*rand(K,1);
%Prepare arrays to store e_k, u_k, and d_k in Algoirthm 7.2
eee = zeros(K,1);
uuu = zeros(K,1);
ddd = zeros(K,1);
%Intialize the difference between current and previous objective values 
diff = 100;

%Initizalize the iteration counter to zero
iterr = 0;
%Go through the algorithm steps if the objective function is improved
%more than 0.0001 or not improved at all
while (diff>0.0001) || (diff<0) 
    %Increase iteration counter by one
    iterr = iterr+1;
    %Update the previous objective value by the current objective value
    objec_old = objec_new;
    
    %Go through all UEs
    for k = 1:K
        %Compute the numerator and denominator in Line 4 of Algorithm 7.2
        numm = sqrt(eta(k,1)*bk(k,1));
        denomm = ck(:,k)'*eta+sigma2(k,1)+eta(k,1)*bk(k,1);
        
        %Update u_k according to Line 4
        uuu(k,1) = numm/denomm;
        %Update e_k and d_k as in Line 5
        eee(k,1) = 1-abs(numm)^2/denomm;
        ddd(k,1) = 1/eee(k,1);
        
    end
    
    %Go through all UEs
    for k = 1:K
        %Compute p_k as in Line 6
        numm = bk(k,1)*ddd(k,1)*uuu(k,1)^2;
        denomm = numm;
        for i = 1:K
            denomm = denomm+ddd(i,1)*abs(uuu(i,1))^2*ck(k,i);
        end
        
        eta(k,1) = min(pmax, ddd(k,1)*numm/(denomm^2));
    end
    
    %Update the current objective value
    objec_new = sum(ddd.*eee-log(ddd));
    %Obtain the difference between current and previous objective values
    diff = objec_old - objec_new;
    
end

%Compute SEs
SE = preLogFactor*log2(ddd);