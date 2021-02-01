function SE = functionUplinkSE_maxmin(bk,ck,sigma2,preLogFactor,K,pmax)
%Compute uplink SE according to Theorem 5.2 with max-min power control 
%in Algorithm 7.1 
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
%SE             = SEs achieved with max-min power control algorithm
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

%Implement Algorithm 7.1

%Initialize iteration counter to zero
iter = 0;
%Initialize the power control coefficients to full power
eta = pmax*ones(K,1);
%Compute the denominator in (7.1) for all UEs
denominator = ck'*eta+sigma2;
%Compute SINRs in (7.1) for all UEs
SINR = eta.*bk./denominator;

while max(SINR)-min(SINR)>0.01 %the condition in Line 2 of Algorithm 7.1 with solution accuracy 0.01
    %Increase iteration counter by one
    iter = iter+1;
    
    eta = denominator./bk; %Line 3 of Algorithm 7.1
    eta = eta*pmax/max(eta); %Line 4 of Algorithm 7.1
    
    %Update the denominator in (7.1) for all UEs
    denominator = ck'*eta+sigma2;
    %Update SINRs in (7.1) for all UEs
    SINR = eta.*bk./denominator;
   
end
%Compute SEs
SE = preLogFactor*log2(1+SINR);

