function SE = functionDownlinkSE_maxmin(bk,ck,wk,preLogFactor,K,rhomax)
%Compute downlink SE according to Theorem 6.1 with max-min power allocation
%in Algorithm 7.3 
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

%Implement Algorithm 7.3

%Initialize iteration counter to zero
iter = 0;
%Initialize the power coefficients
eta = rhomax/K*ones(K,1);
%Compute the denominator in (7.12) for all UEs
denominator = ck'*eta+1;
%Compute SINRs in (7.12) for all UEs
SINR = eta.*bk./denominator;

while max(SINR)-min(SINR)>0.01  %the condition in Line 2 of Algorithm 7.3 with solution accuracy 0.01
    %Increase iteration counter by one
    iter = iter+1;
    
    
    eta = denominator./bk; %Line 3 of Algorithm 7.3
    maxx = max(wk*eta); %Denominator in Line 4 of Algorithm 7.3
    eta = eta*rhomax/maxx; %Line 4 of Algorithö 7.3
    
    %Update the denominator in (7.12) for all UEs
    denominator = ck'*eta+1;
    %Update SINRs in (7.12) for all UEs
    SINR = eta.*bk./denominator;
   
end
%Compute SEs
SE = preLogFactor*log2(1+SINR);

