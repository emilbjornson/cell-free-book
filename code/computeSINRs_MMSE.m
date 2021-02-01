function SINRs = computeSINRs_MMSE(channel)
%Compute the SINR when using MMSE combining defined in (1.18) for Cellular Massive MIMO 
%and in (1.25) for cell-free setup.
%
%INPUT:
%channel        = M x K channel matrix where the kth column corresponds to
%                 kth UE's channel, i.e., g_k in (1.18) for Cellular Massive MIMO 
%                 or h_k in (1.25) for cell-free setup. 
%
%OUTPUT:
%SINRs           = K x 1 UE SINRs in (1.17) for Cellular Massive MIMO or
%                  in (1.24) for cell-free setup.
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


%Extract number of antennas
M = size(channel,1);

%Extract number of UEs
K = size(channel,2);

SINRs = zeros(K,1);

for k = 1:K
    
    %Compute the SINR of the kth UE with the MMSE combining vector
    desiredChannel = channel(:,k);
    
    SINRs(k) = desiredChannel'*((channel*channel' - desiredChannel*desiredChannel' + eye(M) )\desiredChannel);
    
end