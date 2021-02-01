function SE_MMSE = functionComputeSE_cellular_uplink(Hhat,C,BSassignment,tau_c,tau_p,nbrOfRealizations,M,K,nbrBSs,p)
%This function computes the uplink SE using L-MMSE combining for the
%cellular Massive MIMO system.
%
%INPUT:
%Hhat              = Matrix with dimension M*nbrBSs x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel from
%                    all cellular BSs to UE k at channel realization n.
%C                 = Matrix with dimension M x M x nbrBSs x K where (:,:,l,k) is the
%                    spatial correlation matrix of the channel estimation error
%                    between cellular BS l and UE k, normalized by noise variance
%BSassignment      = Matrix with dimension K x 1 containing the
%                    index of the BS that serves a particular UE
%tau_c             = Length of the coherence block
%tau_p             = Length of pilot sequences and number of UEs per cell
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Total number of UEs
%nbrBSs            = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%SE_MMSE           = K x 1 matrix where the kth element is the uplink SE of UE k
%                    achieved with L-MMSE combining 
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

%Store identity matrix of size M
eyeM = eye(M);

%Compute the pre-log factor (normalized with number of channel realizations)
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c*nbrOfRealizations);

%Prepare to store simulation results
SE_MMSE = zeros(K,1);

%Compute sum of the estimation error correlation matrices at every BS
C_tot = sum(C,4);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all BSs
    for l = 1:nbrBSs
        
        %Extract channel estimate realizations from all UEs to BS l
        Hhatallj = reshape(Hhat(1+(l-1)*M:l*M,n,:),[M K]);
        
        %Extract which UEs are served by the BS
        servedUEs = find(BSassignment==l);
        
        
        %Compute MR combining
        V_MR = Hhatallj(:,servedUEs);
        
        
        %Compute L-MMSE combining in (5.29)
        V_MMMSE = p*((p*(Hhatallj*Hhatallj')+p*C_tot(:,:,l)+eyeM)\V_MR);
        
        
        
        %Go through all UEs in cell l
        for ind = 1:tau_p
            
            %Extract UE index
            k = servedUEs(ind);
            
            
            %L-MMSE combining
            v = V_MMMSE(:,ind); %Extract combining vector
            
            %Compute numerator and denominator of the effective SINR in
            %(5.45)
            numerator = p*abs(v'*Hhatallj(:,k))^2;
            denominator = p*norm(v'*Hhatallj)^2 + v'*(p*C_tot(:,:,l)+eyeM)*v - numerator;
            
            %Compute instantaneous SE for one channel realization and
            %update SE by Monte-Carlo estimation according to (5.44)
            SE_MMSE(k) = SE_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator));
            
            
            
        end
        
    end
    
end
