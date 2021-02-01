function [SE_MMSE, SE_P_MMSE, SE_P_RZF,  ...
    SE_L_MMSE, SE_LP_MMSE, SE_MR, ...
    Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR] ...
    = functionComputeSE_downlink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex,rho_dist,gainOverNoisedB,rho_tot)
%Compute downlink SE for different transmit precoding schemes using the capacity
%bound in Theorem 6.1 for the centralized schemes and the capacity bound
%in Corollary 6.3 for the distributed schemes. Compute the genie-aided
%downlink SE from Corollary 6.6 for the centralized and the distributed operations. 
%
%INPUT:
%Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
%H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                    with the true channel realizations. The matrix is
%                    organized in the same way as Hhat.
%D                 = DCC matrix for cell-free setup with dimension L x K 
%                    where (l,k) is one if AP l serves UE k and zero otherwise
%B                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel estimate 
%                    between AP l and UE k, normalized by noise variance
%C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel
%                    estimation error between AP l and UE k,
%                    normalized by noise variance
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Number of UEs 
%L                 = Number of APs
%p                 = Uplink transmit power per UE (same for everyone)
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k,
%                    normalized by noise
%pilotIndex        = Vector containing the pilot assigned to each UE
%rho_dist          = Matrix with dimension L x K where (l,k) is the power
%                    allocated to UE k by AP l in the distributed downlink
%                    operation
%gainOverNoisedB   = Matrix with dimension L x K where (l,k) is the channel
%                    gain (normalized by the noise variance) between AP l
%                    and UE k
%rho_tot           = Maximum allowed transmit power for each AP 
%
%OUTPUT:
%SE_MMSE           = SEs achieved with MMSE precoding in (6.16)
%SE_P_MMSE         = SEs achieved with P-MMSE precoding in (6.17)
%SE_P_RZF          = SEs achieved with P-RZF precoding in (6.18)
%SE_L_MMSE         = SEs achieved with L-MMSE precoding in (6.25)
%SE_LP_MMSE        = SEs achieved with LP-MMSE precoding in (6.33)
%SE_MR             = SEs achieved with MR precoding in (6.26)
%Gen_SE_P_MMSE     = Genie-aided SEs achieved with P-MMSE precoding in (6.17)
%Gen_SE_P_RZF      = Genie-aided SEs achieved with P-RZF precoding in (6.18)
%Gen_SE_LP_MMSE    = Genie-aided SEs achieved with LP-MMSE precoding in (6.33)
%Gen_SE_MR         = Genie-aided SEs achieved with MR precoding in (6.26)
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

%Store the N x N identity matrix
eyeN = eye(N);

%Compute the prelog factor assuming only downlink data transmission
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
Gen_SE_P_MMSE = zeros(K,1);
Gen_SE_P_RZF = zeros(K,1);
Gen_SE_LP_MMSE = zeros(K,1);
Gen_SE_MR = zeros(K,1);


%Prepare to store the terms that appear in SEs
signal_MR = zeros(K,1);
interf_MR = zeros(K,1);
cont_MR  = zeros(K,K);
scaling_MR = zeros(L,K);
interUserGains_MR = zeros(K,K,nbrOfRealizations);


signal_MMSE = zeros(K,1);
interf_MMSE = zeros(K,1);
scaling_MMSE = zeros(K,1);
portionScaling_MMSE = zeros(L,K);
interUserGains_MMSE = zeros(K,K,nbrOfRealizations);

signal_P_MMSE = zeros(K,1);
interf_P_MMSE = zeros(K,1);
scaling_P_MMSE = zeros(K,1);
portionScaling_PMMSE = zeros(L,K);
interUserGains_P_MMSE = zeros(K,K,nbrOfRealizations);

signal_P_RZF = zeros(K,1);
interf_P_RZF = zeros(K,1);
scaling_P_RZF = zeros(K,1);
portionScaling_PRZF = zeros(L,K);
interUserGains_P_RZF = zeros(K,K,nbrOfRealizations);

signal_L_MMSE = zeros(K,1);
interf_L_MMSE = zeros(K,1);
scaling_L_MMSE = zeros(L,K);
interUserGains_L_MMSE = zeros(K,K,nbrOfRealizations);

signal_LP_MMSE = zeros(K,1);
interf_LP_MMSE = zeros(K,1);
scaling_LP_MMSE = zeros(L,K);
interUserGains_LP_MMSE = zeros(K,K,nbrOfRealizations);


%% Compute scaling factors for precoding
%Computation for MR precoding
for l = 1:L
    
    %Extract which UEs are served by AP l
    servedUEs = find(D(l,:)==1);
    
    for ind = 1:length(servedUEs)
        
        %Compute scaling factor using the spatial correlation matrix of the
        %channel estimate
        scaling_MR(l,servedUEs(ind)) = trace(B(:,:,l,servedUEs(ind)));
        
    end
    
end

%Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %L-MMSE, LP-MMSE precoding
    for l = 1:L
        
        %Extract which UEs are served by the AP
        servedUEs = find(D(l,:)==1);
        
        %Compute sum of estimation error covariance matrices of the UEs
        %served by AP l
        Cserved = sum(C(:,:,l,servedUEs),4);
        
        %Compute L-MMSE and LP-MMSE precoding
        V_MR = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);
        V_L_MMSE = p*((p*(V_MR*V_MR')+p*sum(C(:,:,l,:),4)+eyeN)\V_MR(:,servedUEs));
        V_LP_MMSE = p*((p*(V_MR(:,servedUEs)*V_MR(:,servedUEs)')+p*Cserved+eyeN)\V_MR(:,servedUEs));
        
        %Compute scaling factor by Monte Carlo methods
        scaling_L_MMSE(l,servedUEs) = scaling_L_MMSE(l,servedUEs) + sum(abs(V_L_MMSE).^2,1)/nbrOfRealizations;
        scaling_LP_MMSE(l,servedUEs) = scaling_LP_MMSE(l,servedUEs) + sum(abs(V_LP_MMSE).^2,1)/nbrOfRealizations;
        
    end
    
    %MMSE, P-MMSE, and P-RZF precoding
    
    
    %Go through all UEs
    for k = 1:K
        
        %Determine the set of serving APs
        servingAPs = find(D(:,k)==1);
        La = length(servingAPs);
        
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;
        
        
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);
        
        for l = 1:La
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        
        %Compute MMSE, P-MMSE, and P-RZF precoding
        V_MMSE = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
        V_P_MMSE = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));
        V_P_RZF = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));
        
        %Compute scaling factor by Monte Carlo methods
        scaling_MMSE(k) = scaling_MMSE(k) + sum(abs(V_MMSE).^2,1)/nbrOfRealizations;
        scaling_P_MMSE(k) = scaling_P_MMSE(k) + sum(abs(V_P_MMSE).^2,1)/nbrOfRealizations;
        scaling_P_RZF(k) = scaling_P_RZF(k) + sum(abs(V_P_RZF).^2,1)/nbrOfRealizations;
        
        %Go through all the serving APs
        for l=1:La
            
            %Extract the portions of the centralized precoding vectors 
            V_MMSE2 = V_MMSE((l-1)*N+1:l*N,:);
            V_P_MMSE2 = V_P_MMSE((l-1)*N+1:l*N,:);
            V_P_RZF2 = V_P_RZF((l-1)*N+1:l*N,:);
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in Theorem 6.1
            
            portionScaling_MMSE(servingAPs(l),k) = portionScaling_MMSE(servingAPs(l),k) ...
                + sum(abs(V_MMSE2).^2,1)/nbrOfRealizations;
            
            portionScaling_PMMSE(servingAPs(l),k) = portionScaling_PMMSE(servingAPs(l),k) ...
                + sum(abs(V_P_MMSE2).^2,1)/nbrOfRealizations;
            
            portionScaling_PRZF(servingAPs(l),k) = portionScaling_PRZF(servingAPs(l),k) ...
                + sum(abs(V_P_RZF2).^2,1)/nbrOfRealizations;
        end
    end
    
end

%Normalize the norm squares of the portions for the normalized centralized precoders
portionScaling_MMSE = portionScaling_MMSE./repmat(scaling_MMSE.',[L 1]);

portionScaling_PMMSE = portionScaling_PMMSE./repmat(scaling_P_MMSE.',[L 1]);

portionScaling_PRZF = portionScaling_PRZF./repmat(scaling_P_RZF.',[L 1]);


%% Compute MR closed-form expectations

%Go through all APs
for l = 1:L
    
    %Extract which UEs are served by the AP
    servedUEs = find(D(l,:)==1);
    
    %Go through all UEs served by the AP
    for ind = 1:length(servedUEs)
        
        %Extract UE index
        k = servedUEs(ind);
        
        %Desired signal term in (6.27)
        signal_MR(k) = signal_MR(k) + sqrt(rho_dist(l,k)*real(trace(B(:,:,l,k))));
        
        
        for i = 1:K
            
            
            %Non-coherent interference from UE k to UE i (the first term of
            %(6.28))
            interf_MR(i) = interf_MR(i) + rho_dist(l,k)*real(trace(B(:,:,l,k)*R(:,:,l,i)))/real(trace(B(:,:,l,k)));
            
            if pilotIndex(k) == pilotIndex(i)
                
                %Coherent interference from UE k to UE i (the second term
                %of (6.28))
                cont_MR(i,k) = cont_MR(i,k) + sqrt(rho_dist(l,k))*real(trace((B(:,:,l,k)/R(:,:,l,k))*R(:,:,l,i)))/sqrt(real(trace(B(:,:,l,k))));
                
            end
            
        end
        
    end
    
end

%The parameters for the scalable centralized downlink power allocation in
%(7.43)
upsilon = -0.5;
kappa = 0.5;

%Compute the power allocation coefficients for centralized precoding
%according to (7.43)

rho_MMSE = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_MMSE,upsilon,kappa);

rho_PMMSE = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_PMMSE,upsilon,kappa);

rho_PRZF = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_PRZF,upsilon,kappa);

%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    
    
    %Matrix to store Monte-Carlo results for this realization
    interf_MMSE_n = zeros(K,K);
    interf_P_MMSE_n = zeros(K,K);
    interf_P_RZF_n = zeros(K,K);
    interf_L_MMSE_n = zeros(K,K);
    interf_LP_MMSE_n = zeros(K,K);
    
    %Go through all APs
    for l = 1:L
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H((l-1)*N+1:l*N,n,:),[N K]);
        
        %Extract channel estimates from all UEs to AP l
        Hhatallj = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute sum of estimation error covariance matrices of the UEs
        %served by AP l
        Cserved = sum(C(:,:,l,servedUEs),4);
        
        %Compute MR combining
        V_MR = Hhatallj(:,servedUEs);
        
        %Compute L-MMSE combining
        V_L_MMSE = p*((p*(Hhatallj*Hhatallj')+p*sum(C(:,:,l,:),4)+eyeN)\V_MR);
        
        %Compute LP-MMSE combining
        V_LP_MMSE = p*((p*(V_MR*V_MR')+p*Cserved+eyeN)\V_MR);
        
        
        %Go through all UEs served by the AP
        for ind = 1:length(servedUEs)
            
            %Extract UE index
            k = servedUEs(ind);
            
            %Normalize MR precoding
            w = V_MR(:,ind)*sqrt(rho_dist(l,k)/scaling_MR(l,k));
            
            %Compute gain of the signal from UE that arrives at other UEs
            interUserGains_MR(:,k,n) = interUserGains_MR(:,k,n) + Hallj'*w;
            %Normalize L-MMSE precoding
            w = V_L_MMSE(:,ind)*sqrt(rho_dist(l,k)/scaling_L_MMSE(l,k));
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in Corollary 6.3
            signal_L_MMSE(k) = signal_L_MMSE(k) + (Hallj(:,k)'*w)/nbrOfRealizations;
            interf_L_MMSE_n(:,k) = interf_L_MMSE_n(:,k) + Hallj'*w;
            
            %Compute gain of the signal from UE that arrives at other UEs
            interUserGains_L_MMSE(:,k,n) = interUserGains_L_MMSE(:,k,n) + Hallj'*w;
            
            %Normalize LP-MMSE precoding
            w = V_LP_MMSE(:,ind)*sqrt(rho_dist(l,k)/scaling_LP_MMSE(l,k));
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in Corollary 6.3  
            signal_LP_MMSE(k) = signal_LP_MMSE(k) + (Hallj(:,k)'*w)/nbrOfRealizations;
            interf_LP_MMSE_n(:,k) = interf_LP_MMSE_n(:,k) + Hallj'*w;
            
            %Compute gain of the signal from UE that arrives at other UEs
            interUserGains_LP_MMSE(:,k,n) = interUserGains_LP_MMSE(:,k,n) + Hallj'*w;
            
            
            
            
            
        end
        
    end
    
    
    
    %Consider the centralized schemes
    
    
    %Go through all UEs
    for k = 1:K
        
        
        %Determine the set of serving APs
        servingAPs = find(D(:,k)==1);
        
        La = length(servingAPs);
        
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;
        
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hallj_active = zeros(N*La,K);
        
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);
        
        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        %Compute P-MMSE precoding
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));

        %Apply power allocation
        w = w*sqrt(rho_PMMSE(k)/scaling_P_MMSE(k));
        
        
        %Compute realizations of the terms inside the expectations
        %of the signal and interference terms in Theorem 6.1
        signal_P_MMSE(k) = signal_P_MMSE(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_P_MMSE_n(:,k) = interf_P_MMSE_n(:,k) + Hallj_active'*w;
        
        %Compute gain of the signal from UE that arrives at other UEs
        interUserGains_P_MMSE(:,k,n) = interUserGains_P_MMSE(:,k,n) + Hallj_active'*w;
        
               
        %Compute P-RZF combining
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));
            
        %Apply power allocation
        w = w*sqrt(rho_PRZF(k)/scaling_P_RZF(k));
        
        %Compute realizations of the terms inside the expectations
        %of the signal and interference terms in Theorem 6.1
        signal_P_RZF(k) = signal_P_RZF(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_P_RZF_n(:,k) = interf_P_RZF_n(:,k) + Hallj_active'*w;
        
        %Compute gain of the signal from UE that arrives at other UEs
        interUserGains_P_RZF(:,k,n) = interUserGains_P_RZF(:,k,n) + Hallj_active'*w;
        
                
        %Compute MMSE combining
        w = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
            
        %Apply power allocation
        w = w*sqrt(rho_MMSE(k)/scaling_MMSE(k));
        
        %Compute realizations of the terms inside the expectations
        %of the signal and interference terms in Theorem 6.1
        signal_MMSE(k) = signal_MMSE(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_MMSE_n(:,k) = interf_MMSE_n(:,k) + Hallj_active'*w;
        
        %Compute gain of the signal from UE that arrives at other UEs
        interUserGains_MMSE(:,k,n) = interUserGains_MMSE(:,k,n) + Hallj_active'*w;
        
        
        
    end
    
    %Compute interference power in one realization
    interf_MMSE = interf_MMSE + sum(abs(interf_MMSE_n).^2,2)/nbrOfRealizations;
    interf_P_MMSE = interf_P_MMSE + sum(abs(interf_P_MMSE_n).^2,2)/nbrOfRealizations;
    interf_P_RZF = interf_P_RZF + sum(abs(interf_P_RZF_n).^2,2)/nbrOfRealizations;
    
    interf_L_MMSE = interf_L_MMSE + sum(abs(interf_L_MMSE_n).^2,2)/nbrOfRealizations;
    interf_LP_MMSE = interf_LP_MMSE + sum(abs(interf_LP_MMSE_n).^2,2)/nbrOfRealizations;
    
    
    
end



%% Compute the SEs
%Compute SE in Corollary 6.3 with MR  using the closed-form expressions in Corollary 6.4
SE_MR = prelogFactor*real(log2(1+(abs(signal_MR).^2) ./ (interf_MR + sum(abs(cont_MR).^2,2) - abs(signal_MR).^2 + 1)));

%Compute SE in Corollary 6.3 with L-MMSE
SE_L_MMSE = prelogFactor*real(log2(1+(abs(signal_L_MMSE).^2) ./ (interf_L_MMSE - abs(signal_L_MMSE).^2 + 1)));

%Compute SE  in Corollary 6.3 with LP-MMSE
SE_LP_MMSE = prelogFactor*real(log2(1+(abs(signal_LP_MMSE).^2) ./ (interf_LP_MMSE - abs(signal_LP_MMSE).^2 + 1)));

%Compute SE in Theorem 6.1 with MMSE
SE_MMSE = prelogFactor*real(log2(1+(abs(signal_MMSE).^2) ./ (interf_MMSE - abs(signal_MMSE).^2 + 1)));

%Compute SE in Theorem 6.1 with P-MMSE
SE_P_MMSE = prelogFactor*real(log2(1+(abs(signal_P_MMSE).^2) ./ (interf_P_MMSE - abs(signal_P_MMSE).^2 + 1)));

%Compute SE in Theorem 6.1 with P-RZF
SE_P_RZF = prelogFactor*real(log2(1+(abs(signal_P_RZF).^2) ./ (interf_P_RZF - abs(signal_P_RZF).^2 + 1)));




%Compute SEs with perfect CSI at the UEs, i.e., the genie-aided SEs in
%Corollary 6.6




%Go throufh all UEs
for k = 1:K
    
    %Compute SE with MR precoding, assuming perfect CSI at the UE
    
    Gen_SE_MR(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_MR(k,k,:)).^2 ./ ( sum(abs(interUserGains_MR(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
    
    
    %Compute SE with LP-MMSE precoding, assuming perfect CSI at the UE
    
    Gen_SE_LP_MMSE(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_LP_MMSE(k,k,:)).^2 ./ ( sum(abs(interUserGains_LP_MMSE(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
    
    
    %Compute SE with P-MMSE precoding, assuming perfect CSI at the UE
    
    Gen_SE_P_MMSE(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_P_MMSE(k,k,:)).^2 ./ ( sum(abs(interUserGains_P_MMSE(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
    
    %Compute SE with P-RZF precoding, assuming perfect CSI at the UE
    
    Gen_SE_P_RZF(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_P_RZF(k,k,:)).^2 ./ ( sum(abs(interUserGains_P_RZF(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
    
end
%Remove unused large matrices
clear interUserGains_MR interUserGains_MMSE interUserGains_P_MMSE interUserGains_P_RZF interUserGains_L_MMSE interUserGains_LP_MMSE;
