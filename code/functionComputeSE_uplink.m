function [SE_MMSE, SE_P_MMSE, SE_P_RZF, SE_MR_cent, ...
    SE_opt_L_MMSE,SE_nopt_LP_MMSE, SE_nopt_MR, ...
    SE_L_MMSE, SE_LP_MMSE, SE_MR_dist, ...
    Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR_dist, ...
    SE_small_MMSE, Gen_SE_small_MMSE] ...
    = functionComputeSE_uplink(Hhat,H,D,D_small,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex)
%Compute uplink SE for different receive combining schemes using the capacity
%bound in Theorem 5.1 for the centralized schemes and the capacity bound
%in Theorem 5.4 for the distributed schemes. Compute the genie-aided uplink
%SE from Corollary 5.9 for the centralized operation and from Corollary 5.10 
%for the distributed operation. 
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
%D_small           = DCC matrix for small-cell setup with dimension L x K
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
%
%OUTPUT:
%SE_MMSE           = SEs achieved with MMSE combining in (5.11)
%SE_P_MMSE         = SEs achieved with P-MMSE combining in (5.16)
%SE_P_RZF          = SEs achieved with P-RZF combining in (5.18)
%SE_MR_cent        = SEs achieved with centralized MR combining in (5.14)
%SE_opt_L_MMSE     = SEs achieved with opt LSFD in (5.30) and 
%                    L-MMSE combining in (5.29)
%SE_nopt_LP_MMSE   = SEs achieved with n-opt LSFD in (5.41) and 
%                    LP-MMSE combining in (5.39)
%SE_nopt_MR        = SEs achieved with n-opt LSFD in (5.41) and 
%                    local MR combining in (5.32)
%SE_L_MMSE         = SEs achieved with L-MMSE combining in (5.29),
%                    without LSFD   
%SE_LP_MMSE        = SEs achieved with LP-MMSE combining in (5.39),
%                    without LSFD 
%SE_MR_dist        = SEs achieved with local MR combining in (5.32),
%                    without LSFD
%Gen_SE_P_MMSE     = Genie-aided SEs achieved with P-MMSE combining in (5.16)
%Gen_SE_P_RZF      = Genie-aided SEs achieved with P-RZF combining in (5.18)
%Gen_SE_LP_MMSE    = Genie-aided SEs achieved n-opt LSFD in (5.41) and 
%                    LP-MMSE combining in (5.39)
%Gen_SE_MR_dist    = Genie-aided SEs achieved n-opt LSFD in (5.41) and 
%                    local MR combining in (5.32)
%SE_small_MMSE     = SEs achieved with L-MMSE combining for small-cell setup
%Gen_SE_small_MMSE = Genie-aided SEs achieved with L-MMSE combining for 
%                    small-cell setup
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

%Compute the prelog factor assuming only uplink data transmission
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_MMSE =  zeros(K,1);
SE_P_MMSE = zeros(K,1);
SE_P_RZF = zeros(K,1);
SE_MR_cent = zeros(K,1);
SE_opt_L_MMSE = zeros(K,1);
SE_nopt_LP_MMSE = zeros(K,1);
SE_nopt_MR = zeros(K,1);
SE_L_MMSE = zeros(K,1);
SE_LP_MMSE = zeros(K,1);
SE_MR_dist = zeros(K,1);
Gen_SE_P_MMSE = zeros(K,1);
Gen_SE_P_RZF = zeros(K,1);
Gen_SE_LP_MMSE = zeros(K,1);
Gen_SE_MR_dist = zeros(K,1);
SE_small_MMSE = zeros(K,1);
Gen_SE_small_MMSE = zeros(K,1);


%Prepare to store the terms that appear in SEs
gki_MR = zeros(K,L,K);
gki2_MR = zeros(K,L,K);
Fk_MR = zeros(L,K);

gki_L_MMSE = zeros(K,K,L);
gki2_L_MMSE = zeros(K,K,L);
Fk_L_MMSE = zeros(L,K);

gki_LP_MMSE = zeros(K,K,L);
gki2_LP_MMSE = zeros(K,K,L);
Fk_LP_MMSE = zeros(L,K);

gen_gki_MR = zeros(K,L,K,nbrOfRealizations);
gen_Fk_MR = zeros(L,K,nbrOfRealizations);

gen_gki_L_MMSE = zeros(K,L,K,nbrOfRealizations);
gen_Fk_L_MMSE = zeros(L,K,nbrOfRealizations);

gen_gki_LP_MMSE = zeros(K,L,K,nbrOfRealizations);
gen_Fk_LP_MMSE = zeros(L,K,nbrOfRealizations);


%% Compute MR closed-form expectations according to Corollary 5.6

%Go through each AP
for l = 1:L
    %Extract which UEs are served by the AP l
    servedUEs = find(D(l,:)==1);
    
    %Go through all UEs served by the AP l
    for ind = 1:length(servedUEs)
        
        %Extract UE index
        k = servedUEs(ind);
        
        
        %Noise scaling according to (5.35)
        Fk_MR(l,k) = trace(B(:,:,l,k));
        
     
        for i = 1:K
            
            %Compute the first term in (5.34)
            gki2_MR(i,l,k) = real(trace(B(:,:,l,k)*R(:,:,l,i))); 
            
            %If UE i shares the same pilot with UE k
            if pilotIndex(k) == pilotIndex(i)
                
                %The term in (5.33)
                gki_MR(i,l,k) = real(trace((B(:,:,l,k)/R(:,:,l,k))*R(:,:,l,i)));
                %The second term in (5.34)
                gki2_MR(i,l,k) = gki2_MR(i,l,k) + (real(trace((B(:,:,l,k)/R(:,:,l,k))*R(:,:,l,i))))^2;
            end
            
        end
        
    end
    
end


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Consider the distributed schemes

    
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
        
        %Compute MR combining according to (5.32)
        V_MR = Hhatallj(:,servedUEs);
        
        %Compute L-MMSE combining according to (5.29)
        V_L_MMSE = p*((p*(Hhatallj*Hhatallj'+sum(C(:,:,l,:),4))+eyeN)\V_MR);
        
        %Compute LP-MMSE combining according to (5.39)
        V_LP_MMSE = p*((p*(V_MR*V_MR'+Cserved)+eyeN)\V_MR);
        
        %Compute the conjugates of the vectors g_{ki} in (5.23) for three
        %combining schemes above for the considered channel realization
        TemporMatr_MR = Hallj'*V_MR;
        TemporMatr_L_MMSE = Hallj'*V_L_MMSE;
        TemporMatr_LP_MMSE = Hallj'*V_LP_MMSE;
        
        %Update the sample mean estimates of the expectations in (5.27) 
        Fk_L_MMSE(l,servedUEs) = Fk_L_MMSE(l,servedUEs) + vecnorm(V_L_MMSE).^2/nbrOfRealizations;
        Fk_LP_MMSE(l,servedUEs) = Fk_LP_MMSE(l,servedUEs) + vecnorm(V_LP_MMSE).^2/nbrOfRealizations;
        
        %Store the instantaneous combining vector norms for the channel
        %realization n to be used later
        gen_Fk_MR(l,servedUEs,n) = vecnorm(V_MR).^2;
        gen_Fk_L_MMSE(l,servedUEs,n) = vecnorm(V_L_MMSE).^2;
        gen_Fk_LP_MMSE(l,servedUEs,n) = vecnorm(V_LP_MMSE).^2;
        
        %Update the sample mean estimates of the expectations related to g_{ki} in
        %(5.23) to be used in the SE expression of Theorem 5.4
        gki_L_MMSE(:,servedUEs,l) = gki_L_MMSE(:,servedUEs,l) + TemporMatr_L_MMSE/nbrOfRealizations;
        gki_LP_MMSE(:,servedUEs,l) = gki_LP_MMSE(:,servedUEs,l) + TemporMatr_LP_MMSE/nbrOfRealizations;
        
        gki2_L_MMSE(:,servedUEs,l) = gki2_L_MMSE(:,servedUEs,l) + abs(TemporMatr_L_MMSE).^2/nbrOfRealizations;
        gki2_LP_MMSE(:,servedUEs,l) = gki2_LP_MMSE(:,servedUEs,l) + abs(TemporMatr_LP_MMSE).^2/nbrOfRealizations;
        
        %Store the instantaneous entries of g_{ki} in (5.23) for the channel
        %realization n to be used later
        gen_gki_MR(:,l,servedUEs,n) = TemporMatr_MR;
        gen_gki_L_MMSE(:,l,servedUEs,n) = TemporMatr_L_MMSE;
        gen_gki_LP_MMSE(:,l,servedUEs,n) = TemporMatr_LP_MMSE;
        
        
        
    end
    
    
    
    %Consider the centralized schemes
    
    
    %Go through all UEs
    for k = 1:K
        
        
        %Determine the set of serving APs for UE k
        servingAPs = find(D(:,k)==1); %cell-free setup
        servingAP_small = find(D_small(:,k)==1); %small-cell setup
        
        %Compute the number of APs that serve UE k
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
        
        %Compute P-MMSE combining according to (5.16)
        v = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));
        
        %Compute numerator and denominator of instantaneous SINR in (5.5)
        numerator = p*abs(v'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
        
        %Compute numerator and denominator of instantaneous SINR in (5.50)
        numerator_gen = p*abs(v'*Hallj_active(:,k))^2;
        denominator_gen = p*norm(v'*Hallj_active)^2 + v'*v - numerator_gen;
        
        %Update the SE by computing the instantaneous SE for one
        %channel realization according to (5.4)
        SE_P_MMSE(k) = SE_P_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
        %Update the SE by computing the instantaneous SE for one
        %channel realization according to (5.49)        
        Gen_SE_P_MMSE(k) = Gen_SE_P_MMSE(k) + prelogFactor*real(log2(1+numerator_gen/denominator_gen))/nbrOfRealizations;
        
        
        
        %Compute P-RZF combining according to (5.18)
        v = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));
        
        %Compute numerator and denominator of instantaneous SINR in (5.5)
        numerator = p*abs(v'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
        
        %Compute numerator and denominator of instantaneous SINR in (5.50)
        numerator_gen = p*abs(v'*Hallj_active(:,k))^2;
        denominator_gen = p*norm(v'*Hallj_active)^2 + v'*v - numerator_gen;
        
        %Update the SE by computing the instantaneous SE for one
        %channel realization according to (5.4)
        SE_P_RZF(k) = SE_P_RZF(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
        %Update the SE by computing the instantaneous SE for one
        %channel realization according to (5.49)  
        Gen_SE_P_RZF(k) = Gen_SE_P_RZF(k) + prelogFactor*real(log2(1+numerator_gen/denominator_gen))/nbrOfRealizations;
        
        
        %Compute MMSE combining according to (5.11)
        v = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
        
        %Extract the required channel estimates and estimation corelation
        %matrices for the SE computation of UE k in small-cell setup
        Hhatallj_active_small = reshape(Hhat((servingAP_small-1)*N+1:servingAP_small*N,n,:),[N K]);
        C_tot_blk_small = reshape(sum(C(:,:,servingAP_small,:),4),[N,N]);
        %Compute L-MMSE combining for small-cell setup according to (5.29)
        v_small = p*((p*(Hhatallj_active_small*Hhatallj_active_small')+p*C_tot_blk_small+eye(N))\Hhatallj_active_small(:,k));
        
        %Compute numerator and denominator of instantaneous SINR in (5.5)
        %for MMSE combining
        numerator = p*abs(v'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
        
        %Update the SE by computing the instantaneous SE for one
        %channel realization according to (5.4)        
        SE_MMSE(k) = SE_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
        
        %Compute numerator and denominator of instantaneous SINR in (5.5)
        %for L-MMSE combining in small-cell setup
        numerator = p*abs(v_small'*Hhatallj_active_small(:,k))^2;
        denominator = p*norm(v_small'*Hhatallj_active_small)^2 + v_small'*(p*C_tot_blk_small+eye(N))*v_small - numerator;
        
        %Update the SE by computing the instantaneous SE for one
        %channel realization according to (5.4) in small-cell setup
        SE_small_MMSE(k) = SE_small_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
        
        %Compute centralized MR combining according to (5.14)
        v = Hhatallj_active(:,k);
        
        %Compute numerator and denominator of instantaneous SINR in (5.5)     
        numerator = p*abs(v'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
        
        %Update the SE by computing the instantaneous SE for one
        %channel realization according to (5.4)         
        SE_MR_cent(k) = SE_MR_cent(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
        
    end
    
end

%Permute the arrays that consist of the expectations that appear in Theorem
%5.4 to compute LSFD vectors and the corresponding SEs
gki_L_MMSE = permute(gki_L_MMSE,[1 3 2]);
gki_LP_MMSE = permute(gki_LP_MMSE,[1 3 2]);
gki2_L_MMSE = permute(gki2_L_MMSE,[1 3 2]);
gki2_LP_MMSE = permute(gki2_LP_MMSE,[1 3 2]);



%Prepare to store n-opt LSFD vectors to be used later
a_nopt1 = zeros(L,K);
a_nopt2 = zeros(L,K);

%% Compute the SEs for Distributed Case
for k = 1:K
    
    %Determine the set of serving APs for UE k
    servingAPs = find(D(:,k)==1);
    %The number of APs that serve UE k 
    La = length(servingAPs);
    
    %Determine which UEs that are served by partially the same set
    %of APs as UE k, i.e., the set in (5.15)
    servedUEs = find(sum(D(servingAPs,:),1)>=1);
    
 
    %Expected value of g_{kk}, scaled by \sqrt{p} for L-MMSE combining
    num_vector = conj(vec(sqrt(p)*gki_L_MMSE(k,servingAPs,k)));
    %Compute the matrix whose inverse is taken in (5.30) using the first-
    %and second-order moments of the entries in the vectors g_{ki}
    temporMatr = gki_L_MMSE(:,servingAPs,k)'*gki_L_MMSE(:,servingAPs,k);
    denom_matrix = p*(diag(sum(gki2_L_MMSE(:,servingAPs,k),1))...
        +temporMatr-diag(diag(temporMatr)))...
        -num_vector*num_vector'+diag(Fk_L_MMSE(servingAPs,k));
    
    %Compute the opt LSFD according to (5.30)
    a_opt = denom_matrix\num_vector;
    
    %Compute the corresponding weights for the case without LSFD
    a_dist = ones(La,1);
    
    %Compute the SE achieved with opt LSFD and L-MMSE combining according to
    %(5.25)
    SE_opt_L_MMSE(k) = prelogFactor*real(log2(1+abs(a_opt'*num_vector)^2/(a_opt'*denom_matrix*a_opt)));
    
    %Compute the SE achieved with L-MMSE combining and without LSFD according to
    %(5.25)
    SE_L_MMSE(k) = prelogFactor*real(log2(1+abs(a_dist'*num_vector)^2/(a_dist'*denom_matrix*a_dist)));
    
    
    

    %Expected value of g_{kk}, scaled by \sqrt{p} for LP-MMSE combining
    num_vector = conj(vec(sqrt(p)*gki_LP_MMSE(k,servingAPs,k)));
    %Compute the denominator matrix to compute SE in Theorem 5.4 using the first-
    %and second-order moments of the entries in the vectors g_{ki}
    temporMatr = gki_LP_MMSE(:,servingAPs,k)'*gki_LP_MMSE(:,servingAPs,k);
    denom_matrix = p*(diag(sum(gki2_LP_MMSE(:,servingAPs,k),1))...
        +temporMatr-diag(diag(temporMatr)))...
        -num_vector*num_vector'+diag(Fk_LP_MMSE(servingAPs,k));
    
    %Compute the matrix whose inverse is taken in (5.41) using the first-
    %and second-order moments of the entries in the vectors g_{ki}
    temporMatr = gki_LP_MMSE(servedUEs,servingAPs,k)'*gki_LP_MMSE(servedUEs,servingAPs,k);
    
    denom_matrix_partial =  p*(diag(sum(gki2_LP_MMSE(servedUEs,servingAPs,k),1))...
        +temporMatr-diag(diag(temporMatr)))...
        -num_vector*num_vector'+diag(Fk_LP_MMSE(servingAPs,k));
    
    
    %Compute the n-opt LSFD according to (5.41) for LP-MMSE combining
    a_nopt = denom_matrix_partial\num_vector;
    
    %Compute the SE achieved with n-opt LSFD and LP-MMSE combining according to
    %(5.25)
    SE_nopt_LP_MMSE(k) = prelogFactor*real(log2(1+abs(a_nopt'*num_vector)^2/(a_nopt'*denom_matrix*a_nopt)));
    %Compute the SE achieved with LP-MMSE combining and without LSFD according to
    %(5.25)
    SE_LP_MMSE(k) = prelogFactor*real(log2(1+abs(a_dist'*num_vector)^2/(a_dist'*denom_matrix*a_dist)));
    
    %Store the n-opt LSFD vector for LP-MMSE combining to be used later
    a_nopt1(servingAPs,k) = a_nopt;
    
    
 
    %Expected value of g_{kk}, scaled by \sqrt{p} for local MR combining
    num_vector = vec(sqrt(p)*gki_MR(k,servingAPs,k));
    %Compute the denominator matrix to compute SE in Theorem 5.4 using the first-
    %and second-order moments of the entries in the vectors g_{ki}
    temporMatrrr =  gki_MR(:,servingAPs,k).'*conj(gki_MR(:,servingAPs,k));
    denom_matrix = p*(diag(sum(gki2_MR(:,servingAPs,k),1))...
        +temporMatrrr-diag(diag(temporMatrrr)))...
        -num_vector*num_vector'...
        +diag(Fk_MR(servingAPs,k));
    
    %Compute the matrix whose inverse is taken in (5.41) using the first-
    %and second-order moments of the entries in the vectors g_{ki}
    temporMatrrr =  gki_MR(servedUEs,servingAPs,k).'*conj(gki_MR(servedUEs,servingAPs,k));
    
    denom_matrix_partial =  p*(diag(sum(gki2_MR(servedUEs,servingAPs,k),1))...
        +temporMatrrr-diag(diag(temporMatrrr)))...
        -num_vector*num_vector'...
        +diag(Fk_MR(servingAPs,k));
    
    
    %Compute the n-opt LSFD according to (5.41) for local MR combining
    a_nopt = denom_matrix_partial\num_vector;
    
    %Compute the SE achieved with n-opt LSFD and local MR combining according to
    %(5.25)
    SE_nopt_MR(k) = prelogFactor*real(log2(1+abs(a_nopt'*num_vector)^2/(a_nopt'*denom_matrix*a_nopt)));
    %Compute the SE achieved with local MR combining and without LSFD according to
    %(5.25)
    SE_MR_dist(k) = prelogFactor*real(log2(1+abs(a_dist'*num_vector)^2/(a_dist'*denom_matrix*a_dist)));
    
    %Store the n-opt LSFD vector for local MR combining to be used later
    a_nopt2(servingAPs,k) = a_nopt;
    
    
end



%Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all UEs
    for k = 1:K
        
        %Determine the set of serving APs
        servingAPs = find(D(:,k)==1);  %cell-free setup
        servingAP_small = find(D_small(:,k)==1); %small-cell setup
        
        %Compute the numerator and the denominator in (5.53) with a single
        %serving AP in a small-cell setup
        numerator_gen = p*abs(gen_gki_L_MMSE(k,servingAP_small,k,n))^2;
        denominator_gen = p*gen_gki_L_MMSE(:,servingAP_small,k,n)'*...
            gen_gki_L_MMSE(:,servingAP_small,k,n)...
            + gen_Fk_L_MMSE(servingAP_small,k,n)...
            - numerator_gen;
        
        %Update the genie-aided SE by computing the instantaneous SE for one
        %channel realization according to (5.52) in small-cell setup
        Gen_SE_small_MMSE(k) = Gen_SE_small_MMSE(k) + prelogFactor*real(log2(1+numerator_gen/denominator_gen))/nbrOfRealizations;
        
         
        %Compute the numerator and the denominator in (5.53) for n-opt LSFD
        %and LP-MMSE combining
        numerator_gen = p*abs(a_nopt1(servingAPs,k)'*vec(conj(gen_gki_LP_MMSE(k,servingAPs,k,n))))^2;
        temporMatrr =   gen_gki_LP_MMSE(:,servingAPs,k,n)'*gen_gki_LP_MMSE(:,servingAPs,k,n);
        denominator_gen = p*a_nopt1(servingAPs,k)'*...
            temporMatrr*a_nopt1(servingAPs,k)+...
            + a_nopt1(servingAPs,k)'*diag(gen_Fk_LP_MMSE(servingAPs,k,n))*a_nopt1(servingAPs,k)...
            - numerator_gen;
        
        %Update the genie-aided SE by computing the instantaneous SE for one
        %channel realization according to (5.52) 
        Gen_SE_LP_MMSE(k) = Gen_SE_LP_MMSE(k) + prelogFactor*real(log2(1+numerator_gen/denominator_gen))/nbrOfRealizations;
        
        
        %Compute the numerator and the denominator in (5.53) for n-opt LSFD
        %and local MR combining
        numerator_gen = p*abs(a_nopt2(servingAPs,k)'*vec(conj(gen_gki_MR(k,servingAPs,k,n))))^2;
        temporMatrr =   gen_gki_MR(:,servingAPs,k,n)'*gen_gki_MR(:,servingAPs,k,n);
        denominator_gen = p*a_nopt2(servingAPs,k)'*...
            temporMatrr*a_nopt2(servingAPs,k)+...
            + a_nopt2(servingAPs,k)'*diag(gen_Fk_MR(servingAPs,k,n))*a_nopt2(servingAPs,k)...
            - numerator_gen;
        
        %Update the genie-aided SE by computing the instantaneous SE for one
        %channel realization according to (5.52) 
        Gen_SE_MR_dist(k) = Gen_SE_MR_dist(k) + prelogFactor*real(log2(1+numerator_gen/denominator_gen))/nbrOfRealizations;
        
        
    end
end

%Remove unused large arrays
clear gki_MR gki2_MR gki_L_MMSE gki2_L_MMSE gki_LP_MMSE gki2_LP_MMSE;
clear gen_gki_MR gen_Fk_MR gen_gki_L_MMSE gen_Fk_L_MMSE gen_gki_LP_MMSE gen_Fk_LP_MMSE;
