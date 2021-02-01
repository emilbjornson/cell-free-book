%This Matlab script can be used to reproduce Figures 7.1(a) and 7.1(b) in the monograph:
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

%Empty workspace and close figures
close all;
clear;

%% Define simulation setup

%Number of setups with random UE locations
nbrOfSetups = 1400;

%Number of channel realizations per setup
nbrOfRealizations = 1000;

%Number of APs in the cell-free network
L = 100;

%Number of antennas per AP
N = 4;

%Number of UEs in the network
K = 40;

%Length of the coherence block
tau_c = 200;

%Compute number of pilots per coherence block
tau_p = 10;

%Compute the prelog factor assuming only uplink data transmission
preLogFactor = (tau_c-tau_p)/tau_c;

%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(15);  %azimuth angle
ASD_theta = deg2rad(15);   %elevation angle

%Total uplink transmit power per UE (mW)
p = 100;


%Prepare to save simulation results for centralized uplink operation with
%P-MMSE combining
SE_UL_PMMSE_full = zeros(K,nbrOfSetups); %Full
SE_UL_PMMSE_fractional = zeros(K,nbrOfSetups); %FPC, \upsilon = 0.5
SE_UL_PMMSE_fractional2 = zeros(K,nbrOfSetups); %FPC, \upsilon= -0.5
SE_UL_PMMSE_maxmin = zeros(K,nbrOfSetups); %MMF
SE_UL_PMMSE_sumSE = zeros(K,nbrOfSetups); %SumSE

%Prepare to save simulation results for dsitributed uplink operation with
%n-opt LSFD and LP-MMSE combining
SE_UL_nopt_LPMMSE_full = zeros(K,nbrOfSetups); %Full
SE_UL_nopt_LPMMSE_fractional = zeros(K,nbrOfSetups); %FPC, \upsilon = 0.5
SE_UL_nopt_LPMMSE_fractional2 = zeros(K,nbrOfSetups); %FPC, \upsilon= -0.5
SE_UL_nopt_LPMMSE_maxmin = zeros(K,nbrOfSetups); %MMF
SE_UL_nopt_LPMMSE_sumSE = zeros(K,nbrOfSetups); %SumSE




%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs at random locations
    [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K,N,tau_p,1,0,ASD_varphi,ASD_theta);
    
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the cell-free APs
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    
    % Full power
    p_full = p*ones(K,1);
    % Prepare to store FPC coefficients in (7.35) with exponent \upsilon = 0.5
    p_fractional = zeros(K,1);
    % Prepare to store FPC coefficients in (7.35) with exponent \upsilon = -0.5
    p_fractional2 = zeros(K,1);
    
    %Convert channel gains from dB to linear domain
    gainOverNoise = db2pow(gainOverNoisedB);
    
    %Go through all UEs
    for k = 1:K
        
        %Extract which APs that serve UE k
        servingAPs = find(D(:,k)==1);
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = find(sum(D(servingAPs,:),1)>=1);
        %Prepare to compute denominator in (7.35) for two exponent values, i.e., 0.5
        %and -0.5, respectively
        normalization_den = 0;
        normalization_den2 = 0;
        %Go through all UEs that are partially served by the same APs as UE
        %k
        for i = servedUEs
            %Extract APs that serve UE i
            servingAPsi = find(D(:,i)==1);
            %Compute the denominator of (7.35) with exponents 0.5 and -0.5,
            %respectively
            normalization_den = max(normalization_den,  sqrt(sum(gainOverNoise(servingAPsi,i))));
            normalization_den2 = max(normalization_den2,  1/sqrt(sum(gainOverNoise(servingAPsi,i))));
        end
        %Compute p_k according to (7.35) for the two exponents
        p_fractional(k) = p*sqrt(sum(gainOverNoise(servingAPs,k)))/normalization_den;
        p_fractional2(k) = p/sqrt(sum(gainOverNoise(servingAPs,k)))/normalization_den2;
        
        
    end
    
    
    %Obtain the expectations for the computation of the terms in
    %(7.2)-(7.5) for full power control and two FPC schemes
    [signal_P_MMSE_full, signal2_P_MMSE_full, scaling_P_MMSE_full,...
        ~, ~, ~,...
        signal_LP_MMSE_full,signal2_LP_MMSE_full, scaling_LP_MMSE_full] = ...
        functionComputeExpectations(Hhat,H,D,C,nbrOfRealizations,N,K,L,p_full);
    
    [signal_P_MMSE_frac, signal2_P_MMSE_frac, scaling_P_MMSE_frac,...
        ~, ~, ~,...
        signal_LP_MMSE_frac,signal2_LP_MMSE_frac, scaling_LP_MMSE_frac] = ...
        functionComputeExpectations(Hhat,H,D,C,nbrOfRealizations,N,K,L,p_fractional);
    
    [signal_P_MMSE_frac2, signal2_P_MMSE_frac2, scaling_P_MMSE_frac2,...
        ~, ~, ~,...
        signal_LP_MMSE_frac2,signal2_LP_MMSE_frac2, scaling_LP_MMSE_frac2] = ...
        functionComputeExpectations(Hhat,H,D,C,nbrOfRealizations,N,K,L,p_fractional2);
    
    
    
    %Compute (7.2) for the centralized operation with P-MMSE combining
    bk_cent_full = abs(diag(signal_P_MMSE_full)).^2;
    bk_cent_frac = abs(diag(signal_P_MMSE_frac)).^2;
    bk_cent_frac2 = abs(diag(signal_P_MMSE_frac2)).^2;
    
    %Compute (7.3)-(7.4) for the centralized operation with P-MMSE combining
    ck_cent_full = signal2_P_MMSE_full;
    ck_cent_frac = signal2_P_MMSE_frac;
    ck_cent_frac2 = signal2_P_MMSE_frac2;
    
    for k = 1:K
        
        ck_cent_full(k,k) = ck_cent_full(k,k) - bk_cent_full(k);
        ck_cent_frac(k,k) = ck_cent_frac(k,k) - bk_cent_frac(k);
        ck_cent_frac2(k,k) = ck_cent_frac2(k,k) - bk_cent_frac2(k);
    end
    
    %Compute (7.5) for the centralized operation with P-MMSE combining
    sigma2_cent_full = sum(scaling_P_MMSE_full,1).';
    sigma2_cent_frac = sum(scaling_P_MMSE_frac,1).';
    sigma2_cent_frac2 = sum(scaling_P_MMSE_frac2,1).';
    
    %Prepare to store arrays for the compuataion of the n-opt LSFD vectors
    gki_full = zeros(L,K,K);
    gki_frac = zeros(L,K,K);
    gki_frac2 = zeros(L,K,K);

    Gki2_full = zeros(L,L,K,K);
    Gki2_frac = zeros(L,L,K,K);
    Gki2_frac2 = zeros(L,L,K,K);
    
    Gki2p_full = zeros(L,L,K,K);
    Gki2p_frac = zeros(L,L,K,K);
    Gki2p_frac2 = zeros(L,L,K,K);
    
    a_LSFD_full = zeros(L,K);
    a_LSFD_frac = zeros(L,K);
    a_LSFD_frac2 = zeros(L,K);
    
    %Go through all UEs
    for k = 1:K
        
        %Extract which APs that serve UE k
        servingAPs = find(D(:,k)==1);
        
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = find(sum(D(servingAPs,:),1)>=1);
        
        %Go through all UEs
        for i = 1:K
            
            %Compute the expecatation terms in (5.41) for the computation of n-opt LSFD 
            gki_full(:,k,i) = conj(vec(signal_LP_MMSE_full(i,k,:)));
            gki_frac(:,k,i) = conj(vec(signal_LP_MMSE_frac(i,k,:)));
            gki_frac2(:,k,i) = conj(vec(signal_LP_MMSE_frac2(i,k,:)));

            Gki2_full(:,:,k,i) = gki_full(:,k,i)*gki_full(:,k,i)';
            Gki2_frac(:,:,k,i) = gki_frac(:,k,i)*gki_frac(:,k,i)';
            Gki2_frac2(:,:,k,i) = gki_frac2(:,k,i)*gki_frac2(:,k,i)';
            
            for ell = 1:L
                Gki2_full(ell,ell,k,i) = signal2_LP_MMSE_full(i,k,ell);
                Gki2_frac(ell,ell,k,i) = signal2_LP_MMSE_frac(i,k,ell);
                Gki2_frac2(ell,ell,k,i) = signal2_LP_MMSE_frac2(i,k,ell);

            end
            
            Gki2p_full(:,:,k,i) = p_full(i)*Gki2_full(:,:,k,i);
            Gki2p_frac(:,:,k,i) = p_fractional(i)*Gki2_frac(:,:,k,i);
            Gki2p_frac2(:,:,k,i) = p_fractional2(i)*Gki2_frac2(:,:,k,i);

        end
       
        %Compute n-opt LSFD vectors for UE k and for different power
        %control schemes
        a_LSFD_full(servingAPs,k) = p_full(k)*((sum(Gki2p_full(servingAPs,servingAPs,k,servedUEs),4)+diag(scaling_LP_MMSE_full(servingAPs,k)))\gki_full(servingAPs,k,k));
        a_LSFD_frac(servingAPs,k) = p_fractional(k)*((sum(Gki2p_frac(servingAPs,servingAPs,k,servedUEs),4)+diag(scaling_LP_MMSE_frac(servingAPs,k)))\gki_frac(servingAPs,k,k));
        a_LSFD_frac2(servingAPs,k) = p_fractional2(k)*((sum(Gki2p_frac2(servingAPs,servingAPs,k,servedUEs),4)+diag(scaling_LP_MMSE_frac2(servingAPs,k)))\gki_frac2(servingAPs,k,k));
    
    end
    
    %Prepare to store arrays for the computaion of the terms in (7.2)-(7.5)
    %for the distributed opearation
    bk_dist_full = zeros(K,1);
    bk_dist_frac = zeros(K,1);
    bk_dist_frac2 = zeros(K,1);
    
    ck_dist_full = zeros(K,K);
    ck_dist_frac = zeros(K,K);
    ck_dist_frac2 = zeros(K,K);
    
    %Compute (7.5) for distributed operation
    sigma2_dist_full = sum(abs(a_LSFD_full).^2.*scaling_LP_MMSE_full,1).';
    sigma2_dist_frac = sum(abs(a_LSFD_frac).^2.*scaling_LP_MMSE_frac,1).';
    sigma2_dist_frac2 = sum(abs(a_LSFD_frac2).^2.*scaling_LP_MMSE_frac2,1).';
    
    %Go through all UEs
    for k = 1:K
        %Compute (7.2) for distributed operation
        bk_dist_full(k) =  abs(a_LSFD_full(:,k)'*gki_full(:,k,k))^2;
        bk_dist_frac(k) =  abs(a_LSFD_frac(:,k)'*gki_frac(:,k,k))^2;
        bk_dist_frac2(k) =  abs(a_LSFD_frac2(:,k)'*gki_frac2(:,k,k))^2;
        
        %Compute (7.3)-(7.4) for distributed operation

        for i = 1:K
            ck_dist_full(i,k) = real(a_LSFD_full(:,k)'*Gki2_full(:,:,k,i)*a_LSFD_full(:,k));
            ck_dist_frac(i,k) = real(a_LSFD_frac(:,k)'*Gki2_frac(:,:,k,i)*a_LSFD_frac(:,k));
            ck_dist_frac2(i,k) = real(a_LSFD_frac2(:,k)'*Gki2_frac2(:,:,k,i)*a_LSFD_frac2(:,k));
            
        end
        
        ck_dist_full(k,k) = ck_dist_full(k,k) - bk_dist_full(k);
        ck_dist_frac(k,k) = ck_dist_frac(k,k) - bk_dist_frac(k);
        ck_dist_frac2(k,k) = ck_dist_frac2(k,k) - bk_dist_frac2(k);
        
    end
    
    
    %Compute uplink SEs for full and fractional power control schemes using
    %Theorem 5.2
    %Centralized operation 
    SE_UL_PMMSE_full(:,n) = preLogFactor*log2(1+bk_cent_full.*p_full./(ck_cent_full'*p_full+sigma2_cent_full));   
    SE_UL_PMMSE_fractional(:,n) = preLogFactor*log2(1+bk_cent_frac.*p_fractional./(ck_cent_frac'*p_fractional+sigma2_cent_frac));    
    SE_UL_PMMSE_fractional2(:,n) = preLogFactor*log2(1+bk_cent_frac2.*p_fractional2./(ck_cent_frac2'*p_fractional2+sigma2_cent_frac2));
    %Distributed operation
    SE_UL_nopt_LPMMSE_full(:,n) = preLogFactor*log2(1+bk_dist_full.*p_full./(ck_dist_full'*p_full+sigma2_dist_full));   
    SE_UL_nopt_LPMMSE_fractional(:,n) = preLogFactor*log2(1+bk_dist_frac.*p_fractional./(ck_dist_frac'*p_fractional+sigma2_dist_frac));   
    SE_UL_nopt_LPMMSE_fractional2(:,n) = preLogFactor*log2(1+bk_dist_frac2.*p_fractional2./(ck_dist_frac2'*p_fractional2+sigma2_dist_frac2));
    
    %Compute SE with max-min fairness power control in Algorithm 7.1
    %Centralized operation
    SE_UL_PMMSE_maxmin(:,n) = functionUplinkSE_maxmin(bk_cent_full, ck_cent_full, sigma2_cent_full,preLogFactor,K,p);
    %Distributed operation
    SE_UL_nopt_LPMMSE_maxmin(:,n) = functionUplinkSE_maxmin(bk_dist_full, ck_dist_full, sigma2_dist_full,preLogFactor,K,p);
    
    %Compute SE with sum SE maximizing power control in Algorithm 7.2
    %Centralized operation
    SE_UL_PMMSE_sumSE(:,n) = functionUplinkSE_sumSE(bk_cent_full, ck_cent_full, sigma2_cent_full,preLogFactor,K,p);
    %Distributed operation
    SE_UL_nopt_LPMMSE_sumSE(:,n)= functionUplinkSE_sumSE(bk_dist_full, ck_dist_full, sigma2_dist_full,preLogFactor,K,p);

   
end


%% Plot simulation results
% Plot Figure 7.1(a)
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_UL_PMMSE_full(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_UL_PMMSE_fractional(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
ppp = plot(sort(SE_UL_PMMSE_fractional2(:)),linspace(0,1,K*nbrOfSetups),'k:o','LineWidth',2);
ppp.MarkerSize = 6;
ppp.MarkerIndices = 1:ceil(K*nbrOfSetups/7):K*nbrOfSetups;
plot(sort(SE_UL_PMMSE_maxmin(:)),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);
plot(sort(SE_UL_PMMSE_sumSE(:)),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Full','FPC, $\upsilon=0.5$','FPC, $\upsilon=-0.5$','MMF','SumSE' },'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);


% Plot Figure 7.1(b)
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_UL_nopt_LPMMSE_full(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_UL_nopt_LPMMSE_fractional(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
ppp2 = plot(sort(SE_UL_nopt_LPMMSE_fractional2(:)),linspace(0,1,K*nbrOfSetups),'k:o','LineWidth',2);
ppp2.MarkerSize = 6;
ppp2.MarkerIndices = 1:ceil(K*nbrOfSetups/7):K*nbrOfSetups;
plot(sort(SE_UL_nopt_LPMMSE_maxmin(:)),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);
plot(sort(SE_UL_nopt_LPMMSE_sumSE(:)),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Full','FPC, $\upsilon=0.5$','FPC, $\upsilon=-0.5$','MMF','SumSE' },'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);