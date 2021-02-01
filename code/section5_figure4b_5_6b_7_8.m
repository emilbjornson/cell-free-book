%This Matlab script can be used to reproduce Figures 5.4(b), 5.5, 5.6(b), 5.7, and 5.8 in the monograph:
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

%Number of Monte-Carlo setups
nbrOfSetups = 193;

%Number of channel realizations per setup
nbrOfRealizations = 1000;

    
%Number of APs 
L = 100;
    
%Number of antennas per AP
N = 4;
    
%Number of UEs in the network
K = 40;

%Length of coherence block
tau_c = 200;

%Length of pilot sequences
tau_p = 10;

%Compute the prelog factor assuming only uplink data transmission
prelogFactor = (1-tau_p/tau_c);


%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(15);  %azimuth angle
ASD_theta = deg2rad(15);   %elevation angle

%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;

%Prepare to save simulation results
SE_MMSE_original = zeros(K,nbrOfSetups); %MMSE (All)
SE_MMSE_DCC = zeros(K,nbrOfSetups); %MMSE (DCC)
SE_PMMSE_DCC = zeros(K,nbrOfSetups); %P-MMSE(DCC)
SE_PRZF_DCC = zeros(K,nbrOfSetups); %P-RZF (DCC)
SE_MR_DCC = zeros(K,nbrOfSetups); %Centralized MR (DCC)

SE_Genie_PMMSE_DCC = zeros(K,nbrOfSetups); %P-MMSE (Genie-aided)
SE_UatF_PMMSE_DCC = zeros(K,nbrOfSetups); %P-MMSE (UatF)

SE_opt_LMMSE_original = zeros(K,nbrOfSetups); %opt LSFD, L-MMSE (All)
SE_opt_LMMSE_DCC = zeros(K,nbrOfSetups); %opt LSFD, L-MMSE (DCC)
SE_nopt_LPMMSE_DCC = zeros(K,nbrOfSetups); %n-opt LSFD, LP-MMSE (DCC)
SE_nopt_MR_DCC = zeros(K,nbrOfSetups); %n-opt LSFD, MR (DCC)

SE_Genie_nopt_LPMMSE_DCC = zeros(K,nbrOfSetups); %n-opt LSFD, LP-MMSE (Genie-aided)
SE_Genie_nopt_MR_DCC = zeros(K,nbrOfSetups); %n-opt LSFD, MR (Genie-aided)

SE_LMMSE_original = zeros(K,nbrOfSetups); %L-MMSE (All)
SE_LPMMSE_DCC = zeros(K,nbrOfSetups); %LP-MMSE (DCC)
SE_Dist_MR_original = zeros(K,nbrOfSetups); %Distributed MR (All)
SE_Dist_MR_DCC = zeros(K,nbrOfSetups); %Distributed MR (DCC)

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs and APs at random locations
    [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K,N,tau_p,1,0,ASD_varphi,ASD_theta);
    
    %Generate channel realizations with estimates and estimation
    %error correlation matrices
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    %% Original Cell-Free Massive MIMO
    
    %Define the case when all APs serve all UEs
    D_all = ones(L,K);
    
    %Compute SE using combiners and results in Section 5 for centralized
    %and distributed uplink operations for the case when all APs serve all UEs
    [SE_MMSE_all, SE_P_MMSE_all, SE_P_RZF_all, SE_MR_cent_all, ...
          SE_opt_L_MMSE_all,SE_nopt_LP_MMSE_all, SE_nopt_MR_all, ...
          SE_L_MMSE_all, SE_LP_MMSE_all, SE_MR_dist_all, ...
          Gen_SE_P_MMSE_all, Gen_SE_P_RZF_all, Gen_SE_LP_MMSE_all, Gen_SE_MR_dist_all,...
          SE_small_MMSE_all, Gen_SE_small_MMSE_all] ...
          = functionComputeSE_uplink(Hhat,H,D_all,D_small,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex);   
    
    %Save SE values
    SE_MMSE_original(:,n) = SE_MMSE_all;
    SE_opt_LMMSE_original(:,n) = SE_opt_L_MMSE_all;
    SE_LMMSE_original(:,n) = SE_L_MMSE_all;
    SE_Dist_MR_original(:,n) = SE_MR_dist_all;
    
    
    %% Cell-Free Massive MIMO with DCC
    
    %Compute SE using combiners and results in Section 5 for centralized
    %and distributed uplink operations for DCC
    [SE_MMSE, SE_P_MMSE, SE_P_RZF, SE_MR_cent, ...
          SE_opt_L_MMSE,SE_nopt_LP_MMSE, SE_nopt_MR, ...
          SE_L_MMSE, SE_LP_MMSE, SE_MR_dist, ...
          Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR_dist,...
          SE_small_MMSE, Gen_SE_small_MMSE] ...
          = functionComputeSE_uplink(Hhat,H,D,D_small,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex);    
    
    %Save SE values
    SE_MMSE_DCC(:,n) =  SE_MMSE;
    SE_PMMSE_DCC(:,n) = SE_P_MMSE;
    SE_PRZF_DCC(:,n) = SE_P_RZF;
    SE_MR_DCC(:,n) =  SE_MR_cent;
    SE_opt_LMMSE_DCC(:,n) =  SE_opt_L_MMSE;
    SE_nopt_LPMMSE_DCC(:,n) =  SE_nopt_LP_MMSE;
    SE_nopt_MR_DCC(:,n) =  SE_nopt_MR;
    SE_Genie_PMMSE_DCC(:,n) = Gen_SE_P_MMSE;
    SE_Genie_nopt_LPMMSE_DCC(:,n) = Gen_SE_LP_MMSE;
    SE_Genie_nopt_MR_DCC(:,n) = Gen_SE_MR_dist;
    SE_LPMMSE_DCC(:,n) = SE_LP_MMSE;
    SE_Dist_MR_DCC(:,n) = SE_MR_dist;

    
    % Obtain full power vector
    p_full = p*ones(K,1);
    
    %Obtain the expectations for the computation of UatF capacity bound in
    %Theorem 5.2
    [signal_P_MMSE, signal2_P_MMSE, scaling_P_MMSE,...
    signal_P_RZF, signal2_P_RZF, scaling_P_RZF, ...
    signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
     functionComputeExpectations(Hhat,H,D,C,nbrOfRealizations,N,K,L,p_full);
 
    %Prepare to store arrays for the terms in (5.9)
    b_P_MMSE = zeros(K,1);   
    C_P_MMSE = signal2_P_MMSE;

    
    for k = 1:K
        %Compute the square root of the numerator term without p_k in (5.9)
        b_P_MMSE(k) = abs(signal_P_MMSE(k,k));
        %Compute the denominator term without power terms p_i and the noise in (5.9)
        C_P_MMSE(k,k) = C_P_MMSE(k,k) - abs(signal_P_MMSE(k,k))^2;
    end
    %Compute the numerator term without p_k in (5.9) for all UEs
    bk_P_MMSE = b_P_MMSE.^2;    

    %Compute the effective noise variance in (5.9) for all UEs
    sigma2_P_MMSE = sum(scaling_P_MMSE,1).';

    %Compute the SE in (5.8) for P-MMSE combining
    SE_UatF_PMMSE_DCC(:,n) = prelogFactor*log2(1+bk_P_MMSE.*p_full./(ck_P_MMSE'*p_full+sigma2_P_MMSE));

    %Remove large matrices at the end of analyzing this setup
    clear Hhat H B C R;
    
end


%% Plot simulation results
% Plot Figure 5.4(b)
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_MMSE_original(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_MMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_PMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
plot(sort(SE_PRZF_DCC(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(SE_MR_DCC(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',3);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'MMSE (All)','MMSE (DCC)','P-MMSE (DCC)','P-RZF (DCC)','MR (DCC)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);


% Plot Figure 5.5
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_PMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'b','LineWidth',2);
plot(sort(SE_Genie_PMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_UatF_PMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'P-MMSE','P-MMSE (Genie-aided)','P-MMSE (UatF)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);

% Plot Figure 5.6(b)
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_opt_LMMSE_original(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_opt_LMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_nopt_LPMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
plot(sort(SE_nopt_MR_DCC(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'opt LSFD, L-MMSE (All)','opt LSFD, L-MMSE (DCC)','n-opt LSFD, LP-MMSE (DCC)','n-opt LSFD, MR (DCC)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);

% Plot Figure 5.7
figure;
hold on; box on;
set(gca,'fontsize',16);


plot(sort(SE_nopt_LPMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_Genie_nopt_LPMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_nopt_MR_DCC(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
plot(sort(SE_Genie_nopt_MR_DCC(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'n-opt LSFD, LP-MMSE','n-opt LSFD, LP-MMSE (Genie-aided)','n-opt LSFD, MR','n-opt LSFD, MR (Genie-aided)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);

% Plot Figure 5.8
figure;
hold on; box on;
set(gca,'fontsize',16);


plot(sort(SE_nopt_LPMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(SE_LMMSE_original(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_LPMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_Dist_MR_original(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
plot(sort(SE_Dist_MR_DCC(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'n-opt LSFD, LP-MMSE','L-MMSE (All)','LP-MMSE (DCC)','MR (All)','MR (DCC)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);
