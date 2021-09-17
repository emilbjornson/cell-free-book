%This Matlab script can be used to reproduce Figure 5.9(b) in the monograph:
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

%Do not display warnings in the compuation of spatial correlation matrix
%using numerical integration in Cellular Massive MIMO
warning('off','all')

%% Define simulation setup

%Number of Monte-Carlo setups
nbrOfSetups = 126;

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


%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(15);  %azimuth angle
ASD_theta = deg2rad(15);   %elevation angle

%Number of BSs in the cellular network (cannot be changed)
nbrBSs = 4;

%Number of antennas at the 4 BSs
M = 100;

%Length of pilot sequences
tau_p = K/nbrBSs;

%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;

%Prepare to save simulation results
SE_BS_MMSE = zeros(K,nbrOfSetups); %Cellular Massive MIMO, L-MMSE
SE_Genie_Small_MMSE = zeros(K,nbrOfSetups); %Small-cell, L-MMSE
SE_PMMSE_DCC = zeros(K,nbrOfSetups); %Centralized, P-MMSE
SE_nopt_LPMMSE_DCC = zeros(K,nbrOfSetups); %Distributed, n-opt LSFD, LP-MMSE

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs and APs at random locations
    [gainOverNoisedB_AP,gainOverNoisedB_BS,R_AP,R_BS,pilotIndex,BSassignment,D,D_small] = generateSetup_with_cellular(L,K,N,M,1,0,ASD_varphi,ASD_theta);
    
    %Generate channel realizations with estimates and estimation
    %error correlation matrices for Cell-free Massive MIMO
    [Hhat_AP,H_AP,B_AP,C_AP] = functionChannelEstimates(R_AP,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for Cellular Massive MIMO
    [Hhat_BS,~,~,C_BS] = functionChannelEstimates(R_BS,nbrOfRealizations,nbrBSs,K,M,tau_p,pilotIndex,p);   
    
   
    %% Cellular Setup
    
    %Compute SE for the Cellular Massive MIMO system 
    SE_BS_MMSE(:,n) = functionComputeSE_cellular_uplink(Hhat_BS,C_BS,BSassignment,tau_c,tau_p,nbrOfRealizations,M,K,nbrBSs,p);

    %% Cell-Free Massive MIMO with DCC and Small-Cell System
    
    %Compute SEs for the Cell-free Massive MIMO and small-cell system 
    [SE_MMSE, SE_P_MMSE, SE_P_RZF, SE_MR_cent, ...
          SE_opt_L_MMSE,SE_nopt_LP_MMSE, SE_nopt_MR, ...
          SE_L_MMSE, SE_LP_MMSE, SE_MR_dist, ...
          Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR_dist,...
          SE_small_MMSE, Gen_SE_small_MMSE] ...
          = functionComputeSE_uplink(Hhat_AP,H_AP,D,D_small,B_AP,C_AP,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R_AP,pilotIndex);    
    
    %Save SE values
    SE_Genie_Small_MMSE(:,n) = Gen_SE_small_MMSE;
    SE_PMMSE_DCC(:,n) = SE_P_MMSE;
    SE_nopt_LPMMSE_DCC(:,n) = SE_nopt_LP_MMSE;
    
    %Remove large matrices at the end of analyzing this setup
    clear Hhat_AP Hhat_BS H_AP H_BS B_AP C_AP C_BS R_BS R_AP;
    
end


%% Plot simulation results
% Plot Figure 5.9(b)
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_BS_MMSE(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_Genie_Small_MMSE(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_PMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(SE_nopt_LPMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Cellular Massive MIMO, L-MMSE','Small-cell, L-MMSE','Centralized, P-MMSE','Distributed, n-opt LSFD, LP-MMSE'},'Interpreter','Latex','Location','NorthWest');
xlim([0 12]);
