%This Matlab script can be used to reproduce Figure 5.11 in the monograph:
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
nbrOfSetups = 407;

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
tau_p = 5:5:40;

%Angular standard deviation in the local scattering model (in radians)
ASD = deg2rad(15);


%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;

%Prepare to save simulation results
SE_PMMSE = zeros(K,length(tau_p),nbrOfSetups); %P-MMSE
SE_nopt_LPMMSE = zeros(K,length(tau_p),nbrOfSetups); %n-opt LSFD, LP-MMSE
SE_Genie_Small_MMSE = zeros(K,length(tau_p),nbrOfSetups); %Small-cell, L-MMSE

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    
    %Go through all pilot sequence lengths
    for n2=1:length(tau_p)
        
        %Generate one setup with UEs and APs at random locations
        [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K,N,tau_p(n2),1,0,ASD,ASD);
    
        %Generate channel realizations with estimates and estimation
        %error correlation matrices
        [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p(n2),pilotIndex,p);
    
   
        %Compute SEs
        [SE_MMSE, SE_P_MMSE, SE_P_RZF, SE_MR_cent, ...
          SE_opt_L_MMSE,SE_nopt_LP_MMSE, SE_nopt_MR, ...
          SE_L_MMSE, SE_LP_MMSE, SE_MR_dist, ...
          Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR_dist,...
          SE_small_MMSE, Gen_SE_small_MMSE] ...
          = functionComputeSE_uplink(Hhat,H,D,D_small,B,C,tau_c,tau_p(n2),nbrOfRealizations,N,K,L,p,R,pilotIndex);    
        
        %Save SE values
        SE_PMMSE(:,n2,n) = SE_P_MMSE;
        SE_nopt_LPMMSE(:,n2,n) = SE_nopt_LP_MMSE;
        SE_Genie_Small_MMSE(:,n2,n) = Gen_SE_small_MMSE;

  

    end
    
    %Remove large matrices at the end of analyzing this setup
    clear Hhat H B C R;
    
end


%% Plot simulation results
% Plot Figure 5.11
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(tau_p, mean(SE_PMMSE,[1 3]),'kd-','LineWidth',2);
plot(tau_p, mean(SE_nopt_LPMMSE,[1 3]),'rh--','LineWidth',2);
plot(tau_p, mean(SE_Genie_Small_MMSE,[1 3]),'bs-.','LineWidth',2);

xlabel('Pilot sequence length, $\tau_p$','Interpreter','Latex');
ylabel('Average spectral efficiency [bit/s/Hz]','Interpreter','Latex');
legend({'P-MMSE','n-opt LSFD, LP-MMSE','Small-cell, L-MMSE'},'Interpreter','Latex','Location','SouthEast');
