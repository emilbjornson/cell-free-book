%This Matlab script can be used to reproduce Figures 7.2(a) and 7.2(b) in the monograph:
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

%Compute the prelog factor assuming only downlink data transmission
preLogFactor = (tau_c-tau_p)/tau_c;

%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(15);  %azimuth angle
ASD_theta = deg2rad(15);   %elevation angle

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per AP (mW)
rho_tot = 200;


%Prepare to save simulation results for centralized downlink operation with
%P-MMSE precoding
SE_DL_PMMSE_equal = zeros(K,nbrOfSetups); %Equal
SE_DL_PMMSE_fractional1a = zeros(K,nbrOfSetups); %FPA, \upsilon = 0.5, \kappa = 1
SE_DL_PMMSE_fractional1b = zeros(K,nbrOfSetups); %FPA, \upsilon = -0.5, \kappa = 1
SE_DL_PMMSE_fractional2a = zeros(K,nbrOfSetups); %FPA, \upsilon = 0.5, \kappa = 0.5
SE_DL_PMMSE_fractional2b = zeros(K,nbrOfSetups); %FPA, \upsilon = -0.5, \kappa = 0.5
SE_DL_PMMSE_fractional2c = zeros(K,nbrOfSetups); %FPA, \upsilon = -0.5, \kappa = 0
SE_DL_PMMSE_maxmin = zeros(K,nbrOfSetups); %MMF
SE_DL_PMMSE_sumSE = zeros(K,nbrOfSetups); %SumSE

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs at random locations
    [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K,N,tau_p,1,0,ASD_varphi,ASD_theta);
    
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the cell-free APs
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    
    % Full uplink power for the computation of precoding vectors using
    % virtual uplink-downlink duality
    p_full = p*ones(K,1);
   

    %Compute the equal power allocation for centralized precoding
    rho_central = (rho_tot/tau_p)*ones(K,1);
    
    %Obtain the expectations for the computation of the terms in
    %(7.13)-(7.15) 
    [signal_P_MMSE, signal2_P_MMSE, scaling_P_MMSE,...
        signal_P_RZF, signal2_P_RZF, scaling_P_RZF,...
        signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
        functionComputeExpectations(Hhat,H,D,C,nbrOfRealizations,N,K,L,p_full);
    
    %Compute the terms in (7.13)-(7.15)
    bk_PMMSE = zeros(K,1);
    ck_PMMSE = signal2_P_MMSE.';

    for k = 1:K
        
        bk_PMMSE(k) = abs(signal_P_MMSE(k,k))^2;
      
        ck_PMMSE(k,k) = ck_PMMSE(k,k) - bk_PMMSE(k);

    end
    
    %Obtain the scaling factors for the precoding vectors and scale the
    %terms in (7.13)-(7.15) accordingly
    sigma2_PMMSE = sum(scaling_P_MMSE,1).';
   
    for k = 1:K
        
        bk_PMMSE(k) = bk_PMMSE(k)/sigma2_PMMSE(k);
        
        ck_PMMSE(k,:) = ck_PMMSE(k,:)/sigma2_PMMSE(k);
    
    end
    
    %Scale the expected values of the norm squares of the portions of the centralized precoding in
    %accordance with the normalized precoding vectors in (7.16)
    portionScaling_P_MMSE = scaling_P_MMSE./repmat(sum(scaling_P_MMSE,1),[L 1]);

    %Compute the fractional power allocation for centralized precoding
    %according to (7.43) with different \upsilon and \kappa parameters
    rho_1a = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_P_MMSE,0.5,1);
    rho_1b = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_P_MMSE,-0.5,1);
    rho_2a = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_P_MMSE,0.5,0.5);
    rho_2b = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_P_MMSE,-0.5,0.5);
    rho_2c = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_P_MMSE,-0.5,0);

    %Compute SEs according to Theorem 6.1 with several scalable power
    %allocation schemes
    SE_DL_PMMSE_equal(:,n) = preLogFactor*log2(1+bk_PMMSE.*rho_central./(ck_PMMSE'*rho_central+1));       
    SE_DL_PMMSE_fractional1a(:,n) = preLogFactor*log2(1+bk_PMMSE.*rho_1a./(ck_PMMSE'*rho_1a+1));
    SE_DL_PMMSE_fractional1b(:,n) = preLogFactor*log2(1+bk_PMMSE.*rho_1b./(ck_PMMSE'*rho_1b+1));
    SE_DL_PMMSE_fractional2a(:,n) = preLogFactor*log2(1+bk_PMMSE.*rho_2a./(ck_PMMSE'*rho_2a+1));
    SE_DL_PMMSE_fractional2b(:,n) = preLogFactor*log2(1+bk_PMMSE.*rho_2b./(ck_PMMSE'*rho_2b+1));
    SE_DL_PMMSE_fractional2c(:,n) = preLogFactor*log2(1+bk_PMMSE.*rho_2c./(ck_PMMSE'*rho_2c+1));

    %Compute SEs according to Theorem 6.1 with max-min fair power allocation in Algorithm 7.3 
    SE_DL_PMMSE_maxmin(:,n) = functionDownlinkSE_maxmin(bk_PMMSE, ck_PMMSE, portionScaling_P_MMSE,preLogFactor,K,rho_tot);
    %Compute SEs according to Theorem 6.1 with sum SE maximizing power allocation in Algorithm 7.4 
    SE_DL_PMMSE_sumSE(:,n) =  functionDownlinkSE_sumSE(bk_PMMSE, ck_PMMSE, portionScaling_P_MMSE,preLogFactor,L,K,rho_tot,tau_p);
    
end

%% Plot simulation results
% Plot Figure 7.2(a)
figure;
hold on; box on;
set(gca,'fontsize',16);

ppp3 = plot(sort(SE_DL_PMMSE_fractional1a(:)),linspace(0,1,K*nbrOfSetups),'b-.o','LineWidth',2);
plot(sort(SE_DL_PMMSE_fractional2a(:)),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);
ppp3.MarkerSize = 6;
ppp3.MarkerIndices = 1:ceil(K*nbrOfSetups/7):K*nbrOfSetups;

ppp4 = plot(sort(SE_DL_PMMSE_fractional1b(:)),linspace(0,1,K*nbrOfSetups),'k:o','LineWidth',2);
plot(sort(SE_DL_PMMSE_fractional2b(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
ppp4.MarkerSize = 6;
ppp4.MarkerIndices = 1:ceil(K*nbrOfSetups/7):K*nbrOfSetups;

plot(sort(SE_DL_PMMSE_fractional2c(:)),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);



xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'$\upsilon=0.5$, $\kappa=1$','$\upsilon=0.5$, $\kappa=0.5$', '$\upsilon=-0.5$, $\kappa=1$','$\upsilon=-0.5$, $\kappa=0.5$','$\upsilon=-0.5$, $\kappa=0$' },'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);


% Plot Figure 7.2(b)
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_DL_PMMSE_equal(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_DL_PMMSE_fractional2b(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
plot(sort(SE_DL_PMMSE_maxmin(:)),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);
plot(sort(SE_DL_PMMSE_sumSE(:)),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Equal','FPA','MMF','SumSE' },'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);