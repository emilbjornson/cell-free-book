%This Matlab script can be used to reproduce Figures 6.3(a) and 6.5(a) in the monograph:
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
nbrOfSetups = 133;

%Number of channel realizations per setup
nbrOfRealizations = 1000;


%Number of APs 
L = 400;

%Number of antennas per AP
N = 1;

%Number of UEs in the network
K = 40;

%Length of coherence block
tau_c = 200;

%Length of pilot sequences
tau_p = 10;

%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(15);  %azimuth angle
ASD_theta = deg2rad(15);   %elevation angle

%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per AP (mW)
rho_tot = 200;

%Prepare to save simulation results
SE_MMSE_original = zeros(K,nbrOfSetups); %MMSE (All)
SE_MMSE_DCC = zeros(K,nbrOfSetups); %MMSE (DCC)
SE_PMMSE_DCC = zeros(K,nbrOfSetups); %P-MMSE(DCC)
SE_PRZF_DCC = zeros(K,nbrOfSetups); %P-RZF (DCC)

SE_LMMSE_original = zeros(K,nbrOfSetups); %L-MMSE (All)
SE_LMMSE_DCC = zeros(K,nbrOfSetups); %L-MMSE (DCC)
SE_LPMMSE_DCC = zeros(K,nbrOfSetups); %LP-MMSE (DCC)
SE_MR_DCC = zeros(K,nbrOfSetups); %MR (DCC)


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
    
    
    %Compute the power allocation in (6.36) for distributed precoding
    rho_dist = zeros(L,K);
    
    gainOverNoise = db2pow(gainOverNoisedB);
    
    for l = 1:L
        
        %Extract which UEs are served by AP l
        servedUEs = find(D_all(l,:)==1);
        
        %Compute denominator in (6.36)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
        
        for ind = 1:length(servedUEs)
            
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
            
        end
        
    end
    
    
    %Compute SEs for the case where all APs serve all UEs (All)
    [SE_MMSE_all, SE_P_MMSE_all, SE_P_RZF_all,  ...
        SE_L_MMSE_all, SE_LP_MMSE_all, SE_MR_all, ...
        Gen_SE_P_MMSE_all, Gen_SE_P_RZF_all, Gen_SE_LP_MMSE_all, Gen_SE_MR_all] ...
        = functionComputeSE_downlink(Hhat,H,D_all,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex,rho_dist,gainOverNoisedB,rho_tot);
    
    
    
    %Save SE values
    SE_MMSE_original(:,n) = SE_MMSE_all;
    SE_LMMSE_original(:,n) = SE_L_MMSE_all;
    
    
    
    %% Cell-Free Massive MIMO with DCC
    
    %Compute the power allocation in (6.36) for distributed precoding
    rho_dist = zeros(L,K);
    
    for l = 1:L
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute denominator in (6.36)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
        
        for ind = 1:length(servedUEs)
            
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
            

        end
        
    end
    
    %Compute SEs for DCC case
    [SE_MMSE, SE_P_MMSE, SE_P_RZF,  ...
        SE_L_MMSE, SE_LP_MMSE, SE_MR, ...
        Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR] ...
        = functionComputeSE_downlink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex,rho_dist,gainOverNoisedB,rho_tot);
    
    %Save SE values
    SE_MMSE_DCC(:,n) =  SE_MMSE;
    SE_PMMSE_DCC(:,n) = SE_P_MMSE;
    SE_PRZF_DCC(:,n) = SE_P_RZF;
    SE_LMMSE_DCC(:,n) = SE_L_MMSE;
    SE_LPMMSE_DCC(:,n) = SE_LP_MMSE;
    SE_MR_DCC(:,n) = SE_MR;
    
    
    %Remove large matrices at the end of analyzing this setup
    clear Hhat H B C R;
    
end


%% Plot simulation results
% Plot Figure 6.3(a)
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_MMSE_original(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_MMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_PMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
plot(sort(SE_PRZF_DCC(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'MMSE (All)','MMSE (DCC)','P-MMSE (DCC)','P-RZF (DCC)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);


% Plot Figure 6.5(a)
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_LMMSE_original(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_LMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_LPMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
plot(sort(SE_MR_DCC(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'L-MMSE (All)','L-MMSE (DCC)','LP-MMSE (DCC)','MR (DCC)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);

