%This Matlab script can be used to reproduce Figure 6.9 in the monograph:
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
nbrOfSetups = 550;

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

%Angular standard deviation in the local scattering model (in radians)
ASD = deg2rad(5:5:30);


%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;


%Total downlink transmit power per AP (mW)
rho_tot = 200;


%Prepare to save simulation results
SE_PMMSE_UnCorr = zeros(K,nbrOfSetups); %P-MMSE (Uncorrelated)
SE_PMMSE_Corr = zeros(K,length(ASD),nbrOfSetups); %P-MMSE (Correlated)

SE_PRZF_UnCorr = zeros(K,nbrOfSetups); %P-RZF (Uncorrelated)
SE_PRZF_Corr = zeros(K,length(ASD),nbrOfSetups); %P-RZF (Correlated)

SE_LPMMSE_UnCorr = zeros(K,nbrOfSetups); %LP-MMSE (Uncorrelated)
SE_LPMMSE_Corr = zeros(K,length(ASD),nbrOfSetups); %LP-MMSE (Correlated)



%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs and APs at random locations with
    %uncorrelated fading
    [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K,N,tau_p,1,0);
    
    %Generate channel realizations with estimates and estimation
    %error correlation matrices
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    
    %% Cell-Free Massive MIMO with DCC
    
    %Compute the power allocation in (6.36) for distributed precoding
    rho_dist = zeros(L,K);
    
    gainOverNoise = db2pow(gainOverNoisedB);
    
    for l = 1:L
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute denominator in (6.36)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
        
        for ind = 1:length(servedUEs)
            
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
            
        end
        
    end
    
    %Compute SEs for DCC case with uncorrelated fading
    [SE_MMSE, SE_P_MMSE, SE_P_RZF,  ...
        SE_L_MMSE, SE_LP_MMSE, SE_MR, ...
        Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR] ...
        = functionComputeSE_downlink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex,rho_dist,gainOverNoisedB,rho_tot);
    
    %Save SE values
    SE_PMMSE_UnCorr(:,n) = SE_P_MMSE;
    SE_PRZF_UnCorr(:,n) = SE_P_RZF;
    SE_LPMMSE_UnCorr(:,n) = SE_LP_MMSE;
    
    %Go through all angular standard deviations
    for n2=1:length(ASD)
        
        %Generate one setup with UEs and APs at random locations with
        %correlated fading
        [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K,N,tau_p,1,0,ASD(n2),ASD(n2));
        
        %Generate channel realizations with estimates and estimation
        %error correlation matrices
        [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
        
        
        %Compute the power allocation in (6.36) for distributed precoding
        rho_dist = zeros(L,K);
        
        gainOverNoise = db2pow(gainOverNoisedB);
        
        for l = 1:L
            
            %Extract which UEs are served by AP l
            servedUEs = find(D(l,:)==1);
            
            %Compute denominator in (6.36)
            normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
            
            for ind = 1:length(servedUEs)
                
                rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
                
            end
            
        end
        
        %Compute SEs for DCC case with correlated fading
        [SE_MMSE, SE_P_MMSE, SE_P_RZF,  ...
            SE_L_MMSE, SE_LP_MMSE, SE_MR, ...
            Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR] ...
            = functionComputeSE_downlink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex,rho_dist,gainOverNoisedB,rho_tot);
        
        
        
        %Save SE values
        SE_PMMSE_Corr(:,n2,n) = SE_P_MMSE;
        SE_PRZF_Corr(:,n2,n) = SE_P_RZF;
        SE_LPMMSE_Corr(:,n2,n) = SE_LP_MMSE;
        
        
    end
    
    %Remove large matrices at the end of analyzing this setup
    clear Hhat H B C R;
    
end


%% Plot simulation results
% Plot Figure 6.9
figure;
hold on; box on;
set(gca,'fontsize',16);


plot(rad2deg(ASD), mean(SE_PMMSE_Corr,[1 3]),'kd-','LineWidth',2);
plot(rad2deg(ASD), mean(SE_PRZF_Corr,[1 3]),'rh--','LineWidth',2);
plot(rad2deg(ASD), mean(SE_LPMMSE_Corr,[1 3]),'bs-.','LineWidth',2);
plot(rad2deg(ASD), repmat(mean(SE_PMMSE_UnCorr,'all'), 1, length(ASD)),'k:','LineWidth',2);
plot(rad2deg(ASD), repmat(mean(SE_PRZF_UnCorr,'all'),1,length(ASD)),'r:','LineWidth',2);
plot(rad2deg(ASD), repmat(mean(SE_LPMMSE_UnCorr,'all'), 1, length(ASD)),'b:','LineWidth',2);

xlabel('ASD, $\sigma_{\varphi}=\sigma_{\theta}$','Interpreter','Latex');
ylabel('Average spectral efficiency [bit/s/Hz]','Interpreter','Latex');
legend({'P-MMSE','P-RZF','LP-MMSE'},'Interpreter','Latex','Location','SouthEast');
