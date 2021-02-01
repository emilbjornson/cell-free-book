%This Matlab script can be used to reproduce Figures 6.8(a) and 6.8(b) in the monograph:
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
nbrOfSetups = 100;

%Number of channel realizations per setup
nbrOfRealizations = 1000;


%Number of APs 
L = 100;

%Number of antennas per AP
N = 4;

%Number of UEs in the network
K = 20:20:100;

%Length of coherence block
tau_c = 200;

%Length of pilot sequences
tau_p = 10;

%Angular standard deviation in the local scattering model (in radians)
ASD = deg2rad(15);


%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per AP (mW)
rho_tot = 200;


%Prepare to save simulation results (average SEs)
Avg_SE_PMMSE = zeros(length(K),1); %P-MMSE
Avg_SE_PRZF = zeros(length(K),1); %P-RZF
Avg_SE_LPMMSE = zeros(length(K),1); %LP-MMSE
Avg_SE_MR = zeros(length(K),1); %MR

%Prepare to store DCC matrices
D_tot = zeros(L,max(K),length(K),nbrOfSetups); 

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Go through all UE number
    for n2=1:length(K)
        
        %Generate one setup with UEs and APs at random locations
        [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K(n2),N,tau_p,1,0,ASD,ASD);
        %Save the DCC matrix
        D_tot(:,1:K(n2),n2,n) = D;
        %Generate channel realizations with estimates and estimation
        %error correlation matrices
        [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K(n2),N,tau_p,pilotIndex,p);
        
        %Compute the power allocation in (6.36) for distributed precoding
        rho_dist = zeros(L,K(n2));
        
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
        
        %% Cell-Free Massive MIMO with DCC
        
        %Compute SEs for DCC case
        [SE_MMSE, SE_P_MMSE, SE_P_RZF,  ...
            SE_L_MMSE, SE_LP_MMSE, SE_MR, ...
            Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR] ...
            = functionComputeSE_downlink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K(n2),L,p,R,pilotIndex,rho_dist,gainOverNoisedB,rho_tot);
        
        
        %Update average SE values
        Avg_SE_PMMSE(n2,1) = Avg_SE_PMMSE(n2,1) + sum(SE_P_MMSE,'all')/(nbrOfSetups*K(n2));
        Avg_SE_PRZF(n2,1) =  Avg_SE_PRZF(n2,1) +  sum(SE_P_RZF,'all')/(nbrOfSetups*K(n2));
        Avg_SE_LPMMSE(n2,1) = Avg_SE_LPMMSE(n2,1) +  sum(SE_LP_MMSE,'all')/(nbrOfSetups*K(n2));
        Avg_SE_MR(n2,1) = Avg_SE_MR(n2,1) + sum(SE_MR,'all')/(nbrOfSetups*K(n2));
        
        
        
        %Remove large matrices at the end of analyzing this setup
        clear Hhat H B C R;
    end
    
    
    
end


%% Plot simulation results
% Plot Figure 6.8(a)
figure;
hold on; box on;
set(gca,'fontsize',16);


plot(K, Avg_SE_PMMSE','kd-','LineWidth',2);
plot(K, Avg_SE_PRZF','bh--','LineWidth',2);
plot(K, Avg_SE_LPMMSE','rs:','LineWidth',2);
plot(K, Avg_SE_MR','ko--','LineWidth',2);

xlabel('Number of UEs, $K$','Interpreter','Latex');
ylabel('Average spectral efficiency [bit/s/Hz]','Interpreter','Latex');
legend({'P-MMSE','P-RZF','LP-MMSE','MR'},'Interpreter','Latex','Location','NorthEast');


% Plot Figure 6.8(b)
figure;
hold on; box on;
set(gca,'fontsize',16);


plot(K, Avg_SE_PMMSE.*K','kd-','LineWidth',2);
plot(K, Avg_SE_PRZF.*K','bh--','LineWidth',2);
plot(K, Avg_SE_LPMMSE.*K','rs:','LineWidth',2);
plot(K, Avg_SE_MR.*K','ko--','LineWidth',2);

xlabel('Number of UEs, $K$','Interpreter','Latex');
ylabel('Sum spectral efficiency [bit/s/Hz]','Interpreter','Latex');
legend({'P-MMSE','P-RZF','LP-MMSE','MR'},'Interpreter','Latex','Location','NorthWest');


