%This Matlab script can be used to reproduce Figure 7.3 in the monograph:
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
ASD_varphi = deg2rad(15); %azimuth angle
ASD_theta = deg2rad(15);  %elevation angle

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per AP (mW)
rho_tot = 200;


%Prepare to save simulation results

SE_DL_LPMMSE_equal = zeros(K,nbrOfSetups); %Equal
SE_DL_LPMMSE_fractional = zeros(K,nbrOfSetups); %FPA, \upsilon = 0.5
SE_DL_LPMMSE_fractional2 = zeros(K,nbrOfSetups); %FPA, \upsion = -0.5
SE_DL_LPMMSE_maxmin = zeros(K,nbrOfSetups); %MMF
SE_DL_LPMMSE_sumSE = zeros(K,nbrOfSetups); %SumSE


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
   
    gainOverNoise = db2pow(gainOverNoisedB);

    %Equal power allocation
    rho_dist_equal = (rho_tot/tau_p)*ones(L,K);

    %Compute the power allocation in (7.47) for distributed precoding
    rho_dist = zeros(L,K); % with exponent 0.5
    rho_dist2 = zeros(L,K); % with exponent -0.5

    
    for l = 1:L
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute denominator in (7.47)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
        normalizationAPl2 = sum(1./sqrt(gainOverNoise(l,servedUEs)));

        for ind = 1:length(servedUEs)
            
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
            rho_dist2(l,servedUEs(ind)) = rho_tot/sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl2;

        end
        
    end
    
    %Obtain the expectations for the computation of the terms in
    %(7.25)-(7.26)
    [signal_P_MMSE, signal2_P_MMSE, scaling_P_MMSE,...
        signal_P_RZF, signal2_P_RZF, scaling_P_RZF,...
        signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
     functionComputeExpectations(Hhat,H,D,C,nbrOfRealizations,N,K,L,p_full);
    
    %Prepare arrays to store the vectors \tilde{b}_k in (7.25) and matrices
    %\tilde{C}_{ki} in (7.26)
    bk = zeros(L,K);
    Ck = zeros(L,L,K,K);
  
    %Go through all UEs
    for k = 1:K
        %Find the APs that serve UE k
        servingAPs = find(D(:,k)==1);
        %The number of APs that serve UE k
        La = length(servingAPs);
        %Compute the vector in (7.25) for UE k (only the non-zero indices correspondig to 
        %serving APs are considered)
        bk(1:La,k) = real(vec(signal_LP_MMSE(k,k,servingAPs)))./sqrt(scaling_LP_MMSE(servingAPs,k));
        
        %Go through all UEs
        for i = 1:K
            %Find the APs that serve UE i
            servingAPs = find(D(:,i)==1);
            %The number of APs that serve UE i
            La = length(servingAPs);
            %Compute the matrices in (7.26) (only the non-zero indices are
            %considered)
            if i==k
               Ck(1:La,1:La,k,k) = bk(1:La,k)*bk(1:La,k)';
            else
               Ck(1:La,1:La,k,i) = diag(1./sqrt(scaling_LP_MMSE(servingAPs,i)))...
                   *(vec(signal_LP_MMSE(k,i,servingAPs))...
                   *vec(signal_LP_MMSE(k,i,servingAPs))')...
                   *diag(1./sqrt(scaling_LP_MMSE(servingAPs,i)));

            
            end
            
            for j = 1:La
                Ck(j,j,k,i) = signal2_LP_MMSE(k,i,servingAPs(j))/scaling_LP_MMSE(servingAPs(j),i);
            end
        end
     

    end
    
   
    
    %Take the real part (in the SINR expression,the imaginary terms cancel
    %each other)
    Ck = real(Ck);
    
    %Compute hte square roots of the power allocation coefficients
    %corresponding to (7.24)
    tilrho = sqrt(rho_dist_equal);
    tilrho1 = sqrt(rho_dist);
    tilrho2 = sqrt(rho_dist2);

    %Go through all UEs
    for k = 1:K
        %Find APs that serve UE k
        servingAPs = find(D(:,k)==1);
        %The number of APs that serve UE k
        La = length(servingAPs);
        
        %Compute the numerator and denominator of (7.23) for equal and FPA
        %schemes with two different exponents
        numm = abs(bk(1:La,k)'*tilrho(servingAPs,k))^2;
        denomm = 1-numm;
        
        numm1 = abs(bk(1:La,k)'*tilrho1(servingAPs,k))^2;
        denomm1 = 1-numm1;
        
        numm2 = abs(bk(1:La,k)'*tilrho2(servingAPs,k))^2;
        denomm2 = 1-numm2;
        
        for i = 1:K
            servingAPs = find(D(:,i)==1);
            La = length(servingAPs);
            denomm = denomm+tilrho(servingAPs,i)'*Ck(1:La,1:La,k,i)*tilrho(servingAPs,i);
            denomm1 = denomm1+tilrho1(servingAPs,i)'*Ck(1:La,1:La,k,i)*tilrho1(servingAPs,i);
            denomm2 = denomm2+tilrho2(servingAPs,i)'*Ck(1:La,1:La,k,i)*tilrho2(servingAPs,i);


        end
        %Compute SEs using SINRs in (7.23) and Corollary 6.3 for equal and
        %FPA schemes with two different exponents
        SE_DL_LPMMSE_equal(k,n) = preLogFactor*log2(1+numm/denomm);
        SE_DL_LPMMSE_fractional(k,n) = preLogFactor*log2(1+numm1/denomm1);
        SE_DL_LPMMSE_fractional2(k,n) = preLogFactor*log2(1+numm2/denomm2);
    end
    
    %Compute SE according to Corollary 6.3 with max-min fair power
    %allocation in Algorithm 7.5
    SE_DL_LPMMSE_maxmin(:,n) = functionDownlinkSE_maxmin_dist(bk,Ck,preLogFactor,L,K,D,rho_tot);
    
    
    %Compute SE according to Corollary 6.3 with sum SE maximizing power
    %allocation in Algorithm 7.6
    SE_DL_LPMMSE_sumSE(:,n) =  functionDownlinkSE_sumSE_dist(bk,Ck,preLogFactor,L,K,D,rho_tot,tau_p);

    
end

% Plot Figure 7.3
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_DL_LPMMSE_equal(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_DL_LPMMSE_fractional(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
ppp5 = plot(sort(SE_DL_LPMMSE_fractional2(:)),linspace(0,1,K*nbrOfSetups),'k:o','LineWidth',2);
ppp5.MarkerSize = 6;
ppp5.MarkerIndices = 1:ceil(K*nbrOfSetups/7):K*nbrOfSetups;
plot(sort(SE_DL_LPMMSE_maxmin(:)),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);
plot(sort(SE_DL_LPMMSE_sumSE(:)),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);


xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Equal', 'FPA, $\upsilon=0.5$', 'FPA, $\upsilon=-0.5$','MMF','SumSE' },'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);