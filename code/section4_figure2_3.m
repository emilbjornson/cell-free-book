%This Matlab script can be used to reproduce Figures 4.2 and 4.3 in the monograph:
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

%Set the side length of the simulation area
squareLength = 400;

%Total number of antennas in all setups
nbrOfAntennas = 64;

%Set the AP location of the Cellular Massive MIMO setup
APcellular = squareLength/2 + 1i*squareLength/2;

%Set the AP locations for the small-cell and cell-free setups
APperdim = sqrt(nbrOfAntennas);
APcellfree = linspace(squareLength/APperdim,squareLength,APperdim)-squareLength/APperdim/2;
APcellfree = repmat(APcellfree,[APperdim 1]) + 1i*repmat(APcellfree,[APperdim 1])';

%Number of realizations of the random UE locations
nbrOfSetups = 100000;

%Number of UEs in the simulation setup
K = 1;

%Set the pilot length and coherence block length
tau_p = 10;

%Number of APs serving the UE in cell-free setup
Mk = 64;

%Generate the random UE locations for all setups
UElocations = (rand(nbrOfSetups,K)+1i*rand(nbrOfSetups,K))*squareLength;

%Range of transmit powers (in mW)
p = logspace(0,3,100);
pValue = 10; %Value for CDF curves

%Define a function to compute the SNR as function of the horizontal distance
%measured in meter. The AP is 10 meter above the UE. The noise variance is
%-96 dBm and the large-scale fading coefficient is computed according to
%(1.1).
SNR = @(hor_dist,power) db2pow(pow2db(power)+96-30.5-36.7*log10(sqrt(hor_dist.^2+10^2)));


%Prepare to store simulation results
NMSE_cellular_mMIMO = zeros(nbrOfSetups,length(p));
NMSE_cellular_small = zeros(nbrOfSetups,length(p));
NMSE_cellfree = zeros(nbrOfSetups,length(p));

SNR_cellular_mMIMO = zeros(nbrOfSetups,1);
SNR_cellular_small = zeros(nbrOfSetups,1);
SNR_cellfree = zeros(nbrOfSetups,1);


%% Go through all random realizations of the UE locations
for n = 1:nbrOfSetups
    
    for m = 1:length(p)
        
        
        %Compute the SNR and effective SNR in (4.15) in the Cellular Massive MIMO setup
        SNRCellular = SNR(abs(APcellular-UElocations(n,1)),p(m));
        SNRCellular_effective = SNRCellular*tau_p;
        
        %Compute the NMSE for i.i.d. fading according to (4.21)
        NMSE_cellular_mMIMO(n,m) = 1/(SNRCellular_effective+1);
        
        
        %Compute the SNRs and effective SNRs in the cell-free and small-cell setup
        SNRCellfree = SNR(abs(APcellfree(:) - UElocations(n,1)),p(m));
        SNRCellfree = sort(SNRCellfree,'descend');
        SNRCellfree_effective = SNRCellfree*tau_p;
        
        
        
        %Compute the NMSE in (4.23) for the cell-free setup when served by only
        %eight APs
        NMSE_cellfree(n,m) = 1 - sum(SNRCellfree_effective(1:Mk)./(1+1./SNRCellfree_effective(1:Mk)))/sum(SNRCellfree_effective(1:Mk));
        

        %Compute the NMSE in the small-cell setup when only served by one AP
        %according to (4.23)
        NMSE_cellular_small(n,m) = 1/(SNRCellfree_effective(1)+1);
        
        
        if p(m) == pValue
        
            %Compute SNR in the Cellular Massive MIMO setup using (4.26)
            %with i.i.d. fading
            SNRCellular = SNR(abs(APcellular-UElocations(n,1)),pValue);
            SNR_cellular_mMIMO(n,1) = (nbrOfAntennas+1)*SNRCellular*(1-NMSE_cellular_mMIMO(n,m));
            
            %Compute the SNRs in the cell-free and small-cell setups using
            %(4.26) with single antenna case
            SNRCellfree = SNR(abs(APcellfree(:) - UElocations(n,1)),pValue);
            SNRCellfree = sort(SNRCellfree,'descend');
            
            RCterm = SNRCellfree(1:Mk).*SNRCellfree_effective(1:Mk)./(1+SNRCellfree_effective(1:Mk));
            
            SNR_cellfree(n,1) = sum(RCterm.^2)/sum(RCterm) +sum(RCterm);
            
            SNR_cellular_small(n,1) = 2*SNRCellfree(1)*(1-NMSE_cellular_small(n,m));
            
        end
    
    end
    
end




%% Plot simulation results
figure;
hold on; box on;
plot(p,mean(NMSE_cellular_mMIMO,1),'k','LineWidth',2);
plot(p,mean(NMSE_cellular_small,1),'r-.','LineWidth',2);
plot(p,mean(NMSE_cellfree,1),'b--','LineWidth',2);
xlabel('Transmit power [mW]','Interpreter','latex');
ylabel('Average NMSE','Interpreter','latex');
legend({'Cellular: Massive MIMO','Cellular: Small cells','Cell-free'},'Interpreter','latex','Location','SouthWest');
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'fontsize',16);
ylim([1e-5 1]);
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1]);

figure;
hold on; box on;
plot(pow2db(sort(SNR_cellfree(:),'ascend')),linspace(0,1,nbrOfSetups),'b--','LineWidth',2);
plot(pow2db(sort(SNR_cellular_small(:),'ascend')),linspace(0,1,nbrOfSetups),'r-.','LineWidth',2);
plot(pow2db(sort(SNR_cellular_mMIMO(:),'ascend')),linspace(0,1,nbrOfSetups),'k','LineWidth',2);
xlabel('SNR of UE $k$ [dB]','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend({'Cell-free','Cellular: Small cells','Cellular: Massive MIMO'},'Interpreter','latex','Location','SouthEast');
set(gca,'fontsize',16);
xlim([0 60]);
