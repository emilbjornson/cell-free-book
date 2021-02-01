%This Matlab script can be used to reproduce Figure 4.1 in the monograph:
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

% SNR of AP 1
SNR_1dB = 10;
SNR_1 = 10^(SNR_1dB/10);

% SNR of AP 2
SNR_2dB = -30:1:30;
SNR_2 = 10.^(SNR_2dB/10);

% NMSE of UE when served by both APs according to (4.23)
NMSE = 1 - (SNR_1/(1 + 1/SNR_1) + SNR_2./(1 + 1./SNR_2))./(SNR_1 + SNR_2);

% NMSE of UE when served by only AP 1 according to (4.23)
NMSE_onlyAP1 = 1 - (SNR_1/(1 + 1/SNR_1))./(SNR_1);

% NMSE of UE when served by only AP 2 according to (4.23)
NMSE_onlyAP2 = 1 - (SNR_2./(1 + 1./SNR_2))./(SNR_2);


%% Plot simulation results
figure;
hold on; box on;
set(gca,'fontsize',16);
plot(SNR_2dB,NMSE,'k','LineWidth',2);
plot(SNR_2dB,NMSE_onlyAP1*ones(size(SNR_2dB)),'b-.','LineWidth',2);
plot(SNR_2dB,NMSE_onlyAP2,'r--','LineWidth',2);
xlabel('Effective SNR of AP 2 [dB]','Interpreter','Latex');
ylabel('NMSE','Interpreter','Latex');
set(gca,'YScale','lin');
set(gca,'XScale','lin');
legend({'Both APs','Only AP 1','Only AP 2'},'Interpreter','Latex','Location','SouthWest');
ylim([0 0.3])
