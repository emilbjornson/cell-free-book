%This Matlab script can be used to reproduce Figure 4.5 in the monograph:
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

%Distance of interfering UE
d = 10:0.01:500;


%Define a function to compute the SNR as function of the horizontal distance 
%measured in meter. The AP is 10 meter above the UE. The UE transmit power
%is 10 dBm (10 mW). The noise variance is -96 dBm and the large-scale fading coefficient
%is computed according to (1.1).
SNR = @(hor_dist) db2pow(10+96-30.5-36.7*log10(sqrt(hor_dist.^2+10^2)));


%Compute SNRs at different distances
SNR1 = pow2db(SNR(d)); %tau_p = 1
SNR10 = 10+SNR1; %tau_p = 10


%% Plot the simulation results
figure;
hold on; box on;
set(gca,'fontsize',16);
plot(d, SNR10,'b-.','LineWidth',2);
plot(d, SNR1,'k-','LineWidth',2);
plot(d, -10*ones(size(d)),'k:','LineWidth',2);
xlabel('Distance to the interfering UE [m]','Interpreter','Latex');
ylabel('Effective SNR [dB]','Interpreter','Latex');
legend({'$\tau_p = 10$','$\tau_p = 1$'},'Interpreter','Latex','Location','NorthEast');
set(gca,'YScale','lin');
ylim([-30 40]);
