%This Matlab script can be used to reproduce Figure 4.4 in the monograph:
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

%Number of AP antennas
N = 1;

% Define the SNR values of the desired UE
SNR1dB = [0, 10, 20];
SNR1 = 10.^(SNR1dB/10);

% Define the SNR values of the interfering UE
SNR2dB = 20:-0.05:-30;
SNR2 = 10.^(SNR2dB/10);
        

%Prepare to save the NMMSE results
NMSE = zeros(length(SNR1dB), length(SNR2dB));
NMSE_uncorr = zeros(length(SNR1dB), length(SNR2dB));


%% Compute the NMSE in the presence of pilot contamination

%Spatial correlation matrix of the desired UE
R1 = 1;

%Spatial correlation matrix of the interfering UE
R2 = 1;

for n1 = 1:length(SNR1)
    
    for n2 = 1:length(SNR2)
        
        %Compute the NMSE using (4.27), i.e., the simplified version of (4.12)
        NMSE(n1,n2) = real(trace(R1 - SNR1(n1)*R1*((SNR1(n1)*R1+SNR2(n2)*R2+eye(N))\R1)))/trace(R1);
        
    end
    
end



%% Consider the case without pilot contamination
NMSE_no_contamination = zeros(length(SNR1),1);

for n1 = 1:length(SNR1)
    
    %Compute the NMSE using (4.27)
    NMSE_no_contamination(n1) = 1 - SNR1(n1)/(SNR1(n1)+1);
    
end


%% Plot the simulation results
figure;
hold on; box on;
set(gca,'fontsize',16);
plot(SNR2dB, NMSE(1,:),'k-','LineWidth',2);
plot(SNR2dB, NMSE(2,:),'r--','LineWidth',2);
plot(SNR2dB, NMSE(3,:),'b-.','LineWidth',2);


plot(SNR2dB, NMSE_no_contamination(1)*ones(size(SNR2dB)),'k:','LineWidth',2);
plot(SNR2dB, NMSE_no_contamination(2)*ones(size(SNR2dB)),'r:','LineWidth',2);
plot(SNR2dB, NMSE_no_contamination(3)*ones(size(SNR2dB)),'b:','LineWidth',2);

xlabel('SNR of interfering signal [dB]','Interpreter','latex');
ylabel('NMSE','Interpreter','latex');
set(gca,'YScale','log');
legend({'SNR = 0 dB','SNR = 10 dB','SNR = 20 dB'},'Interpreter','latex','Location','SouthEast');
ylim([1e-3 1e0]);
