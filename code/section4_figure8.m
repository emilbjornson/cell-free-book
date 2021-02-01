%This Matlab script can be used to reproduce Figure 4.8 in the monograph:
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
N = 4;

%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(10);
ASD_theta = 0;

%Range of nominal angles of the desired UE
varphiDesiredDegrees = 30;
varphiDesiredRadians = deg2rad(varphiDesiredDegrees);

%Range of nominal angles of the interfering UE
varphiInterfererDegrees = linspace(-60,60,73);
varphiInterfererRadians = deg2rad(varphiInterfererDegrees);

%Elevation angle
theta = -pi/12;


% Define the SNR values of the desidered UE
SNR1dB = 10;
SNR1 = 10.^(SNR1dB/10);

% Define the SNR values of the interfering UE
SNR2dB = [10,-10];
SNR2 = 10.^(SNR2dB/10);

%Prepare to save the NMMSE results
NMSE_corr = zeros(length(varphiInterfererRadians),length(SNR2dB));
NMSE_uncorr = zeros(2,length(SNR2dB));

%% Go through all setups

%Compute the spatial correlation matrix of the interfering UE 
%using local scattering model in (2.18) with Gaussian angular distribution
R1 = functionRlocalscattering(N,varphiDesiredRadians,theta,ASD_varphi,ASD_theta);

for s2 = 1:length(varphiInterfererRadians)
    
    %Compute the spatial correlation matrix of the interfering UE
    R2 = functionRlocalscattering(N,varphiInterfererRadians(s2),theta,ASD_varphi,ASD_theta);
        
        for n2 = 1:length(SNR2)
            
            %Compute the NMSE according to (4.31)
            NMSE_corr(s2,n2) = real(trace(R1 - SNR1*R1*((SNR1*R1+SNR2(n2)*R2+eye(N))\R1)))/trace(R1);
            
        end

end

%Go through all number of antennas
for n2 = 1:length(SNR2)
    
    %Compute the NMSE when having uncorrelated fading according to (4.31)
    NMSE_uncorr(:,n2) = 1 - SNR1/(SNR1+SNR2(n2)+1);
    
end


%% Plot the simulation results
figure;
hold on; box on;
set(gca,'fontsize',16);
plot(varphiInterfererDegrees, NMSE_corr(:,1),'b--','LineWidth',2);
plot(varphiInterfererDegrees, NMSE_corr(:,2),'r-','LineWidth',2);

plot([varphiInterfererDegrees(1); varphiInterfererDegrees(end)],NMSE_uncorr(:,1),'k:','LineWidth',2);
plot([varphiInterfererDegrees(1); varphiInterfererDegrees(end)],NMSE_uncorr(:,2),'k:','LineWidth',2);

xlabel('Angle of interfering UE [degree]','Interpreter','latex');
ylabel('NMSE','Interpreter','latex');
ylabel('NMSE');
set(gca,'YScale','log');
legend({'Same SNR','20 dB weaker'},'Interpreter','latex','Location','SouthWest');
xlim([-60 60]);
ylim([1e-2 1e0]);
