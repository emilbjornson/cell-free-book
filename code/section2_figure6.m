%This Matlab script can be used to reproduce Figure 2.6 in the monograph:
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


%Number of AP antennas
N = 8;

%Set the azimuth angle of the UE in radians
varphi = pi/6;

%Set the elevation angle of the UE in radians
theta = -pi/12;

%Set the ASD for azimuth angle deviation in radians
ASD_varphi = deg2rad([5 10 20]);

%Set the ASD for elevation angle deviation in radians
ASD_theta = deg2rad([5 10 20]);



%Compute spatial correlation matrix with local scattering model and
%different Gaussian angular distributions 
R_ULA_5 = functionRlocalscattering(N,varphi,theta,ASD_varphi(1), ASD_theta(1));
R_ULA_10 = functionRlocalscattering(N,varphi,theta,ASD_varphi(2), ASD_theta(2));
R_ULA_20 = functionRlocalscattering(N,varphi,theta,ASD_varphi(3), ASD_theta(3));


%Channel correlation matrix with uncorrelated fading
R_uncorrelated = eye(N);


%Extract the eigenvalues and place them in decreasing order
eigenvalues_ULA_5 = flipud(eig(R_ULA_5));
eigenvalues_ULA_10 = flipud(eig(R_ULA_10));
eigenvalues_ULA_20 = flipud(eig(R_ULA_20));
eigenvalues_uncorr = flipud(eig(R_uncorrelated));



%% Plot the simulation results
figure;
hold on; box on;

plot(1:N,10*log10(eigenvalues_ULA_5),'rs--','LineWidth',2);
plot(1:N,10*log10(eigenvalues_ULA_10),'bd-.','LineWidth',2);
plot(1:N,10*log10(eigenvalues_ULA_20),'ko-','LineWidth',2);
plot(1:N,10*log10(eigenvalues_uncorr),'k:','LineWidth',2);

xlabel('Eigenvalue number','Interpreter','latex');
ylabel('Normalized eigenvalue [dB]','Interpreter','latex');
legend({'ASD $=5^\circ$','ASD $=10^\circ$','ASD $=20^\circ$','Uncorrelated'},'Location','SouthWest','Interpreter','latex');
ylim([-50 10]);
set(gca,'fontsize',16);
