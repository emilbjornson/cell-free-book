%This Matlab script can be used to reproduce Figure 1.12 in the monograph:
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
clear

%Set the wavelength (cm)
wavelength = 10;

%Side length of the simulation area
squareSide = 10*wavelength;

%Generate the grid of points where the normalized SNR is computed
x = 4:1:squareSide-4;
y = 4:1:squareSide-4;


%Put out the antennas
pointsperdim = 10; %Number of antennas per side in the square
antenna_locations= squareSide*[linspace(0,1,pointsperdim+1) 1i*linspace(1/pointsperdim,1,pointsperdim) 1i+linspace(1/pointsperdim,1,pointsperdim) 1+1i*linspace(1/pointsperdim,1-1/pointsperdim,pointsperdim-1)];


%Select the target location of the transmission
targetlocation = round(squareSide/4)+1i*round(squareSide/2);


%Prepare to save simulation results
phaseshifts = zeros(length(x),length(y),length(antenna_locations));
distances = zeros(length(x),length(y),length(antenna_locations));
angletarget = zeros(length(x),length(y),length(antenna_locations));


%% Go through all antennas
for n = 1:length(antenna_locations)
    
    %Go through all spatial sample points (potential UE locations)
    for k = 1:length(x)
        
        for j = 1:length(y)
            
            %Compute the distance from the transmit antenna to the sample point
            distances(k,j,n) = abs(x(k)+1i*y(j) - antenna_locations(n));
            
            %Compute the phase-shift from the antenna to the sample point
            phaseshifts(k,j,n) = 2*pi*distances(k,j,n)/wavelength;
            
        end
        
    end
    
    %Compute the phase-shift from the antenna to the intended target point
    angletarget(:,:,n) = 2*pi*abs(targetlocation - antenna_locations(n))/wavelength;
    
end

%Compute the distances from all antennas to the target point
distances_target = abs(antenna_locations-targetlocation);


%Compute the combined channel gain from all antennas to the different
%sample points. We are considering a free-space propagation pathloss
channel_gain = abs(sum(exp(1i*(phaseshifts-angletarget))./sqrt(4*pi*distances.^2),3)).^2;

%Compute the combined channel gain from all antennas to the target point
channel_gain_target = abs(sum(exp(1i*0)./sqrt(4*pi*distances_target.^2))).^2;

%Set a normalization factor for the simulations
normalization = channel_gain_target;


%% Plot simulation results

%Part a
figure;

surf(y/wavelength,x/wavelength,channel_gain/normalization);
colormap(flipud(parula));
colorbar;
shading interp;
hold on;

for n = 1:length(antenna_locations)
    plot3(imag(antenna_locations(n))/wavelength,real(antenna_locations(n))/wavelength,0,'k*','MarkerSize',5);
end

plot3(imag(targetlocation)/wavelength,real(targetlocation)/wavelength,(channel_gain_target)/normalization,'ko','MarkerSize',5,'MarkerFaceColor','k');
hold off;

xlabel('Distance (wavelengths)','Interpreter','latex');
ylabel('Distance (wavelengths)','Interpreter','latex');
zlabel('Normalized SNR','Interpreter','latex');
set(gca,'fontsize',16);
xlim([0 squareSide/wavelength]);
ylim([0 squareSide/wavelength]);
view([-20 26]);


%Part b
figure;

surf(y/wavelength,x/wavelength,channel_gain/normalization);
colormap(flipud(parula));
colorbar;
shading interp;
hold on;

for n = 1:length(antenna_locations)
    plot3(imag(antenna_locations(n))/wavelength,real(antenna_locations(n))/wavelength,0,'k*','MarkerSize',5);
end

plot3(imag(targetlocation)/wavelength,real(targetlocation)/wavelength,(channel_gain_target)/normalization,'ko','MarkerSize',5,'MarkerFaceColor','k');
hold off;

xlabel('Distance (wavelengths)','Interpreter','latex');
ylabel('Distance (wavelengths)','Interpreter','latex');
zlabel('Normalized SNR','Interpreter','latex');
set(gca,'fontsize',16);
view(2);
xlim([0 squareSide/wavelength]);
ylim([0 squareSide/wavelength]);
axis equal;
