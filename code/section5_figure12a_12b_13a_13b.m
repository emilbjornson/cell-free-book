%This Matlab script can be used to reproduce Figures 5.12(a), 5.12(b), 5.13(a), and 5.13(b) in the monograph:
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

%Prepare to store DCC matrices 
D_tot = zeros(L,max(K),length(K),nbrOfSetups);

%% Go through all setups
for n = 1:nbrOfSetups

    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    %Go through all the UE numbers
    for n2=1:length(K)
        
        %Generate one setup with UEs and APs at random locations
        [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K(n2),N,tau_p,1,0,ASD,ASD);
        %Save DCC matrix
        D_tot(:,1:K(n2),n2,n) = D;

        %Remove large matrices at the end of analyzing this setup
        clear R;
    end



end

%Compute or prepare to store the total number of complex scalars sent from 
%the APs to the CPU over the fronthaul per coherence block according to
%Table 5.2
Cent_Front = tau_c*N*L*ones(1,length(K)); %Centralized
Dist_Front_all = (tau_c-tau_p)*L*K; %Distributed (All)
Dist_Front_DCC = zeros(1,length(K)); %Distributed (DCC)


%Compute or prepare to store the total number of complex scalars sent from 
%the APs to the CPU over the fronthaul for each realization channel statistics 
%according to Table 5.2
opt_LSFD_Front_all = (3*K+1)/2.*L.*K; %opt LSFD (All)
opt_LSFD_Front_DCC = zeros(1,length(K)); %opt LSFD (DCC)
nopt_LSFD_Front_DCC = zeros(1,length(K)); %n-opt LSFD (DCC)

%Compute or prepare to store the number of complex multiplications required 
%per coherence block to compute the combining vectors of all UEs for the
%centralized operation according to Table 5.1. The common operations will
%be counted once by modifying the terms in Table 5.1 accordingly.
MMSE_all_Multip = (N*tau_p+N^2)*K*L+((N*L)^2+N*L)/2*K+(N*L)^2*K+((N*L)^3-N*L)/3; %MMSE (All)
MMSE_DCC_Multip = zeros(1,length(K)); %MMSE (DCC)
PMMSE_DCC_Multip = zeros(1,length(K)); %P-MMSE (DCC)
PRZF_DCC_Multip = zeros(1,length(K)); %P-RZF (DCC)
MR_all_Multip = (N*tau_p+N^2)*K*L; %Centralized MR (All)
MR_DCC_Multip = zeros(1,length(K)); %Centralized MR (DCC)

%Compute or prepare to store the number of complex multiplications required 
%per coherence block to compute the combining vectors of all UEs for the
%dsitributed operation according to Table 5.3. The common operations will
%be counted once by modifying the terms in Table 5.3 accordingly.
LMMSE_all_Multip = (N*tau_p+N^2)*K*L+(N^2+N)/2*K*L+N^2*K*L+(N^3-N)/3*L; %L-MMSE (All)
LMMSE_DCC_Multip = zeros(1,length(K)); %L-MMSE (DCC)
LPMMSE_DCC_Multip = zeros(1,length(K)); %LP-MMSE (DCC)
MR_Dist_DCC_Multip = zeros(1,length(K)); %Distributed MR (DCC)

%Go through all the UE numbers
for n2 = 1:length(K)
    
    %For UE number K(n2), obtain the mean of \sum_{k=1}^K|M_k| 
    %(the mean of \sum_{l=1}^L|D_l|) over all APs and UEs.
    sumDl = mean(D_tot(:,1:K(n2),n2,:),'all')*L*K(n2);
    
    %Compute the fronthaul signaling load per coherence block load for 
    %"Distributed (DCC)" according to Table 5.2.
    Dist_Front_DCC(1,n2) = (tau_c-tau_p)*sumDl;
    
    %Compute the fronthaul signaling load per statistics for 
    %"opt LSFD (DCC)" according to Table 5.2.
    opt_LSFD_Front_DCC(1,n2) = (3*K(n2)+1)/2*sumDl;
    
    %Compute the number of complex multiplications for 
    %centralized MR combining according to Table 5.1
    MR_DCC_Multip(1,n2) = (N*tau_p+N^2)*sumDl;
    
    %Compute the number of complex multiplications for 
    %distributed MR combining according to Table 5.3
    MR_Dist_DCC_Multip(1,n2) = (N*tau_p+N^2)*sumDl;
    
    
    
    %Go through all the setups
    for n = 1:nbrOfSetups
        
        %Extract the DCC matrix for setup n and UE number K(n2)
        Dn = reshape(D_tot(:,1:K(n2),n2,n), [L, K(n2)]);
        
        %Number of APs that serve at least one UE
        L_used = length(find(sum(Dn,2)>=1));
        
        %Update the number of complex multiplication for MMSE combining
        %according to Table 5.1 (the common operations for the computation
        %of different UEs' combining vectors are counted once)
        MMSE_DCC_Multip(1,n2) = MMSE_DCC_Multip(1,n2) + ...
            (N*tau_p+N^2)*K(n2)*L_used + ...
            ((N*L_used)^2+N*L_used)/2*K(n2);
        
        %Update the number of complex multiplication for L-MMSE and LP-MMSE
        %combining according to Table 5.3 (the common operations for the
        %computation of different UEs' combining vectors are counted once)
        LMMSE_DCC_Multip(1,n2) = LMMSE_DCC_Multip(1,n2) +...
            (N*tau_p+N^2)*K(n2)*L_used + ...
            (N^2+N)/2*K(n2)*L_used+N^2*sumDl+(N^3-N)/3*L_used;
        
        LPMMSE_DCC_Multip(1,n2) = LPMMSE_DCC_Multip(1,n2) + ...
            (N*tau_p+N^2)*sumDl + ...
            (N^2+N)/2*sumDl+N^2*sumDl+(N^3-N)/3*L_used;
        
        %Go through all UEs
        for k = 1:K(n2)
            
            %Find the APs that serve UE k
            servingAPs = find(Dn(:,k)==1);
            %Find the UEs that are served partially by the same set of APs
            %as UE k, i.e., the set in (5.15)
            servedUEs = find(sum(Dn(servingAPs,:),1)>=1);
            %Compute the corresping number of UEs and APs for the above
            %sets
            Sk = length(servedUEs);
            La = length(servingAPs);

            %This is computed to count the number of operations once
            Bk = length(find(sum(Dn(:,servedUEs),2)>=1));

            %Update the fronthaul signaling load per statistics for 
            %"n-opt LSFD (DCC)" according to Table 5.2. 
            nopt_LSFD_Front_DCC(1,n2) = nopt_LSFD_Front_DCC(1,n2)+(3*Sk+1)/2*La/nbrOfSetups;
            
            %Update the number of complex multiplication for MMSE, P-MMSE,
            %and P-RZF combining with DCC in the centralized operation
            %according to Table 5.1 (the common operations for the computation
            %of different UEs' combining vectors are counted once)
            MMSE_DCC_Multip(1,n2) = MMSE_DCC_Multip(1,n2) + ...
                (N*La)^2+((N*La)^3-N*La)/3;
            
            PMMSE_DCC_Multip(1,n2) = PMMSE_DCC_Multip(1,n2) + ...
                ((N*tau_p+N^2)*Bk+...
                ((N*Bk)^2+N*Bk)/2+(N*La)^2+((N*La)^3-N*La)/3);
            
            PRZF_DCC_Multip(1,n2) = PRZF_DCC_Multip(1,n2) + ...
                Sk^2 + Sk*N*La + (Sk^3-Sk)/3;
        end
        
        %Go through all APs
        for l = 1:L
            
            %Find the UEs that are served by AP l
            servedUEs0 = find(Dn(l,:)==1);
            %The following is done to find the union of the sets in (5.15)
            %to count the common operation for different UEs once
            tempVec = [];
            for ind = 1:length(servedUEs0)
                k = servedUEs0(ind);
                servingAPs = find(Dn(:,k)==1);
                servedUEs = find(sum(Dn(servingAPs,:),1)>=1);
                tempVec = [tempVec servedUEs];
            end
            %The cardinality of the above set
            Cl = 0;
            for k = 1:K(n2)
                if sum(tempVec==k)>= 1
                    Cl = Cl + 1;
                end
            end
            
            %Update the number of complex multiplication for P-RZF
            %combining with DCC in the centralized operation
            %according to Table 5.1 (the common operations for the computation
            %of different UEs' combining vectors are counted once)
            PRZF_DCC_Multip(1,n2) = PRZF_DCC_Multip(1,n2)+ ...
                (N*tau_p+N^2)*Cl + (Cl^2+Cl)/2*N; ...
        end
    
    end
end

%Obtain the averages over different random setups
MMSE_DCC_Multip = MMSE_DCC_Multip/nbrOfSetups;
PMMSE_DCC_Multip =  PMMSE_DCC_Multip/nbrOfSetups;
PRZF_DCC_Multip = PRZF_DCC_Multip/nbrOfSetups;
LMMSE_DCC_Multip = LMMSE_DCC_Multip/nbrOfSetups;
LPMMSE_DCC_Multip = LPMMSE_DCC_Multip/nbrOfSetups;

%% Plot simulation results
% Plot Figure 5.12(a)
figure;
hold on; box on;
set(gca,'fontsize',16);


plot(K, Cent_Front/L,'kd-','LineWidth',2);
plot(K, Dist_Front_all/L,'bh--','LineWidth',2);
plot(K, Dist_Front_DCC/L,'rs:','LineWidth',2);

xlabel('Number of UEs, $K$','Interpreter','Latex');
ylabel('Number of complex scalars (data+pilot)','Interpreter','Latex');
legend({'Centralized','Distributed (All)','Distributed (DCC)'},'Interpreter','Latex','Location','SouthEast');

% Plot Figure 5.12(b)
figure;
hold on; box on;
set(gca,'fontsize',16);


plot(K, opt_LSFD_Front_all/L,'kd-','LineWidth',2);
plot(K, opt_LSFD_Front_DCC/L,'bh--','LineWidth',2);
plot(K, nopt_LSFD_Front_DCC/L,'rs:','LineWidth',2);

xlabel('Number of UEs, $K$','Interpreter','Latex');
ylabel('Number of complex scalars (statistics)','Interpreter','Latex');
legend({'opt LSFD (All)','opt LSFD (DCC)','n-opt LSFD (DCC)'},'Interpreter','Latex','Location','SouthEast');

% Plot Figure 5.13(a)
figure;
hold on; box on;
set(gca,'fontsize',16);


plot(K, MMSE_all_Multip/L,'kd-','LineWidth',2);
plot(K, MMSE_DCC_Multip/L,'rh--','LineWidth',2);
plot(K, PMMSE_DCC_Multip/L,'bs--','LineWidth',2);
plot(K, PRZF_DCC_Multip/L,'bs-','LineWidth',2);
plot(K, MR_all_Multip/L,'kd:','LineWidth',2);
plot(K, MR_DCC_Multip/L,'ko--','LineWidth',2);

xlabel('Number of UEs, $K$','Interpreter','Latex');
ylabel('Number of complex multiplications','Interpreter','Latex');
legend({'MMSE (All)','MMSE (DCC)','P-MMSE (DCC)','P-RZF (DCC)', 'MR (All)','MR (DCC)'},'Interpreter','Latex','Location','SouthEast');


% Plot Figure 5.13(b)
figure;
hold on; box on;
set(gca,'fontsize',16);


plot(K, LMMSE_all_Multip/L,'kd-','LineWidth',2);
plot(K, LMMSE_DCC_Multip/L,'bh--','LineWidth',2);
plot(K, LPMMSE_DCC_Multip/L,'rs--','LineWidth',2);
plot(K, MR_Dist_DCC_Multip/L,'ko:','LineWidth',2);

xlabel('Number of UEs, $K$','Interpreter','Latex');
ylabel('Number of complex multiplications','Interpreter','Latex');
legend({'L-MMSE (All)','L-MMSE (DCC)','LP-MMSE (DCC)','MR (DCC)'},'Interpreter','Latex','Location','SouthEast');