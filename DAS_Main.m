%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main File for Implementing the DAS RF beamforming                      %
%                                                                        %                          
% Authors: Gayathri Malamal (121814001@smail.iitpkd.ac.in)               %
% Center for Computational Imaging, IIT Palakkad                         %
% G. Malamal and M. R. Panicker, “On the Physics of Ultrasound           %
% Transmission for In-Plane Needle Tracking in Guided Interventions,”    %
% in Biomedical Physics and Engineering Express.                         %
% https://doi.org/10.1088/2057-1976/acc338.                              %
% Date: 31-03-2023                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load dataset and variables

clc;
clearvars;

% Add the dataset to be beamformed to the present working directory

% Set the type of dataset
acq = 'PW'; % Uncomment for PW datasets
%acq = 'STA'; % Uncomment for STA datasets

% Load the dataset
Files=dir(fullfile(pwd,'*.mat'));
load(Files(1).name);

probe_geometry = dataset.probe_geometry(:,1);
Fnum = 1.5;
rawDataRF = double(dataset.rawData);
timeVector=(0:size(rawDataRF,1)-1)'/(dataset.Fs*1e6);
z_axis=0.5*1540*timeVector;
[x_grid,z_grid]=meshgrid(probe_geometry,z_axis);

%% DAS RF Beamforming

for txIdx=1:size(rawDataRF,3)
    clc
    disp(strcat('DAS Beamforming...',num2str(round(txIdx/size(rawDataRF,3)*100,2)),'%'))
  
    if(strcmp(acq, 'PW'))
       % tx = -dataset.angles(txIdx); % Uncomment for k-Wave simulated data
         tx = dataset.angles(txIdx); % Uncomment for experimental data
          
    elseif(strcmp(acq, 'STA'))
         tx = probe_geometry(txIdx);
    end
    
    beamformedDataDAS(:,:,txIdx)=DAS_RF(acq,rawDataRF(:,:,txIdx), timeVector, x_grid, z_grid, probe_geometry, tx,1540*ones(size(x_grid)), Fnum);
    
end

clc
disp('******DAS Beamforming Completed******')

%% Display the image
beamformedDataDASSum= single(sum(beamformedDataDAS, 3));

envelopeDAS = abs(hilbert(beamformedDataDASSum));
beamformedDataDASImage = envelopeDAS./max(envelopeDAS(:));

figure,imagesc(probe_geometry.*100, z_axis.*100, 20*log10(beamformedDataDASImage));

colormap(gray);
colorbar;
vrange = [-60 0];
caxis(vrange);
shading('interp');
set(gca,'TickLabelInterpreter','latex')
colorbar('TickLabelInterpreter','latex');
hold on;
xlabel('\bf{x (cm)}','interpreter','latex');
ylabel('\bf{z (cm)}','interpreter','latex');
title('\bf{DAS IMAGE}');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',24);
set(gcf, 'Position',  [100, 100, 600, 600])
hold off;


