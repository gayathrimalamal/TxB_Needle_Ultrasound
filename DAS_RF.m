%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements the DAS RF beamforming with pixel wise addition             %
%                                                                        %
% This is inspired from DAS in USTB                                      %
% UltraSound ToolBox (https://www.ustb.no/)                              %
% Authors: Gayathri Malamal and Mahesh Raveendranatha Panicker           %
% Date: 31-03-2023                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function beamformedData = DAS_RF(acq,rawData,timeVector,x_grid,z_grid,probe_geometry,tx,c_map, Fnum)

N_c=size(rawData,2);
rows=size(x_grid,1);
columns=size(x_grid,2);
beamformedData = zeros(rows,columns);
delay_compensation=zeros(rows*columns,N_c);

%% depth based apodization
rx_f_number=Fnum;
rx_aperture=z_grid(:)/rx_f_number;
aperture=rx_aperture*ones(1,N_c);
rx_aperture_distance=abs(x_grid(:)*ones(1,N_c)-ones(rows*columns,1)*probe_geometry(:,1).'); 

%hanning apodization
receive_apodization=double(rx_aperture_distance<=aperture/2).*(0.5 + 0.5*cos(2*pi*rx_aperture_distance./aperture)); 

%% delay compensation and beamforming
if(strcmp(acq, 'PW'))
    transmit_distance = z_grid(:)*cos(tx)+x_grid(:)*sin(tx);
elseif (strcmp(acq, 'STA'))
    transmit_distance = sqrt((x_grid(:)-tx).^2+z_grid(:).^2); %ang will be the transmit center in SA transmission
end

for nrx=1:N_c
    
    receive_distance=sqrt((probe_geometry(nrx,1)-x_grid(:)).^2+z_grid(:).^2);
    delay=(transmit_distance+receive_distance)./c_map(:);
       
    delay_compensation(:,nrx)=interp1(timeVector,rawData(:,nrx),delay,'spline',0);
end

beamformedData(:)=sum(receive_apodization.*delay_compensation,2);
beamformedData(isnan(beamformedData))=0;
beamformedData=reshape(beamformedData,size(x_grid));