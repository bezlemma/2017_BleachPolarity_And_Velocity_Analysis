%clear all;

function [TPH, smoothData, data1D_HOLD] = func_TPH_read(pth_sdt, tph_name, frames,FrameStart)

tph_dir = [pth_sdt tph_name]; 
global yBin
%Height of data
global yBotEnd

D_tph = dir([tph_dir '\*.tif']);

%Get size info from the first frame
TPH_Temp = double(imread([tph_dir D_tph(1).name]));
yBotEnd = length(TPH_Temp(:,1));
xRgt = length(TPH_Temp(1,:));
FrameStart = FrameStart - 1; %Fix some numbering problems

%% ---------------- Fit gaussians to the given data
for fr = frames(end):-1:frames(1)
    %Read in a frame
    TPH(:,:,fr) = double(imread([tph_dir D_tph(fr+FrameStart).name]));
    %[tph_dir D_tph(fr+FrameStart).name]
    for yTop = 1:(yBotEnd-yBin)
        %Set up the y counter
        yy = yTop;  
        %Crop the frame
        data(:,:,fr)  = imcrop(TPH(:,:,fr),[1 yTop xRgt yBin]);  %[xmin ymin width height]
        %Flatten the data into 1D
        data1D(:,fr) = mean(data(:,:,fr),1);
        %Flip the data
        data1D(:,fr) = abs(max(data1D(:,fr)) - data1D(:,fr) );
        %Store the data
        data1D_HOLD(:,yy,fr) = data1D(:,fr);
        %Create a smooth version of the data if necessary
        smoothData(:,yy,fr) = smooth(data1D(:,fr),4);
    end
end  %End Fr Loop
 

end