%This is a macro to run func_BLEACHPOLvsVEL for many different data sets
%and pool them together

%time between frames in seconds, you'll want to replace this with actual
%timing data from the microscope
global timeStep;
timeStep = 1;

%Size per Pixel, again you'll want to replace this with real data that may
%change from sample to sample
global umPerPixel;
umPerPixel = 2/3;

%pick a root folder
pth_root = 'C:\Users\USERNAME\LOCATION';
paths = dir(pth_root);

%Fitting for all expected images
trim = 1; %How much to trim off the front of the picture. Keep this as 1, it's not really an  implemented feature yet.
set_size = 10; %number of frames per data set, this needn't bee a single number.
sensitivity = 50; %Roughly the count distance from the non-bleach signal to the bleach signal


%Here's the loop to run through the data. You can break this into seperate
%loops to change individual global fitting and timing parameters
for jj = 3:(length(paths)) %ingnore the first two locations, they are always "." and ".."
    
    pth_sdt = [pth_root paths(jj).name '\'];
   % paths(jj).name
    
    %How many images are in this folder?
    numel_data = numel(   dir( [ pth_sdt  '*.tif']  )  );
    if numel_data > 0
    for ii = 1:floor(numel_data/set_size)
       
        FrameStart = 1+ii*set_size-set_size;
        FrameEnd = FrameStart + set_size;
        try
           func_BLEACHPOLvsVEL(trim,sensitivity,FrameStart,FrameEnd,pth_sdt)
              ii %just print out where we are in the data set.
        catch
            warning('BleachPolvsVel Failed, skipping');
       end
    end
    end
    
end