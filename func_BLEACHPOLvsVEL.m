function[] = func_BLEACHPOLvsVEL(trim,sensitivity,FrameStart,FrameEnd,pth_sdt)

%time between frames in seconds
global timeStep;

%Size per Pixel
global umPerPixel; 

%Bin size for y-data
global yBin; 
yBin =10;

%Global for how wide the data is
global xRgt;

%Global for how tall the data is
global yBotEnd;

%Global to hold raw images in RAM
global TPH;

Pol_list = [];
Vel_list = [];
tph_name = '\';

%% Set up frames list with sane numbering
frames = (FrameStart:FrameEnd) - (FrameStart -1);

%% Read in TPH data and flatten into 1D strips
[TPH, smoothData, data1D] = func_TPH_read(pth_sdt, tph_name, frames,FrameStart);
xRgt = length(smoothData(:,1,1));
yBotEnd = length(TPH(:,1,1));

%% Determine a guesss set of Peaks
graph = 0;
[SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,avgwid] =...
    func_guesspeaks(data1D,TPH,frames,sensitivity,graph,trim);

%% Fit double gaussians to 1D strips
graph = 0;
[lftloc,rgtloc,SkipList] =  func_doubleline_fit( ...
    data1D,frames,TPH,graph,SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,avgwid);


%% Load Time Stamp information
%This is to be used if you have a pre-defined list of times for each frame
%time_list_path = [pth_sdt '\time_list'];
%load(time_list_path,'time_list');
%If you don't, and just want a rough idea of velocities, you could always
%use something like:
time_list = frames/2;


%% Calculate velocities
graph = 0;
[V_list_L,V_list_R] = func_vel_fit_list(frames, lftloc,rgtloc,SkipList,time_list,graph);
Vel_list = (V_list_R - V_list_L);


%% Calculate Polarity
Pol_list= func_amp_fit_list_first(frames(end),lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,SkipList,0 );

%%Save the data 
save([pth_sdt '\data' num2str(FrameStart)],'Pol_list','Vel_list');