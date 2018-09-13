function [loc1,amp1,wid1] = ...
    func_singlegauss_fit(slice1D,idxValid,amp1,loc1,wid1,minloc,maxloc,graph)

global xRgt

bckgrndNoise = 5;
trimStart = 1;
trimEnd = xRgt;

xList = (trimStart:trimEnd)';

minExpectedWidth =  1; %The min half-width we expect for a bleach line
maxExpectedWidth = 6; %The max half-width we expect for a bleach line


%fitting parameters go in alphabetical order so
% a1 a2 b1 b2 c1 c2 d
f = fit(xList(idxValid),slice1D(idxValid), ...
    'a1*exp(-((x-b1)/c1)^2) + d1 + e1*x',...
    'start',[amp1,loc1,wid1,bckgrndNoise,0], ...
    'Lower',[amp1/4,minloc,minExpectedWidth,-amp1,-amp1], ...
    'Upper',[amp1*4,maxloc,maxExpectedWidth,amp1,amp1] );

%Save the data
amp1 = f.a1; 
loc1 = f.b1; 
wid1 = f.c1;
lin = f.d1; slope = f.e1;
 
end