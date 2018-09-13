function [lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,lftamp_err,rgtamp_err] = ...
    func_doublegauss_fit(slice1D,idxValid,lftamp, ...
                         rgtamp,lftloc,rgtloc,lftwid,rgtwid,...
                         minlftloc,minrgtloc,maxlftloc,maxrgtloc,graph)

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
    'a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + d1 + e1*x',...
    'start',[lftamp,rgtamp,lftloc,rgtloc,lftwid,rgtwid,bckgrndNoise,0], ...
    'Lower',[lftamp/4,rgtamp/4, ...
    minlftloc,minrgtloc, ...
    minExpectedWidth,minExpectedWidth, -lftamp,-lftamp], ...
    'Upper',[lftamp*4,rgtamp*4,...
    maxlftloc,maxrgtloc,...
      maxExpectedWidth, maxExpectedWidth, lftamp,lftamp] );

%Save the data
lftamp = f.a1; rgtamp = f.a2;
lftloc = f.b1; rgtloc = f.b2;
lftwid = f.c1; rgtwid = f.c2;
ConfInt_95 = confint(f,0.95);

%To Get SE from ConfInt is silly, but so is Matlab
%RMSE = SE = (upper limit – lower limit) / 3.92.
stupidRMSE_left = (   ConfInt_95(2,1) -ConfInt_95(1,1)   )  / 3.92;
stupidRMSE_right = (   ConfInt_95(2,2) -ConfInt_95(1,2)   )  / 3.92;
lftamp_err = stupidRMSE_left;
rgtamp_err = stupidRMSE_right;
lin = f.d1; slope = f.e1;
 
end