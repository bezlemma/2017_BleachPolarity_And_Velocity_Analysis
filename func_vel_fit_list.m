function [vFitsL, vFitsR, vErrL, vErrR] = func_vel_fit_list( frames, lftloc, rgtloc,SkipList,time_list,graph)
global yBin
global yBotEnd
global TPH
global umPerPixel
%This function simply creates a list of velocity fits based off the polyfit
%function

%Uses R-squared > 0.50 as a threshold for keeping a linear fit.
%See https://www.mathworks.com/help/curvefit/evaluating-goodness-of-fit.html
%for details
yIntList = 1:(yBotEnd-yBin) ;

%% Hampel Filter the Guesses
%Hampel can show the sigma of the data, we could delete or
%bend data depending on the sigma


%Compare frames before and after
  for  yy = 1:(yBotEnd-yBin)
      [y,i,rgtloc(frames,yy),xsigma] = hampel(rgtloc(frames,yy),1);
      [y,i,lftloc(frames,yy),xsigma] = hampel(lftloc(frames,yy),1);
  end
  
minFit = 80/100; %(value from 0 to 1 that tells how close a linear fit neads to be to be accepted)

%Fit a velocity to the motion
for yy = 1:length(rgtloc(1,:))
    %Remove outliers
    k = yIntList(  isoutlier(rgtloc(frames,yy))  );
    rgtloc(k,yy) = NaN;
    k = yIntList( isoutlier(lftloc(frames,yy)) );
    lftloc(k,yy) = NaN;
    %Make Fit
    LinFitsR = fitlm( time_list(frames)'   ,  rgtloc(frames,yy)*umPerPixel );
    LinFitsL = fitlm(   time_list(frames)'  ,  lftloc(frames,yy)*umPerPixel );  
    
    %Only save the data if the fit was good
    if  LinFitsR.Rsquared.Adjusted > minFit %&& LinFitsR.Coefficients.Estimate(2) > 0
        vFitsR(yy) = LinFitsR.Coefficients.Estimate(2);
         vErrR(yy) = LinFitsR.RMSE;%LinFitsR.Coefficients.SE(2);
    end
    if LinFitsL.Rsquared.Adjusted > minFit %&& LinFitsL.Coefficients.Estimate(2) < 0
        vFitsL(yy) = LinFitsL.Coefficients.Estimate(2);
        vErrL(yy) = LinFitsL.RMSE;%LinFitsL.Coefficients.SE(2);
    end
    
end


%Now go through again, and save the data if one fit is good, using the mean
%of the other. Otherwise, throw them both out.
for yy = 1:length(rgtloc(1,:))
    LinFitsR = fitlm( time_list(frames)'   ,  rgtloc(frames,yy)*umPerPixel );
    LinFitsL = fitlm(   time_list(frames)'  ,  lftloc(frames,yy)*umPerPixel );
    
    if LinFitsR.Rsquared.Adjusted > minFit && LinFitsL.Rsquared.Adjusted > minFit
            %do nothing
    elseif LinFitsR.Rsquared.Adjusted > minFit && LinFitsL.Rsquared.Adjusted < minFit
        vFitsR(yy) = LinFitsR.Coefficients.Estimate(2);
        vFitsL(yy) = nanmean(vFitsL);
    elseif LinFitsR.Rsquared.Adjusted < minFit && LinFitsL.Rsquared.Adjusted > minFit
        vFitsL(yy) = LinFitsL.Coefficients.Estimate(2);
        vFitsR(yy) = nanmean(vFitsR);
    else
        vFitsL(yy) = NaN;
        vFitsR(yy) = NaN;
    end
    
end

end

