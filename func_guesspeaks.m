function [SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,avgwid]  = ...
    func_guesspeaks(data1D,TPH,frames,sensitivity,graph,trim)

global yBin
%Height of data
global yBotEnd
global xRgt

%% ---------------- Fitting parameters
trimStart = trim;
trimEnd = xRgt-trim;
clear rgtloc;
clear lftloc;


%Intialize/Reset Variables

loc   = NaN*zeros(length(frames),yBotEnd-yBin,2);
lftloc   = NaN*zeros(length(frames),yBotEnd-yBin);
rgtloc   = NaN*zeros(length(frames),yBotEnd-yBin);
wid   = NaN*zeros(length(frames),yBotEnd-yBin,2);
lftwid   = NaN*zeros(length(frames),yBotEnd-yBin);
rgtwid   = NaN*zeros(length(frames),yBotEnd-yBin);
amp  = NaN*zeros(length(frames),yBotEnd-yBin,2);
lftamp   = NaN*zeros(length(frames),yBotEnd-yBin);
rgtamp   = NaN*zeros(length(frames),yBotEnd-yBin);
SkipList =  zeros(yBotEnd-yBin,length(frames));
ylist = (1/2):(yBotEnd-yBin);
yIntList = 1:(yBotEnd-yBin) ;


%% Make an intial set of guesses with little trim
for fr = (frames(end)):-1:frames(1)
    for yy = yIntList
        
        FrameData = data1D(trimStart:trimEnd,yy,fr) ;
        
        %Find the peaks
        [amp_pks, locs_pks, wid_pks, prom_pks] = ...
            findpeaks(FrameData,'MinPeakProminence',max(FrameData)/(sensitivity));
        
        [prom_pks, a_order] = sort(prom_pks,'descend');
        locs_pks =  locs_pks(a_order,:);
        if length(amp_pks) > 1
            loc(fr,yy,2) = locs_pks(2)+trimStart-1;
            loc(fr,yy,1) = locs_pks(1)+trimStart-1;
        elseif length(amp_pks) == 1
            loc(fr,yy,1) = locs_pks(1)+trimStart-1;
        end
        
        
    end
end



%% Sort out our data into a left column and a right column

LeftMaster  = (loc(:,:,1) > loc(:,:,2)) .* (loc(:,:,1) > 0);
RightMaster = (loc(:,:,1) < loc(:,:,2)) .* (loc(:,:,2) > 0);

lftloc =    LeftMaster.*loc(:,:,2) + RightMaster.*loc(:,:,1);
rgtloc =   LeftMaster.*loc(:,:,1) + RightMaster.*loc(:,:,2);

%% Possible Graphing Step
if graph == 1
    
    for fr = frames(end):-1:frames(1)
        figure('Name',['Frame: ' num2str(fr)],'NumberTitle','off');
        imagesc(TPH(:,:,fr))
        colormap('gray')
        hold on;
        plot(rgtloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0),'r.','MarkerSize',15)
        plot(lftloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0),'b.','MarkerSize',15)
        
        title(['Raw Guess ' num2str(fr) ' Fit'],'FontSize',18,'interpreter','latex');
        ylabel('Bleach Axis [um]','FontSize',18,'interpreter','latex');
        xlabel('Ordered Axis[um]','FontSize',18,'interpreter','latex');
        set(gca,'fontsize',18)
        
    end
end


% %Compare y under and over
for fr = frames
    [y,i,rgtloc(fr,:),xsigma] = hampel(rgtloc(fr,:),yBin*10);
    [y,i,lftloc(fr,:),xsigma] = hampel(lftloc(fr,:),yBin*10);
    [y,i,rgtloc(fr,:),xsigma] = hampel(rgtloc(fr,:),yBin*4);
    [y,i,lftloc(fr,:),xsigma] = hampel(lftloc(fr,:),yBin*4);
end

%Fit slopes to the lines
for fr = frames
    %BleachSlopeR = fitlm( ylist'   ,  rgtloc(fr,:) );
    % BleachSlopeL = fitlm( ylist'   ,  lftloc(fr,:) );
    
    % predictR(fr,:) = BleachSlopeR.Coefficients.Estimate(1) + BleachSlopeR.Coefficients.Estimate(2)*yIntList;
    
    BleachSlopeL = polyfit( ylist(~isnan(lftloc(fr,:))),lftloc(fr,~isnan(lftloc(fr,:))),3);
    BleachSlopeR = polyfit( ylist(~isnan(rgtloc(fr,:))),rgtloc(fr,~isnan(rgtloc(fr,:))),3);
    
    predictR(fr,:) = BleachSlopeR(4) + BleachSlopeR(3)*yIntList + BleachSlopeR(2)*yIntList.^2 + BleachSlopeR(1)*yIntList.^3;
    predictL(fr,:) = BleachSlopeL(4) + BleachSlopeL(3)*yIntList + BleachSlopeL(2)*yIntList.^2 + BleachSlopeL(1)*yIntList.^3;
    
end


%% Make a second set of guesses with trim calculated from the first set

for fr = (frames(end)):-1:frames(1)
    for yy = yIntList
        SkipList(yy,fr) = 1;
        lftloc(fr,yy) = NaN;
        rgtloc(fr,yy) = NaN;
        
        buffer = 3;
        
        trimStart = floor(predictL(fr,yy) - buffer);
        idxValidL = trimStart:ceil(predictL(fr,yy) + buffer);
        FrameData = data1D(idxValidL,yy,fr) ;
        
        %Find the peaks
        [amp_pks, locs_pks, wid_pks, prom_pks] = ...
            findpeaks(FrameData,'MinPeakProminence',max(FrameData)/sensitivity);
        if length(amp_pks) > 0
            [prom_pks, a_order] = sort(prom_pks,'descend');
            locs_pks =  locs_pks(a_order,:);
            wid_pks =  wid_pks(a_order,:);
            amp_pks =  amp_pks(a_order,:);
            
            lftwid(fr,yy) = wid_pks(1);
            lftamp(fr,yy) = amp_pks(1);
            lftloc(fr,yy) = locs_pks(1)+trimStart-1; %Check about removing the 1 to see if this is wrong
            SkipList(yy,fr) = 0;
        end
              
        trimStart = floor(predictR(fr,yy) - buffer);
        idxValidR = trimStart:ceil(predictR(fr,yy) + buffer);
        FrameData = data1D(idxValidR,yy,fr) ;
        
        %Find the peaks
        [amp_pks, locs_pks, wid_pks, prom_pks] = ...
            findpeaks(FrameData,'MinPeakProminence',max(FrameData)/sensitivity);
        if length(amp_pks) > 0
        [prom_pks, a_order] = sort(prom_pks,'descend');
        locs_pks =  locs_pks(a_order,:);
        wid_pks =  wid_pks(a_order,:);
        amp_pks =  amp_pks(a_order,:);
        
        rgtwid(fr,yy) =  wid_pks(1);
        rgtamp(fr,yy) =  amp_pks(1);
        rgtloc(fr,yy) = locs_pks(1)+trimStart-1; %Check about removing the 1 to see if this is wrong
        SkipList(yy,fr) = 0;
        end
        
    end
end

avgwid = nanmean(nanmean(lftwid+rgtwid))/2;


%% Check to see if velocities are sane

%Fit a velocity to the motion
for yy = yIntList
    LinFitsR = fitlm( frames , rgtloc(frames,yy) );
    vFitsR(yy) = LinFitsR.Coefficients.Estimate(2);
    LinFitsL = fitlm( frames , lftloc(frames,yy) );
    vFitsL(yy) = LinFitsL.Coefficients.Estimate(2);
end

%% Grab the insane velocities and reset the points
for fr = frames
    
    BleachSlopeL = polyfit( ylist(~isnan(lftloc(fr,:))),lftloc(fr,~isnan(lftloc(fr,:))),3);
    BleachSlopeR = polyfit( ylist(~isnan(rgtloc(fr,:))),rgtloc(fr,~isnan(rgtloc(fr,:))),3);
    
    predictR(fr,:) = BleachSlopeR(4) + BleachSlopeR(3)*yIntList + BleachSlopeR(2)*yIntList.^2 + BleachSlopeR(1)*yIntList.^3;
    predictL(fr,:) = BleachSlopeL(4) + BleachSlopeL(3)*yIntList + BleachSlopeL(2)*yIntList.^2 + BleachSlopeL(1)*yIntList.^3;

    for yy = yIntList( abs(vFitsL) > 1.5*mean(abs(vFitsL)) )
        lftloc(fr,yy) = predictL(fr,yy);
        lftwid(fr,yy) = avgwid;
        lftamp(fr,yy) = data1D(ceil(lftloc(fr,yy)),yy,fr);
    end
    for yy = yIntList( abs(vFitsR) > 1.5*mean(abs(vFitsR)) )
        rgtloc(fr,yy) = predictR(fr,yy);
        rgtwid(fr,yy) = avgwid;
        rgtamp(fr,yy) = data1D(ceil(rgtloc(fr,yy)),yy,fr);
    end
    % Grab outliers from prediction and reset the points
%     predictR = BleachSlopeR.Coefficients.Estimate(1) + BleachSlopeR.Coefficients.Estimate(2)*yIntList;
%     predictL = BleachSlopeL.Coefficients.Estimate(1)   + BleachSlopeL.Coefficients.Estimate(2)*yIntList;
%
%          kr = yIntList( abs(rgtloc(fr,:) - predictR) > 3  );
%          rgtloc(fr,kr) = BleachSlopeR.Coefficients.Estimate(1) + BleachSlopeR.Coefficients.Estimate(2)*kr;
%          kl = yIntList( abs(lftloc(fr,:) - predictL) > 3  );
%          lftloc(fr,kl) = BleachSlopeL.Coefficients.Estimate(1) + BleachSlopeL.Coefficients.Estimate(2)*kl;
end

% %% Hampel Filter the Guesses
% %Hampel can show the sigma of the data, we could delete or
% %bend data depending on the sigma
%
% %Compare y under and over
% for fr = frames
%     [y,i,rgtloc(fr,:),xsigma] = hampel(rgtloc(fr,:),yBin);
%     [y,i,lftloc(fr,:),xsigma] = hampel(lftloc(fr,:),yBin);
% end

%Compare frames before and after, BE CAREFUL, THIS MAKES NANS DISAPPEAR
% for  yy = yIntList
%     [y,i,rgtloc(frames,yy),xsigma] = hampel(rgtloc(frames,yy),1);
%     [y,i,lftloc(frames,yy),xsigma] = hampel(lftloc(frames,yy),1);
% end



%% Possible Graphing Step
if graph == 1
    
    for fr = frames(end):-1:frames(1)
        figure('Name',['Frame: ' num2str(fr)],'NumberTitle','off');
        imagesc(TPH(:,:,fr))
        colormap('gray')
        hold on;
        plot(rgtloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0),'r.','MarkerSize',15)
        plot(lftloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0),'b.','MarkerSize',15)
        plot(predictR(fr,:),ylist,'g-','MarkerSize',15)
        plot(predictL(fr,:),ylist,'g-','MarkerSize',15)
        
        title(['Frame ' num2str(fr) ' Fit'],'FontSize',18,'interpreter','latex');
        ylabel('Bleach Axis [um]','FontSize',18,'interpreter','latex');
        xlabel('Ordered Axis[um]','FontSize',18,'interpreter','latex');
        set(gca,'fontsize',18)
        
    end
end


disp([ 'Using ' num2str(yBotEnd-sum(SkipList)) ' out of '  num2str(yBotEnd) ' rows'])



end