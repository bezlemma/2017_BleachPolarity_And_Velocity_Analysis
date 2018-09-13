function [lftloc,rgtloc,lftamp,rgtamp,lftamp_err,rgtamp_err,SkipList] = ...
    func_doubleline_fit(data1D,frames,TPH,graph,SkipList,lftloc,rgtloc,...
    lftamp,rgtamp,lftwid,rgtwid,avgwid)

global yBin
global yBotEnd
global xRgt

%% ---------------- Fitting parameters
trimStart = 1; %THERE IS AN ERROR HERE BETWEEN HOW idxValid works and how xList works! Don't try to trim the data unless you fix it or your indexing will be off!
trimEnd = xRgt;
ylist = (1/2):(yBotEnd-yBin);
yIntList =  1:(yBotEnd-yBin);
NARROW = 3.5;
graphgraph = 0;

%% ---------------- Fit gaussians to the last frame
fr = frames(end);
disp(['Now analayzing frame: ' num2str(fr)] )
    
    for yy = yIntList
        
        if SkipList(yy,fr) == 1
            continue
        end
        
        %Narrow the window to NARROW pixels centered around each
        %gaussian (a minimum of 8 points is needed to fit, we have NARROW*2*2)
        removeList2 = NaN*zeros(length(1:trimEnd),1);
        for ii = (trimStart):(trimEnd)
            if  (ii < (lftloc(fr,yy)+ NARROW)) &&  (ii > (lftloc(fr,yy) - NARROW))
                removeList2(ii) = 1;
            end
            if  (ii < (rgtloc(fr,yy)+ NARROW)) &&  (ii > (rgtloc(fr,yy) - NARROW))
                removeList2(ii) = 1;
            end
        end
        
        idxValid2 = ~isnan(removeList2);
        
        if ~isnan(lftloc(fr,yy)) && ~isnan(rgtloc(fr,yy))
                [lftloc(fr,yy), rgtloc(fr,yy),lftamp(fr,yy), rgtamp(fr,yy),~,~,~,~] = ...
            func_doublegauss_fit(data1D(trimStart:trimEnd,yy,fr),...
            idxValid2,  lftamp(fr,yy),rgtamp(fr,yy),...
            lftloc(fr,yy),rgtloc(fr,yy),...
            lftwid(fr,yy),rgtwid(fr,yy),...  %Guesses at Widths 
            (lftloc(fr,yy)-NARROW/2),(rgtloc(fr,yy)-NARROW/2),... %Min Locations
            (lftloc(fr,yy)+NARROW/2),(rgtloc(fr,yy)+NARROW/2),graphgraph); %Max Locations
        elseif ~isnan(lftloc(fr,yy))
             [lftloc(fr,yy),lftamp(fr,yy),~] = ...
            func_singlegauss_fit(data1D(trimStart:trimEnd,yy,fr),...
            idxValid2,  lftamp(fr,yy),lftloc(fr,yy),...
            lftwid(fr,yy),(lftloc(fr,yy)-NARROW/2),... 
            (lftloc(fr,yy)+NARROW/2),graphgraph); %Max Locations
        elseif ~isnan(rgtloc(fr,yy))
              [rgtloc(fr,yy),rgtamp(fr,yy),~] = ...
            func_singlegauss_fit(data1D(trimStart:trimEnd,yy,fr),...
            idxValid2,  rgtamp(fr,yy),rgtloc(fr,yy),...
            rgtwid(fr,yy),(rgtloc(fr,yy)-NARROW/2),... 
            (rgtloc(fr,yy)+NARROW/2),graphgraph); %Max Locations
        else
            fuck
        end
        
        
    end   %End Y Loop


%% ---------------- Fit gaussians to the following frames
for fr = (frames(end)-1):-1:frames(1)
    disp(['Now analayzing frame: ' num2str(fr)] )
    
    for yy = yIntList
        
        if SkipList(yy,fr) == 1 
            continue
        end
        
        %Narrow the window to NARROW pixels centered around each
        %gaussian (a minimum of 8 points is needed to fit, we have NARROW*2*2)
        removeList2 = NaN*zeros(length(1:trimEnd),1);
        for ii = (trimStart+3):(trimEnd-3)
            if  (ii < (lftloc(fr,yy)+ NARROW)) &&  (ii > (lftloc(fr,yy) - NARROW))
                removeList2(ii) = 1;
            end
            if  (ii < (rgtloc(fr,yy)+ NARROW)) &&  (ii > (rgtloc(fr,yy) - NARROW))
                removeList2(ii) = 1;
            end
        end
        
        idxValid2 = ~isnan(removeList2);

    %later frames have spread further out
    %which means you can always bound a rgt gaussian to be further away after each frame 
    %if going backwards, each rgt gaussian needs to be closer than the
    %previous one.
    
        %if the two frames are too far apart, screw it
         if lftloc(fr+1,yy) > (lftloc(fr,yy)+NARROW-1) 
            SkipList(yy,fr) = 1;
            lftloc(fr,yy) = NaN;
            rgtloc(fr,yy) = NaN;
            continue
         end
        
          if (rgtloc(fr,yy)-NARROW+1) > rgtloc(fr+1,yy)
            SkipList(yy,fr) = 1;
            lftloc(fr,yy) = NaN;
            rgtloc(fr,yy) = NaN;
            continue
          end
            
            %By default say the errors are 0, they will get real errors if
            %a real fit occurs
            lftamp_err(fr,yy) = 0;
            rgtamp_err(fr,yy) = 0;
            
        if ~isnan(lftloc(fr,yy)) && ~isnan(rgtloc(fr,yy))
            [lftloc(fr,yy), rgtloc(fr,yy), lftamp(fr,yy), rgtamp(fr,yy),~,~,lftamp_err(fr,yy),rgtamp_err(fr,yy)] = ...
            func_doublegauss_fit(data1D(trimStart:trimEnd,yy,fr),...
            idxValid2,  lftamp(fr,yy),rgtamp(fr,yy),...
            lftloc(fr,yy),rgtloc(fr,yy),...
            lftwid(fr,yy),rgtwid(fr,yy),...  %Guesses at Widths 
            lftloc(fr+1,yy)-0.1,(rgtloc(fr,yy)-NARROW),... %Min Locations
            (lftloc(fr,yy)+NARROW),rgtloc(fr+1,yy)+0.1,graphgraph); %Max Locations
        
        
        elseif ~isnan(lftloc(fr,yy))
             [lftloc(fr,yy),lftamp(fr,yy),~] = ...
            func_singlegauss_fit(data1D(trimStart:trimEnd,yy,fr),...
            idxValid2,  lftamp(fr,yy),lftloc(fr,yy),...
            lftwid(fr,yy),(lftloc(fr,yy)-0.1),... 
            (lftloc(fr,yy)+NARROW/2),graphgraph); %Max Locations
        elseif ~isnan(rgtloc(fr,yy))
              [rgtloc(fr,yy),rgtamp(fr,yy),~] = ...
            func_singlegauss_fit(data1D(trimStart:trimEnd,yy,fr),...
            idxValid2,  rgtamp(fr,yy),rgtloc(fr,yy),...
            rgtwid(fr,yy),(rgtloc(fr,yy)-NARROW/2),... 
            (rgtloc(fr,yy)+0.1),graphgraph); %Max Locations
        else
            fuck
        end
            
            

        
    end   %End Y Loop
end  %End Fr Loop




end