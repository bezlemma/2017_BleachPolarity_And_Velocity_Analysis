function [Pol_list] = func_amp_fit_list_first( fr, lftloc,rgtloc,lftamp, rgtamp,lftamp_err, rgtamp_err,lftwid,rgtwid,SkipList,graph)
global yBin
global yBotEnd
%This function simply creates a list of ratios of amplitudes of the one
%image (fr)

for yy = 1:length(rgtamp(1,:))
    
    maxxx =  max(lftamp(fr,yy),rgtamp(fr,yy)) + rgtamp(fr,yy) + lftamp(fr,yy) -  rgtamp(fr,yy) - lftamp(fr,yy); 
    minnn =  min(lftamp(fr,yy),rgtamp(fr,yy)) + rgtamp(fr,yy) + lftamp(fr,yy) -  rgtamp(fr,yy) - lftamp(fr,yy); 
    Pol_list(yy) = (maxxx - minnn) / (maxxx + minnn);
end

end

