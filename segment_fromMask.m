function [Objects] = segment_fromMask(mask, par)

%SegmentCells from a mask

    minSize=par(1);
    maxSize=par(2);
    
    mask = imfill(mask,'holes'); % fill holes
    ThrImg = bwareaopen(mask,minSize);
    
    DistImg = -bwdist(~ThrImg);  %lower values in the center of the objects
    dummy = imextendedmin(DistImg,2); %Search the minima and somehow expands them
    DistImg = imimposemin(DistImg,dummy); %Modifies the image so the local minima are only where dummy is nonzero
    WSHD = watershed(DistImg); %Watershedregions are the basin for the minima, all having the same value
    ThrImg(WSHD == 0) = 0; %This is to identify the boundaries between basins for the watershed.

    
    ThrImg = xor(bwareaopen(ThrImg,minSize),bwareaopen(ThrImg,maxSize));
            
    Objects = bwconncomp(ThrImg);
    
end

