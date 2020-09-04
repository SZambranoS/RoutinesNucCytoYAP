function [Objects] = segment_fromMask_NOWATERSHED(mask, par)

%SegmentCells from a mask

    minSize=par(1);
    maxSize=par(2);
    
    mask = imfill(mask,'holes'); % fill holes
    ThrImg = bwareaopen(mask,minSize);
   
    ThrImg = xor(bwareaopen(ThrImg,minSize),bwareaopen(ThrImg,maxSize));
            
    Objects = bwconncomp(ThrImg);
    
end

