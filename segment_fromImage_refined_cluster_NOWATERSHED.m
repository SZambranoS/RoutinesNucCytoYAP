function [Objects,roughmask,refinedmask] = segment_fromImage_refined_cluster_NOWATERSHED(img, par)

sigma = par(1);    %sigma of the gaussian blurring of the image
minSize = par(2);  %minimum size of the objects for the BLUE image
nclusters = par(3);   %clusters for the first kmeans.
factor_fine = par(4);   %for the fine thresholding.
maxSize=par(5);    %maximum size of the objects for the BLUE image

% Blur roughly the  image
    FilterImg = imfilter(img,fspecial('gaussian', 20,  1),'same');

    
    fmeantiles = @(block_struct) mean2(block_struct.data);
    meantiles = blockproc(FilterImg,[32 32],fmeantiles);
    
    meantiles(:);
    
    
% repeat the clustering 3 times to avoid local minima
[~, cluster_center] = kmeans(meantiles(:),3,'distance','sqEuclidean', ...
                                      'Replicates',3);
                                  
                                 
    
    [valuesclustersorted,indexessorted]=sort(cluster_center);
    
    levelTr=mean(valuesclustersorted(1:nclusters-1));
    
    
    
    
% Threshold the image 
   
     ThrImg = (FilterImg>levelTr);
  
% Exclude the small objects in the  thresholded image
    ThrImg = imfill(ThrImg,'holes'); % fill holes
    ThrImg = bwareaopen(ThrImg,minSize);
    
    
% Get the objects thresholded 

    ObjectsProv=segment_fromMask_NOWATERSHED(ThrImg, [0.5*minSize,1.5*maxSize]);    %Notice that for the provisory objects we are more loose about the cells considered
    LabelsMapProv=labelmatrix(ObjectsProv);   
    BoundingBox = regionprops(LabelsMapProv,'BoundingBox');  
    
    roughmask=(LabelsMapProv>0);
    
%    figure(999)
%       imagesc(FilterImg)
%       colorbar
% % %     
%     figure(1000)
%      imagesc(ThrImg);
%   
%      figure(1001)
%      imagesc(LabelsMapProv);
  
    refinedmask=zeros(size(LabelsMapProv));   
    doubleimg=double(img);
    doubleimg=imfilter(doubleimg,fspecial('gaussian', 20,sigma),'same');

%We apply now a threshold applied for the mean in each object. 
    
for n=1:size(BoundingBox)
        Boundingcell=BoundingBox(n).BoundingBox;
        smallmat=doubleimg(ceil(Boundingcell(2)): floor(Boundingcell(2)+Boundingcell(4)), ceil(Boundingcell(1)): floor(Boundingcell(1)+Boundingcell(3)));
        matrixprov=zeros(size(LabelsMapProv));
        matrixprov(ceil(Boundingcell(2)): floor(Boundingcell(2)+Boundingcell(4)), ceil(Boundingcell(1)): floor(Boundingcell(1)+Boundingcell(3)))=1;            
        threshprov=factor_fine*mean(smallmat(:));
        refinedmask=or(refinedmask,((matrixprov.*doubleimg)>threshprov));
             
 end;
       

refinedmask=refinedmask.*(LabelsMapProv>0); 
refinedmask=imfill(refinedmask,'holes');
refinedmask = xor(bwareaopen(refinedmask,minSize),bwareaopen(refinedmask,maxSize));
Objects = bwconncomp(refinedmask);
% 
%    figure(1002)
%      imagesc(labelmatrix(Objects));
% %     
%     pause(3)
%   

end

