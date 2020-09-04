function [Objects,roughmask,refinedmask] = segment_fromImage_refined_cluster(img, par)

sigma = par(1);    %sigma of the gaussian blurring of the image
minSize = par(2);  %minimum size of the objects for the BLUE image
nclusters = par(3);   %clusters for the first kmeans.
factor_fine = par(4);   %for the fine thresholding.
maxSize=par(5);    %maximum size of the objects for the BLUE image

% Blur roughly the  image
    FilterImg = imfilter(img,fspecial('gaussian', 20, 2),'same');

    
    fmeantiles = @(block_struct) mean2(block_struct.data);
    meantiles = blockproc(FilterImg,[32 32],fmeantiles);
    
    meantiles(:);
    
    
% repeat the clustering 3 times to avoid local minima
[~, cluster_center] = kmeans(meantiles(:),nclusters,'distance','sqEuclidean', ...
                                      'Replicates',3);
                                  
                                 
    
    [valuesclustersorted,indexessorted]=sort(cluster_center);
    
   
    if nclusters>=2
    levelTr=mean(valuesclustersorted(2:nclusters))
    else
    levelTr=mean(valuesclustersorted); 
    end; 
    %levelTr=valuesclustersorted(nclusters)

    
    
    
% Threshold the image 
   
     %ThrImg = (FilterImg>1.5*levelTr);
     ThrImg = (FilterImg>0.8*levelTr);
   
  
% Exclude the small objects in the  thresholded image
    ThrImg = imfill(ThrImg,'holes'); % fill holes

     ThrImg = bwareaopen(ThrImg,minSize);
     
     imagesc(ThrImg);
    
    % figure(800)
     %imagesc(ThrImg);
% Get the objects thresholded 

    ObjectsProv=segment_fromMask(ThrImg, [0.5*minSize,1.5*maxSize]);    %Notice that for the provisory objects we are more loose about the cells considered
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
    if sigma>0
    doubleimg=imfilter(doubleimg,fspecial('gaussian', 20,sigma),'same');
    end
        
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


Objects=segment_fromMask(refinedmask, [minSize,maxSize]);


end

