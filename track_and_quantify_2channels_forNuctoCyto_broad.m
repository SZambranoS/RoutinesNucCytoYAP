function [LabelsMap, LabelsRing, OUT, nCells,TotalintensityTrack,TotalintensityQuant, TotalintensityQuantSignal, ObjectsPerFrame,AverageBGQuant] = track_and_quantify_2channels_forNuctoCyto_broad(cellROIs, StackTrack, StackQuant, nFrames, TolArea,WidthRing,FactorBG,sigmafilter)

% Uses nearest neighbours to track cells in one movie. Cells
% need to be previously identified as "closesd components" of a binary
% image. This routine generates a labels map with the tracked ROIs, the
% list of the baricenter of the cells and of their bounding boxes.

%(SZ, Dec 2015) we keep track of the areas, which might be useful to
%discard weird cells authomatically.

%(SZ, April 2017) we add the possibility of adding a second channel to quantify.

%(SZ, April 2017) it also works to quantify single frames.

%(SZ, June 2017). It computes the background BG of the image, the  average intensity in a ring around the object
%detected of size WidthRing but for pixels above FactorBG*BG 

%(SZ, August 2017) A number of corrections were applied to compute properly
%the ring oc cytosol around each nucleus.

%(SZ, June 2018): It now fills the holes of the cytoplasms. 

LabelsMap(1).data =labelmatrix(cellROIs(1)); %ROIs are specified as connected components

[M,N]=size(LabelsMap(1).data);

BaricenterPre = regionprops(LabelsMap(1).data,'Centroid');
Baricenter{1} = cell2mat((struct2cell(BaricenterPre))');
nCells = length(Baricenter{1}(:,1));




if nFrames==1

    BoundingBoxPre = regionprops(LabelsMap(1).data,'BoundingBox');
    BoundingBox{1}=cell2mat((struct2cell(BoundingBoxPre))'); 
    
else
    
    
    
    for j = 1:nFrames - 1;
        % Compute baricenter of cells at frame j + 1
        BaricenterPost = regionprops(cellROIs(j+1),'Centroid');
        BaricenterPost = cell2mat((struct2cell(BaricenterPost))');
        AreaPost=regionprops(cellROIs(j+1),'Area');
        AreaPost=cell2mat((struct2cell(AreaPost))');
        
        % Compute baricenter of cells at frame j
        BaricenterPre = regionprops(LabelsMap(j).data,'Centroid');
        Baricenter{j} = cell2mat((struct2cell(BaricenterPre))');
        AreaPre0=regionprops(LabelsMap(j).data,'Area');
        AreaPre0=cell2mat((struct2cell(AreaPre0))');
        
        
        %%%If in one frame, by some disaster, things are not working properly,
            %%%you copy last image's cells.
            
        
        if isempty(BaricenterPost)
            disp('WARNING, THERE WAS A WEIRD FRAME')
            LabelsMap(j+1).data = LabelsMap(j).data;
            Baricenter{j+1}=Baricenter{j};
            BoundingBox{j}=cell2mat((struct2cell(BoundingBoxPre))');
       
        elseif isempty(Baricenter{j})
            
            Baricenter{j}=NaN*ones(nCells,2); 
        
            for k=j:nFrames-1
                
                Baricenter{k+1}=NaN*ones(nCells,2);
                LabelsMap(k+1).data=zeros(M,N);
                
                
                
            end;
        
            
            break; 
            
        else
            
            
      
            if length(Baricenter{j}(:,1)) < nCells
                Baricenter{j}(end+1:nCells,:) =NaN;
            end
            BaricenterPre0=Baricenter{j};
            
            
            
            
            %Vector with areas of the previous frame.
            AreaPre=zeros(nCells,1);
            AreaPre(find(~isnan(Baricenter{j}(:,1))))=AreaPre0(find(AreaPre0));
            
            
            
            
            
            
            % Compute Bounding box of cells at frame j
            BoundingBoxPre = regionprops(LabelsMap(j).data,'BoundingBox');
            BoundingBox{j} = cell2mat((struct2cell(BoundingBoxPre))');
            
            
            
            % Nearest Neighbour search
            [indexPre,Distance] = knnsearch(Baricenter{j},BaricenterPost);
            
            % Initiate labels map
            LabelsMapTemp = zeros(size(LabelsMap(j).data));
            
            % Link cells from frame j to frame j + 1;
            
            
            for i = 1:max(indexPre);
                
                multInd = find(indexPre ==i);
                
                if length(multInd)> 1       % if the cell i at frame j is linked
                    indexPre(multInd) = 0;  % to more than one cell at frame j + 1
                    [~,idMin] = min(Distance(multInd)); % only keep the link
                    indexPre(multInd(idMin)) = i;           % that minimize the
                    % the distance
                end
                
                
                cellpre=i;
                cellpost=find(indexPre==i);
                
                %We check if the area changed a lot.
                if ~isempty(cellpost)
                    areapre=AreaPre(i);
                    areapost=AreaPost(cellpost);
                    
                    areamean=mean([areapre,areapost]);
                    
                    if (abs(areapre-areapost)/areamean)>TolArea
                        indexPre(cellpost)=0;
                        disp('Change area!!!!')%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        cellpre
                    end;
                    
                end;
                
                %We also check if the baricentres are too far appart
                
                if ~isempty(cellpost)
                    BariPre=BaricenterPre0(i,:);
                    BariPost=BaricenterPost(cellpost,:);
                    
                    if (norm(BariPre-BariPost))>50
                        indexPre(cellpost)=0;
                       % disp('Change position!!!!')
                    end;
                    
                end;
                
                
                
            end
            
            % Generate Label map for frame j + 1.
            for i = 1:cellROIs(j+1).NumObjects
                LabelsMapTemp(cellROIs(j+1).PixelIdxList{i}) = indexPre(i);
            end
            LabelsMap(j+1).data = LabelsMapTemp;
            
        end;
        
        
        
        
    end
    
    % Calculate Baricenter and Bounding Box for last Frame
    
    BaricenterPre = regionprops(LabelsMap(j+1).data,'Centroid');
   
    
    
    
    if isempty(BaricenterPre)
       
    
        Baricenter{j+1}=NaN*ones(size(Baricenter{j}));
        
    else
 
    Baricenter{j+1} = cell2mat((struct2cell(BaricenterPre))');
    
    if length(Baricenter{j+1}(:,1)) < nCells
                Baricenter{j+1}(end+1:nCells,:) =NaN;
    end
    
    
    BoundingBoxPre = regionprops(LabelsMap(j+1).data,'BoundingBox');
    BoundingBox{j+1} = cell2mat((struct2cell(BoundingBoxPre))');
    
    end;
    
    
    % Rearrange Data cell-wise

    
   
 
end;


fmeantiles = @(block_struct) mean2(block_struct.data);
fstdtiles = @(block_struct) std2(block_struct.data);

  
    




for iFrame = 1:nFrames
    
    TotalintensityTrack(iFrame)=sum(sum(double(StackTrack(iFrame).data)));
    TotalintensityQuant(iFrame)=sum(sum(double(StackQuant(iFrame).data)));
    ObjectsPerFrame(iFrame)=length(unique(LabelsMap(1).data)-1);
    meantiles = blockproc(StackQuant(iFrame).data,[32 32],fmeantiles);
    stdtiles= blockproc(StackQuant(iFrame).data,[32 32],fstdtiles);
    
    [valueminmean,idminmean]=min(meantiles(:));
    
    vstd=stdtiles(:);
    
    stdbg=vstd(idminmean);
    
    threshold(iFrame)=valueminmean+0.5*(stdbg);
    
    if sigmafilter>0
    FilterImgQuant(iFrame).data = imfilter(StackQuant(iFrame).data,fspecial('gaussian', 20,  sigmafilter),'same');
    else
    FilterImgQuant(iFrame).data = StackQuant(iFrame).data;
    end;
    
    QUANTsignalimage=(FilterImgQuant(iFrame).data)>threshold(iFrame);
    
    BGSignalimage=valueminmean*QUANTsignalimage; 
    
    
    AverageBGQuant(iFrame)=valueminmean; 

%     figure(901)
%     imagesc(FilterImgQuant); 
%     colorbar
%     figure(900)
%     imagesc(QUANTsignalimage);
%     colorbar
%     figure(899)
%     imagesc(BGSignalimage);
%     colorbar
%     
    TotalintensityQuantSignal(iFrame)=sum(sum(QUANTsignalimage.*(double(StackQuant(iFrame).data))))-sum(sum(BGSignalimage)); 
    
%     pause(0.5); 


    LabelsRing(iFrame).data=zeros(M,N);
    
end







% Cycle on cells
for iCell = 1:nCells        %cycle on cell
    for iFrame = 1:nFrames; %cycle on frame
        
        
        
       
        %iFrame
        %size(Baricenter{iFrame})
        if ~isnan(Baricenter{iFrame}(iCell,1))      % if the cell exist at iFrame
            
            OUT{iCell}.BoundingBox(iFrame,:) = BoundingBox{iFrame}(iCell,:);   % Reararange Bounding box
            OUT{iCell}.Baricenter(iFrame,:) = Baricenter{iFrame}(iCell,:);     % Rearrange Baricenter
            
            MaskiCell=(LabelsMap(iFrame).data == iCell);
            
            OUT{iCell}.TotalIntensityTrack(iFrame,1) =...                           % Compute total Intensity of Stack for Tracking
                sum(sum(double(StackTrack(iFrame).data).*(MaskiCell)));
            
            OUT{iCell}.TotalIntensityQuant(iFrame,1) =...                           % Compute total Intensity of Stack for Quantification
                sum(sum(double(StackQuant(iFrame).data).*(MaskiCell)));
         
            
          se = strel('square',WidthRing);
          dilatedMask = imdilate(MaskiCell,se);
            
         
           RingCell=dilatedMask-MaskiCell; 
           
          
           
           
           RingCell=double(RingCell);
           
           MatrixRing=(double(FilterImgQuant(iFrame).data)).*(RingCell); 
           
           MaskRing=(MatrixRing>(FactorBG*AverageBGQuant(iFrame)));
           
                      
           matrixRingvalues=MatrixRing.*MaskRing;
          
           MaskRing(MaskiCell)=1;
           
           MaskRing = imfill(MaskRing,'holes'); % fill holes
           
           %Get the bigger one
           
           ObjectsSingleRing = bwconncomp(MaskRing);
           
           SingleRings=labelmatrix(ObjectsSingleRing);
                      
           AreaObjectsSingleRing=regionprops(ObjectsSingleRing,'Area');
           
           AreaFinal=cell2mat((struct2cell(AreaObjectsSingleRing))');
           
           [valuearea,indexarea]=max(AreaFinal);
           
           MaskRing=(SingleRings==indexarea);
           
           MaskRing(MaskiCell)=0; 
          
            
            
            OUT{iCell}.Area(iFrame,1) =length(find(MaskiCell));
            
            
             [Boundary,L] = bwboundaries(MaskiCell);
            
            OUT{iCell}.LengthBoundary(iFrame,1)=length(Boundary{1});
            
        else                                   % If the cell doesn't exist
            OUT{iCell}.maxFrame = iFrame - 1;   % Set last frame for iCell
            break;                              % Exit from the for loop on iFrame
        end
    OUT{iCell}.maxFrame = iFrame;
%         
%       size(LabelsRing(iFrame).data) 
%       size(MaskRing)
    LabelsRing(iFrame).data=LabelsRing(iFrame).data + iCell*MaskRing; 
    
    %LabelsRing(iFrame).data(find(MaskRing))=0;
    LabelsRing(iFrame).data((MaskRing>0))=0;
    LabelsRing(iFrame).data=LabelsRing(iFrame).data + iCell*MaskRing; 
        
    end
    % identify the largest bounding box and store its length and its width
    OUT{iCell}.maxBB = [max(OUT{iCell}.BoundingBox(:,3)), max(OUT{iCell}.BoundingBox(:,4))];
    
    
end



for iFrame = 1:nFrames; %cycle on frame
      %LabelsRing(iFrame).data(find(LabelsMap(iFrame).data))=0; 
      LabelsRing(iFrame).data((LabelsMap(iFrame).data>0))=0; 
end;

for iCell = 1:nCells        %cycle on cell
    for iFrame = 1:nFrames; %cycle on frame

          if ~isnan(Baricenter{iFrame}(iCell,1))      % if the cell exist at iFrame

              
           MaskRing=(LabelsRing(iFrame).data==iCell);
           
%            imagesc(MaskRing+2*(LabelsMap(iFrame).data == iCell));
%            title(num2str(iCell)); 
%            pause(0.5);
           
           MaskRing=double(MaskRing);
           
           MatrixRing=MaskRing.*(double(StackQuant(iFrame).data)); 

           vRing=MatrixRing(:);
           
           %vRing=vRing(find(vRing));
           vRing=vRing(vRing>0);
           
           LRing=length(vRing>0);
             
           if LRing>0           
          OUT{iCell}.RingIntQuant(iFrame,1) = mean(vRing);      %%%%%%%%%%%
           else
          OUT{iCell}.RingIntQuant(iFrame,1) = 0;     
           end;
           
           
         OUT{iCell}.RingArea(iFrame,1) = LRing;    
            
            
          end;
          
    end;
    
end; 
























