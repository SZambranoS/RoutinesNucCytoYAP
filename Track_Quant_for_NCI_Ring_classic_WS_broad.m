function [LabelsMap, cellROIs_refined,refinedmask,roughmask] = Track_Quant_for_NCI_Ring_classic_WS_broad(StackTrack,StackQuant, parCellTrack, filenameprint, filenamesave, shallweplot,pathoutput)



parCellSeg=parCellTrack(1:5);
TolArea=parCellTrack(6);
RingWidth=parCellTrack(7);
FactorBG=parCellTrack(8);
sigmaforring=parCellTrack(9);


nFrames = length(StackTrack);

for i = 1:nFrames;
    [cellROIs_refined(i),roughmask(i).data,refinedmask(i).data] = segment_fromImage_refined_cluster(StackTrack(i).data, parCellSeg);
    disp(i)
end
clear i


%% Track Cells based on Nearest Neighbor






[LabelsMap, LabelsRing, OUT, nCells,TotalintensityTrack,TotalintensityQuant,TotalintensitySignalQuant,ObjectsPerFrame,AverageBGQUANT] = track_and_quantify_2channels_forNuctoCyto_broad(cellROIs_refined, StackTrack, StackQuant, nFrames, TolArea,RingWidth,FactorBG,sigmaforring);





% for n=1:nFrames
%     figure(101)
%     imagesc(LabelsMap(n).data)
%     figure(102)
%     imagesc(LabelsRing(n).data)
% 
% end;


%%%Extract information of interest in matrices, each column corresponds to
%%%a cell.

matrixareas=-1*ones(nFrames,nCells);
matrixQUANT=-1*ones(nFrames,nCells);
matrixTRACK=-1*ones(nFrames,nCells);
matrixLengthboundaries=-1*ones(nFrames,nCells);


for n=1:nCells
    nframes=OUT{n}.maxFrame;
    matrixareas(1:nframes,n)=OUT{n}.Area;
    matrixQUANT(1:nframes,n)=OUT{n}.TotalIntensityQuant;
    matrixTRACK(1:nframes,n)=OUT{n}.TotalIntensityTrack;
    matrixINTRING(1:nframes,n)=OUT{n}.RingIntQuant;
    matrixAREARING(1:nframes,n)=OUT{n}.RingArea;
    matrixLengthboundaries(1:nframes,n)=OUT{n}.LengthBoundary;
end;





%%%This is to print jpegs of the tracking





cd(pathoutput);


for i=1:nFrames
    
    [B,L] = bwboundaries(LabelsMap(i).data>0);
    
    
    Bnuc{i}=B;
    
    [Bring,Lring] = bwboundaries(LabelsRing(i).data);
    
    Bcyto{i}=Bring;
    
    
    
    
    
    if (strcmp(shallweplot,'QUANT')||(strcmp(shallweplot,'BOTH'))||(strcmp(shallweplot,'BOTHANDRINGS'))||(strcmp(shallweplot,'QUANTANDRINGS')))
        
        
        
        
        
        figure(1)
        %    FilterImg = imfilter(StackQuant(i).data,fspecial('gaussian', 20,  0),'same');
        FilterImg=StackQuant(i).data;
        imagesc(FilterImg);
        colormap('hot');
        %colormap('jet');
        hold on;
        
        
        if strcmp(shallweplot,'BOTH')|| strcmp(shallweplot,'BOTHANDRINGS')
            
            figure(2)
            %       FilterImg = imfilter(StackTrack(i).data,fspecial('gaussian', 20,  0),'same');
            FilterImg=StackTrack(i).data;
            imagesc(FilterImg);
            colormap('hot');
            %colormap('jet');
            hold on;
            
        end;
        
        
        
        
        
        for k = 1:length(B)
            figure(1)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1);
            if strcmp(shallweplot,'BOTH')||strcmp(shallweplot,'BOTHANDRINGS')
                figure(2)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1);
                
            end;
        end
        
        %[Bring,Lring] = bwboundaries(LabelsRing(i).data>0);
        
        
        
        
        if strcmp(shallweplot,'QUANTANDRINGS')||strcmp(shallweplot,'BOTHANDRINGS')
        
        for k = 1:length(Bring)
            figure(1)
            boundary = Bring{k};
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1);
        end
        
        end;
      
        
        
        
        for k=1:nCells
            
            if i<=OUT{k}.maxFrame
                figure(1)
                plot(OUT{k}.Baricenter(i,1), OUT{k}.Baricenter(i,2),'w+');
                text(OUT{k}.Baricenter(i,1), OUT{k}.Baricenter(i,2),num2str(k),'color','w');
                
                
                if strcmp(shallweplot,'BOTH')||strcmp(shallweplot,'BOTHANDRINGS')
                    figure(2)
                    plot(OUT{k}.Baricenter(i,1), OUT{k}.Baricenter(i,2),'w+');
                    text(OUT{k}.Baricenter(i,1), OUT{k}.Baricenter(i,2),num2str(k),'color','w');
                end;
                
            end;
            
        end;
        
        
        
        
        
        filetoprint=strcat('QUANT',filenameprint,'F',num2str(i),'.jpg');
        figure(1)
        title(num2str(i));
        print(filetoprint,'-djpeg');
        
        if strcmp(shallweplot,'BOTH') ||strcmp(shallweplot,'BOTHANDRINGS') 
            filetoprint=strcat('TRACK',filenameprint,'F',num2str(i),'.jpg');
            figure(2)
            title(num2str(i));
            print(filetoprint,'-djpeg');
        end;
        
        hold off;
        %close all;
        
        
    end;
    
end;

Bnuc;

Bcyto;


save(filenamesave,'OUT', 'matrixareas', 'matrixQUANT','matrixINTRING','matrixAREARING', 'matrixLengthboundaries', 'matrixTRACK','TotalintensityTrack','TotalintensityQuant','ObjectsPerFrame','AverageBGQUANT','Bnuc','Bcyto');


end

