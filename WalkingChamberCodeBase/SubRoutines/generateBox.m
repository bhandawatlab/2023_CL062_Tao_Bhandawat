function [] = generateBox(params,extensions,directory,folders)

d_side = directory{1};
d_bot = directory{2};
id_side = extensions{1};
id_bot = extensions{2};

boxFold = folders.boxFolder;
videoFolder = folders.vidFolder;

DLT_params = params.DLT_params;

assert(numel(d_side)==numel(d_bot));
nVid = numel(d_side);

%for currGen = 1:numel(allGen)%3
    close all
    
    % for selecting the bounding box
    vidFile = cell(nVid,1);
    for i = 1:nVid
        fNameSide = d_side(i).name;
        C = strsplit(d_side(i).name,'_');
        vidFile{i} = strjoin(C(1:2),'_');
        if exist([boxFold strjoin(C(1:3),'_') '_BoundingBox.mat'],'file')~=2
            
            v = VideoReader([videoFolder strjoin(C(1:3),'_') '_side_init2.avi']);
            v2 = VideoReader([videoFolder strjoin(C(1:3),'_') '_bot_init2.avi']);
%             v = VideoReader([videoFolder strjoin(C(1:3),'_') '_side.avi']);
%             v2 = VideoReader([videoFolder strjoin(C(1:3),'_') '_bot.avi']);
            I_side = v.read(1);
            I_bot = v2.read(1);
            
            fig = figure;set(gcf,'Position',[2 42 958 958])
            suptitle([num2str(i) '/' num2str(nVid)])
            subplot(2,1,1);hold on;
            imshow(I_side);
            for input = 1:4
                [x_side(input),y_side(input)] = ginputWhite(1);
                scatter(x_side(input),y_side(input))
            end
            
            subplot(2,1,2);hold on;
            imagesc(I_bot);
            for input = 1:4
                [x_bot(input),y_bot(input)] = ginputWhite(1);
                scatter(x_bot(input),y_bot(input))
            end
            
            xy_side = [x_side',y_side'];
            xy_sideBin = [x_side'>mean(x_side),y_side'>mean(y_side)];
            
            xy_bot = [x_bot',y_bot'];
            xy_botBin = [x_bot'>mean(x_bot),y_bot'>mean(y_bot)];
            
            k = 1;
            for j = 0:1
                for jj = 0:1
                    sideCorner = xy_side(xy_sideBin(:,1)==j & xy_sideBin(:,2)==jj,:);
                    for jjj = 0:1
                        for jjjj = 0:1
                            botCorner = xy_bot(xy_botBin(:,1)==jjj & xy_botBin(:,2)==jjjj,:);
                            cornerUV(k,:) = [sideCorner,botCorner];
                            k = k+1;
                        end
                    end
                end
            end
            err = abs(cornerUV(:,1)-cornerUV(:,3));
            cornerUV(err>mean(err),:) = [];
            
            [boxXYZ,rmse] = dlt_reconstructbol(cornerUV,DLT_params);
            boxXYZ = boxXYZ./1000;
            
            % plot the box
            figure;hold on;
            k = 1;suptitle([num2str(i) '/' num2str(nVid)])
            for j = 1:size(boxXYZ,1)-1
                for jj = j+1:size(boxXYZ,1)
                    if sum(abs(boxXYZ(j,:)-boxXYZ(jj,:))<3)==2
                        edge = [boxXYZ(j,:);boxXYZ(jj,:)];
                        plot3(edge(:,1),edge(:,2),edge(:,3),'k')
                        allEdges{k} = edge;
                        k = k+1;
                    end
                end
            end
            view([45 45])
            save([boxFold strjoin(C(1:3),'_') '_BoundingBox.mat'],'cornerUV','boxXYZ','allEdges')
        end
    end
%end
    
end