function [genotype] = interp_shortUnconfidentFrames(genotype,conf_thresh,dur_thresh)
nFlies = numel(genotype.trialID);
nBodyPart = size(genotype.bodyPartUV_sideConf{1},2);
nTrial = 15;
nPts = size(genotype.bodyPartUV_sideConf{1},1)./nTrial;
lowConfFrames = cellfun(@(x,y) x<conf_thresh | y<conf_thresh, genotype.bodyPartUV_sideConf,genotype.bodyPartUV_botConf,'UniformOutput',false);

for fly = 1:nFlies
    for bodyPart = 1:nBodyPart
        nAll{fly,bodyPart} = [];
        spdAll{fly,bodyPart} = [];
        tmp = lowConfFrames{fly}(:,bodyPart)';
        tmp = reshape(tmp,[],nTrial);
        bodyPartXYZ_confidenceAdjusted = cell(nTrial,1);
        for trial = 1:nTrial
            [startNdx,endNdx,type] = startEndSeq(tmp(:,trial)');
            startNdx(type == 0) = [];
            endNdx(type == 0) = [];
            bodyPart_tmp = genotype.bodyPartXYZ{fly,bodyPart}((trial-1).*nPts+[1:nPts],:);
            bodyPart_tmp2 = bodyPart_tmp;
            n = endNdx-startNdx+1;
            spd = nan(size(n));
            for i = 1:numel(startNdx)
                if n(i)>dur_thresh
                    bodyPart_tmp2(startNdx(i):endNdx(i),:) = nan;
                else
                    if startNdx(i) == 1 && endNdx(i) == nPts
                        bodyPart_tmp2(:,:) = nan;
                    elseif startNdx(i) == 1
                        bodyPart_tmp2(startNdx(i):endNdx(i),:) = repmat(bodyPart_tmp(endNdx(i)+1,:),n(i),1);
                    elseif endNdx(i) == nPts
                        bodyPart_tmp2(startNdx(i):endNdx(i),:) = repmat(bodyPart_tmp(startNdx(i)-1,:),n(i),1);
                    else
                        for j = 1:3
                            bodyPart_tmp2(startNdx(i):endNdx(i),j) = ...
                                interp1([startNdx(i)-1 endNdx(i)+1],...
                                bodyPart_tmp([startNdx(i)-1 endNdx(i)+1],j),...
                                [startNdx(i):endNdx(i)]);
                        end
                        spd(i) = norm(diff(bodyPart_tmp2([startNdx(i)-1 endNdx(i)+1],:)))./(n(i)+1).*100;% 100 hz hardcoded
                    end
                end
            end
            bodyPartXYZ_confidenceAdjusted{trial,1} = bodyPart_tmp2;
            nAll{fly,bodyPart} = [nAll{fly,bodyPart} n];
            spdAll{fly,bodyPart} = [spdAll{fly,bodyPart} spd];
        end
        genotype.bodyPartXYZ_confidenceAdjusted{fly,bodyPart} = cell2mat(bodyPartXYZ_confidenceAdjusted);
    end
end
figure;histogram(cell2mat(nAll(:,4)'),[0:2:100]);hold on;histogram(cell2mat(nAll(:,5)'),[0:2:100]);
xlabel('duration (frames)');ylabel('count');legend('left wing unconfident','right wing unconfident')
title([genotype.id ' confThresh = ' num2str(conf_thresh)])
end