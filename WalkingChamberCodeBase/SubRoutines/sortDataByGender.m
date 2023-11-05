function [data,cellLabels,gender] = sortDataByGender(data,cellLabels,gender)

[gender,ndx] = sort(gender);
data = data(ndx,:);
cellLabels = cellLabels(ndx,:);
end