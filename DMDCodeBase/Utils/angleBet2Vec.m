function [theta] = angleBet2Vec(u,v)
nPt = size(u,1);
nPt2 = size(v,1);
if nPt==1
    u = repmat(u,nPt2,1);
end
if nPt2==1
    v = repmat(v,nPt,1);
end
theta = nan(nPt,1);
for i = 1:nPt
    theta(i,1) = atan2(norm(cross(u(i,:),v(i,:))),dot(u(i,:),v(i,:))).*180./pi;
end
end