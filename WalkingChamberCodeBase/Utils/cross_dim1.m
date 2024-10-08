function c = cross_dim1(a,b)
% c = cross_dim1(a,b)
% Calculate cross product along the first dimension
% NOTE: auto expansion allowed
c = zeros(max(size(a),size(b)));
c(1,:) = a(2,:).*b(3,:)-a(3,:).*b(2,:);
c(2,:) = a(3,:).*b(1,:)-a(1,:).*b(3,:);
c(3,:) = a(1,:).*b(2,:)-a(2,:).*b(1,:);
end