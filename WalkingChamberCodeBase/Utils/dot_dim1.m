function d = dot_dim1(a,b)
% d = dot_dim1(a,b)
% Calculate dot product along the first dimension
% NOTE: auto expansion allowed
d = sum(a.*b,1);
end % dot_dim1