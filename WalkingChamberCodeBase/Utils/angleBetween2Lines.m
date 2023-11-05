function [angle] = angleBetween2Lines(Va,Vb,Vn)
Va = Va./norm(Va);
Vb = Vb/norm(Vb);
%Vn = cross(Va, Vb);
%Vn = Vn / norm( Vn ); % just to make it unit length
%angle = atan2d(cross(Va,Vb), dot(Va,Vb));%atan2((Vb x Va) . Vn, Va . Vb)
angle = atan2d(dot(cross(Va,Vb),Vn), dot(Va,Vb));%atan2((Vb x Va) . Vn, Va . Vb)

end