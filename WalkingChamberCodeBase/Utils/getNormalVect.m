function [normal_vector] = getNormalVect(P0,P1,P2,plotFig)

normal_vector = cross(P0-P1, P0-P2);
normal_vector = normal_vector / norm( normal_vector );
P_all = [P1;P0;P2];

if isempty(plotFig)
    plotFig = false;
end

if plotFig
    figure;
    h=patch('Faces',1:3,'Vertices',P_all);hold on;
    set(h,'FaceColor','r','EdgeColor','k','LineWidth',2,'FaceAlpha',0.5)
    scatter3(P_all(:,1),P_all(:,2),P_all(:,3));
    quiver3(P0(1), P0(2), P0(3), normal_vector(1), normal_vector(2), normal_vector(3));
    axis equal
    disp(dot((P0 - P1), normal_vector));
    disp(dot((P0 - P2), normal_vector));
end
end