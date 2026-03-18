function h = plot_solution(p,e,t,u,fig_title)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

gcf; hold on
% Plots the solution on the interfaces and the bulk areas
h = [];

x = reshape(p(1,e(1:2,:)),2,[]);
y = reshape(p(2,e(1:2,:)),2,[]);
z = reshape(u(e(1:2,:)),2,[]);
h_temp = plot3(x,y,z,'k','linewidth',2);
h = [h; h_temp];
x = reshape(p(1,t(1:3,:)),3,[]);
y = reshape(p(2,t(1:3,:)),3,[]);
z = reshape(u(t(1:3,:)),3,[]);
h_temp = fill3(x,y,z,'r','FaceAlpha',.5,'LineWidth',.25);
h = [h; h_temp];
view(30,30)

if nargin == 5
    title(fig_title)
end


end