function h = plot_solution_comp(p,e,t,e_coarse,t_coarse,u1,u2,title1,title2)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

% Plot solution comparison
figure()
h = [];
subplot(2,2,1)
h_temp = plot_solution(p,e,t,u1); h = [h; h_temp];
view(30,30)
try title(title1), end
subplot(2,2,2)
h_temp = plot_solution(p,e,t,u2); h = [h; h_temp];
view(30,30)
try title(title2), end
subplot(2,2,3)
h_temp = plot_network(p,e,t,e_coarse,t_coarse); h = [h; h_temp];
subplot(2,2,4)
h_temp = plot_solution(p,e,t,u1 - u2); h = [h; h_temp];
view(30,30)
set(gcf, 'Position', [100, 100, 1000, 900])