function pde_problem(h,H,p_tol)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

% ------------------- Problem definition ------------------
rng(10)

A0 = @(x,y) 1 + 0*x; %sin(x.*y);
A1 = @(x,y) 1 + 0*x; %sin(30*pi*x).*sin(30*pi*y) + 2;
A1_cell = {A1};

% % A0 = p.w. constant on fine mesh, random between [0.01,1]
% grid = 0:1/2^h/6:1; rand_mat = .99*rand(length(grid(1:end-1))) + 0.01;
% A0 = @(x,y) rand_mat(find((x >= grid(1:end-1)).*(x < grid(2:end))),find((y >= grid(1:end-1)).*(y < grid(2:end))));

% rand_v = .99*rand(1,150000) + 0.01;
% A0 = .99*rand(1,150000) + 0.01; % ones(1,150000); %
% %A1_v = .99*rand(1,150000) + 0.01; % ones(1,150000); %

B1 = @(x,y) 1;
f_bulk = @(x,y) sin(pi*x).*sin(pi*y); % 1 + 0*x; % sin(30*pi*x).*sin(pi*y);
f_interface = @(x,y) x + 2*y; % sin(30*pi*x).*sin(pi*y); %1 + 0*x;

% % Non-convex coarse triangles-experiment
% [e,e_coarse,t,t_coarse] = nonconvex_coarse_triangles(p,e,t,p_inter_logical,p_boundary_logical,e_inter,e_coarse,t_coarse);

max_ell = 4;
energy_error_fem_lod = zeros(H+1,max_ell);

for i_cell = 1:length(A1_cell)
parfor j = 0:H % 0:H
    % ------------------- Data structure ------------------
    A1 = A1_cell{i_cell};

    [p,e,t,p_inter_logical,p_boundary_logical,e_inter,e_coarse,t_coarse] = ...
        read_variables(h,j,p_tol);

    % for i = 1:size(e_inter,2)
    %     current_edge = e_inter(1:2,i);
    %     A1_v(i) = A1(mean(p(1,current_edge)),mean(p(2,current_edge)));
    % end

    % p = [------ x-coord --------------------------
    %      ------ y-coord --------------------------
    %      ------ subdomain index ------------------
    %      ------ coincides w/ interf point no -----]

    % -------------------- Solution ---------------------------

    quad_p = 2; % quad_p = 2 means 2 Legendre points
    % FEM solution
    [A,b,M,u_fem] = fem_fcn(A0,A1,B1,f_interface,f_bulk,...
        p,p_inter_logical,p_boundary_logical,...
        e_inter,t,quad_p);

    lod_logical = 0;
    lod_enh_loc_logical = 1;
    if lod_logical, u_lod = zeros(size(A,1),max_ell); end
    if lod_enh_loc_logical, u_lod_enh_loc = zeros(size(A,1),max_ell); end

    tmp = zeros(1,max_ell);
    for i = 1:max_ell % for i = j
        ell = i;
        disp(['Layers = ',num2str(ell)])
        [~,u_lod_enh_loc(:,i)] = patch_fcn(A,b,A0,A1, ... % u_lod(:,i), u_lod_enh_loc(:,i)
            p,p_inter_logical,p_boundary_logical,e,e_inter,e_coarse, ...
            t,t_coarse,lod_logical,lod_enh_loc_logical,ell);
        tmp(i) = norm_diff(u_fem, u_lod_enh_loc(:,i), A, M);
        %disp([i,j,energy_error_fem_lod])
    end
    energy_error_fem_lod(j+1,:) = tmp;
    
    % plot_solution_comp(p,e,t,e_coarse,t_coarse,u_fem,u_lod_enh_loc(:,level_patches),'FEM solution','LOD solution, enhanced localisation');

end
end
energy_error_fem_lod

end



