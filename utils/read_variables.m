function [p,e,t,p_inter_logical,p_boundary_logical,e_inter,e_coarse,t_coarse] = ...
    read_variables(h,H,p_tol)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

if mod(h,1), error('Mesh level must be integer'), end
% h = 2^(h-1);
% ------------------- Create mesh -------------------------
if ~isempty(H)
    if mod(H,1), error('Mesh level must be integer'), end

    run(['utils/data_mesh/mesh_h', num2str(h), '_H', num2str(H),'.m']) % Load fine mesh
    p = msh.POS';
    e = msh.LINES';
    t = msh.TRIANGLES';

    %H = 2^(H-1);
    run(['utils/data_mesh/mesh_h', num2str(H),'.m']) % Load coarse mesh
%     p_coarse = msh.POS(:,:)';
    e_coarse = msh.LINES(:,1:2)';
    t_coarse = msh.TRIANGLES(:,1:3)';
else
    run(['utils/data_mesh/mesh_h', num2str(h),'.m']) % Load fine mesh
    p = msh.POS';
    e = msh.LINES';
    t = msh.TRIANGLES';
    e_coarse = [];
    t_coarse = [];
end

nonconvex_case = 0;
if nonconvex_case
    [e,e_coarse,t,t_coarse] = nonconvex_coarse_triangles(p,e,t,e_coarse,t_coarse);%,p_inter_logical,p_boundary_logical,e_inter,e_coarse,t_coarse);
end

[e_inter,p_inter_logical,p_boundary_logical] = select_inter(e,p,p_tol);


% check_fine_coarse_grid(1,0,p,e_inter,t,p_coarse,e_coarse,t_coarse)

% Enumerate subdomains and add interface points to bulk subdomains
[p,p_inter_logical,p_boundary_logical,t,t_coarse] = ...
    enumerate_subdomains(p,p_inter_logical,p_boundary_logical,e_inter(1:2,:),t,t_coarse);

% Reorder points to get interface points last (corresponding to 
% [A00 A01; A10 A11])
[p,p_inter_logical,p_boundary_logical,e,e_inter,e_coarse,t,t_coarse] = reorder_points(p,p_inter_logical,p_boundary_logical,e,e_inter,e_coarse,t,t_coarse);

% figure('Position',[0 0 1500 1500])
% plot_network(p,e,t,e_coarse,t_coarse);
% exportgraphics(gca,'figures/nonconvex_domain.png')

% p = [------ x-coord --------------------------
%      ------ y-coord --------------------------
%      ------ subdomain index ------------------
%      ------ coincides w/ interf point no -----]

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [e_inter,p_inter_logical,p_boundary_logical] = select_inter(e,p,p_tol)

e_inter_logical = zeros([1,length(e(1,:))]);
% Find all edges that are not part of the boundary
for i = 1:size(e,2)
    if ~isOnSameBoundary(e(1,i),e(2,i),p,p_tol)
        e_inter_logical(i) = 1;
    end
end
e_inter = e(:,logical(e_inter_logical));
p_inter_logical = ismember(1:length(p(1,:)), e_inter(1:2,:));
% p_inter = p(:,p_inter_logical);
p_boundary_logical = logical((abs(p(1,:))<p_tol) + (abs(p(2,:))<p_tol) + ...
    (abs(p(1,:)-1)<p_tol) + (abs(p(2,:)-1)<p_tol));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [on_same_boundary] = isOnSameBoundary(p1,p2,p,p_tol)
% Går att effektivisera genom att bara undersöka de punkter som faktiskt är
% på boundary?!
    xmax = max(p(1,:)); ymax = max(p(2,:));
    p1 = p(:,p1);
    p2 = p(:,p2);
    on_same_boundary = 0;
    if ((abs(p1(1)) < p_tol) && (abs(p2(1)) < p_tol)) ||...
        ((abs(p1(2)) < p_tol) && (abs(p2(2)) < p_tol)) ||...
        ((abs(p1(1)-xmax) < p_tol) && (abs(p2(1)-xmax) < p_tol)) ||...
        ((abs(p1(2)-ymax) < p_tol) && (abs(p2(2)-ymax) < p_tol))
        on_same_boundary = 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,p_inter_logical,p_boundary_logical,e,e_inter,e_coarse,t,t_coarse] = ...
    reorder_points(p,p_inter_logical,p_boundary_logical,e,e_inter,e_coarse,t,t_coarse)
% Reorder p in order bulk - interface. Change indices in e and t to match
% this order

% Change fine scale mesh indices
ind_bulk = find(p_inter_logical == 0); ind_inter = find(p_inter_logical == 1);
ind_from = [ind_bulk, ind_inter];
ind_to = 1:length(ind_from);
p = p(:,[ind_bulk,ind_inter]); % Reorder p
p(4,:) = changem(p(4,:),ind_to,ind_from); % Change references to interface points
p_inter_logical = p_inter_logical([ind_bulk, ind_inter]); % Reorder p_inter_logical
p_boundary_logical = p_boundary_logical([ind_bulk, ind_inter]); % Reorder p_boundary_logical

t(1:3,:) = changem(t(1:3,:), ind_to, ind_from); % Change references to points
try
    t_coarse(1:3,:) = changem(t_coarse(1:3,:), ind_to, ind_from); % Change references to points
end
e(1:2,:) = changem(e(1:2,:), ind_to, ind_from); % Change references to points
e_inter(1:2,:) = changem(e_inter(1:2,:), ind_to, ind_from); % Change references to points
try 
    e_coarse(1:2,:) = changem(e_coarse(1:2,:), ind_to, ind_from); % Change references to points
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_new,p_inter_logical,p_boundary_logical,t,t_coarse] = enumerate_subdomains(...
    p,p_inter_logical,p_boundary_logical,e_inter,t,t_coarse)

% Find subdomains from edges
subdomains_p = readmatrix(...
    'utils/data_domain/subdomains_p.txt');
    % , delimitedTextImportOptions('DataLines')) %,[1 Inf]
points = readmatrix(...
    'utils/data_domain/points.txt');
[~,domainp2meshp] = ismember(points,p(1:2,:)','rows');
p_new = [p;zeros(size(p_inter_logical))];

t_coarse_orig = t_coarse;
t_orig = t;

on_interface = zeros(1,size(p,2));
on_interface(reshape(e_inter,1,numel(e_inter))) = 1;

for i = 1:size(subdomains_p,1)
    current_subdomain = subdomains_p(i,subdomains_p(i,:)>0);
    current_subdomain = domainp2meshp(current_subdomain);

    % Find all points in subdomain
    t_coarse_in_subd = all(inpolygon(reshape(p(1,t_coarse_orig),3,[]),reshape(p(2,t_coarse_orig),3,[]),p(1,current_subdomain),p(2,current_subdomain)));
    t_in_subd = ismember(t(4,:),find(t_coarse_in_subd));
    in = zeros(1,size(p,2)); in(t_orig(1:3,t_in_subd)) = 1; % points in subdomain
    % Mark all points INSIDE subdomain with subdomain number
    p_new(3,logical(max((in-on_interface),0))) = i;

    % Find all points on subdomain boundary
    subd_boundary = find(in);
    subd_boundary = subd_boundary(ismember(subd_boundary,e_inter));

    % Add new subdomain boundary/interface points to the end of p
    subd_boundary_new = (size(p_new,2) + 1):(size(p_new,2) + length(subd_boundary));
    p_new = [p_new [p(1:2,subd_boundary); i*ones(1,length(subd_boundary)); subd_boundary]];
    p_boundary_logical = [p_boundary_logical, p_boundary_logical(subd_boundary)];

    % Change interface indices in triangles
    t(1:3,t_in_subd) = changem(t(1:3,t_in_subd),subd_boundary_new,subd_boundary);
    try
    t_coarse(1:3,t_coarse_in_subd) = changem(t_coarse(1:3,t_coarse_in_subd),subd_boundary_new,subd_boundary);
    end

end %for

p_inter_logical = logical([p_inter_logical, ...
    zeros(1,size(p_new,2)-size(p,2))]);

end %function




function check_fine_coarse_grid(check_triangles,check_edges,p,e,t,p_coarse,e_coarse,t_coarse)

% Visual check that t and t_coarse correspond to each other
if check_triangles
    plot_network(p,e,t_coarse,t)
    for i = 1:length(t_coarse(1,:))
        x = [p_coarse(1,t_coarse(1,i)) p_coarse(1,t_coarse(2,i)) p_coarse(1,t_coarse(3,i))];
        y = [p_coarse(2,t_coarse(1,i)) p_coarse(2,t_coarse(2,i)) p_coarse(2,t_coarse(3,i))];
        h1 = fill(x,y,'r');
        pause(.2)
        corr_small_triangles = t(4,:) == i;
        x = [p(1,t(1,corr_small_triangles))' p(1,t(2,corr_small_triangles))' p(1,t(3,corr_small_triangles))']';
        y = [p(2,t(1,corr_small_triangles))' p(2,t(2,corr_small_triangles))' p(2,t(3,corr_small_triangles))']';
        h2 = fill(x,y,'b', 'EdgeColor','k');
        pause(.2)
    end
end

% Visual check that e and e_coarse correspond to each other
if check_edges
    plot_network(p,e,t_coarse,t)
    for i = 1:length(e_coarse(1,:))
        x = [p_coarse(1,e_coarse(1,i)), p_coarse(1,e_coarse(2,i))];
        y = [p_coarse(2,e_coarse(1,i)), p_coarse(2,e_coarse(2,i))];
        h1 = plot(x,y,'r','LineWidth',3);
        pause(1)
        corr_small_edges = e(3,:) == i;
        x = [p(1,e(1,corr_small_edges)), p(1,e(2,corr_small_edges))];
        y = [p(2,e(1,corr_small_edges)), p(2,e(2,corr_small_edges))]; 
        h2 = plot(x,y,'b','LineWidth',2);
        pause(1)
    end
end

end















