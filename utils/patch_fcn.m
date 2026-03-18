function [u_global,u_global_enh_loc] = patch_fcn(A_global,b_global,A0,A1,...
    p,p_inter_logical,p_boundary_logical,...
    e,e_inter,e_coarse,t,t_coarse,lod_logical,lod_enh_loc_logical,ell)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

% Initiate (global) solution vectors in order to not omit an output argument
u_global = []; u_global_enh_loc = [];

n_tri = size(t_coarse,2); n_edg = size(e_coarse,2);
n = n_tri + n_edg;

C_global = C_fcn(p,e_inter,t,e_coarse,t_coarse);

I = sparse(n,size(p,2));
elements_subdomains = zeros(n,1);
for j = 1:n_edg
    edge = e_coarse(:,j);
    bulk_points = ismember(p(4,:),edge);
    I(j,edge) = 1;
    I(j,bulk_points) = 1;
    % elements_subdomains(j) = 0;
end
for i = 1:n_tri
    triangle = t_coarse(:,i);
    interface_points = setdiff(p(4,triangle),0);
    bulk_points = ismember(p(4,:),interface_points);
    I(n_edg+i,triangle) = 1;
    I(n_edg+i,interface_points) = 1;
    I(n_edg+i,bulk_points) = 1;
    elements_subdomains(n_edg+i) = unique(p(3,triangle));
end

if ell < Inf
    patches = logical(I*I');
    % Number of levels the patch should spread
    patches = logical(patches^ell);
end


if lod_enh_loc_logical
    % -------------------- LOD enhanced localisation ----------------------
    % Number of neighbouring elements (in subdomain) for every coarse point
    no_of_neighbours = zeros(1,size(p,2));
    % plot_network(p,e,t,e_coarse,t_coarse);
    for current_subdomain = 0:max(p(3,:))
        indices = find(p(3,:) == current_subdomain);
        no_of_neighbours(indices) = sum(I(elements_subdomains == current_subdomain,indices));
    end
    % Values
    val = no_of_neighbours; val(val > 0) = 1./val(val > 0);
    val(p_boundary_logical) = 0;
    
    % Interpolation matrix for triangles
    I_H = sparse(size(p,2),size(p,2));
    for i = 1:n_tri
        current_triangle = t_coarse(:,i);
        points = unique(t(1:3,t(4,:) == i));
        xa = p(1,current_triangle(1)); ya = p(2,current_triangle(1));
        xb = p(1,current_triangle(2)); yb = p(2,current_triangle(2));
        xc = p(1,current_triangle(3)); yc = p(2,current_triangle(3));
        for k = 1:length(points)
            x = p(1,points(k)); y = p(2,points(k));
            s = [xb-xa xc-xa; yb-ya yc-ya]\[x-xa; y-ya];
            I_H(points(k),current_triangle) = [1-s(1)-s(2);s(1);s(2)];
        end
    end
    % Interpolation matrix for edges
    for i = 1:n_edg
        current_edge = e_coarse(:,i);
        points = unique(e(1:2,e(3,:) == i));
        s = sqrt(sum((p(1:2,points) - p(1:2,current_edge(1))).^2)./sum((p(1:2,current_edge(2)) - p(1:2,current_edge(1))).^2));
        I_H(points,current_edge) = [1-s',s'];
    end
end

if lod_logical, Phi = zeros(size(A_global,1),0); Phi_time = 0; end
if lod_enh_loc_logical
    Phi_enh_loc = sparse([],[],[],size(A_global,1),n_edg+n_tri);
    Kv = sparse([],[],[],size(A_global,1),n_edg+n_tri);
    Phi_enh_loc_time = 0;
end

if ell < Inf
for j = 1:size(patches,1)
    % Create data structure corresponding to patch
    % Need the following:
    % p,p_inter_logical,p_boundary_logical,
    % e,e_inter,e_coarse,
    % t,t_coarse

    % If edge is on boundary, skip this edge
    if j <= n_edg
        if all(p_boundary_logical(e_coarse(:,j)))
            continue
        % else, disp([num2str(j),', edge'])
        end
    % else, disp([num2str(j),', triangle'])
    end


    % Edges belonging to patch
        e_coarse_indices = find(patches(j,1:n_edg));
        e_coarse_patch = e_coarse(:,e_coarse_indices);
        e_patch = e(:,ismember(e(3,:),e_coarse_indices));
        e_inter_patch = e_inter(:,ismember(e_inter(3,:),e_coarse_indices));
        e_patch(3,:) = changem(e_patch(3,:),1:length(e_coarse_indices),e_coarse_indices);
        e_inter_patch(3,:) = changem(e_inter_patch(3,:),1:length(e_coarse_indices),e_coarse_indices);

    % Triangles belonging to patch
        t_coarse_indices = find(patches(j,(n_edg+1):end));
        t_coarse_patch = t_coarse(:,t_coarse_indices);
        t_patch = t(:,ismember(t(4,:),t_coarse_indices));
        t_patch(4,:) = changem(t_patch(4,:),1:length(t_coarse_indices),t_coarse_indices);

    % Points belonging to patch
        % Spara de punkter som återfinns i någon triangel eller någon edge.
        occurring_points = union(t_patch(1:3,:),e_patch(1:2,:));
        occurring_points = union(occurring_points,union(e_coarse_patch,t_coarse_patch));
        p_patch = p(:,occurring_points);
        p_inter_logical_patch = p_inter_logical(:,occurring_points);
        p_patch(4,:) = changem(p_patch(4,:),1:length(occurring_points),occurring_points);

        % Add additional boundary to p_boundary_logical_patch
        outside_t = t(1:3,ismember(t(4,:),find(~patches(j,(1+n_edg):end))));
        outside_e = e(1:2,ismember(e(3,:),find(~patches(j,1:n_edg))));
        % Add points that are equivalent to interface points
        outside_e = union(outside_e,find(ismember(p(4,:),outside_e)));
        % All points that belong to the outside of patch incl border
        outside_p = union(outside_e,outside_t);
        % Border
        boundary_p = intersect(outside_p,occurring_points);
        p_boundary_logical_patch = p_boundary_logical;
        p_boundary_logical_patch(boundary_p) = 1;
        p_boundary_logical_patch = p_boundary_logical_patch(:,occurring_points);

        % Global boundary
        p_global_boundary_logical = p_boundary_logical(occurring_points);

    % Change indices
        new = 1:length(occurring_points);
        t_patch(1:3,:) = changem(t_patch(1:3,:),new,occurring_points);
        t_coarse_patch = changem(t_coarse_patch,new,occurring_points);
        e_coarse_patch = changem(e_coarse_patch,new,occurring_points);
        e_inter_patch(1:2,:) = changem(e_inter_patch(1:2,:),new,occurring_points);
        e_patch(1:2,:) = changem(e_patch(1:2,:),new,occurring_points);

    % Incidence matrix
    I_patch = I(patches(j,:),occurring_points);

    % Retrieve A on patch
    A_patch = A_global(occurring_points,occurring_points);
    A_patch(p_boundary_logical_patch,:) = 0;
    A_patch(:,p_boundary_logical_patch) = 0;
    A_patch(p_boundary_logical_patch,p_boundary_logical_patch) = ...
        spdiags(ones(sum(p_boundary_logical_patch),1),0,sum(p_boundary_logical_patch),sum(p_boundary_logical_patch));
    % Retrieve C on patch
    C_patch = C_global(patches(j,:),occurring_points);

    if isa(A0,'double')
        A0_patch = A0(ismember(t(4,:),t_coarse_indices));
    else % function handle
        A0_patch = A0;
    end
    if isa(A1,'double')
        A1_patch = A1(ismember(e_inter(3,:),e_coarse_indices));
    else % function handle
        A1_patch = A1;
    end

    if j <= n_edg
         middle_ind = find(e_coarse_indices == j);
    else
         middle_ind = find(t_coarse_indices == j - n_edg) + length(e_coarse_indices);
    end

    % --------------------- LOD classical -----------------------
    if lod_logical
    Phi(occurring_points,end+1) = lod_basis_fcn(A_patch,C_patch,...
        p_boundary_logical_patch,e_inter_patch,e_coarse_patch,middle_ind);
    end
    % --------------- LOD enhanced localisation -----------------
    if lod_enh_loc_logical
        if j > n_edg
        [Kv_temp,~,~] = ...
            lod_basis_fcn_enh_loc(A_patch,C_patch,A0_patch,A1_patch,...
            p_patch,p_inter_logical_patch,p_boundary_logical_patch,p_global_boundary_logical,...
            e_patch,e_inter_patch,e_coarse_patch,t_patch,t_coarse_patch,...
            middle_ind,I_patch,val(occurring_points),I_H(occurring_points,occurring_points));
        Kv(occurring_points,patches(j,:)) = Kv(occurring_points,patches(j,:)) + Kv_temp;
        %A_lod_enh_loc_total(occurring_points,occurring_points) = A_lod_enh_loc_total(occurring_points,occurring_points) + A_T;
        end
    end

end

else % ell = Inf: Global solution
   % --------------------- LOD classical -----------------------
    if lod_logical
        for j = 1:n
            Phi(:,j) = lod_basis_fcn(A_global,C_global,...
                p_boundary_logical,e_inter,e_coarse,j);
        end
    
        Phi = Phi(:,sum(Phi) > 0);
    end
    % --------------- LOD enhanced localisation -----------------
    if lod_enh_loc_logical
        for j = 1:n
            if j > n_edg
            [Kv_temp,~,~] = ...
                lod_basis_fcn_enh_loc(A_global,C_global,A0,A1,...
                p,p_inter_logical,p_boundary_logical,p_boundary_logical,...
                e,e_inter,e_coarse,t,t_coarse,j,I,val,I_H);
            Kv = Kv + Kv_temp;
            end
        end
    end
end

% -------------------------- Solve the LOD system -------------------------

if lod_logical
    % Solution LOD
    A_lod = Phi'*(A_global*Phi);
    b_lod = Phi'*b_global;
    u_lod = A_lod\b_lod;
    u_global = Phi*u_lod;
end
%save('lod_basis_h6_H2_B1_ell4','Phi','p','e','t','p_boundary_logical','e_coarse','t_coarse','-v7.3')

if lod_enh_loc_logical
    % Add I_H 1_S to every column
    % 1_S-matrix
    % Edges
    one_S_edg = zeros(size(p,2),n_edg);
    linInd = sub2ind(size(one_S_edg),e_coarse,repmat(1:size(e_coarse,2),size(e_coarse,1),1));
    one_S_edg(linInd) = val(e_coarse);
    % Triangles
    one_S_tri = zeros(size(p,2),n_tri);
    linInd = sub2ind(size(one_S_tri),t_coarse,repmat(1:size(t_coarse,2),size(t_coarse,1),1));
    one_S_tri(linInd) = val(t_coarse);
    one_S_ = [one_S_edg, one_S_tri];
    % Add I_H 1_S to every column
    I_H_one_S_ = I_H*one_S_;
    Phi_enh_loc = Kv + I_H_one_S_;

    Phi_enh_loc = Phi_enh_loc(:,sum(Phi_enh_loc) > 0);

    % % Mean values are checked by multiplying C and phi
    % (C_global*Phi_enh_loc) %(:,1:n_edg)

    % Solution, LOD with enhanced localisation
    A_lod_enh_loc = Phi_enh_loc'*(A_global*Phi_enh_loc);
    b_enh_loc = Phi_enh_loc'*b_global;
    u_enh_loc = A_lod_enh_loc\b_enh_loc;
    u_global_enh_loc = Phi_enh_loc*u_enh_loc;
end % if

% figure()
% h = plot_solution(p,e,t,Phi_enh_loc(:,100));

end
