function C = C_fcn(p,e_inter,t,e_coarse,t_coarse)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

% e_coarse must be original e_coarse because e(3,:) refers to columns in
% e_coarse => cannot be removed columns in e_coarse. Same goes for t_coarse
C = spalloc(size(e_coarse,2) + size(t_coarse,2),size(p,2),2*size(e_inter,2)+3*size(t,2)); % Change to spalloc, compute max non-zero entries
n_tri = size(t_coarse,2);
n_edg = size(e_coarse,2);

% Add values to coarse edges
for i = 1:n_edg %unique(e_inter(3,:)) % References to coarse edges are in e_inter(3,:)
    E = e_coarse(:,i);
    if ~all(ismember(E,e_inter(1:2,:)))
        continue
    end
    x_E = p(1,E); y_E = p(2,E);
    length_E = sqrt(diff(x_E).^2 + diff(y_E)^2);
    edges_logical = e_inter(3,:) == i;
    n_e_edg = sum(edges_logical);
    edges = e_inter(1:2,edges_logical);
    p1 = p(1:2,edges(1,:)); p2 = p(1:2,edges(2,:));
    x = [p1(1,:); p2(1,:)]; y = [p1(2,:);p2(2,:)];
    length_e = sqrt(diff(x).^2 + diff(y).^2);
    for current_edge = 1:size(edges,2)
        C(i,edges(:,current_edge)) = ...
            C(i,edges(:,current_edge)) + 1/length_E*length_e(current_edge)*.5; %n_edg*   n_edg*1/n_e_edg
    end
end

% Add values to coarse triangles
for i = 1:n_tri
    T = t_coarse(:,i);
    area_T = polyarea(p(1,T),p(2,T));
    triangles_logical = (t(4,:) == i);
    triangles = t(1:3,triangles_logical);
    % n_p = length(unique(triangles));
    p1 = p(1:2,triangles(1,:));
    p2 = p(1:2,triangles(2,:));
    p3 = p(1:2,triangles(3,:));
    x = [p1(1,:);p2(1,:);p3(1,:)];
    y = [p1(2,:);p2(2,:);p3(2,:)];
    area_t = polyarea(x,y);
    % subd_tri = unique(p(3,t_coarse(:,i)));
    % n_tri_subd = subd_c(subd == subd_tri);
    n_t_tri = sum(triangles_logical);
    for current_triangle = 1:size(triangles,2)
        C(i + size(e_coarse,2),triangles(:,current_triangle)) = ...
            C(i + size(e_coarse,2),triangles(:,current_triangle)) + 1/area_T*1/3*area_t(current_triangle); %n_tri*n_tri_subd*1/n_t_tri
    end
    

end

% area_t = polyarea(reshape(p(1,t(1:3,:)),3,[]),reshape(p(2,t(1:3,:)),3,[]))
% area_T = polyarea(reshape(p(1,t_coarse(1:3,:)),3,[]),reshape(p(2,t_coarse(1:3,:)),3,[]))

end