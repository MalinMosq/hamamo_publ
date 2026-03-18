function refine_mesh(h,H, plotta)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

if mod(h,1) > 0, error('(Fine) mesh level must be integer'), end
if ~isempty(H)
    if mod(H,1) > 0, error('(Coarse) mesh level must be integer'), end
    if H >= h
        error('The coarse mesh level must be lower than the fine mesh level')
    end
end
% ------------------- Create mesh -------------------------
tic
run('utils/data_mesh/mesh_h0.m') % Läs in meshet som finns sparat i filen mesh.m
TH.p = msh.POS(:,1:2);
TH.e = msh.LINES(:,1:2);
TH.t = msh.TRIANGLES(:,1:3);

if plotta
    figure()
    subplot(ceil((h+1)/2),2,1)
    plot_network(TH.p',TH.e',TH.t');
end

level = 0;
track_old_structure = false;
nref = 1;
for i = 1:h
    level = level + nref;
    if level > H,track_old_structure = true; end
%     if level == coarse_mesh_level, t_before = TH.t'; end
    [TH] = refineMesh(TH,track_old_structure,nref);
    msh.POS = [TH.p(:,1:2), zeros(size(TH.p,1),1)];
    msh.LINES = TH.e;
    msh.TRIANGLES = TH.t;
    if level == H
        matlab.io.saveVariablesToScript(['utils/data_mesh/mesh_h', num2str(i)], {'msh'})
    end
    % disp(['Points appended afterwards: ', num2str(all(all(TH.p(1:n,1:2) == p(1:2,:)')))])
    if plotta
        subplot(ceil((h+1)/2),2,i+1)
        plot_network(TH.p',TH.e',TH.t');
        title(num2str(i))
    end
end
matlab.io.saveVariablesToScript(['utils/data_mesh/mesh_h',num2str(i),'_H', num2str(H)], {'msh'})
% t_before
% t_after = TH.t'

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function refineMesh 
% refines a coarse triangulation TH nref times by applying the 
% function REFINE nref times. Suitable restriction and prolongation
% operator between corresponding mesh functions are computed
function [Th] = refineMesh(TH, track_old_structure, nref)
    % default value for number of refinements
    if nargin<3
        nref = 1;
    end
    if ~isfield(TH,'e')
        TH.e = computeEdges(TH);
    end
    for k=1:nref
        [Th] = refine(TH,track_old_structure);
    end

end
% End of function refineMesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function REFINE
% refines the current triangulation by dividing
% each marked element into 2^dim congruent elements (red refinement)
function [T] = refine(T,track_old_structure)
    % Construct data structure
    [np,d] = size(T.p);
    nt = size(T.t,1);
    [e,~] = computeEdges(T);
    ne = size(e,1);
    d2p = sparse(e,e(:,[2,1]),(1:ne)'*[1 1],np,np);

    % New nodes from the mid points of each edge
    newnode = 0.5.*(T.p(e(:,1),:)+T.p(e(:,2),:));
    T.p = [T.p; newnode];
    emarker = (np+1:np+ne)';
    p = zeros(nt,d+1+sum(1:d));
    p(:,1:3) = T.t(:,1:3);
    p(:,4) = emarker(d2p(T.t(:,1)+np*(T.t(:,2)-1)));
    p(:,5) = emarker(d2p(T.t(:,2)+np*(T.t(:,3)-1)));
    p(:,6) = emarker(d2p(T.t(:,3)+np*(T.t(:,1)-1)));
    if track_old_structure
        if size(T.t,2) == 3
            p(:,7) = (1:size(T.t,1))';
        elseif size(T.t,2) == 4
            p(:,7) = T.t(:,4);
        else
            error('Tracking past triangles isnt built for this no of dimensions')
        end
        T.t = [p(:,[1,4,6,7])
              p(:,[4,2,5,7])
              p(:,[6,5,3,7])
              p(:,[4,5,6,7])];
    else
        T.t = [p(:,[1,4,6])
              p(:,[4,2,5])
              p(:,[6,5,3])
              p(:,[4,5,6])];
    end

    
    remove_columns = [];
    for i = 1:ne
        j = ne+1-i;
        if ~any(all(ismember(T.e(:,1:2),e(j,:)),2))
            %e(j,:) = [];
            remove_columns = [remove_columns j];
            continue
        end
        temp = [e(j,1), np+j; np+j e(j,2)];
        remove_columns = [remove_columns j];
        %e(j,:) = [];
        e = [e; temp];
    end
    e(remove_columns,:) = [];
    if track_old_structure
        T.e = trackEdges(T.e,e);
    else
        T.e = e;
    end

end
% End of function REFINE



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function computeEdges 
% computes the edges of the triangulation T
function [e,nmbe] = computeEdges(T)
    np = size(T.p,1);
    e = [T.t(:,[1,2]); T.t(:,[1,3]); T.t(:,[2,3])];
    e = sort(e,2);
    d2p = sparse(e(:,1),e(:,2),1,np,np);
%     full(d2p)
    [e1,e2,nmbe] = find(d2p);
    e = [e1,e2];
end
% End of function computeEdges


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function trackEdges 

function e = trackEdges(e_old,e_new)
% Takes edges variable from previous step and the new edges in order to
% equip new edges with references to old edges
%     disp('e_old'),disp(e_old')

    if all(e_new(1:2:end-1,2) == e_new(2:2:end,1))
        e_new(:,3) = 0;
        if size(e_old,2) == 2
            for i = 1:size(e_old,1)
                e_new(2*i-1:2*i,3) = find(all(ismember(e_old, [e_new(2*i-1,1); e_new(2*i,2)]),2));
            end
        else % already exists tracked edges
            for i = 1:size(e_old,1)
                %[e_new(2*i-1,1); e_new(2*i,2)]
                old_edge = all(ismember(e_old(:,1:2), [e_new(2*i-1,1); e_new(2*i,2)]),2);
                e_new(2*i-1:2*i,3) = e_old(old_edge,3);
            end
        end
    else
        error('Track edges function failed')
    end
    e = e_new;
%     disp('e_new'),disp(e')
end

