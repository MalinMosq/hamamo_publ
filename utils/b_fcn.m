function [b0,b1] = b_fcn(p,p_inter_logical,p_boundary_logical,e_inter,f_interface,f_bulk,t,quad_p)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

e_inter = e_inter(1:2,:);
t = t(1:3,:);

n = sum(p_inter_logical);
b1 = sparse(n,1);

p_interboundary_logical = p_boundary_logical(p_inter_logical);
p_inter_indices = find(p_inter_logical);
for i_ = 1:n
    if p_interboundary_logical(i_)
        continue
    end
    i = p_inter_indices(i_);
    nbrs_index = (e_inter == i);
    nbrs_index = nbrs_index([2 1],:);
    x = [p(1,i)*ones(sum(sum(nbrs_index)),1),p(1,e_inter(nbrs_index))'];
    y = [p(2,i)*ones(sum(sum(nbrs_index)),1),p(2,e_inter(nbrs_index))'];
    f = @(t) f_interface(x(:,1) + t.*(x(:,2)-x(:,1)), y(:,1) + t.*(y(:,2)-y(:,1))).*...
            (1-t).*sqrt((x(:,2) - x(:,1)).^2 + (y(:,2) - y(:,1)).^2);
    %b1(i_) = sum(integral(f,0,1,'ArrayValued',1));
    b1(i_) = sum(quad(f,quad_p));
end

% b0: For every triangle, compute
% (f,phi_i^0)_{\Omega_i^0}

n = length(p_inter_logical);
b0 = sparse(n,1);

% Area of triangles
x = reshape(p(1,t),3,size(t,2));
y = reshape(p(2,t),3,size(t,2));
area = .5*abs(x(1,:).*(y(2,:)-y(3,:)) + ...
    x(2,:).*(y(3,:)-y(1,:)) + x(3,:).*(y(1,:)-y(2,:)));

for i = 1:size(t,2)
    triangle = t(:,i);
    x = p(1,triangle)';
    y = p(2,triangle)';
    % Ändpunktskvadratur % *[1 0 0] för phi_1 - använd hela I
    I = area(i)/3*f_bulk(x,y);
    % Only add values if the basis function is not on the border
    I = I(~p_boundary_logical(triangle));
    triangle = triangle(~p_boundary_logical(triangle));
    b0(triangle) = b0(triangle) + I;
end

b0 = b0(~p_inter_logical);

end
