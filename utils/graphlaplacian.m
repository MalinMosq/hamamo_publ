function [L,M] = graphlaplacian(p,e,Aj,quad_p) % [L,M]=graphlaplacian(p,e,w)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

% Assemble reciprocal edge length weighted graph Laplacian
n=size(p,2);le=size(e,2);
el=sqrt((p(1,e(2,:))-p(1,e(1,:))).^2+(p(2,e(2,:))-p(2,e(1,:))).^2);
% off diagonal
% L = sparse(e(1,:),e(2,:),-w'./el,n,n,n+2+length(e)); % L = L + L';
L = sparse([],[],[],n,n,n+2+length(e));
L = L;
M = sparse(n,n);
% Lvec=zeros(n,1); Mvec=zeros(n,1);
for i=1:le
    current_edge = [e(1,i),e(2,i)];
    x = p(1,current_edge); y = p(2,current_edge);
    if isa(Aj,'function_handle')
        %I = quad(@(t) Aj(x(1)+t*diff(x), y(1)+t*diff(y)),quad_p)/el(i);
        I = Aj(x(1) + .5*diff(x), y(1) + .5*diff(y))/el(i);
    else % Aj is a vector
        I = Aj(i)/el(i);
    end
    L(current_edge,current_edge) = L(current_edge,current_edge) + I*[1 -1; -1 1];
    M(current_edge,current_edge) = M(current_edge,current_edge) + el(i)/6*[2 1; 1 2];
end


% L = L+spdiags(Lvec,0,n,n);
% M=spdiags(M,0,n,n);
end