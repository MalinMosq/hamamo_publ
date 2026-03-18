function [A,b,M,u_fem] = fem_fcn(Ai,Aj,Bj,f_interface,f_bulk,...
    p,p_inter_logical,p_boundary_logical,e_inter,t,quad_p)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS


% -------------------- Create A ---------------------------

[A00,A01,A10,A11,M00,M11] = A_fcn(Ai,Aj,Bj,p,p_inter_logical,p_boundary_logical,e_inter,t,quad_p);
A = [A00, A01; A10, A11];
M01 = zeros(size(M00,1),size(M11,2));
M = [M00, M01; M01', M11];

% -------------------- Create b ---------------------------

[b0,b1] = b_fcn(p,p_inter_logical,p_boundary_logical,e_inter,f_interface,f_bulk,t,quad_p);
%[b0,b1] = b_fcn_pw_const(p,p_inter_logical,p_boundary_logical,e_inter,t);

b = [b0;b1];
u_fem = A\b; % Not too big problems

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A00,A01,A10,A11,M00,M11] = A_fcn(Ai,Aj,Bj,p,p_inter_logical,p_boundary_logical,e_inter,t,quad_p)
e_inter = e_inter(1:2,:);
t = t(1:3,:);

% Initiate stiffness matrix A00
%A00 = sparse(length(p_inter_logical),length(p_inter_logical));
A00 = spalloc(length(p_inter_logical),length(p_inter_logical),sum(~p_inter_logical)*6);

% Initiate mass matrix M00
M00 = spalloc(length(p_inter_logical),length(p_inter_logical),sum(~p_inter_logical)*6);

% Area of triangles
x = reshape(p(1,t),3,size(t,2));
y = reshape(p(2,t),3,size(t,2));
% area = .5*abs(det([x(:,1), y(:,1), ones(3,1)]))
area = .5*abs(x(1,:).*(y(2,:)-y(3,:)) + ...
    x(2,:).*(y(3,:)-y(1,:)) + x(3,:).*(y(1,:)-y(2,:)));

% A00: For every triangle, compute
% (A_i \grad u_i^0, \grad w_i^0)_{\Omega_i^0}
for i = 1:size(t,2)
    triangle = t(:,i);
    x = p(1,triangle)';
    y = p(2,triangle)';
    basis = [-1 1 0;-1 0 1].*(~p_boundary_logical(triangle));
    gradphi = [x(2)-x(1),y(2)-y(1);...
            x(3)-x(1),y(3)-y(1)]\basis;
    if isa(Ai,'function_handle')
        A00(triangle,triangle) = A00(triangle,triangle) + ...
            quad2d(Ai,x,y)*(gradphi'*gradphi); %Ai*area(i)*(gradphi'*gradphi);
    else % Ai is a vector
        A00(triangle,triangle) = A00(triangle,triangle) + ...
            Ai(i)*area(i)*(gradphi'*gradphi); %Ai*area(i)*(gradphi'*gradphi);
    end
    % Mass matrix
    M00(triangle,triangle) = M00(triangle,triangle) + ...
        area(i)*(eye(3,3) + ones(3,3))/12;
end


% Dirichlet condition
A00(p_boundary_logical,p_boundary_logical) = ...
    spdiags(ones(sum(p_boundary_logical),1), 0, ...
    sum(p_boundary_logical),sum(p_boundary_logical));
M00(p_boundary_logical,p_boundary_logical) = ...
    spdiags(ones(sum(p_boundary_logical),1), 0, ...
    sum(p_boundary_logical),sum(p_boundary_logical));

% Initiate A01
%A01 = sparse(length(p_inter_logical),length(p_inter_logical));
A01 = spalloc(length(p_inter_logical),length(p_inter_logical),sum(p_inter_logical)*3);
% Initiate A11
A11 = spalloc(length(p_inter_logical),length(p_inter_logical),sum(~p_inter_logical)*6);
% % Initiate M01
% M01 = spalloc(length(p_inter_logical),length(p_inter_logical),sum(p_inter_logical)*3);


el = sqrt((p(1,e_inter(2,:))-p(1,e_inter(1,:))).^2+(p(2,e_inter(2,:))-p(2,e_inter(1,:))).^2);
% For every edge segment, compute
% b(v,w) = (Bj(v_i^0 - v_j^1), w_i^0 - w_j^1)_{\Omega_j^1}
for i = 1:size(e_inter,2)
    ip = e_inter(1:2,i); % interface point
    ip1 = ip(1); ip2 = ip(2);
    x = p(1,ip); y = p(2,ip);
    %diag_add = 1/3*Bj*sqrt(diff(x)^2 + diff(y)^2);
    %I = k*integral(@(t) Bj2(x(1)+t*diff(x), y(1)+t*diff(y)).*(1-t).^2,0,1);
    diag_add = el(i)*quad(@(t) Bj(x(1)+t*diff(x),y(1)+t*diff(y)).*(1-t).^2,quad_p);
    nondiag_add = el(i)*quad(@(t) Bj(x(1)+t*diff(x),y(1)+t*diff(y)).*(1-t).*t,quad_p);

    % Find corresponding point indices belonging to bulk and remove
    % boundary points due to homog Dirichlet condition
    bp1 = find(ip1 == p(4,:)); bp2 = find(ip2 == p(4,:));
    [~,ind1,ind2] = intersect(p(3,bp1),p(3,bp2));
    bp1 = bp1(ind1); bp1 = bp1(~p_boundary_logical(bp1));
    bp2 = bp2(ind2); bp2 = bp2(~p_boundary_logical(bp2));

    % A00: Bulk entries as per (B_j v_i^0, w_i^0)_{\Omega_j^1}
    A00([bp1 bp2],[bp1 bp2]) = A00([bp1 bp2], [bp1 bp2]) + ...
        diag_add * eye(length([bp1 bp2]),length([bp1 bp2]));
    A00(bp1,bp2) = A00(bp1,bp2) + ...
        nondiag_add * eye(length(bp1),length(bp2));
    A00(bp2,bp1) = A00(bp2,bp1) + ...
        nondiag_add * eye(length(bp2),length(bp1));

    % A01, A10: Interface-bulk coupling as per -(B_j v_i^0, w_j^0)_{\Omega_j^1}
    ip1 = ip1(~p_boundary_logical(ip1));
    ip2 = ip2(~p_boundary_logical(ip2));
    A01(bp1,ip1) = A01(bp1,ip1) - diag_add;
    A01(bp2,ip2) = A01(bp2,ip2) - diag_add;
    A01(bp2,ip1) = A01(bp2,ip1) - nondiag_add;
    A01(bp1,ip2) = A01(bp1,ip2) - nondiag_add;

    % % M01, M10: Mass matrix
    % M01(bp1,ip1) = M01(bp1,ip1) + el(i)/3; % Diagonal: same point
    % M01(bp2,ip2) = M01(bp2,ip2) + el(i)/3; % Diagonal: same point
    % M01(bp1,ip2) = M01(bp1,ip2) + el(i)/6; % Non-diagonal: two different basis functions
    % M01(bp2,ip1) = M01(bp2,ip1) + el(i)/6; % Non-diagonal: two different basis functions

    % A11: Interface entries as per (B_j v_j^1, w_j^1)_{\Omega_j^1}
    A11([ip1 ip2],[ip1 ip2]) = A11([ip1 ip2],[ip1 ip2]) + ...
        2*diag_add * eye(length([ip1 ip2]),length([ip1 ip2]));
    A11(ip1,ip2) = A11(ip1,ip2) + 2*nondiag_add;
    A11(ip2,ip1) = A11(ip2,ip1) + 2*nondiag_add;


end

A00 = A00(~p_inter_logical,~p_inter_logical);
A01 = A01(~p_inter_logical,p_inter_logical);
A10 = A01';
A11 = A11(p_inter_logical,p_inter_logical);

M00 = M00(~p_inter_logical,~p_inter_logical);
% M01 = M01(~p_inter_logical,p_inter_logical);
% M10 = M01';

% A11: Compute (A_j \grad u_j^1, \grad w_j^1)_{\Omega_j^1}
new_ind = cumsum(p_inter_logical);
new_e = new_ind(e_inter);
[L,M11] = graphlaplacian(p(1:2,p_inter_logical),new_e,Aj,quad_p); % [L,~] = graphlaplacian(p(1:2,p_inter_logical),new_e,Aj);
% eig(L)
% Dirichlet condition on L
logical_boundary_inter = logical(p_boundary_logical(p_inter_logical));
L(logical_boundary_inter,:) = 0;
L(:,logical_boundary_inter) = 0;
L(logical_boundary_inter,logical_boundary_inter) = ...
    spdiags(ones(sum(logical_boundary_inter),1), 0, ...
    sum(logical_boundary_inter),sum(logical_boundary_inter));
A11 = A11 + L;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = quad2d(f,x,y) % Change to midpoint rule
%I = sum(f(x,y))/3.*polyarea(x,y);
I = f(sum(x)/length(x),sum(y)/length(y)).*polyarea(x,y); % Midpoint rule
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b0,b1] = b_fcn_pw_const(p,p_inter_logical,p_boundary_logical,e_inter,t)
rng(33)
e_constants = rand(1,max(e_inter(3,:)));
t_constants = rand(1,max(t(4,:)));

f_interface = e_constants(e_inter(3,:)); % Divided into value per edge
f_bulk = t_constants(t(4,:)); % Divided into value per triangle

e_inter = e_inter(1:2,:);
t = t(1:3,:);

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
    triangle = t(1:3,i);
    x = p(1,triangle)';
    y = p(2,triangle)';
    % Ändpunktskvadratur % *[1 0 0] för phi_1 - använd hela I
    I = area(i)/3*f_bulk(i); I = [I I I]';
    % Only add values if the basis function is not on the border
    I = I(~p_boundary_logical(triangle));
    triangle = triangle(~p_boundary_logical(triangle));
    b0(triangle) = b0(triangle) + I;
end

b0 = b0(~p_inter_logical);

% b1: For every edge, compute (f,phi_j^1)_{\Omega_j^1}
n = sum(p_inter_logical);
b1 = sparse(n,1);

% Lengths of edges
x = reshape(p(1,e_inter),2,[]);
y = reshape(p(2,e_inter),2,[]);
lengths = sqrt(diff(x).^2 + diff(y).^2);
p_inter_indices = find(p_inter_logical);
%p_interboundary_logical = p_boundary_logical(p_inter_logical);
for i = 1:size(e_inter,2)
    edge = e_inter(1:2,i);
    x = p(1,edge); y = p(2,edge);
    I = lengths(i)/2*f_interface(i); I = [I I]';
    % Only add values if the basis function is not on the border
    I = I(~p_boundary_logical(edge));
    edge = edge(~p_boundary_logical(edge));
    edge = changem(edge,1:length(p_inter_indices),p_inter_indices);
    b1(edge) = b1(edge) + I;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function x = coloncatrld(start, stop)
% % COLONCAT Concatenate colon expressions
% %    X = COLONCAT(START,STOP) returns a vector containing the values
% %    [START(1):STOP(1) START(2):STOP(2) START(END):STOP(END)].
% 
% % Based on
% % https://blogs.mathworks.com/loren/2008/10/13/vectorizing-the-notion-of-colon/
% lengths = stop - start + 1;
% 
% endlocs = cumsum(lengths);
% incr = ones(1, endlocs(end));
% incr(1) = start(1);
% jumps = start(2:end) - stop(1:end-1);
% incr(endlocs(1:end-1)+1) = jumps;
% x = cumsum(incr);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = isPointOnLine(P1, P2, Q, p_tol)

    % Is point Q=[x3,y3] on line through P1=[x1,y1] and P2=[x2,y2]
    % From:
    % https://se.mathworks.com/matlabcentral/answers/351581-points-lying-within-line
    % Normal along the line:
    P12 = P2 - P1;
    L12 = sqrt(P12' * P12);
    N   = P12 / L12;
    % Line from P1 to Q:
    PQ = Q - P1;
    % Norm of distance vector: LPQ = N x PQ
    Dist = abs(N(1) * PQ(2,:) - N(2) * PQ(1,:));
    % Consider rounding errors:
    Limit = 1e-5; %10 * eps(max(abs(cat(1, P1(:), P2(:), Q(:)))))
    R = (Dist < Limit);
end