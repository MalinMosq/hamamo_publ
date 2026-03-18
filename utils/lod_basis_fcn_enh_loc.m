function [phi,A_T,C_T] = lod_basis_fcn_enh_loc(A,C,Ai,Aj,p,p_inter_logical,p_boundary_logical,p_global_boundary_logical,...
    e,e_inter,e_coarse,t,t_coarse,j,I,val,I_H)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

n_edg = size(e_coarse,2);
n_tri = size(t_coarse,2);
phi = sparse(size(A,1),n_edg+n_tri);
A_T = sparse(size(A,1),size(A,2));
C_T = sparse(n_edg+n_tri,size(A,1));

if j <= n_edg
    return
end

% -------------------- LOD matrix ---------------------------

% [ A   C^T] * [K_T 1_S ] = [-A_T I_H 1_S      ]
% [ C   0  ]   [lambda_S]   [C_T(b_S - I_H 1_S)]
% phi_T = K_T 1_S + I_H_one_S

% Dirichlet conditions
% Remove coarse boundary edges
e_coarse_int = unique(e_inter(3,:));
C_int = C([e_coarse_int, (n_edg+1):end],:);

% Remove fine boundary points
C_int(:,p_boundary_logical) = [];
A_int = A(~p_boundary_logical,~p_boundary_logical);

% AC matrix
AC = [A_int C_int'; C_int zeros(size(C_int,1))];

% -------------------- Basis functions ---------------------------

% [Q,R] = qr(AC);
% AC_inv = inv(R)*Q';
% disp([time_qr, time_inv])

one_T = zeros(size(A,1),1);

T = t_coarse(:,j - n_edg);


% A_T
if isa(Ai,'function_handle')
    Ai_T = @(x,y) Ai(x,y).*inpolygon(x,y,p(1,T),p(2,T));
else
    Ai_T = Ai.*(t(4,:) == j - n_edg);
end
if isa(Aj,'function_handle')
    Aj_T = @(x,y) Aj(x,y)*1/2.*inpolygon(x,y,p(1,T),p(2,T));
else
    e_T = (sum(ismember(e_coarse(1:2,e_inter(3,:)),p(4,T))) == 2);
    Aj_T = Aj*1/2.*e_T;
end
Bj = @(x,y) 1/2*inpolygon(x,y,p(1,T),p(2,T));
[A00_T,A01_T,A10_T,A11_T] = A_fcn(Ai_T,Aj_T,Bj,p,p_inter_logical,p_boundary_logical,e_inter,t,2);
A00_T(p_boundary_logical(~p_inter_logical),p_boundary_logical(~p_inter_logical)) = 0;
A_T = [A00_T,A01_T;A10_T,A11_T];

% C_T
neighbouring_edges = logical(sum(ismember(e_coarse, p(4,T))) == 2);
current_elem = zeros(n_edg+n_tri,1);
current_elem(j) = 1; current_elem(neighbouring_edges) = .5;
C_T = C.*current_elem;

% RHS using 1_S
rhs_A = sparse(size(p,2),1);%size(C,1));
rhs_C = sparse(size(I,1),1);%,size(C,1));
I = I*I';

neighbours = find((I(j,:)));
for k = 1:length(neighbours)

    % % I_H 1_S
    if neighbours(k) > n_edg % If neighbouring elmnt is a triangle
        S_coarse_index = neighbours(k) - n_edg;
        S = t_coarse(:,S_coarse_index);
        one_S = zeros(size(p,2),1); one_S(S) = val(S);
        I_H_1_S = I_H*one_S;
    else % Neighbouring elmnt is an edge
        S = e_coarse(:,neighbours(k));
        if ~ismember(neighbours(k),e_coarse_int)
            continue
        end
        one_S = zeros(size(p,2),1); one_S(S) = val(S);
        I_H_1_S = I_H*one_S;
    end

    C_T_b_S = zeros(n_edg+n_tri,1);

    % C_T*bubble_S = 0 or 1 (triangle) or .5 (edge)
    if neighbours(k) == j
        C_T_b_S(neighbours(k)) = 1;
    elseif neighbours(k) <= n_edg && neighbouring_edges(neighbours(k))
        C_T_b_S(neighbours(k)) = .5;
    end

    % RHS
    rhs_A = sparse(-A_T*I_H_1_S);
    rhs_A(p_boundary_logical) = [];
    rhs_C = C_T_b_S - C_T*I_H_1_S;
    rhs_C = rhs_C([e_coarse_int, (n_edg+1):end]);

    rhs = [rhs_A; rhs_C];
    
    phi_int = AC\rhs;
    phi(~p_boundary_logical,neighbours(k)) = phi_int(1:sum(~p_boundary_logical));

end

% % Check that all boundary points are 0 (hard-coded, just extra check)
% any(any(phi(p_boundary_logical,:)))

end


% Need another A_fcn with different quadrature and diff Dirichlet cond
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
            quadd2d(Ai,x,y)*(gradphi'*gradphi); %Ai*area(i)*(gradphi'*gradphi);
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
% L(logical_boundary_inter,logical_boundary_inter) = ...
%     spdiags(ones(sum(logical_boundary_inter),1), 0, ...
%     sum(logical_boundary_inter),sum(logical_boundary_inter));
A11 = A11 + L;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = quadd2d(f,x,y)
% I = f(mean(x),mean(y)).*polyarea(x,y); % sum(f(x,y))/3.*polyarea(x,y);
I = f(sum(x)/length(x),sum(y)/length(y)).*polyarea(x,y); % Faster, midpoint
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = quad(f,quad_p)
% Gaussian quadrature:
% int_{-1}^1 f(x) dx ~= sum_1^n w_i f(x_i)
% This code is written for integral from 0 to 1 instead

if quad_p == 2 % 2 points
    w = [1 1];
    p = [-1/sqrt(3), 1/sqrt(3)];
elseif quad_p == 3 % 3 points
    w = [5/9 8/9 5/9];
    p = [-sqrt(3/5) 0 sqrt(3/5)];
elseif quad_p == 4 % 4 points
    w = [(18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36 (18-sqrt(30))/36];
    p = [sqrt(3/7 - 2/7*sqrt(6/5)) -sqrt(3/7 - 2/7*sqrt(6/5)) ...
        sqrt(3/7 + 2/7*sqrt(6/5)) -sqrt(3/7 + 2/7*sqrt(6/5))];
elseif quad_p == 5 % 5 points
    w = [128/225 (322+13*sqrt(70))/900 (322+13*sqrt(70))/900 ...
        (322-13*sqrt(70))/900 (322-13*sqrt(70))/900];
    p = [0 1/3*sqrt(5-2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) ...
        1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5+2*sqrt(10/7))];
end
p = (p+1)/2;
I = f(p)*w';
I = I/2; % Multiply by 1/2 due to interval [0 1] instead of [-1 1]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = quad2d(f,x,y)
x = reshape(x,3,[]); y = reshape(y,3,[]);
% Check: area = sum(polyarea(x,y))
I = sum(sum(f(x,y))/3.*polyarea(x,y));
end