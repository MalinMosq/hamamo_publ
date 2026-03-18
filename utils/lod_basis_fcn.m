function phi = lod_basis_fcn(A,C,p_boundary_logical,e_inter,e_coarse,j)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

n_edg = size(e_coarse,2);
phi = zeros(size(A,1),1);
if j <= n_edg && all(p_boundary_logical(e_coarse(:,j)))
    % No basis on boundary of domain. Since j is the middle element there
    % is no need to check if it is part of the global boundary, it will
    % never only be part of the local patch boundary
    return
end

% -------------------- LOD matrix ---------------------------

% [ A   C^T] * [ phi_T] = [      0      ]
% [ C   0  ]   [lambda]   [0 ... 1 ... 0]^T

% Dirichlet conditions
% Remove coarse boundary edges
e_coarse_int = unique(e_inter(3,:));
C_int = C([e_coarse_int, (n_edg+1):end],:);
% Remove fine boundary points
C_int(:,p_boundary_logical) = [];
A_int = A(~p_boundary_logical,~p_boundary_logical);

AC = [A_int C_int'; C_int zeros(size(C_int,1))];

% -------------------- Basis functions ---------------------------

% [Q,R] = qr(AC);
% AC_inv = inv(R)*Q';

rhs_b = zeros(size(C,1),1);
rhs_b(j) = 1;

rhs_check = [zeros(size(A,1),1); rhs_b];
rhs_b = rhs_b([e_coarse_int, (n_edg+1):end]);
rhs_A = zeros(size(A_int,1),1);
rhs = [rhs_A; rhs_b];
if sum(rhs) == 0, error('One basis in LOD is wrong'), end

% phi_lambda_int = AC_inv*rhs;
% phi_lambda_int = R\(Q'*rhs);
phi_lambda_int = AC\rhs;
phi(~p_boundary_logical) = phi_lambda_int(1:sum(~p_boundary_logical));

% % Mean values are checked by multiplying C and phi
% % Should be 1 in one place and 0 in the rest
% C*phi

end