function [norm_e_diff,norm_ltwo_diff] = norm_diff(u1, u2, A, M)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

format long

norm_e_sq_1 = full(u1'*A*u1);
norm_e_sq_2 = u2'*(A*u2);
% norm_e_sq_lod = u_lod'*A_lod*u_lod
% diff = norm_e_sq_lod - norm_e_sq_lod_alt

norm_e_diff = sqrt(abs(norm_e_sq_1 - norm_e_sq_2));
norm_ltwo_diff = sqrt((u1 - u2)'*M*(u1 - u2));

end