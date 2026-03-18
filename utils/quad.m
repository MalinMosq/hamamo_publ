function I = quad(f,quad_p)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS
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
