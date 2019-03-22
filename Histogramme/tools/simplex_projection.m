function y = projection_simplexe(x)
% perform_simplex_projection - compute the projection on the simplex of size rho in the diagonal metric ||.||
%
%   y = perform_simplex_projection_weighted(x,rho);
%
%   x is the projection of y on the set {a \ a >= 0 sum_i a_i = rho },
%
%   Copyright (c) 2014 Jalal Fadili
    

[n,d] = size(x);


[xs,I] = sort(x,2,'descend');
xs = reshape(x(sub2ind([n d],reshape(repmat([1:n]',[1 d]),n*d,1),I(:))),[n d]);

xtmp = (cumsum(xs,2)-1)./cumsum(1./ones(n,d),2);
y = (max(bsxfun(@minus,x,xtmp(sub2ind([n,d],(1:n)',sum(xs>xtmp,2)))),0));



