function [y] = itrnd(dim,x,p,meth)
%  ITRND  Random number generator using inverse transform method
%
%  Usage: [y] = ITRND(dim,x,p,meth)
%
%  where dim  - output array dimension(s)
%        x    - pdf x-values
%        p    - pdf p-values
%        meth - generation method (exact or interp)
%        y    - output array
%
%  Do not call this function directly unless you know what you are doing!
%  Valentin Wyart <valentin.wyart@ens.fr>

if nargin < 4
    error('meth should be ''exact'' or ''interp''.');
end
if nargin < 3
    error('missing input arguments.');
end
if size(x) ~= size(p)
    error('x and p should have the same size.');
end
if numel(dim) == 1
    dim = [dim,1];
end

switch meth
    case 'exact'
        p = p/sum(p);
        c = cumsum(p);
    case 'interp'
        p = p/trapz(x,p);
        c = cumtrapz(x,p);
end

imin = find(c == 0,1,'last');
if isempty(imin)
    imin = 1;
end
imax = find(c == 1,1,'first');
if isempty(imax)
    imax = length(c);
end
c = c(imin:imax);
x = x(imin:imax);

switch meth
    case 'exact'
        [~,i] = max(bsxfun(@gt,c(:),reshape(rand(dim),1,[])),[],1);
        y = reshape(x(i),dim);
    case 'interp'
        y = interp1(c,x,rand(dim));
end

end