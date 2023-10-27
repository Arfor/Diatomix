function vecC = SphericalToCartesian(vec)
%vecC = SphericalToCartesian(vecS)
% works for both cell array or numeric arrays
assert(length(vec)==3);
if isnumeric(vec)
    vec = reshape(vec,[],1);
    U = [1/sqrt(2), 0,  -1/sqrt(2);
         1i/sqrt(2),0,  1i/sqrt(2);
         0,         1,      0];
    vecC = U*vec;
else
    assert(iscell(vec));
    [dm,d0,dp] = vec{:};
    x = (dm-dp)/sqrt(2);
    y = 1i*(dm+dp)/sqrt(2);
    z = d0;
    vecC = {x,y,z};
end
end

