function vecS = CartesianToSpherical(vec)
%vecS = CartesianToSpherical(vec)
assert(length(vec)==3);
if isnumeric(vec)
    vec = reshape(vec,[],1);
    U = [1/sqrt(2), -1i/sqrt(2), 0;
         0,          0,          1;
         -1/sqrt(2), -1i/sqrt(2),0];
    vecS = U*vec;
else
    assert(iscell(vec));
    [x,y,z] = vec{:};
    dm = (x-1i*y)/sqrt(2);
    d0 = z;
    dp = -(x+1i*y)/sqrt(2);
    vecS = {dm,d0,dp};
end
end

