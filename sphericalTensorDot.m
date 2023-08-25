function x = sphericalTensorDot(T1,T2)
%x = tensorDot(T1,T2)
%   Calculates dot product between spherical tensors: See Brown and Carrington 5.111
    assert(length(T1)==length(T2))
    [m,n] = size(T1{1});
    x = sparse(m,n);
    rank = (length(T1)-1)/2;
    p = -rank:rank;

    for ip = 1:length(T1)
        x = x + (-1)^p(ip) * T1{ip} * T2{end+1-ip};
    end
end

