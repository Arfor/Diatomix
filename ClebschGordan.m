function cg = ClebschGordan(j1,j2,j,m1,m2,m, opts)
arguments
    j1
    j2
    j
    m1
    m2
    m
    opts.warnings = 1
end
% ClebschGordan(j1,j2,j,m1,m2,m)
% Describes how to get to j through a coupling of |j1,m1> with |j2,m2>.An erbium atom in the GS
% coupled with a photon with sigma- polarisation would be for example 
% |j1,m1>=|6,-6>,|j2,m2>=|1,-1> going to th exited state |j,m>=|7,-7>.
% Example:
% J = 6
% for j=1:13
%     sigm(j)=misc.ClebschGordan(1,6,7,-1,gManif(j));
%     pip(j) = misc.ClebschGordan(1,6,7,0,gManif(j));
%     sigp(j)=misc.ClebschGordan(1,6,7,1,gManif(j));
% end
% j=1:13;
% figure(1); clf;
% hold on;
% plot(j,sigm, DisplayName="\sigma^-")
% plot(j,pip, DisplayName="\pi")
% plot(j,sigp, DisplayName="\sigma^+")
% legend()
% xlim([1,13])
% m = m1+m2;
warning('on','cg:warn')
if ~opts.warnings
    warning('off','cg:warn')
end
% error checking
if ( 2*j1 ~= floor(2*j1) || 2*j2 ~= floor(2*j2) || 2*j ~= floor(2*j) ...
        || 2*m1 ~= floor(2*m1) || 2*m2 ~= floor(2*m2) || 2*m ~= floor(2*m) )
    error('All arguments must be integers or half-integers.');
    return;
end

if m1 + m2 ~= m
    warning('cg:warn','m1 + m2 must equal m.');
    cg = 0;
    return;
end

if ( j1 - m1 ~= floor ( j1 - m1 ) )
    warning('cg:warn','2*j1 and 2*m1 must have the same parity');
    cg = 0;
    return;
end

if ( j2 - m2 ~= floor ( j2 - m2 ) )
    warning('cg:warn','2*j2 and 2*m2 must have the same parity');
    cg = 0;
    return;
end

if ( j - m ~= floor ( j - m ) )
    warning('cg:warn','2*j and 2*m must have the same parity');
    cg = 0;
    return;
end

if j > j1 + j2 || j < abs(j1 - j2)
    warning('cg:warn','j is out of bounds.');
    cg = 0;
    return;
end

if abs(m1) > j1
    warning('cg:warn','m1 is out of bounds.');
    cg = 0;
    return;
end

if abs(m2) > j2
    warning('cg:warn','m2 is out of bounds.');
    cg = 0;
    return;
end

if abs(m) > j
    warning('cg:warn','m is out of bounds.');
    cg = 0;
    return;
end

cg = (-1)^(j1-j2+m) * sqrt(2*j + 1) * Wigner3j([j1 j2 j],[m1 m2 -m]);


% Reference: Clebsch-Gordan Coefficient entry of Eric Weinstein's Mathworld: http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html
