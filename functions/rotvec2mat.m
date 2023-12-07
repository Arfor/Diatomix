function rotMat = rotvec2mat(k)
%ROTVEC2MAT Basically a copy of rotvec2mat3d, but just this function such
%that you dont have to download the entire package from MATLAB
% Uses the Rodrigues formula en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation

assert(length(k)==3); %make sure its a 3D vector
t = norm(k); %get norm/angle of rotation
r = k/t; %normalise vector

K = [0, -r(3), r(2);
     r(3), 0, -r(1);
     -r(2), r(1), 0];
K2 = K*K;

rotMat = eye(3) + sin(t)*K + (1-cos(t))*K2;
end

