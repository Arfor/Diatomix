function a = struct2array(s)
%a = struct2array(s)
% converts structure to array via cell
    temp = struct2cell(s);
    a = [temp{:}];
end

