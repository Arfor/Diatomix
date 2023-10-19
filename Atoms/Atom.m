classdef Atom
    %ATOM SuperClass 
    properties 

    end
    
    methods
        function info = AtomInfoString(obj)
            %info = AtomInfoString(obj)
            %returns a string with all info of the Atom
            ps = properties(obj);
            info = [];
            for k = 1:length(ps)
                p = ps{k};
                if strcmp(p,'name')
                    info = [info; sprintf("%s = %s",p,obj.(p))];
                else
                    info = [info; sprintf("%s = %.4g",p,obj.(p))];
                end
            end
        end
    end
end