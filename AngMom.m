classdef AngMom
    %ANGMOM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name string
        js (1,:)
        NStates
        constituents
        qnumbers
        coupled
        restrictedJs
        restrictedProjection (1,:)
        states
        storeStates
    end
    
    methods
        function obj = AngMom(js,name, opts)
            arguments
                js
                name
                opts.constituents = [];
                opts.states = [];
                opts.qnumbers = [];
                opts.coupled = [];
                opts.restrictedProjection = [];
                opts.restrictedJs = [];
                opts.storeStates = 1;
            end
            %ANGMOM Construct an instance of this class
            %   Detailed explanation goes here
            obj.name = name;
            obj.js = js;
            obj.NStates = sum(2*(js)+1);
            %set all optional fields if given
            fieldNames = fields(opts);
            for i = 1:length(fieldNames) 
                obj.(fieldNames{i})=opts.(fieldNames{i});
            end
            if isempty(obj.qnumbers)
                obj.qnumbers = [obj.name,"m"+obj.name];
            end
            if isempty(obj.constituents)
                obj.constituents = obj.name;
            end
            if isempty(obj.restrictedProjection)
                obj.restrictedProjection = -max(js):max(js);
            end
            if isempty(obj.restrictedJs)
                obj.restrictedJs = min(js):max(js);
            end
            %if states not given, calculate them. If states given (from
            %addition), do nothing
            if obj.storeStates
                if isempty(obj.states)
                    states = [];
                    for j=js
                        for mj=-j:j
                            if any(obj.restrictedProjection==mj)
                                states=[states; j,mj];
                            end
                        end
                    end
                    obj.states=states;
                end
            end
        end
      
        function coupledAngMom = couple(j1,j2, name, opts) 
            arguments
                j1
                j2
                name
                opts.storeStates = 1
            end
            if any(matches(name,[j1.name,j2.name]))
                error("Names of momenta must be unique!")
            end
            st= [];
            Js = [];
            for j1s=j1.js
                for j2s=j2.js
                    Js=[Js, min(abs(j1s-j2s)):max(j1s+j2s)];
                end
            end
            qs = [j1.constituents,j2.constituents, name,"m"+name];
            j1q = j1.qnumbers(1:end-1);
            j2q = j1.qnumbers(1:end-1);

            
            %horrible implementation, sorry --> actually limiting
            %performance, think of clever way to do it, like in basis
            if all([j1.storeStates,j2.storeStates])&&opts.storeStates
                storeStates = 1;
                nStates = sum(2*Js+1);
                st = nan(nStates, size(j1.states,2)+size(j2.states,2));
               %Now we have to make sure to include the other quantum numbers
                %to know where it came from. Removes degeneracy from
                %projection, should be as long a length(js)
                %using tables adds a lot of time to the process
                sCount = 1;
                uniquej1States = unique(j1.states(:,1:end-1),'rows');
                uniquej2States = unique(j2.states(:,1:end-1),'rows');
                for a=1:height(uniquej1States)
                    j1s = j1.js(a);
                    j1v = uniquej1States(a,:);
                    for b=1:height(uniquej2States)
                        j2s = j2.js(b);
                        j2v = uniquej2States(b,:);
                        for J=abs(j1s-j2s):(j1s+j2s)
                            %construct whole mj part in one go
                            mJ = -J:J;
                            expand1 = repmat([j1v, j2v,J],length(mJ),1);
                            st(sCount:sCount+length(mJ)-1,:) = [expand1,mJ'];
                            sCount = sCount+length(mJ);
                        end
                    end
                end
            else
                st = [];
                storeStates = [];
            end
%             cTab = array2table(st,"VariableNames",qs);
% 
            coupledAngMom = AngMom(Js,name,states=st,constituents=[j1.constituents,j2.constituents, name], qnumbers=qs, coupled=[j1,j2], storeStates=storeStates);   
        end

        function obj = restrictJ(obj,rJs)
            obj.restrictedJs = rJs;
            obj.js = obj.js(ismember(obj.js,rJs));
            obj.states = obj.states(ismember(obj.states(:,matches(obj.qnumbers,obj.name)),rJs),:);
            obj.NStates = height(obj.states);
        end

%         function constructCG(obj)
%             j1 = obj.j1;
%             j2 = obj.j2;
% %             o
%             CGs = cell(1,length(obj.restrictedJs));
%             for k = 1:length(obj.restrictedJs)
%                 J = obj.restrictedJs(k);
%                 CG = CGs{k}; %select CG cell
%                 for mJ = -J:J
%                     
%                 end
%             end
%         end

        function stateTable = printStates(obj)
            stateTable = array2table(obj.states,"VariableNames",obj.qnumbers);
        end
    end

end

