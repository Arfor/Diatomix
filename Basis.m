classdef Basis
    %BASIS Summary of this class goes here
    %   Detailed explanation goes here
properties
    states
    momenta
    qnumbers
    NStates
end

methods
    function obj = Basis(varargin)
        %BASIS Construct an instance of this class
        %   Creates a basis out of the given angular momentum objects
        angulars = [varargin{:}];
        obj.states = [];
        obj.qnumbers = [];
        obj.NStates = 1;
        for j=angulars
            obj.momenta.(j.name) = j; 
            obj.NStates = obj.NStates*sum(2*(j.js)+1);
            obj.qnumbers = [obj.qnumbers, j.qnumbers];
        end
        obj = calcStates(obj);
    end

    function [vec,states] = prepStates(varargin)
        %%% [vec,states] = prepState(varargin)
        % returns normalised state as column vector and tabel with labels
        base = varargin{1};
        qs = struct(varargin{2:end});
        qCols = matches(base.qnumbers,string(fieldnames(qs)));
        qnames = base.qnumbers(qCols); 
        qs = orderfields(qs,qnames); %orders the q numbers in the right way
        for q=qnames %check that arguments are valid
            if ~any(matches(base.qnumbers,q))%check that all qnumbers are contained in basis
                error("Basis does not contain '%s'",q{:})
            end
        end

        vec = all(base.states(:,qCols)==struct2array(qs),2);
        if ~any(vec)% check if not empty vector
            error("Empty vector, check qnumber values")
        end
        states = base.states(vec,:);
        vec = vec/norm(vec*1.0);
    end
    
    function U = calcTransform(fromBasis,toBasis,coupledMom)
        %%% U = calcTransform(fromBasis,toBasis,coupledMom)
        % Gives you U such that a state w in toBasis can be calculated from a
        % state v in fromBasis as w = U*v. Or the other way: v=conj(U')*w
        nStates = fromBasis.NStates;        
        U = sparse(nStates,nStates); %predefine sparse matrix
        if ~checkConstituents(fromBasis, toBasis, coupledMom)
            return
        end
        %find the relevant columns to compare qnumbers of states to
        JStates = fromBasis.states;
        JName = coupledMom.name;
        [j1name,j2name]=coupledMom.coupled.name;
        j1idx = find(matches(fromBasis.qnumbers,j1name));
        j2idx = find(matches(fromBasis.qnumbers,j2name));
        Jidx = find(matches(fromBasis.qnumbers,[JName,"m"+JName]));
        JStates = JStates(:,[j1idx,j2idx,Jidx]);

        j12States = toBasis.states;
        j1idx = find(matches(toBasis.qnumbers,[j1name,"m"+j1name]));
        j2idx = find(matches(toBasis.qnumbers,[j2name,"m"+j2name]));
        j12States = j12States(:,[j1idx,j2idx]); %just to be sure the ordering is right

        j1AM = coupledMom.coupled(1);
        j2AM = coupledMom.coupled(2);
        [uj1,uj2]=coupledMom.coupled.restrictedJs;
        uJ = coupledMom.restrictedJs;

        maxCalcs = length(uj1)*length(uj2)*length(uJ)*...
            (2*length(uj1)+1)*(2*length(uj2)+1);
        %for loops all the way, but only over the relevant qnumbers
        calcCount = 1;
        CGs = nan(maxCalcs,1);
        iCol = cell(maxCalcs,1); %column indices for U matrix
        iRow = cell(maxCalcs,1); %row indices for U matrix

        skipCount=0;
        for j1 = uj1 %this doesnt loop through uj1!!, but somehow it works? Accidentally vectorised the problem?
            j1idxTo = find(j12States(:,1)==j1); %definitely a smarter way to do this by saving indices and selecting from them (cuts down computation along the way
            j1idxFrom = find(JStates(:,1)==j1);
            for j2 = uj2
                j2idxTo = find(j12States(j1idxTo,3)==j2);
                j2idxFrom = find(JStates(j1idxFrom,2)==j2);
                for J = abs(j1-j2):(j1+j2)
                    JidxFrom = find(JStates(j1idxFrom(j2idxFrom),3)==J);
                    for mj1 = -j1:j1
                        mj1idxTo = find(j12States(j1idxTo(j2idxTo),2)==mj1);
                        for mj2 = -j2:j2
                            mJ = mj1+mj2;
                            if abs(mJ)>J
                                skipCount = skipCount+1;
                                continue
                            end
                            mj2idxTo = find(j12States(j1idxTo(j2idxTo(mj1idxTo)),4)==mj2);
                            j12x = j1idxTo(j2idxTo(mj1idxTo(mj2idxTo)));

                            mJidxFrom = find(JStates(j1idxFrom(j2idxFrom(JidxFrom)),4)==mJ);
                            Jx = j1idxFrom(j2idxFrom(JidxFrom(mJidxFrom)));
                            %if really always the same indices, would only have to do one?
                            if isempty(Jx)||isempty(j12x)
                                continue
                            end
                            CGs(calcCount) = ClebschGordan(j1,j2,J,mj1,mj2,mJ);
                            iRow{calcCount} = Jx;%always seem to be equal?
                            iCol{calcCount} = j12x;
                            calcCount = calcCount+1;
                        end
                    end
                end
            end
        end
        %prune the arrays
        pruneA = ~isnan(CGs);
        iCol = iCol(pruneA);
        iRow = iRow(pruneA);
        CGs = CGs(pruneA);
        %Expand iCol and iRow
        iCols = cell2mat(iCol);
        iRows = cell2mat(iRow);
        CGexp = nan(size(iCols));
        %construct sparse matrix
        count=1;
        for k=1:length(CGs)
            n = length(iRow{k});
            CGexp(count:count+n-1) = CGs(k);
            count = count+n;
        end
        U=sparse(iRows,iCols,CGexp,nStates,nStates);

        if ~isequal(round(U*U',9),speye(nStates)) %check unitarity of matrix (up to numerical errors)
            warning("Matrix not really unitary, check code");
        end
    end
    
    function [reducedStates,reducedMatrix, groupIdx] = reduceMatrix(obj,matrixToReduce,reduceToQs)
        keepCols = matches(obj.qnumbers,reduceToQs);
        [reducedStates,~,groupIdx]=unique(obj.states(:,keepCols),"rows", "stable");
        NStatesReduced = height(reducedStates);
        reducedMatrix = sparse(obj.NStates,NStatesReduced);
        for k = 1:NStatesReduced %loop over all groupidx
            colIdx = find(groupIdx==k);
            summedSparseMat = sum(abs(matrixToReduce(:,colIdx)).^2,2);
            [rows,~,vals]=find(summedSparseMat);
            cols = ones(1,length(rows))*k;
            reducedMatrix(rows',k) = vals';
        end
    end

    function c = checkConstituents(fromBasis, toBasis, coupledMom)
        try isobject(fromBasis.momenta.(coupledMom.name));
            if sum(matches(fieldnames(toBasis.momenta),[coupledMom.coupled.name]))==2
                c=1;
                return
            else
                error("Constituent momenta not contained in toBasis")
            end
        catch
            error("Coupled momentum not contained in fromBasis")
        end
    end            

    function states = getStates(obj,n)
        arguments
            obj
            n=1
        end
        if isnumeric(n)
            states = (array2table(obj.states(n,:), VariableNames=string(obj.qnumbers)));
        end

        if strcmp(n,"all")
            states = (array2table(obj.states(:,:), VariableNames=string(obj.qnumbers)));
        end

        if any(strcmp(obj.qnumbers,n))
            states = obj.states(:,strcmp(obj.qnumbers,n));
        end

    end

    function obj = addToBasis(obj,j)        
        obj.momenta.(j.name) = j; 
        obj.NStates = obj.NStates*sum(2*(j.js)+1);
        obj.qnumbers = [obj.qnumbers, j.qnumbers];
        if ~isempty(obj.states)
            expand1 = repelem(obj.states,height(j.states),1);
            expand2 = repmat(j.states,height(obj.states),1);
            obj.states = [expand1,expand2];
        end
    end

    function obj = calcStates(obj)
        for j = struct2array(obj.momenta)
            if isempty(obj.states)
                obj.states = j.states;
            else %--> do with repmat instead of looping over states 1000x faster! Effectively krondelta
                expand1 = repelem(obj.states,height(j.states),1);
                expand2 = repmat(j.states,height(obj.states),1);
                obj.states = [expand1,expand2];
            end
        end
    end
    
    function Jplus = raisingOperator(obj, mom)
        assert(length(mom.js)==length(unique(mom.js))) %only works now if angular momentum has unique vales
        J = obj.getStates(mom.name);
        mJ = obj.getStates("m"+mom.name);

        uJ = unique(J); umJ = unique(mJ); 
        nStates = length(mJ);
        maxCalcs = (length( uJ)*length(umJ))^2;
        X = nan(maxCalcs,1);
        iCol = cell(maxCalcs,1); %column indices for matrix
        iRow = cell(maxCalcs,1); %row indices for matrix
        
        calcCount = 1;
        for j = reshape(uJ,1,[]) %hurrah, for loops o_O
            for mj1 = -j:(j-1)
                mj2 = mj1+1;
                jLogic = (J == j);
                colIdx = find(jLogic&(mJ==mj1));
                rowIdx = find(jLogic&(mJ==mj2));

                x = sqrt( j*(j+1) - mj1*(mj1+1) ); 

                X(calcCount) = x;
                iRow{calcCount} = rowIdx;
                iCol{calcCount} = colIdx;
                calcCount = calcCount+1;
            end
        end
        Jplus = makeSparseMatrix(obj, X, iCol, iRow, nStates);
    end

    function [Jx,Jy,Jz,Jmin,Jplus] = AngMomOperators(obj, mom)
        % J = AngMomOperators(obj, mom)
        % returns J = {Jx,Jy,Jz,Jmin,Jplus};
        Jplus = raisingOperator(obj, mom);
        Jmin = Jplus';
        Jx = 0.5*(Jplus + Jmin);
        Jy = -1i*0.5*(Jplus - Jmin);
        Jz = 0.5*(Jplus*Jmin - Jmin*Jplus);
    end

    function T = coupleSpherically(obj,j1,j2)
        % returns spherical tensor of rank 2
        % see Brown and Carrington Eq. 5.113-5.118
        T = cell(5,1); %rank 2 tensor
        [J1x,J1y,J1z,~,~] = obj.AngMomOperators(j1);

        [J2x,J2y,J2z,~,~] =  obj.AngMomOperators(j2);
        
        T{1} = (1/2)*(J1x*J2x - J1y*J2y - 1i*(J1x*J2y + J1y*J2x));%T2_-2
        T{2} = (1/2)*(J1x*J2z + J1z*J2x - 1i*(J1y*J2z + J1z*J2y));%T2_-1
        T{3} = (1/sqrt(6))*(2*J1z*J2z - J1x*J2x - J1y*J2y);%T2_0
        T{4} = -(1/2)*(J1x*J2z + J1z*J2x + 1i*(J1y*J2z + J1z*J2y));%T2_1
        T{5} = (1/2)*(J1x*J2x - J1y*J2y + 1i*(J1x*J2y + J1y*J2x));%T2_2
    end

    function s = makeSparseMatrix(~, X, iCol, iRow, mSize)
        % Prune lists
        pruneA = ~isnan(X);
        iCol = iCol(pruneA);
        iRow = iRow(pruneA);
        X = X(pruneA);
        
        %Expand iCol and iRow
        iCols = cell2mat(iCol);
        iRows = cell2mat(iRow);
        Xm = nan(size(iCols));
    
        %construct sparse matrix
        count=1;
        for k=1:length(X)
            n = length(iRow{k});
            Xm(count:count+n-1) = X(k);
            count = count+n;
        end
    
        %return
        s = sparse(iRows,iCols,Xm,mSize,mSize);
    end

end

end

