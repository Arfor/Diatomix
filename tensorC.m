function T = tensorC(N,mN, rank)
% T = TensorC(N,mN, component, rank)
% Computes the spherical tensors for the spherical harmonics in the
% rotational basis
arguments
    N
    mN
    rank
end
T = cell(2*rank+1,1);

uN = unique(N); umN = unique(mN); 
matrixSize = length(N);
maxCalcs = (length( uN)*length(umN))^2;

p = -rank:rank;
for ip = 1:length(p)
    calcCount = 1;
    X = nan(maxCalcs,1);
    iCol = cell(maxCalcs,1); %column indices for matrix
    iRow = cell(maxCalcs,1); %row indices for matrix
    for n1 = reshape(uN,1,[]) %hurrah, for loops o_O
        for n2 = reshape(uN,1,[]) %make sure, they're row vectors
            for m1 = -n1:n1
                for m2 = -n2:n2
                    colIdx = find((N == n1)&(mN==m1));
                    rowIdx = find((N == n2)&(mN==m2));

                    if (m1-m2+p(ip))~=0; continue; end
                    x = (-1)^m2 * sqrt((2*n2+1)*(2*n1+1))*...
                        Wigner3j([n2, rank, n1],[-m2, p(ip), m1])*Wigner3j([n2, rank, n1],[0,0,0]);

                    X(calcCount) = x;
                    iRow{calcCount} = rowIdx;
                    iCol{calcCount} = colIdx;
                    calcCount = calcCount+1;
                end
            end
        end
    end
    T{ip} = makeSparseMatrix(X, iCol, iRow, matrixSize);
end

function s = makeSparseMatrix(X, iCol, iRow, mSize)
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
