classdef Hamiltonian < handle
    %HAMILTONIAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hyperfine
        zeeman
        dc_stark
        ac_stark
        Basis
        Fields
        Molecule
        ReferenceAxis
        electricGradient % a.k.a quadrupole operator, tensor rank 2
        dipoleOperator  %tensor rank 1
        maxN
    end
    
    methods
        function obj = Hamiltonian(opts)
            arguments
                opts.Molecule = []
                opts.maxN = 1;
                opts.ReferenceAxis = [1,0,0];
                opts.Fields = [];
            end
            %HAMILTONIAN Construct an instance of this class
            %set all optional fields if given
            fieldNames = fields(opts);
            for i = 1:length(fieldNames) 
                obj.(fieldNames{i})=opts.(fieldNames{i});
            end
            if isempty(obj.Fields)
                obj.Fields.E.dir = obj.ReferenceAxis;
                obj.Fields.B.dir = obj.ReferenceAxis;
                obj.Fields.I.dir = obj.ReferenceAxis;
                obj.Fields.I.pol = [1,0,0]; 
            end
            if ~isempty(opts.Molecule)
                C = Constants;
                Mol = opts.Molecule;
                obj.Molecule=Mol;
                Atom1 = Mol.Atom1;
                Atom2 = Mol.Atom2;
                i1 = AngMom(Atom1.spin,"i1");
                i2 = AngMom(Atom2.spin,"i2");
                N = AngMom([0:obj.maxN],"N");
                I = couple(i1,i2,"I");      %couple the two nuclear momenta
                Fi1 = couple(i1,N,"Fi1");   %couple only a single nucleus to the rotation
                Fi2 = couple(i2,N,"Fi2");   %couple only a single nucleus to the rotation
                % F = couple(I,N,"F");        %fully coupled
                UCBasis = Basis(N,i1,i2);
                obj.Basis = UCBasis;
                SCBasis = Basis(N,I); 
                Fi1Basis = Basis(Fi1,i2);
                Fi2Basis = Basis(Fi2,i1);
                % FCBasis = Basis(F);
                nStates = UCBasis.NStates;
                SCtoUC = calcTransform(SCBasis,UCBasis,I); %gives you U such that UCBasis = U*SCBasis (UC is a sparse matrix) (unitarity is checked)
                Fi1toUC = calcTransform(Fi1Basis,UCBasis,Fi1);
                Fi2toUC = calcTransform(Fi2Basis,UCBasis,Fi2);
                %% Set standard Fields
                % Fields.B = struct(value=0,dir=[0,0,1]);
                % Fields.E = struct(value=0,dir=[0,0,1]);
                % Fields.I = struct(value=0,dir=[0,0,1], pol = [0,1,0]);

                %% Hyperfine
                %Rigid Rotor
                st = UCBasis.getStates("all");
                obj.hyperfine.rigidRotor = obj.rotational(st.N,Mol.Brot,0); 

                % Nuclear Electric Quadrupole
                h_electricQuadrupole1 = obj.electricQuadrupole(st.i1,st.mi1,UCBasis);
                h_electricQuadrupole2 = obj.electricQuadrupole(st.i2,st.mi2,UCBasis);
                obj.hyperfine.electricQuadrupole = Mol.Q1*h_electricQuadrupole1 + Mol.Q2*h_electricQuadrupole2;

                %Scalar Spin-spin coupling
                st = SCBasis.getStates("all");
                spinSpinScalar = obj.scalarNuclear(Mol.c4,st.I,st.i2, st.i1);
                obj.hyperfine.spinSpinScalar =  SCtoUC'*spinSpinScalar*SCtoUC; 

                %Scalar Spin-Rotation Coupling 
                st = Fi1Basis.getStates("all");
                h_spinRot1 = obj.scalarNuclear(Mol.c1,st.Fi1,st.N, st.i1);
                st = Fi2Basis.getStates("all");
                h_spinRot2 = obj.scalarNuclear(Mol.c2,st.Fi2,st.N, st.i2);
                obj.hyperfine.spinRotation = Fi1toUC'*h_spinRot1*Fi1toUC + Fi2toUC'*h_spinRot2*Fi2toUC;

                %Tensor Spin-Spin Coupling
                obj.hyperfine.spinSpinTensor = obj.tensorNuclear(Mol.c3,i1,i2, UCBasis); 

                obj.hyperfine.total = obj.makeHyperfine();
                %% Zeeman
                g1 = Atom1.gFactor; sig1 = Atom1.nuclearShielding;
                g2 = Atom2.gFactor; sig2 = Atom2.nuclearShielding;
                st = UCBasis.getStates("all");
                % obj.zeeman = obj.makeZeeman(g1*(1-sig1)*C.muN,st.mi1) + obj.makeZeeman(g2*(1-sig2)*C.muN,st.mi2) + obj.makeZeeman(Mol.gr*C.muN,st.mN); 
                obj.zeeman = C.muN *( g1*(1-sig1)*obj.makeZeeman(UCBasis,i1) + g2*(1-sig2)*obj.makeZeeman(UCBasis,i2) + Mol.gr*obj.makeZeeman(UCBasis,N) ); 
                %% DC Stark
                obj.dc_stark = obj.makeDCStark(Mol.d0,st.N,st.mN, dir=obj.Fields.E.dir);
                %% AC Stark
                obj.ac_stark = obj.makeACStark(Mol.a0,Mol.a2,st.N,st.mN, Field=obj.Fields.I);
            end

        end
        function [H,s] = scalarNuclear(~,c4,jtot,j1,j2)
            n = length(jtot);
            s = 0.5*c4*(jtot.*(jtot+1) - j1.*(j1+1) - j2.*(j2+1));
            H = spdiags(s,0,n,n); %create diagonal sparse matrix
        end
        function H = tensorNuclear(obj,c3,j1,j2, basis)           
            assert(length(j1.js)==1 && length(j2.js)==1); %Only checked for single vectors, do some thorough checking of basis.AngMomOperators first if you want to extend functionality!
            I1 = j1.js;
            I2 = j2.js;
            N = basis.getStates("N");
            mN = basis.getStates("mN");
            
            if isempty(obj.electricGradient)
                T2C = tensorC(N,mN, 2); 
                obj.electricGradient = T2C;
            else
                T2C = obj.electricGradient;
            end
            T2_Coupled = basis.coupleSpherically(j1,j2);
            
            H = sqrt(6)*c3*sphericalTensorDot(T2C,T2_Coupled); 
        end
        function [H,r] = rotational(~, N, Brot, Drot)
            n = length(N);
            N2 = N.*(N+1);
            r = Brot*N2 - Drot*N2.*N2;
            %r = Brot*N.^2 - Drot*(N.^2).*(N.^2);
            H = spdiags(r ,0,n,n); %create diagonal sparse matrix
        end

        function H = electricQuadrupole(obj,I,mI,basis) 
            if isempty(obj.electricGradient)
                N = obj.Basis.getStates("N");
                mN = obj.Basis.getStates("mN");
                obj.electricGradient = tensorC(N,mN, 2);
            end
            QM = obj.quadrupoleMoment(I,mI);
            H = sphericalTensorDot(QM,obj.electricGradient)/4;
        end
        function T = quadrupoleMoment(obj,I,mI)
            %Calculate nuclear electric quadrupole moment and return 2nd
            %order Tensor
            T = cell(5,1);

            uI = unique(I); umI = unique(mI); 
            matrixSize = length(I);
            maxCalcs = (length( uI)*length(umI))^2;
            
            p = -2:2;
            for ip = 1:length(p)
                calcCount = 1;
                X = nan(maxCalcs,1);
                iCol = cell(maxCalcs,1); %column indices for matrix
                iRow = cell(maxCalcs,1); %row indices for matrix
                for I1 = reshape(uI,1,[]) %For a single atom, this will just be one value
                    w2 = Wigner3j([I1, 2, I1],[-I1,0,I1]);
                    for mI1 = -I1:I1
                        for mI2 = -I1:I1
                            colIdx = find((I == I1)&(mI==mI1));
                            rowIdx = find((I == I1)&(mI==mI2));
        
                            if (-mI1+mI2+p(ip))~=0; continue; end
                            x = (-1)^(I1-mI1) * Wigner3j([I1, 2, I1],[-mI1, p(ip), mI2]) / w2;
        
                            X(calcCount) = x;
                            iRow{calcCount} = rowIdx;
                            iCol{calcCount} = colIdx;
                            calcCount = calcCount+1;
                        end
                    end
                end
                T{ip} = obj.makeSparseMatrix(X, iCol, iRow,matrixSize);
            end
        end
        function H = makeHyperfine(obj, opts)
            arguments
                obj
                opts.useRigidRotor = 1
                opts.useSpinSpinScalar = 1
                opts.useSpinSpinTensor = 1
                opts.useSpinRotation = 1
                opts.useNuclearElectric= 1
            end
            h_hyperfine = sparse(obj.Basis.NStates,obj.Basis.NStates);
            if opts.useRigidRotor
                h_hyperfine = h_hyperfine + obj.hyperfine.rigidRotor;
            end
            if opts.useSpinSpinScalar
                h_hyperfine = h_hyperfine + obj.hyperfine.spinSpinScalar;
            end
            if opts.useSpinSpinTensor
                h_hyperfine = h_hyperfine + obj.hyperfine.spinSpinTensor;
            end            
            if opts.useSpinRotation
                h_hyperfine = h_hyperfine + obj.hyperfine.spinRotation;
            end
            if opts.useNuclearElectric
                h_hyperfine = h_hyperfine + obj.hyperfine.electricQuadrupole;
            end
            
            h_hyperfine = Constants.h *h_hyperfine/Constants.h;
            H = (h_hyperfine'+h_hyperfine)/2;  %make sure that it is hermitian
        end
        function H = makeZeeman(obj,basis,AngMom, opts)
            arguments
                obj
                basis
                AngMom
                opts.Bdir = [0,0,1]; %electric field direction
            end
            Bdir = opts.Bdir/norm(opts.Bdir);   
            [Jx,Jy,Jz,~,~] = basis.AngMomOperators(AngMom);
            % H = round( Bdir(1)*Jx + Bdir(2)*Jy + Bdir(3)*Jz, 12);
            H = Bdir(1)*Jx + Bdir(2)*Jy + Bdir(3)*Jz;
            H = (H'+H)/2; %make sure that it is hermitian
        end
        function H = makeDCStark(obj,d0, N,mN, opts)
            arguments
                obj
                d0
                % basis
                N
                mN
                opts.dir = [0,0,1]; %electric field direction
            end
            Edir = opts.dir/norm(opts.dir);

            if isempty(obj.dipoleOperator)
                tC = tensorC(N,mN,1);
                obj.dipoleOperator = {tC{3},tC{2},tC{1}}; %tensorC returns sigma_plus as tC{1}
                % obj.dipoleOperator = {}
            end
            dipOp = obj.dipoleOperator; 
            cart = SphericalToCartesian(dipOp);
            [dX,dY,dZ] = cart{:};
            % H = -d0*round( Edir(1)*dX + Edir(2)*dY + Edir(3)*dZ ,12);
            H = -d0* (Edir(1)*dX + Edir(2)*dY + Edir(3)*dZ);
            H = (H'+H)/2;
        end
        function H = makeACStark(obj,a0,a2, N,mN, opts)
            % follows 10.1103/PhysRevResearch.2.013251  
            % only for 1Sigma molecules as of now
            % Requires multiplication by Intensity, not electric field
            arguments
                obj
                a0
                a2
                N
                mN
                opts.Field struct = struct(dir=[1,0,0], pol=[0,0,1]); % should contain polarisation (Jones vector) and direction
            end
            e0 = 8.8541878128e-12;
            c = 299792458;
            n = length(N);
            Field = opts.Field;
            pol = reshape(Field.pol/norm(Field.pol),[],1); %make sure they're vectors
            dir = reshape(Field.dir/norm(Field.dir),[],1);

            if ~(dot(pol,dir)==0)
                error("Polarisation should be perpendicular to propagation direction"); 
            end
            
            %First find the matrix to rotate the beam/polarisation into the reference axis
            rotAx = cross(dir,obj.ReferenceAxis); %rotation axis
            angle = acos(dot(dir,obj.ReferenceAxis)); %rotation angle in radians
            polRot = rotvec2mat3d(rotAx*angle)*pol;

            %The isotropic part is diagonal
            H0 = speye(n,n);

            %The anisotropic is a bit more involved, but very similar to DC stark
            %The constants in front of the polarisation and polarisability
            %tensors follow the convention from Blackmore et al. https://doi.org/10.1103/PhysRevA.102.053316

            % A0 = tensorC(N,mN, 0); %Is just diagonal
            % A1 = tensorC(N,mN, 1);
            if isempty(obj.electricGradient) %calculated multiple times
                A2 = tensorC(N,mN, 2); %second order tensor, for rotational states
                obj.electricGradient = A2;
            else
                A2 = obj.electricGradient;
            end
            [P0,P1,P2] = coupleCartesianSpherically(polRot,conj(polRot)); % Polarisation Tensor
            P0 = cellfun(@(x) -sqrt(3)*x,P0,'un',0); 
            P1 = cellfun(@(x) -sqrt(2)*x,P1,'un',0); %eugh matlab, sometimes....
            P2 = cellfun(@(x) sqrt(3/2)*x,P2,'un',0);
            % A2 = A2;

            % H0 = sphericalTensorDot(A0,P0);
            % H1 = sphericalTensorDot(A1,P1);
            H2 = sphericalTensorDot(A2,P2);
            H = -(2/(e0*c))*(a0*H0 + a2*H2)/4; %+a1*H1;
            H = (H'+H)/2; %make sure that it is hermitian
        end
        function s = makeSparseMatrix(~,X, iCol, iRow, matrixSize)
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
            s = sparse(iRows,iCols,Xm,matrixSize,matrixSize);
        end
    end
end

