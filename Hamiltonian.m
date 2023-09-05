classdef Hamiltonian
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
    end
    
    methods
        function obj = Hamiltonian(opts)
            arguments
                opts.Molecule = []
                opts.maxN = 1;
            end
            %HAMILTONIAN Construct an instance of this class
            %   Detailed explanation goes here
            if ~isempty(opts.Molecule)
                C = Constants;
                Mol = opts.Molecule;
                obj.Molecule=Mol;
                Atom1 = Mol.Atom1;
                Atom2 = Mol.Atom2;
                i1 = AngMom(Atom1.spin,"i1");
                i2 = AngMom(Atom2.spin,"i2");
                N = AngMom([0:opts.maxN],"N");
                I = couple(i1,i2,"I");      %couple the two nuclear momenta
                Fi1 = couple(i1,N,"Fi1");   %couple only a single nucleus to the rotation
                Fi2 = couple(i2,N,"Fi2");   %couple only a single nucleus to the rotation
                F = couple(I,N,"F");        %fully coupled
                UCBasis = Basis(N,i1,i2);
                SCBasis = Basis(N,I); 
                Fi1Basis = Basis(Fi1,i2);
                Fi2Basis = Basis(Fi2,i1);
                FCBasis = Basis(F);
                nStates = UCBasis.NStates;
                SCtoUC = calcTransform(SCBasis,UCBasis,I); %gives you U such that UCBasis = U*SCBasis (UC is a sparse matrix) (unitarity is checked)
                Fi1toUC = calcTransform(Fi1Basis,UCBasis,Fi1);
                Fi2toUC = calcTransform(Fi2Basis,UCBasis,Fi2);
                %% Set standard Fields
                Fields.B = struct(value=0,dir=[0,0,1]);
                Fields.E = struct(value=0,dir=[0,0,1]);
                Fields.I = struct(value=0,dir=[0,0,1], pol = [0,1,0]);

                %% Zeeman term
                g1 = Atom1.gFactor; sig1 = Atom1.nuclearShielding;
                g2 = Atom2.gFactor; sig2 = Atom2.nuclearShielding;
                st = UCBasis.getStates("all");
                obj.zeeman = obj.makeZeeman(g1*(1-sig1)*C.muN,st.mi1) + obj.makeZeeman(g2*(1-sig2)*C.muN,st.mi2) + obj.makeZeeman(Mol.gr*C.muN,st.mN); 
                %% DC Stark term
                obj.dc_stark = obj.makeDCStark(Mol.d0,st.N,st.mN);
                %% Hyperfine
                %Rigid Rotor
                obj.hyperfine.rigidRotor = obj.rotational(st.N,Mol.Brot,0); 

                %Scalar Spin-spin coupling
                st = SCBasis.getStates("all");
                spinSpinScalar = obj.scalarNuclear(Mol.c4,st.I,st.i2, st.i1);
                obj.hyperfine.spinSpinScalar =  SCtoUC'*spinSpinScalar*SCtoUC; 

                %Scalar Spin-Rotation Coupling & Nuclear Electric Quadrupole
                st = Fi1Basis.getStates("all");
                h_spinRot1 = obj.scalarNuclear(Mol.c1,st.Fi1,st.N, st.i1);
                h_electricQuadrupole1 = obj.electricQuadrupole(Mol.Q1,st.N,st.i1,st.Fi1);

                st = Fi2Basis.getStates("all");
                h_spinRot2 = obj.scalarNuclear(Mol.c2,st.Fi2,st.N, st.i2);
                h_electricQuadrupole2 = obj.electricQuadrupole(Mol.Q2,st.N,st.i2,st.Fi2);

                obj.hyperfine.spinRotation = Fi1toUC'*h_spinRot1*Fi1toUC + Fi2toUC'*h_spinRot2*Fi2toUC;
                obj.hyperfine.electricQuadrupole = Fi1toUC'*h_electricQuadrupole1*Fi1toUC + Fi2toUC'*h_electricQuadrupole2*Fi2toUC;

                %Tensor Spin-Spin Coupling
                obj.hyperfine.spinSpinTensor = obj.tensorNuclear(Mol.c3,i1,i2, UCBasis); 
                obj.Basis = UCBasis;
            end
        end
        function H = makeHamiltonian(obj, opts)
            arguments
                obj
                opts.useRigidRotor = 1
                opts.useSpinSpinScalar = 1
                opts.useSpinSpinTensor = 1
                opts.useSpinRotation = 1
                opts.useNuclearElectric= 1
                opts.B struct = struct(value=[0],dir=[0,0,1]) %Magnetic field in Tesla, maybe should just be a vector?
                opts.E struct = struct(value=[0],dir=[0,0,1]) %Electric field in V/m
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
            DC_stark = squeeze(opts.E.dir(1)*obj.dc_stark(1,:,:)+opts.E.dir(2)*obj.dc_stark(2,:,:)+opts.E.dir(3)*obj.dc_stark(3,:,:));
            H = h_hyperfine + opts.B.value*obj.zeeman + opts.E.value*DC_stark;           
        end
        function [H,s] = scalarNuclear(~,c,jtot,j1,j2)
            n = length(jtot);
            s = 0.5*c*(jtot.*(jtot+1) - j1.*(j1+1) - j2.*(j2+1));
            H = spdiags(s,0,n,n); %create diagonal sparse matrix
        end
        function H = tensorNuclear(~,c,j1,j2, basis)           
            assert(length(j1.js)==1 && length(j2.js)==1); %Only checked for single vectors, do some thorough checking of basis.AngMomOperators first if you want to extend functionality!
            I1 = j1.js;
            I2 = j2.js;
            N = basis.getStates("N");
            mN = basis.getStates("mN");
            
            T2C = tensorC(N,mN, 2); 
            T2_Coupled = basis.coupleSpherically(j1,j2);
            
            H = sqrt(6)*c*sphericalTensorDot(T2C,T2_Coupled);
        end
        function [H,r] = rotational(~, N, Brot, Drot)
            n = length(N);
            N2 = N.*(N+1);
            r = Brot*N2 - Drot*N2.*N2;
            % r = Brot*N.^2 - Drot*(N.^2).*(N.^2);
            H = spdiags(r ,0,n,n); %create diagonal sparse matrix
        end    
        % function [H,q] = electricQuadrupole(obj,Q,N,I,Fi)
        %     n = length(N);            
        %     q = Q*(3/4)*( (I.^2 + N.^2 - Fi.^2).*((I.^2 + N.^2 - Fi.^2)-1) - I.*(I+1).*N.*(N+1)) ./ ...
        %         (2*I.*(2*I-1).*(2*N-1).*(2*N+3));
        %     H = spdiags(q,0,n,n); %create diagonal sparse matrix
        % end 
        function [H,q] = electricQuadrupole(~,Q,N,I,Fi) 
            n = length(N);            
            NdotI = 0.5*(Fi.*(Fi+1) - N.*(N+1) - I.*(I+1));            
            q = -Q*( 3*(NdotI.^2) + 1.5*NdotI - I.*(I+1).*N.*(N+1)) ./ ( 2*I.*(2*I-1).*(2*N-1).*(2*N+3) );
            H = spdiags(q,0,n,n); %create diagonal sparse matrix
        end 
        function [H,z] = makeZeeman(~,mu,m)
            n = length(m);
            z = -mu*m;
            H = spdiags(z,0,n,n); %only works if bfield is in z-direction
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
            dir = opts.dir/norm(opts.dir);
           
            TC = tensorC(N,mN, 1); 
            dcZ = TC{2};
            hm = TC{1}; 
            % hp = TC{3};
            
            % dcX = -(hm-hp)/sqrt(2);%see Tills code
            % dcY = 1i*(hm-hp)/sqrt(2);
            dcX = -(1/sqrt(2))*(hm + conj(hm'));
            dcY = 1i*(1/sqrt(2))*(hm - conj(hm'));

            H = nan(3,size(dcY,1),size(dcY,2));
            H(1,:,:) = d0*dcX; 
            H(2,:,:) = d0*dcY; 
            H(3,:,:) = d0*dcZ; 
        end
        function s = makeSparseMatrix(X, iCol, iRow)
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
            s = sparse(iRows,iCols,Xm);
        end
    end
end

