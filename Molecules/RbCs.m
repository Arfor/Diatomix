classdef RbCs < Molecule
    %RbCs 
%  Most recent Rb87Cs133 Constants are given in the supplementary 
% of Gregory et al., Nat. Phys. 17, 1149-1153 (2021)
% https://www.nature.com/articles/s41567-021-01328-7
%  Polarisabilities are for 1064 nm reported 
% in Blackmore et al., PRA 102, 053316 (2020)
% https://journals.aps.org/pra/abstract/10.1103/PhysRevA.102.053316    
    properties (Constant)
        d0  = 1.225 * 3.33564e-30; %Dipole Moment (V/m), 1Debye = 3.33564e-30 C*m  from Ni et al., Science 322, 231-235 (2008) 0.574
    end
    
    methods
        function obj = RbCs(nRb,nCs)
            arguments
                nRb {mustBeMember(nRb,[87])} = 87;
                nCs {mustBeMember(nCs,[133])} = 133;
            end
            %RbCs Construct an instance of this class            
            if nCs==133&&nRb==87
                obj.Atom1=Rb87;
                obj.Atom2=Cs133;
                idx=1;
            end
            eps0 = 8.8541878128e-12;
            bohr = 5.29177210903e-11;
            h = 6.62607015e-34; 
            %   Constants taken from DOI: 10.1038/s41567-021-01328-7
            c1 = [98.4] * h;%for Rb
            c2 = [194.2] * h;%for Cs
            c3 = [192.4] * h;
            c4 = [19.0189557e3] * h;
            gr = [0.0062];
            Brot = 1e6*[490.173994326310] * h; 
            Drot = [207.3] * h; 
            Q1 = 1e3*[-809.29]*h; %nuclear electric quadrupole constant (at position of Atom 1) 
            Q2 = 1e3*[59.98]*h;%nuclear electric quadrupole constant (at position of Atom 2) 
            a0 = [2020]*4*pi*eps0*bohr^3; 
            a2 = [1997]*4*pi*eps0*bohr^3;
            
            obj.c1 = c1(idx);
            obj.c2 = c2(idx);            
            obj.c3 = c3(idx);
            obj.c4 = c4(idx);
            obj.gr = gr(idx);
            obj.Brot = Brot(idx);
            obj.Drot = Drot(idx);
            obj.Q1 = Q1(idx);
            obj.Q2 = Q2(idx);
            obj.a0 = a0(idx);
            obj.a2 = a2(idx);

            obj.name = obj.Atom1.name + obj.Atom2.name;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

