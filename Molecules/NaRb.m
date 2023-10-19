classdef NaRb < Molecule
%NaRb
% Parameters from Guo et al., PRA 97, 020501(R) (2018)
% https://journals.aps.org/pra/abstract/10.1103/PhysRevA.97.020501
% and from Aldegunde et al., PRA 96, 042506 (2017)
% https://journals.aps.org/pra/abstract/10.1103/PhysRevA.96.042506
    properties (Constant)
        d0  = 3.2 * 3.33564e-30; %Dipole Moment (V/m), 1Debye = 3.33564e-30 C*m
    end
    
    methods
        function obj = NaRb(nNa,nRb)
            arguments
                nNa {mustBeMember(nNa,[23])} = 23;
                nRb {mustBeMember(nRb,[87])} = 87;
            end
            %RbCs Construct an instance of this class            
            if nRb==87&&nNa==23
                obj.Atom1=Na23;
                obj.Atom2=Rb87;
                idx=1;
            end
            %   Constants taken from DOI: 10.1038/s41567-021-01328-7
            h = 6.62606896e-34; 
            c1 = [60.7] * h;%for Rb
            c2 = [983.8] * h;%for Cs
            c3 = [259.3] * h;
            c4 = [6.56e3] * h;
            gr = [0.001];
            Brot = 1e9*[2.0896628] * h; 
            Drot = [0] * h; 
            Q1 = 1e6*[-0.139]*h; %nuclear electric quadrupole constant (at position of Atom 1) 
            Q2 = 1e6*[-3.048]*h;%nuclear electric quadrupole constant (at position of Atom 2) 
            
            obj.c1 = c1(idx);
            obj.c2 = c2(idx);            
            obj.c3 = c3(idx);
            obj.c4 = c4(idx);
            obj.gr = gr(idx);
            obj.Brot = Brot(idx);
            obj.Drot = Drot(idx);
            obj.Q1 = Q1(idx);
            obj.Q2 = Q2(idx);

            obj.name = obj.Atom1.name + obj.Atom2.name;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

