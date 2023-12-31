classdef NaCs < Molecule
    %NaCs 
%   Parameters from Aldegunde & Hutson (2017): 10.1103/PhysRevA.96.042506
    properties (Constant)
        d0  = 4.75 * 3.33564e-30; %Dipole Moment (V/m), 1Debye = 3.33564e-30 C*m 
    end
    
    methods
        function obj = NaCs(nNa,nCs)
            arguments
                nNa {mustBeMember(nNa,[23])} = 23;
                nCs {mustBeMember(nCs,[133])} = 133;
            end
            %RbCs Construct an instance of this class            
            if nCs==133&&nNa==23
                obj.Atom1=Na23;
                obj.Atom2=Cs133;
                idx=1;
            end
            eps0 = 8.8541878128e-12;
            bohr = 5.29177210903e-11;
            h = 6.62607015e-34; 
            %   Constants taken from DOI: 
            c1 = [14.2] * h;%for Na
            c2 = [854.5] * h;%for CS
            c3 = [105.6] * h;
            c4 = [3941.8] * h;
            gr = [0.0062];
            Brot = 1e9*[1.7396] * h; 
            Drot = [0] * h; 
            Q1 = 1e6*[-0.097]*h; %nuclear electric quadrupole constant (at position of Atom 1) 
            Q2 = 1e6*[0.15]*h;%nuclear electric quadrupole constant (at position of Atom 2) 
            a0 = [935.399]*4*pi*eps0*bohr^3; 
            a2 = [936.721]*4*pi*eps0*bohr^3;

            obj.c1 = c1(idx);
            obj.c2 = c2(idx);            
            obj.c3 = c3(idx);
            obj.c4 = c4(idx);
            obj.gr = gr(idx);
            obj.Brot = Brot(idx);
            obj.Drot = Drot(idx);
            obj.Q1 = Q1(idx);
            obj.Q2 = Q2(idx);
            obj.a0 = a0;
            obj.a2 = a2;

            obj.name = obj.Atom1.name + obj.Atom2.name;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

