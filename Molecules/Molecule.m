classdef Molecule
    %Molecule SuperClass
    
    properties 
        name    ="Molecule";%Name Molecule
        Atom1   ="Atom1";   %Name Atom1
        Atom2   ="Atom2";   %Name Atom2
        c1      =0;         %for Atom 1
        c2      =0;         %for Atom 2
        c3      =0;         %
        c4      =0;         %
        gr      =0;         %Molecular Magnetic g-factor
        Brot    =0;         %Rotational constant
        Drot    =0;         %Centrifugal constant
        Q1      =0;         %nuclear electric quadrupole constant (at position of Atom1) 
        Q2      =0;         %nuclear electric quadrupole constant (at position of Atom2) 
        a0      =0;         %h*Hz/(W/cm^2) at 1064nm
        a2      =0;         % tensor polarisability at 1064nm
    end
    
    methods
        function obj = Molecule(inputArg1,inputArg2)
            %MOLECULE Construct an instance of this class
            %   Detailed explanation goes here
            % obj.Property1 = inputArg1 + inputArg2;
        end
        
        function info = infoString(obj)
            %info = moleculeInfoString(obj)
            %returns a string with all info of the molecule
            ps = properties(obj);
            info = [];
            for k = 1:length(ps)
                p = ps{k};
                if any(strcmp(p,{'Atom1','Atom2'}))
                    try
                        info = [info; sprintf("%s = %s",p,obj.(p).name)];
                    catch
                        info = [info; sprintf("%s = %s",p,p)];
                    end
                elseif strcmp(p,'name')
                    info = [info; sprintf("%s = %s",p,obj.(p))];
                else
                    info = [info; sprintf("%s = %.4g",p,obj.(p))];
                end
            end
        end
    end
end