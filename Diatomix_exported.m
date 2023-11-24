classdef Diatomix_exported < matlab.ui.componentcontainer.ComponentContainer

    % Properties that correspond to underlying components
    properties (Access = private, Transient, NonCopyable)
        GridLayout                  matlab.ui.container.GridLayout
        GridLayout2                 matlab.ui.container.GridLayout
        GridLayout3                 matlab.ui.container.GridLayout
        Button3                     matlab.ui.control.Button
        LoadButton                  matlab.ui.control.Button
        SaveButton                  matlab.ui.control.Button
        ToWorkspaceButton           matlab.ui.control.Button
        statusLabel                 matlab.ui.control.Label
        statusLamp                  matlab.ui.control.Lamp
        Label                       matlab.ui.control.Label
        AutoUpdateButton            matlab.ui.control.StateButton
        CalculateButton             matlab.ui.control.Button
        TabGroup                    matlab.ui.container.TabGroup
        FieldsTab                   matlab.ui.container.Tab
        GridLayout4                 matlab.ui.container.GridLayout
        PolarisationLabel           matlab.ui.control.Label
        DirectionLabel_2            matlab.ui.control.Label
        FieldAmplitudeLabel         matlab.ui.control.Label
        Laser1PolField              matlab.ui.control.EditField
        BLabel                      matlab.ui.control.Label
        BField                      matlab.ui.control.EditField
        BDirField                   matlab.ui.control.EditField
        Laser1PolPhaseField         matlab.ui.control.EditField
        aEditFieldLabel             matlab.ui.control.Label
        Laser1PolAngleField         matlab.ui.control.EditField
        degreeLabel                 matlab.ui.control.Label
        ILabel                      matlab.ui.control.Label
        ELabel                      matlab.ui.control.Label
        Laser1DirField              matlab.ui.control.EditField
        EDirField                   matlab.ui.control.EditField
        EField                      matlab.ui.control.EditField
        Laser1Field                 matlab.ui.control.EditField
        MoleculeTab                 matlab.ui.container.Tab
        GridLayout7                 matlab.ui.container.GridLayout
        MoleculeInfo                matlab.ui.control.TextArea
        MoleculeEditField           matlab.ui.control.EditField
        MoleculeLabel               matlab.ui.control.Label
        Tree                        matlab.ui.container.CheckBoxTree
        RigidRotorNode              matlab.ui.container.TreeNode
        SpinSpinScalarNode          matlab.ui.container.TreeNode
        SpinSpinTensorNode          matlab.ui.container.TreeNode
        SpinRotationNode            matlab.ui.container.TreeNode
        NuclearElectricNode         matlab.ui.container.TreeNode
        BasisTab                    matlab.ui.container.Tab
        GridLayout9                 matlab.ui.container.GridLayout
        BasisChoice                 matlab.ui.control.DropDown
        BasisDropDownLabel          matlab.ui.control.Label
        NmaxLabel                   matlab.ui.control.Label
        NmaxEditField               matlab.ui.control.NumericEditField
        PlotTab                     matlab.ui.container.Tab
        GridLayout5                 matlab.ui.container.GridLayout
        LabelsCheckBox              matlab.ui.control.CheckBox
        N2EditField                 matlab.ui.control.EditField
        N2EditFieldLabel            matlab.ui.control.Label
        N1EditField                 matlab.ui.control.EditField
        N1EditFieldLabel            matlab.ui.control.Label
        Y2DropDown                  matlab.ui.control.DropDown
        Y2DropDownLabel             matlab.ui.control.Label
        Y1DropDown                  matlab.ui.control.DropDown
        Y1DropDownLabel             matlab.ui.control.Label
        ShowPointsCheckBox          matlab.ui.control.CheckBox
        PlotAdiabatsCheckBox        matlab.ui.control.CheckBox
        TDMTab                      matlab.ui.container.Tab
        GridLayout8                 matlab.ui.container.GridLayout
        PlotTDMLabel                matlab.ui.control.Label
        plotTDMCheckbox             matlab.ui.control.CheckBox
        TDMVarValue                 matlab.ui.control.NumericEditField
        TDMcurrVarLabel             matlab.ui.control.Label
        TDMTable                    matlab.ui.control.Table
        TDMInitialStateEditField    matlab.ui.control.NumericEditField
        InitialStateEditFieldLabel  matlab.ui.control.Label
        UITable                     matlab.ui.control.Table
        UIAxes                      matlab.ui.control.UIAxes
    end

    properties (Access = private)
        Fields % Container for all fields
        Ham % Hamiltonian
        Molecule
    end
    
    methods (Access = private) %Most calculations done with these functions
        function updatePlot(comp)
            setStatus(comp, 'yellow', 'Update Plot');
            ud = comp.UserData;
            ax = comp.UIAxes;

            delete(findobj(ax,'Type','line'))
            delete(findobj(ax,'Type','ConstantLine'))
            delete(findobj(ax,'Type','Scatter'))
            delete(findobj(ax,'Type','Text'))
            delete(findobj(ax,'Tag','updatePlot'))
            % spec = ud.spectrum(:,ud.selectStates);
            y1 = ud.yVar1.value;
            if sum(ud.selectStates)<size(y1,2)
                y1 = y1(:,ud.selectStates);
            end
            %Part for the first axis
            % yyaxis(ax,'left')
            % ax.YAxis(1).Color = "k";
            switch string(ud.yVar1.name)
                case "E"
                    set(ax, 'YDir','normal')
                case "E-<E>"
                    y1 = y1 - mean(y1,2);
                    set(ax, 'YDir','normal')
                case "µ"
                    set(ax, 'YDir','reverse')
                case "d"
                    set(ax, 'YDir','normal')
            end
            if ~isempty(y1)
                plot(ax,ud.xVar.value,y1,'-', "Color",[0,0,0,0.4], LineWidth=1, Tag='Spec');
            end
            % for p=1:length(ps)
            %     ps(p).DataTipTemplate(1).DataTipRows(2).Label="E(MHz)=";
            % %     ps(p).DataTipTemplate(1).DataTipRows(3).Label="E(MHz)=";
            % end
            xlabel(ax,sprintf("%s (%s)",ud.xVar.name,ud.xVar.unit))
            ylabel(ax,sprintf("%s (%s)",ud.yVar1.name,ud.yVar1.unit))
            title(ax,comp.Molecule.name, Interpreter="tex")
            xlim(ax,[min(ud.xVar.value),max(ud.xVar.value)]);
            ylim(ax,'auto')
            ud.currPlotXlim = ax.XLim;
            ud.currPlotYlim = ax.YLim;
            legend(ax,'off');    

            % turn off HitTest for all children of plot to be able to click them
            for k = 1:length(ax.Children)
                ax.Children(k).HitTest = 'off';
            end
            comp.UserData = ud;
            ShowPointsCheckBoxValueChanged(comp);
            LabelsCheckBoxValueChanged(comp);
            plotTDMCheckboxValueChanged(comp);
            updateTDMTable(comp);
            setStatus(comp, 'green', '');
        end

        function calculateSpectrum(comp)
            C = Constants;
            ud = comp.UserData;
            F = comp.Fields;
            h = comp.Ham;
            nStates = ud.NStates;
            energies = nan(length(ud.xVar.value),nStates);
            states = nan(length(ud.xVar.value),nStates,nStates);
            % E.dir = F.E.dir;
            % B.dir = F.B.dir;
            E = F.E.value(1)*F.E.scale;
            B = F.B.value(1)*F.B.scale;
            I = F.I.value(1)*F.I.scale;
            tmp = namedargs2cell(comp.UserData.hamOpts); %because eugh, matlab...
            H0 = h.makeHyperfine(tmp{:});
            assert(ishermitian(H0 + h.zeeman + h.dc_stark + h.ac_stark));
            for k = 1:length(ud.xVar.value)
                setStatus(comp, 'yellow', sprintf('Calc. %.0f/%.0f',k,length(ud.xVar.value)));
                if strcmp(ud.xVar.name,"B")
                    B = F.B.value(k)*F.B.scale;
                    H = H0 + h.zeeman*B + h.dc_stark*E + h.ac_stark*I;
                elseif strcmp(ud.xVar.name,"E")
                    E = F.E.value(k)*F.E.scale;
                    H = H0 + h.zeeman*B + h.dc_stark*E + h.ac_stark*I;
                elseif strcmp(ud.xVar.name,"I")
                    I = F.I.value(k)*F.I.scale;
                    H = H0 + h.zeeman*B + h.dc_stark*E + h.ac_stark*I;
                elseif strcmp(ud.xVar.name,"Pol. Angle")
                    %remake ACstark for every polarisation
                    Mol = h.Molecule;
                    st = h.Basis.getStates('all');
                    b = ud.xVar.value(k)*pi/180;
                    a = F.I.polPhase(1)*pi/180;%needed to evaluate Laser1PolField
                    F.I.pol = round(eval(comp.Laser1PolField.Value),6); %makes calculation a little bit easier
                    H_ac = h.makeACStark(Mol.a0,Mol.a2,st.N,st.mN,Field=F.I);
                    H = H0 + h.zeeman*B + h.dc_stark*E + H_ac*I;
                elseif strcmp(ud.xVar.name,"Pol. Phase")
                    Mol = h.Molecule;
                    st = h.Basis.getStates('all');
                    b = F.I.polAngle(1)*pi/180;
                    a = ud.xVar.value(k)*pi/180;%needed to evaluate Laser1PolField
                    F.I.pol = round(eval(comp.Laser1PolField.Value),6); %makes calculation a little bit easier
                    H_ac = h.makeACStark(Mol.a0,Mol.a2,st.N,st.mN,Field=F.I);
                    H = H0 + h.zeeman*B + h.dc_stark*E + H_ac*I;
                else
                end
                [V,D] = eig(full(H/C.h),'vector'); %no need to use sparse matrices for matrices smaller than 1000x1000. Divide by Plancks constant to get energy in Hz
                [d,sortIdx] = sort(real(D));
                energies(k,:) = d;        
                states(k,:,:) = V(:,sortIdx);
            end
            ud.spectrum = energies*1e-6; %energies in MHz (xVar,States)
            ud.states = states;
            ud.currStates = states;
            ud.stateLabels = repmat(1:size(ud.spectrum,2),size(ud.spectrum,1),1);
            ud.currLabels = ud.stateLabels; 

            ud.diabatStates=[];
            ud.muEffective = [];
            ud.inducedDipole = [];
            ud.TDM = [];
            comp.UserData = ud;

            updateVars(comp);
            setStatus(comp, 'green', '');
        end
        function updateHamiltonian(comp)
            ud = comp.UserData;
        
            Hamil = Hamiltonian(Molecule = comp.Molecule, Fields=comp.Fields, maxN=ud.Nmax);
            comp.Ham = Hamil; 
            UCBasis = Hamil.Basis;
            ud.NStates = UCBasis.NStates;
            statesUCBasis = UCBasis.getStates('all');
            stateNames = ["N","mN","mi1","mi2"];
            ud.statesUCBasis = statesUCBasis(:,stateNames);
            comp.UITable.ColumnName = ["Comp",stateNames];
            ud.selectStates = ismember(ud.statesUCBasis.N,ud.selectN);
            
            ud.Bases.FC = [];
            ud.Bases.SC = [];
            ud.diabatStates=[];
            ud.muEffective = [];
            ud.inducedDipole = [];
            ud.clickedState = [];
            comp.UserData = ud;
        end
        function sortStates(comp)
            % sortStates(comp) calculates all diabats and returns sorted states
            ud = comp.UserData;
            if isempty(ud.diabatStates) %only do it if it hasnt been done before for this hamiltonian
                y = ud.spectrum;
                states = ud.states;
                x = ud.xVar.value;  
                diabatStates = nan(size(states));
                diabatSorting = nan(size(y), "like",1);
                for stIdx = 1:size(states,3)
                    [~, diabatState, ~, ~, diabatSort] = findAdiabat(x,y,states, length(x),stIdx); 
                    diabatStates(:,:,stIdx) = diabatState;
                    diabatSorting(:,stIdx) = diabatSort;
                end
                ud.diabatStates = diabatStates;
                ud.diabatSorting = diabatSorting;
            end
            comp.UserData = ud;
        end
        function mu = calculateMuEffective(comp)
            sortStates(comp); %first sort states as adiabat
            ud = comp.UserData;
            Hz = comp.Ham.zeeman;
            mu = nan(size(ud.spectrum));
            psis = ud.diabatStates;
            for x = 1:length(ud.xVar.value)
                psi = squeeze(psis(x,:,:));
                mu(x,:) = diag(psi'*Hz*psi); %A'=ctranspose(A), A.' = transpose(A)
            end
            ud.muEffective = -real(mu)/Constants.muN;
            comp.UserData = ud;
        end
        function d = calculateInducedDipole(comp)
            sortStates(comp); %first sort states as diabats
            ud = comp.UserData;
            H_dc = comp.Ham.dc_stark/Constants.D;
            d = zeros(size(ud.spectrum));
            psis = ud.diabatStates;
            for xidx = 1:length(ud.xVar.value)
                psi = squeeze(psis(xidx,:,:)); %A'=ctranspose(A), A.' = transpose(A)
                d(xidx,:) = diag( psi'*(H_dc)*psi ); %expectation value of Stark part of the Hamiltonian
            end

            ud.inducedDipole = -real(d);
            comp.UserData = ud;
        end
        function TDM = calculateTDM(comp)
            C = Constants;
            Mol = comp.Molecule;
            setStatus(comp, 'yellow', 'Update TDM');
            ud = comp.UserData;
            h = comp.Ham;
            initState = squeeze(ud.states(:,:,ud.TDMInitialState)); %always use the statelabels from the energy-spectrum
            TDM = nan(length(ud.xVar.value),ud.NStates,3);
            for x = 1:length(ud.xVar.value)
                for p = 1:3
                    TDM(x,:,p) = Mol.d0/C.D*abs( conj(squeeze(initState(x,:))) * h.dipoleOperator{p} * (squeeze(ud.states(x,:,:))) ); %in Debye, should be real() instead of abs(), don't know the reason why
                end
            end
            ud.TDM = TDM;
            comp.UserData = ud;
            setStatus(comp, 'green', '');
        end

        function updateVars(comp)
            setStatus(comp, 'yellow', sprintf('Update Vars'));
            comp.UserData.updateSuccess = 0;
            updateX(comp);            
            updateY(comp);
            comp.UserData.updateSuccess = 1;
        end
        function updateY(comp)
           ud = comp.UserData;
           C=Constants;
           x = ud.xVar.value * ud.xVar.scale;
           switch string(comp.Y1DropDown.Value)
                case "E"
                    yVar1.name = "E";
                    yVar1.unit = "MHz";
                    yVar1.value = ud.spectrum;
                    comp.UserData.currStates = comp.UserData.states;
                    comp.UserData.currLabels = comp.UserData.stateLabels;
                case "E - <E>"
                    yVar1.name = "E - <E>";
                    yVar1.unit = "MHz";
                    y1 = ud.spectrum;
                    yVar1.value = y1 - mean(comp.UserData.spectrum(:,ud.selectStates),2);
                    comp.UserData.currStates = comp.UserData.states;
                    comp.UserData.currLabels = comp.UserData.stateLabels;
                case "E - a0"
                    yVar1.name = "E - Ia_0/(2*\epsilon_0c)";
                    yVar1.unit = "MHz";
                    y1 = ud.spectrum;
                    yVar1.value = y1 + x'*comp.Molecule.a0/(2*C.h*C.e0*C.c)*1e-6;
                    comp.UserData.currStates = comp.UserData.states;
                    comp.UserData.currLabels = comp.UserData.stateLabels;
                case "E - GS"
                    yVar1.name = "E - GS";
                    yVar1.unit = "MHz";
                    y1 = ud.spectrum;
                    yVar1.value = y1 - ud.spectrum(:,1);
                    comp.UserData.currStates = comp.UserData.states;
                    comp.UserData.currLabels = comp.UserData.stateLabels;
                case "µ"
                    yVar1.name = "µ";
                    yVar1.unit = "µ_N";
                    if isempty(ud.muEffective)
                        calculateMuEffective(comp);
                    end
                    yVar1.value = comp.UserData.muEffective;
                    comp.UserData.currStates = comp.UserData.diabatStates;
                    comp.UserData.currLabels = comp.UserData.diabatSorting;
                case "d"
                    yVar1.name = "d";
                    yVar1.unit = "D";
                    if isempty(ud.inducedDipole)
                        calculateInducedDipole(comp);
                    end
                    yVar1.value = comp.UserData.inducedDipole;
                    comp.UserData.currStates = comp.UserData.diabatStates;
                    comp.UserData.currLabels = comp.UserData.diabatSorting;
           end
           comp.UserData.yVar1 = yVar1;
        end
        function updateX(comp)
           ud = comp.UserData;
           Field = comp.Fields;

            %update xVar
           try 
                Field.B.value = eval(comp.BField.Value);
                Field.E.value = eval(comp.EField.Value);
                Field.I.value = eval(comp.Laser1Field.Value);
                a = Field.I.polPhase; b = Field.I.polAngle; %needed to evaluate Laser1PolField
                if length(Field.B.value)>1
                    ud.xVar = Field.B;
                elseif length(Field.E.value)>1
                    ud.xVar = Field.E;
                elseif length(Field.I.value)>1
                    ud.xVar = Field.I;
                elseif length(Field.I.polPhase)>1
                    ud.xVar.value = a;
                    ud.xVar.name = "Pol. Phase";
                    ud.xVar.unit = "°";
                elseif length(Field.I.polAngle)>1
                    ud.xVar.value = b;
                    ud.xVar.name = "Pol. Angle";
                    ud.xVar.unit = "°";
                else
                    return
                end
            catch
                error("Not a valid entry for fields!")
                return
            end
            NPoints = length(Field.B.value)*length(Field.E.value)*length(Field.I.value);
            if NPoints > 302
                warning("Do you really want to calculate %d points?",NPoints)
                return
            end
            ud.Nmax = comp.NmaxEditField.Value;
            comp.Fields = Field;
            comp.TDMcurrVarLabel.Text = sprintf("%s (%s)",ud.xVar.name, ud.xVar.unit);
            % comp.TDMVarValue.Value = ud.xVar.value(end);
            comp.UserData = ud;
        end

        function bool = checkZoomed(comp)
            ax = comp.UIAxes;
            ud = comp.UserData;
            if any(ax.XLim~=ud.currPlotXlim)&&any(ax.YLim~=ud.currPlotYlim)
                ud.currPlotXlim = ax.XLim;
                ud.currPlotYlim = ax.YLim;
                bool = 1;
                comp.UserData = ud;
                LabelsCheckBoxValueChanged(comp)
            else %didnt change
                bool = 0;
            end
        end
        function setStatus(comp, color, message)
            comp.statusLamp.Color = color;
            comp.statusLabel.Text = message;
            pause(0.001);
        end
        
        function updateTDMTable(comp)
            stateIdx = comp.UserData.TDMInitialState;

            x = comp.TDMVarValue.Value;
            xVals = comp.UserData.xVar.value;
            [~,xIdx] = min(abs(x-xVals));
            comp.TDMVarValue.Value = xVals(xIdx);

            spec = comp.UserData.spectrum;
            % states = comp.UserData.currStates;
            if isempty(comp.UserData.TDM) %only recalculate if empty
                calculateTDM(comp);
            end
            TDM = comp.UserData.TDM;

            tdm = round(squeeze(TDM(xIdx,:,:)),4);
            Ediff = round(reshape(spec(xIdx,:) - spec(xIdx,stateIdx),[],1), 5);
            [~,sortIdx] = sort(max(abs(tdm),[],2),"descend");
            
            E = table(num2str(Ediff(sortIdx),'%.3f'), VariableNames="E (MHz)");
            tdmTable = [E, round(array2table([tdm(sortIdx,:),sortIdx],"VariableNames",["s-","p","s+","StateIdx"]),4)];
            comp.TDMTable.Data = tdmTable(1:min(50,comp.UserData.NStates),:);
        end
        
        function basis = getBasis(comp) %calculates basis if empty, otherwise just returns the right basis
            UCBasis = comp.Ham.Basis;
            i1 = UCBasis.momenta.i1;
            i2 = UCBasis.momenta.i2;
            N = UCBasis.momenta.N;
            switch string(comp.UserData.BasisChoice)
                case "Uncoupled"
                    comp.UserData.qnumbers = ["N","mN","mi1","mi2"];
                    basis = UCBasis;
                    basis.transforms.toUC = speye(basis.NStates);
                case "Spin Coupled"
                    comp.UserData.qnumbers = ["N","mN","I","mI"];    
                    if isempty(comp.UserData.Bases.SC)
                        I = couple(i1,i2,"I");      %couple the two nuclear momenta
                        SCBasis = Basis(N,I); 
                        SCtoUC = calcTransform(SCBasis,UCBasis,I); %gives you U such that UCBasis = U*SCBasis
                        SCBasis.transforms.toUC = SCtoUC';
                        comp.UserData.Bases.SC = SCBasis;
                    end
                    basis = comp.UserData.Bases.SC;
                case "Fully Coupled"
                    comp.UserData.qnumbers = ["N","I","F","mF"];
                    if isempty(comp.UserData.Bases.SC)
                        I = couple(i1,i2,"I");      %couple the two nuclear momenta
                        SCBasis = Basis(N,I); 
                        SCtoUC = calcTransform(SCBasis,UCBasis,I); %gives you U such that UCBasis = U*SCBasis
                        SCBasis.transforms.toUC = SCtoUC;
                        comp.UserData.Bases.SC = SCBasis;
                    else
                        SCBasis = comp.UserData.Bases.SC;
                        I = SCBasis.momenta.I;
                    end
                    if isempty(comp.UserData.Bases.FC)
                        F = couple(I,N,"F");        %fully coupled
                        FCBasis = Basis(F);
                        SCtoUC = SCBasis.transforms.toUC;
                        FCtoSC = calcTransform(FCBasis,SCBasis,F); %function gives you U such that SCBasis = U*FCBasis
                        FCBasis.transforms.toUC = (FCtoSC * SCtoUC); %gives you U such that FCBasis = U*UCBasis
                        comp.UserData.Bases.FC = FCBasis;
                    end
                    basis = comp.UserData.Bases.FC;
            end            
        end
        
        function updateCompTable(comp)
            clickedState = comp.UserData.clickedState;
            if ~isempty(clickedState)
                B = getBasis(comp);
                statesBasis = B.getStates('all');
                statesBasis = statesBasis(:,comp.UserData.qnumbers);
                st = B.transforms.toUC * clickedState;
                stateComp = (round(squeeze(abs(st).^2),6));
                [statesCompTable, ~] = sortrows([array2table(stateComp,"VariableNames","Comp") , statesBasis], "Comp","descend");
                comp.UITable.ColumnName = statesCompTable.Properties.VariableNames;
                comp.UITable.Data = statesCompTable(1:30,:);
            end
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function initialiseApp(comp)
                % fig = comp.Parent;
                addpath("Atoms\")
                addpath("Molecules\")
                addpath("functions\")

                try %sometimes gives error
                comp.UIAxes.InteractionOptions.BrushSupported = 'off';
                comp.UIAxes.InteractionOptions.RotateSupported = 'off';
                comp.UIAxes.InteractionOptions.PanLimitsBounded = 'on';
                end

                ud = comp.UserData;
                ud.autoUpdate = comp.AutoUpdateButton.Value;
                ud.hamChanged = 0;
                Field = comp.Fields;
                Field.B.name = "B"; Field.B.unit = "G"; Field.B.scale = 1e-4; %from Tesla to Gauss
                Field.E.name = "E"; Field.E.unit = "V/cm"; Field.E.scale = 1e2; %from V/m to V/cm
                Field.I.name = "I"; Field.I.unit = "kW/cm^{2}"; Field.I.scale = 1e7;
                comp.Fields = Field;
                Field.B.dir = eval(comp.BDirField.Value);
                Field.E.dir = eval(comp.EDirField.Value);
                Field.I.dir = eval(comp.Laser1DirField.Value);
                Field.I.polAngle = eval(comp.Laser1PolAngleField.Value);
                Field.I.polPhase = eval(comp.Laser1PolPhaseField.Value);
                a = Field.I.polPhase(1)*pi/180; b = Field.I.polAngle(1)*pi/180;
                Field.I.pol = eval(comp.Laser1PolField.Value);
                Field.B.value = eval(comp.BField.Value);
                Field.E.value = eval(comp.EField.Value);
                Field.I.value = eval(comp.Laser1Field.Value);
                ud.selectN = 0;
                ud.xVar = Field.B;
                ud.Nmax = comp.NmaxEditField.Value;
                ud.qnumbers = ["N","mN","mi1","mi2"];
                ud.BasisChoice = comp.BasisChoice.Value;
                ud.Bases.SC = [];
                ud.Bases.FC = [];
                ud.hamOpts.useRigidRotor = 1;
                ud.hamOpts.useSpinSpinScalar = 1;
                ud.hamOpts.useSpinSpinTensor = 1;
                ud.hamOpts.useSpinRotation = 1;
                ud.hamOpts.useNuclearElectric = 1;
                ud.TDMInitialState = comp.TDMInitialStateEditField.Value;
                comp.UserData = ud;
                comp.Fields = Field;

                comp.Molecule = KRb(40,87);
                comp.MoleculeInfo.Value = comp.Molecule.infoString;
                updateHamiltonian(comp);
                calculateSpectrum(comp);
                updateVars(comp);
                calculateTDM(comp);
                % collisionE = 0.1; %MHz
                updatePlot(comp);  
                comp.UserData.currPlotXlim = comp.UIAxes.XLim;
                comp.UserData.currPlotYlim = comp.UIAxes.YLim;

                drawnow;
                comp.Parent.WindowState = 'maximized';
                drawnow; % Force to draw the uifigure first
                event.IntersectionPoint = [comp.UIAxes.XLim(2),comp.UIAxes.YLim(1),0];
                UIAxesButtonDown(comp, event)
        end

        % Button down function: UIAxes
        function UIAxesButtonDown(comp, event)
            zoomed = checkZoomed(comp);
            if zoomed
                return
            end
            ud = comp.UserData;
            ax = comp.UIAxes;
            delete(findobj(ax,'Tag','updatePlot'))

            cStates = ud.currStates(:,:,ud.selectStates);
            cSpec = ud.yVar1.value(:,ud.selectStates);
            stateLabels = ud.currLabels(:,ud.selectStates);

            %snap to nearest point and highlight
            [~,xIdx] = min(abs(event.IntersectionPoint(1)-ud.xVar.value));
            [~,yIdx] = min(abs(event.IntersectionPoint(2)-cSpec(xIdx,:)));
            ax.NextPlot ="add";
            % xline(ax,ud.xVar.value(xIdx),'--',color="k",Tag="updatePlot") %not needed anymore with legend
            % yline(ax,cSpec(xIdx,yIdx),'--',color="k",Tag="updatePlot")
            st = sprintf("State %.0f\n%s = %.3g%s\n%s = %.2g%s",stateLabels(xIdx,yIdx),ud.yVar1.name,cSpec(xIdx,yIdx),ud.yVar1.unit,...
                ud.xVar.name,ud.xVar.value(xIdx),ud.xVar.unit);
            clickedPoint = scatter(ax,ud.xVar.value(xIdx), cSpec(xIdx,yIdx),50,'green','filled',...
                'o', MarkerFaceAlpha=0.5, Tag="updatePlot", DisplayName=st);
            
            %annotate
            % text(ax,min(ud.xVar.value),cSpec(xIdx,yIdx), sprintf("%.3g",cSpec(xIdx,yIdx)), HorizontalAlignment="right", Tag='updatePlot');

            %retrieve state composition of point and update table
            clickedState = reshape(cStates(xIdx,:,yIdx),[],1);
            comp.UserData.clickedState = clickedState;
            updateCompTable(comp);

            %find diabat and adiabat and highlight. Also plots legend
            if comp.PlotAdiabatsCheckBox.Value
                x = ud.xVar.value;
                if any(strcmp(ud.yVar1.name, ["µ","d"]))
                adiabatLine = plot(ax,x, cSpec(:,yIdx), '-', Color="red", LineWidth=1.5, ...
                    DisplayName="Diabat",Tag="updatePlot");   
                legend(ax,[adiabatLine,clickedPoint], Location="northeast", AutoUpdate="off");                
                else
                [diabat, ~, adiabat, ~] = findAdiabat(x,cSpec,cStates, xIdx, yIdx);
                diabatLine = plot(ax,x, diabat, '-', Color="red", LineWidth=1.7, ...
                    DisplayName="Diabat",Tag="updatePlot");
                adiabatLine = plot(ax,x, adiabat, '-', Color="blue",LineWidth=1.2, ...
                    DisplayName="Adiabat",Tag="updatePlot");
                legend(ax,[adiabatLine,diabatLine,clickedPoint], Location="northeast", AutoUpdate="off");
                end
            else
                legend(ax, clickedPoint, Location="northeast", AutoUpdate="off");
            end

            % turn off HitTest of plot to be able to click them
            for k = 1:length(findobj(ax,'Tag','updatePlot'))
                hmm = findobj(ax,'Tag','updatePlot');
                hmm(k).HitTest = 'off';
            end
        end

        % Value changed function: BField, EField, Laser1Field
        function FieldValueChanged(comp, event)
            updateX(comp);
            if comp.UserData.autoUpdate
                calculateSpectrum(comp);
                updateY(comp);
                updatePlot(comp);            
            end
        end

        % Value changed function: BDirField, EDirField
        function FieldDirectionChanged(comp, event)
            Field = comp.Fields;
            Bdir = eval(comp.BDirField.Value);
            Edir = eval(comp.EDirField.Value);
            Idir = eval(comp.Laser1DirField.Value);
            assert(length(Bdir)==3 && sum(Bdir)~=0);
            assert(length(Edir)==3 && sum(Edir)~=0);
            assert(length(Idir)==3 && sum(Idir)~=0);

            Field.B.dir=Bdir/norm(Bdir); 
            Field.E.dir=Edir/norm(Edir); 
            Field.I.dir=Idir/norm(Idir); 
            
            comp.Fields = Field;
            if comp.UserData.autoUpdate
                updateHamiltonian(comp);
                calculateSpectrum(comp);
                updatePlot(comp);   
            else
                comp.UserData.hamChanged = 1;
            end
        end

        % Value changed function: MoleculeEditField
        function MoleculeChanged(comp, event)
            value = comp.MoleculeEditField.Value;
            try
                molecule = eval(value);
            catch
                error("%s is not a valid Molecule Class", value);
            end
            comp.Molecule = molecule;
            comp.MoleculeInfo.Value = comp.Molecule.infoString;
            
            if comp.UserData.autoUpdate
                updateHamiltonian(comp);
                calculateSpectrum(comp);
                updatePlot(comp); 
            else
                comp.UserData.hamChanged = 1;
            end
        end

        % Value changed function: Y1DropDown, Y2DropDown
        function YAxisChanged(comp, event)
            updateVars(comp);
            updatePlot(comp);   
        end

        % Value changed function: NmaxEditField
        function NmaxChanged(comp, event)
            Nmax = comp.NmaxEditField.Value;
            comp.UserData.Nmax = Nmax;
            assert(max(comp.UserData.selectN)<= Nmax);

            % updateVars(comp);
            if comp.UserData.autoUpdate
                updateHamiltonian(comp);
                calculateSpectrum(comp);
                updatePlot(comp);
            else
                comp.UserData.hamChanged = 1;
            end
        end

        % Callback function: Tree
        function HyperfineChanged(comp, event)
            hOpts.useRigidRotor = 0;
            hOpts.useSpinSpinScalar = 0;
            hOpts.useSpinSpinTensor = 0;
            hOpts.useSpinRotation = 0;
            hOpts.useNuclearElectric = 0;
            n = length(comp.Tree.CheckedNodes);
            for k = 1:n
                t = comp.Tree.CheckedNodes(k);
                switch string(t.Text)
                    case "Rigid Rotor"
                        hOpts.useRigidRotor =1;
                    case "Spin-Spin (Scalar)"
                        hOpts.useSpinSpinScalar =1;
                    case "Spin-Spin (Tensor)"
                        hOpts.useSpinSpinTensor =1;
                    case "Spin-Rotation"
                        hOpts.useSpinRotation =1;
                    case "Nuclear Electric"
                        hOpts.useNuclearElectric =1;
                    otherwise
                end
            end
            comp.UserData.hamOpts = hOpts;
            
            if comp.UserData.autoUpdate
                calculateSpectrum(comp);
                updatePlot(comp); 
            else
                comp.UserData.hamChanged = 1;
            end
        end

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(comp, event)
            if comp.UserData.hamChanged
                updateHamiltonian(comp);
                comp.UserData.hamChanged = 0;
            end
            calculateSpectrum(comp);
            updatePlot(comp); 
        end

        % Value changed function: AutoUpdateButton
        function AutoUpdateButtonValueChanged(comp, event)
            value = comp.AutoUpdateButton.Value;
            comp.UserData.autoUpdate = value;
            if value %turn off manual calculate button
                comp.CalculateButton.Enable = "off";
            else
                comp.CalculateButton.Enable = "on";
            end
        end

        % Value changed function: N1EditField
        function plotNChanged(comp, event)
            %update N Manifold to plot
            ud = comp.UserData;
            selectN = int32(eval(comp.N1EditField.Value));
            assert(max(selectN)<=ud.Nmax);
            if ~isequal(ud.selectN,selectN)
            ud.selectN = selectN;
            ud.selectStates = ismember(ud.statesUCBasis.N,selectN);
            end
            comp.UserData = ud;

            updateY(comp);
            updatePlot(comp);
            plotTDMCheckboxValueChanged(comp);
        end

        % Value changed function: TDMInitialStateEditField
        function TDMInitStateChanged(comp, event)
            stateIdx = comp.TDMInitialStateEditField.Value;
            assert(stateIdx <= comp.UserData.NStates);
            comp.UserData.TDMInitialState = stateIdx;

            calculateTDM(comp);
            updateTDMTable(comp);
            plotTDMCheckboxValueChanged(comp);
        end

        % Selection changed function: TDMTable
        function TDMTableSelectionChanged(comp, event)
            ud = comp.UserData;
            yIdx = table2array(comp.TDMTable.Data(comp.TDMTable.Selection,"StateIdx")); %is always only one
            x = comp.TDMVarValue.Value;
            xVals = comp.UserData.xVar.value;
            [~,xIdx] = min(abs(x-xVals)); 
            comp.UserData.clickedState = reshape(ud.states(xIdx,:,yIdx),[],1);
            updateCompTable(comp);
        end

        % Value changed function: plotTDMCheckbox
        function plotTDMCheckboxValueChanged(comp, event)
            value = comp.plotTDMCheckbox.Value;
            ax = comp.UIAxes;
            delete(findobj(ax,'Tag','TDM'));
            specLines = findobj(ax,Type='line',Tag='Spec');
            if value
                setStatus(comp, 'yellow', 'Plotting');
                if isempty(comp.UserData.TDM)
                    calculateTDM(comp);
                end
                for k = 1:length(specLines)
                    specLines(k).Color = [0,0,0,0.1]; %make spectrum more transparent
                end
                % colormap(ax,colorcet('L13'));
                ud = comp.UserData;
                x = ud.xVar.value;
                idx = 1:length(ud.selectStates);
                tdms = rescale(round(abs(ud.TDM),4));
                for k=idx(ud.selectStates)
                    y = ud.yVar1.value(:,k)';
                    z= max(squeeze(tdms(:,k,:)),[],2);
                    if any(z) %do with patchline, just overlay with transparent color
                        patch(ax,XData=[x nan],YData=[y nan], edgecolor='r',Tag='TDM', linewidth=1,...
                            FaceVertexAlphaData=[z; nan], EdgeAlpha='interp', AlphaDataMapping='none', HitTest='off'); 
                    end
                end
                setStatus(comp, 'green', '');
            else
                for k = 1:length(specLines)
                    specLines(k).Color = [0,0,0,0.4]; %set spectral lines back to initial transparency
                end
            end

        end

        % Value changed function: LabelsCheckBox
        function LabelsCheckBoxValueChanged(comp, event)
            ax = comp.UIAxes;
            ud = comp.UserData;
            delete(findobj(ax,Tag='StateLabels'))
            if comp.LabelsCheckBox.Value
                y1 = ud.yVar1.value(:,ud.selectStates);
                idx = find(ud.selectStates);
                for k=1:size(y1,2)
                    text(ax,ud.currPlotXlim(2),y1(end,k), sprintf("%.0f",idx(k)), Tag='StateLabels')
                end
            end    
        end

        % Value changed function: ShowPointsCheckBox
        function ShowPointsCheckBoxValueChanged(comp, event)
            ax = comp.UIAxes;
            ud = comp.UserData;
            delete(findobj(ax,Tag='SpecPoints'));
            if comp.ShowPointsCheckBox.Value
                hold(ax,'on');
                scatter(ax,ud.xVar.value,ud.yVar1.value(:,ud.selectStates), 5,'k','filled',MarkerFaceAlpha=0.5, MarkerEdgeColor='none', Tag='SpecPoints', HitTest='off');
                hold(ax,'off');
            end
        end

        % Value changed function: TDMVarValue
        function TDMVarValueValueChanged(comp, event)
            updateTDMTable(comp);
        end

        % Value changed function: Laser1DirField, Laser1PolAngleField, 
        % ...and 2 other components
        function ACStarkChanged(comp, event)
            Field = comp.Fields;
            Idir = eval(comp.Laser1DirField.Value);
            assert(length(Idir)==3 && sum(Idir)~=0);
            Field.I.dir=Idir/norm(Idir); 

            Field.I.polAngle = (eval(comp.Laser1PolAngleField.Value));
            Field.I.polPhase = (eval(comp.Laser1PolPhaseField.Value));
            a = Field.I.polPhase*pi/180;
            b = Field.I.polAngle*pi/180;%needed to evaluate Laser1PolField
            if (length(a)>1)||(length(b)>1)
                comp.Fields = Field;
                FieldValueChanged(comp);
            else
                Field.I.pol = round(eval(comp.Laser1PolField.Value),10); %makes calculation a little bit easier
                comp.Fields = Field;
                H = comp.Ham;
                Mol = comp.Molecule;
                UCBasis = comp.UserData.statesUCBasis;
                try
                    comp.Ham.ac_stark = H.makeACStark(Mol.a0,Mol.a2,UCBasis.N,UCBasis.mN, Field=Field.I);
                    setStatus(comp, 'green', '');
                catch ME
                    setStatus(comp, 'red', ME.message);
                end
    
                if comp.UserData.autoUpdate
                    calculateSpectrum(comp);
                    updatePlot(comp);   
                else
                    comp.UserData.hamChanged = 1;
                end
            end
        end

        % Button pushed function: ToWorkspaceButton
        function ToWorkspaceButtonPushed(comp, event)
            assignin('base', 'Ham',comp.Ham);
            assignin('base', 'X',comp.UserData.xVar);
            assignin('base', 'Y',comp.UserData.yVar1);
            assignin('base', 'States',comp.UserData.currStates);
        end

        % Value changed function: BasisChoice
        function BasisChoiceValueChanged(comp, event)
            comp.UserData.BasisChoice = comp.BasisChoice.Value;
            updateCompTable(comp);
        end
    end

    methods (Access = protected)
        
        % Code that executes when the value of a public property is changed
        function update(comp)
            % Use this function to update the underlying components
            % BE VERY CAREFUL WITH THIS FUNCTION --> can lead to infinite
            % loops in the debugger, rendering it completely useless
            % I suspect MATLAB calls this function a LOT of times in the
            % background
        end

        % Create the underlying components
        function setup(comp)

            comp.Position = [0 0 1280 720];
            comp.BackgroundColor = [0.94 0.94 0.94];

            % Create GridLayout
            comp.GridLayout = uigridlayout(comp);
            comp.GridLayout.ColumnWidth = {'1x', 350};
            comp.GridLayout.RowHeight = {'2x', '1x'};
            comp.GridLayout.ColumnSpacing = 0;
            comp.GridLayout.RowSpacing = 5;
            comp.GridLayout.Padding = [2 2 15 2];
            comp.GridLayout.Interruptible = 'off';
            comp.GridLayout.BackgroundColor = [0.9804 0.9804 0.9804];

            % Create UIAxes
            comp.UIAxes = uiaxes(comp.GridLayout);
            xlabel(comp.UIAxes, 'X')
            ylabel(comp.UIAxes, 'E (MHz)')
            zlabel(comp.UIAxes, 'Z')
            comp.UIAxes.GridColor = [0 0 0];
            comp.UIAxes.GridAlpha = 0.05;
            comp.UIAxes.Box = 'on';
            comp.UIAxes.XGrid = 'on';
            comp.UIAxes.YGrid = 'on';
            comp.UIAxes.FontSize = 12;
            comp.UIAxes.Layout.Row = 1;
            comp.UIAxes.Layout.Column = [1 2];
            comp.UIAxes.ButtonDownFcn = matlab.apps.createCallbackFcn(comp, @UIAxesButtonDown, true);
            comp.UIAxes.Interruptible = 'off';

            % Create UITable
            comp.UITable = uitable(comp.GridLayout);
            comp.UITable.ColumnName = {'Comp'; 'N'; 'mN'; 'mi1'; 'mi2'};
            comp.UITable.ColumnWidth = {'2x', '1x', '1x', '1x', '1x'};
            comp.UITable.RowName = {};
            comp.UITable.SelectionType = 'row';
            comp.UITable.Interruptible = 'off';
            comp.UITable.Layout.Row = 2;
            comp.UITable.Layout.Column = 2;

            % Create GridLayout2
            comp.GridLayout2 = uigridlayout(comp.GridLayout);
            comp.GridLayout2.ColumnWidth = {'1x'};
            comp.GridLayout2.RowHeight = {'6x', '1x'};
            comp.GridLayout2.ColumnSpacing = 2;
            comp.GridLayout2.RowSpacing = 2;
            comp.GridLayout2.Padding = [2 2 2 0];
            comp.GridLayout2.Interruptible = 'off';
            comp.GridLayout2.Layout.Row = 2;
            comp.GridLayout2.Layout.Column = 1;

            % Create TabGroup
            comp.TabGroup = uitabgroup(comp.GridLayout2);
            comp.TabGroup.Interruptible = 'off';
            comp.TabGroup.Layout.Row = 1;
            comp.TabGroup.Layout.Column = 1;

            % Create FieldsTab
            comp.FieldsTab = uitab(comp.TabGroup);
            comp.FieldsTab.Title = 'Fields';

            % Create GridLayout4
            comp.GridLayout4 = uigridlayout(comp.FieldsTab);
            comp.GridLayout4.ColumnWidth = {20, '3.37x', '1.62x', '4x', 20, '2x', 20, '2x'};
            comp.GridLayout4.RowHeight = {23, 23, 23, 23, 23};
            comp.GridLayout4.Interruptible = 'off';

            % Create Laser1Field
            comp.Laser1Field = uieditfield(comp.GridLayout4, 'text');
            comp.Laser1Field.CharacterLimits = [0 100];
            comp.Laser1Field.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @FieldValueChanged, true);
            comp.Laser1Field.Tooltip = {'Intensity (kW/cm^2)'};
            comp.Laser1Field.Layout.Row = 4;
            comp.Laser1Field.Layout.Column = 2;
            comp.Laser1Field.Value = '0';

            % Create EField
            comp.EField = uieditfield(comp.GridLayout4, 'text');
            comp.EField.CharacterLimits = [0 100];
            comp.EField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @FieldValueChanged, true);
            comp.EField.Tooltip = {'Electric Field (V/cm)'};
            comp.EField.Layout.Row = 3;
            comp.EField.Layout.Column = 2;
            comp.EField.Value = '0';

            % Create EDirField
            comp.EDirField = uieditfield(comp.GridLayout4, 'text');
            comp.EDirField.CharacterLimits = [7 100];
            comp.EDirField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @FieldDirectionChanged, true);
            comp.EDirField.Tooltip = {'Electric Field Direction'};
            comp.EDirField.Layout.Row = 3;
            comp.EDirField.Layout.Column = 3;
            comp.EDirField.Value = '[0,1,0]';

            % Create Laser1DirField
            comp.Laser1DirField = uieditfield(comp.GridLayout4, 'text');
            comp.Laser1DirField.CharacterLimits = [7 100];
            comp.Laser1DirField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @ACStarkChanged, true);
            comp.Laser1DirField.Tooltip = {'Light Field Direction'};
            comp.Laser1DirField.Layout.Row = 4;
            comp.Laser1DirField.Layout.Column = 3;
            comp.Laser1DirField.Value = '[1,0,0]';

            % Create ELabel
            comp.ELabel = uilabel(comp.GridLayout4);
            comp.ELabel.Layout.Row = 3;
            comp.ELabel.Layout.Column = 1;
            comp.ELabel.Text = 'E';

            % Create ILabel
            comp.ILabel = uilabel(comp.GridLayout4);
            comp.ILabel.Layout.Row = 4;
            comp.ILabel.Layout.Column = 1;
            comp.ILabel.Text = 'I';

            % Create degreeLabel
            comp.degreeLabel = uilabel(comp.GridLayout4);
            comp.degreeLabel.HorizontalAlignment = 'center';
            comp.degreeLabel.Layout.Row = 4;
            comp.degreeLabel.Layout.Column = 7;
            comp.degreeLabel.Text = 'b';

            % Create Laser1PolAngleField
            comp.Laser1PolAngleField = uieditfield(comp.GridLayout4, 'text');
            comp.Laser1PolAngleField.CharacterLimits = [0 100];
            comp.Laser1PolAngleField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @ACStarkChanged, true);
            comp.Laser1PolAngleField.Tooltip = {'Polarisation Angle (°)'};
            comp.Laser1PolAngleField.Layout.Row = 4;
            comp.Laser1PolAngleField.Layout.Column = 8;
            comp.Laser1PolAngleField.Value = '0';

            % Create aEditFieldLabel
            comp.aEditFieldLabel = uilabel(comp.GridLayout4);
            comp.aEditFieldLabel.HorizontalAlignment = 'center';
            comp.aEditFieldLabel.Layout.Row = 4;
            comp.aEditFieldLabel.Layout.Column = 5;
            comp.aEditFieldLabel.Text = 'a';

            % Create Laser1PolPhaseField
            comp.Laser1PolPhaseField = uieditfield(comp.GridLayout4, 'text');
            comp.Laser1PolPhaseField.CharacterLimits = [0 100];
            comp.Laser1PolPhaseField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @ACStarkChanged, true);
            comp.Laser1PolPhaseField.Tooltip = {'Ellipticity (°)'};
            comp.Laser1PolPhaseField.Layout.Row = 4;
            comp.Laser1PolPhaseField.Layout.Column = 6;
            comp.Laser1PolPhaseField.Value = '0';

            % Create BDirField
            comp.BDirField = uieditfield(comp.GridLayout4, 'text');
            comp.BDirField.CharacterLimits = [7 100];
            comp.BDirField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @FieldDirectionChanged, true);
            comp.BDirField.Tooltip = {'Magnetic Field Direction'};
            comp.BDirField.Layout.Row = 2;
            comp.BDirField.Layout.Column = 3;
            comp.BDirField.Value = '[0,0,1]';

            % Create BField
            comp.BField = uieditfield(comp.GridLayout4, 'text');
            comp.BField.CharacterLimits = [1 100];
            comp.BField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @FieldValueChanged, true);
            comp.BField.Tooltip = {'Magnetic Field (G)'};
            comp.BField.Layout.Row = 2;
            comp.BField.Layout.Column = 2;
            comp.BField.Value = 'linspace(1,50,50)';

            % Create BLabel
            comp.BLabel = uilabel(comp.GridLayout4);
            comp.BLabel.Layout.Row = 2;
            comp.BLabel.Layout.Column = 1;
            comp.BLabel.Text = 'B';

            % Create Laser1PolField
            comp.Laser1PolField = uieditfield(comp.GridLayout4, 'text');
            comp.Laser1PolField.CharacterLimits = [7 100];
            comp.Laser1PolField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @ACStarkChanged, true);
            comp.Laser1PolField.Tooltip = {'Polarisation Vector. Should be perpendicular to light field.'};
            comp.Laser1PolField.Layout.Row = 4;
            comp.Laser1PolField.Layout.Column = 4;
            comp.Laser1PolField.Value = '[0, sin(b)*exp(1i*a), cos(b)]';

            % Create FieldAmplitudeLabel
            comp.FieldAmplitudeLabel = uilabel(comp.GridLayout4);
            comp.FieldAmplitudeLabel.HorizontalAlignment = 'center';
            comp.FieldAmplitudeLabel.Layout.Row = 1;
            comp.FieldAmplitudeLabel.Layout.Column = 2;
            comp.FieldAmplitudeLabel.Text = 'Field Amplitude';

            % Create DirectionLabel_2
            comp.DirectionLabel_2 = uilabel(comp.GridLayout4);
            comp.DirectionLabel_2.HorizontalAlignment = 'center';
            comp.DirectionLabel_2.Layout.Row = 1;
            comp.DirectionLabel_2.Layout.Column = 3;
            comp.DirectionLabel_2.Text = 'Direction';

            % Create PolarisationLabel
            comp.PolarisationLabel = uilabel(comp.GridLayout4);
            comp.PolarisationLabel.HorizontalAlignment = 'center';
            comp.PolarisationLabel.Layout.Row = 1;
            comp.PolarisationLabel.Layout.Column = 4;
            comp.PolarisationLabel.Text = 'Polarisation';

            % Create MoleculeTab
            comp.MoleculeTab = uitab(comp.TabGroup);
            comp.MoleculeTab.Title = 'Molecule';

            % Create GridLayout7
            comp.GridLayout7 = uigridlayout(comp.MoleculeTab);
            comp.GridLayout7.ColumnWidth = {56, '1.45x', '2.58x', '1x', '1x'};
            comp.GridLayout7.RowHeight = {25, '1x', 40, '1.8x'};

            % Create Tree
            comp.Tree = uitree(comp.GridLayout7, 'checkbox');
            comp.Tree.Tooltip = {'Terms of the hyperfine interaction to include'};
            comp.Tree.Layout.Row = [2 4];
            comp.Tree.Layout.Column = [1 2];

            % Create RigidRotorNode
            comp.RigidRotorNode = uitreenode(comp.Tree);
            comp.RigidRotorNode.Text = 'Rigid Rotor';

            % Create SpinSpinScalarNode
            comp.SpinSpinScalarNode = uitreenode(comp.Tree);
            comp.SpinSpinScalarNode.Text = 'Spin-Spin (Scalar)';

            % Create SpinSpinTensorNode
            comp.SpinSpinTensorNode = uitreenode(comp.Tree);
            comp.SpinSpinTensorNode.Text = 'Spin-Spin (Tensor)';

            % Create SpinRotationNode
            comp.SpinRotationNode = uitreenode(comp.Tree);
            comp.SpinRotationNode.Text = 'Spin-Rotation';

            % Create NuclearElectricNode
            comp.NuclearElectricNode = uitreenode(comp.Tree);
            comp.NuclearElectricNode.Text = 'Nuclear Electric';

            % Assign Checked Nodes
            comp.Tree.CheckedNodes = [comp.RigidRotorNode, comp.SpinSpinScalarNode, comp.SpinSpinTensorNode, comp.SpinRotationNode, comp.NuclearElectricNode];
            % Assign Checked Nodes
            comp.Tree.CheckedNodesChangedFcn = matlab.apps.createCallbackFcn(comp, @HyperfineChanged, true);

            % Create MoleculeLabel
            comp.MoleculeLabel = uilabel(comp.GridLayout7);
            comp.MoleculeLabel.HorizontalAlignment = 'right';
            comp.MoleculeLabel.Layout.Row = 1;
            comp.MoleculeLabel.Layout.Column = 1;
            comp.MoleculeLabel.Text = 'Molecule:';

            % Create MoleculeEditField
            comp.MoleculeEditField = uieditfield(comp.GridLayout7, 'text');
            comp.MoleculeEditField.CharacterLimits = [0 15];
            comp.MoleculeEditField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @MoleculeChanged, true);
            comp.MoleculeEditField.HorizontalAlignment = 'center';
            comp.MoleculeEditField.Tooltip = {'Supported Molecules:'; 'KRb,'; 'NaCs,'; 'NaK,'; 'NaRb,'; 'RbCs.'; 'See Molecule classes for supported and default isotopes.'};
            comp.MoleculeEditField.Placeholder = 'Molecule';
            comp.MoleculeEditField.Layout.Row = 1;
            comp.MoleculeEditField.Layout.Column = 2;
            comp.MoleculeEditField.Value = 'KRb(40,87)';

            % Create MoleculeInfo
            comp.MoleculeInfo = uitextarea(comp.GridLayout7);
            comp.MoleculeInfo.WordWrap = 'off';
            comp.MoleculeInfo.Layout.Row = [1 4];
            comp.MoleculeInfo.Layout.Column = 3;

            % Create BasisTab
            comp.BasisTab = uitab(comp.TabGroup);
            comp.BasisTab.Title = 'Basis';

            % Create GridLayout9
            comp.GridLayout9 = uigridlayout(comp.BasisTab);
            comp.GridLayout9.ColumnWidth = {50, 150, '1x', '1x'};
            comp.GridLayout9.RowHeight = {25, 25, '1x'};

            % Create NmaxEditField
            comp.NmaxEditField = uieditfield(comp.GridLayout9, 'numeric');
            comp.NmaxEditField.Limits = [0 5];
            comp.NmaxEditField.RoundFractionalValues = 'on';
            comp.NmaxEditField.ValueDisplayFormat = '%.0f';
            comp.NmaxEditField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @NmaxChanged, true);
            comp.NmaxEditField.HorizontalAlignment = 'center';
            comp.NmaxEditField.Tooltip = {'Maximum Number of rotational states to include in calculation. N>4 will take significant time to calculate. Currently limited to N=5 by design.'};
            comp.NmaxEditField.Layout.Row = 1;
            comp.NmaxEditField.Layout.Column = 2;
            comp.NmaxEditField.Value = 2;

            % Create NmaxLabel
            comp.NmaxLabel = uilabel(comp.GridLayout9);
            comp.NmaxLabel.HorizontalAlignment = 'right';
            comp.NmaxLabel.Layout.Row = 1;
            comp.NmaxLabel.Layout.Column = 1;
            comp.NmaxLabel.Text = 'Nmax';

            % Create BasisDropDownLabel
            comp.BasisDropDownLabel = uilabel(comp.GridLayout9);
            comp.BasisDropDownLabel.HorizontalAlignment = 'right';
            comp.BasisDropDownLabel.Layout.Row = 2;
            comp.BasisDropDownLabel.Layout.Column = 1;
            comp.BasisDropDownLabel.Text = 'Basis';

            % Create BasisChoice
            comp.BasisChoice = uidropdown(comp.GridLayout9);
            comp.BasisChoice.Items = {'Uncoupled', 'Fully Coupled'};
            comp.BasisChoice.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @BasisChoiceValueChanged, true);
            comp.BasisChoice.Layout.Row = 2;
            comp.BasisChoice.Layout.Column = 2;
            comp.BasisChoice.Value = 'Uncoupled';

            % Create PlotTab
            comp.PlotTab = uitab(comp.TabGroup);
            comp.PlotTab.Title = 'Plot';

            % Create GridLayout5
            comp.GridLayout5 = uigridlayout(comp.PlotTab);
            comp.GridLayout5.ColumnWidth = {25, 65, '1x', 36, '1.39x', '4.02x'};
            comp.GridLayout5.RowHeight = {20, 20, '1x', 20};

            % Create PlotAdiabatsCheckBox
            comp.PlotAdiabatsCheckBox = uicheckbox(comp.GridLayout5);
            comp.PlotAdiabatsCheckBox.Text = 'Plot Adiabats';
            comp.PlotAdiabatsCheckBox.Layout.Row = 4;
            comp.PlotAdiabatsCheckBox.Layout.Column = 5;

            % Create ShowPointsCheckBox
            comp.ShowPointsCheckBox = uicheckbox(comp.GridLayout5);
            comp.ShowPointsCheckBox.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @ShowPointsCheckBoxValueChanged, true);
            comp.ShowPointsCheckBox.Text = 'Show Points';
            comp.ShowPointsCheckBox.Layout.Row = 4;
            comp.ShowPointsCheckBox.Layout.Column = [1 2];

            % Create Y1DropDownLabel
            comp.Y1DropDownLabel = uilabel(comp.GridLayout5);
            comp.Y1DropDownLabel.Layout.Row = 1;
            comp.Y1DropDownLabel.Layout.Column = 1;
            comp.Y1DropDownLabel.Text = 'Y1';

            % Create Y1DropDown
            comp.Y1DropDown = uidropdown(comp.GridLayout5);
            comp.Y1DropDown.Items = {'E', 'µ', 'd', 'E - <E>', 'E - a0', 'E - GS'};
            comp.Y1DropDown.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @YAxisChanged, true);
            comp.Y1DropDown.Layout.Row = 1;
            comp.Y1DropDown.Layout.Column = [2 3];
            comp.Y1DropDown.Value = 'E';

            % Create Y2DropDownLabel
            comp.Y2DropDownLabel = uilabel(comp.GridLayout5);
            comp.Y2DropDownLabel.Enable = 'off';
            comp.Y2DropDownLabel.Visible = 'off';
            comp.Y2DropDownLabel.Layout.Row = 1;
            comp.Y2DropDownLabel.Layout.Column = 5;
            comp.Y2DropDownLabel.Text = 'Y2';

            % Create Y2DropDown
            comp.Y2DropDown = uidropdown(comp.GridLayout5);
            comp.Y2DropDown.Items = {'E', 'E - <E>', 'µ', 'd'};
            comp.Y2DropDown.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @YAxisChanged, true);
            comp.Y2DropDown.Enable = 'off';
            comp.Y2DropDown.Visible = 'off';
            comp.Y2DropDown.Layout.Row = 1;
            comp.Y2DropDown.Layout.Column = [5 6];
            comp.Y2DropDown.Value = 'E';

            % Create N1EditFieldLabel
            comp.N1EditFieldLabel = uilabel(comp.GridLayout5);
            comp.N1EditFieldLabel.Layout.Row = 2;
            comp.N1EditFieldLabel.Layout.Column = 1;
            comp.N1EditFieldLabel.Text = 'N1';

            % Create N1EditField
            comp.N1EditField = uieditfield(comp.GridLayout5, 'text');
            comp.N1EditField.CharacterLimits = [1 10];
            comp.N1EditField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @plotNChanged, true);
            comp.N1EditField.Tooltip = {'Rotational state manifold to plot. Must be smaller than Nmax.'};
            comp.N1EditField.Layout.Row = 2;
            comp.N1EditField.Layout.Column = [2 3];
            comp.N1EditField.Value = '0';

            % Create N2EditFieldLabel
            comp.N2EditFieldLabel = uilabel(comp.GridLayout5);
            comp.N2EditFieldLabel.HorizontalAlignment = 'right';
            comp.N2EditFieldLabel.Enable = 'off';
            comp.N2EditFieldLabel.Layout.Row = 2;
            comp.N2EditFieldLabel.Layout.Column = 4;
            comp.N2EditFieldLabel.Text = 'N2';

            % Create N2EditField
            comp.N2EditField = uieditfield(comp.GridLayout5, 'text');
            comp.N2EditField.Enable = 'off';
            comp.N2EditField.Layout.Row = 2;
            comp.N2EditField.Layout.Column = [5 6];

            % Create LabelsCheckBox
            comp.LabelsCheckBox = uicheckbox(comp.GridLayout5);
            comp.LabelsCheckBox.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @LabelsCheckBoxValueChanged, true);
            comp.LabelsCheckBox.Text = 'Plot Labels';
            comp.LabelsCheckBox.Layout.Row = 4;
            comp.LabelsCheckBox.Layout.Column = 3;

            % Create TDMTab
            comp.TDMTab = uitab(comp.TabGroup);
            comp.TDMTab.Title = 'TDM';

            % Create GridLayout8
            comp.GridLayout8 = uigridlayout(comp.TDMTab);
            comp.GridLayout8.ColumnWidth = {80, 100, 50, '1x', '1x'};
            comp.GridLayout8.RowHeight = {'1x', '1x', '1x', '1x', '1x'};

            % Create InitialStateEditFieldLabel
            comp.InitialStateEditFieldLabel = uilabel(comp.GridLayout8);
            comp.InitialStateEditFieldLabel.HorizontalAlignment = 'center';
            comp.InitialStateEditFieldLabel.Layout.Row = 1;
            comp.InitialStateEditFieldLabel.Layout.Column = 1;
            comp.InitialStateEditFieldLabel.Text = 'Initial State';

            % Create TDMInitialStateEditField
            comp.TDMInitialStateEditField = uieditfield(comp.GridLayout8, 'numeric');
            comp.TDMInitialStateEditField.Limits = [1 100000];
            comp.TDMInitialStateEditField.RoundFractionalValues = 'on';
            comp.TDMInitialStateEditField.ValueDisplayFormat = '%.0f';
            comp.TDMInitialStateEditField.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @TDMInitStateChanged, true);
            comp.TDMInitialStateEditField.Tag = 'TDMinitState';
            comp.TDMInitialStateEditField.HorizontalAlignment = 'center';
            comp.TDMInitialStateEditField.Layout.Row = 1;
            comp.TDMInitialStateEditField.Layout.Column = 2;
            comp.TDMInitialStateEditField.Value = 1;

            % Create TDMTable
            comp.TDMTable = uitable(comp.GridLayout8);
            comp.TDMTable.ColumnName = {'E (MHz)'; 's-'; 'p'; 's+'; 'State'};
            comp.TDMTable.ColumnWidth = {'3x', '1x', '1x', '1x', '1x'};
            comp.TDMTable.RowName = {};
            comp.TDMTable.ColumnSortable = true;
            comp.TDMTable.SelectionType = 'row';
            comp.TDMTable.SelectionChangedFcn = matlab.apps.createCallbackFcn(comp, @TDMTableSelectionChanged, true);
            comp.TDMTable.Multiselect = 'off';
            comp.TDMTable.Layout.Row = [1 5];
            comp.TDMTable.Layout.Column = [4 5];

            % Create TDMcurrVarLabel
            comp.TDMcurrVarLabel = uilabel(comp.GridLayout8);
            comp.TDMcurrVarLabel.HorizontalAlignment = 'center';
            comp.TDMcurrVarLabel.Layout.Row = 2;
            comp.TDMcurrVarLabel.Layout.Column = 1;
            comp.TDMcurrVarLabel.Text = 'xVar';

            % Create TDMVarValue
            comp.TDMVarValue = uieditfield(comp.GridLayout8, 'numeric');
            comp.TDMVarValue.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @TDMVarValueValueChanged, true);
            comp.TDMVarValue.Tag = 'TDMVarValue';
            comp.TDMVarValue.HorizontalAlignment = 'center';
            comp.TDMVarValue.Layout.Row = 2;
            comp.TDMVarValue.Layout.Column = 2;

            % Create plotTDMCheckbox
            comp.plotTDMCheckbox = uicheckbox(comp.GridLayout8);
            comp.plotTDMCheckbox.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @plotTDMCheckboxValueChanged, true);
            comp.plotTDMCheckbox.Text = '';
            comp.plotTDMCheckbox.Layout.Row = 3;
            comp.plotTDMCheckbox.Layout.Column = 2;

            % Create PlotTDMLabel
            comp.PlotTDMLabel = uilabel(comp.GridLayout8);
            comp.PlotTDMLabel.HorizontalAlignment = 'center';
            comp.PlotTDMLabel.Layout.Row = 3;
            comp.PlotTDMLabel.Layout.Column = 1;
            comp.PlotTDMLabel.Text = 'Plot TDM';

            % Create GridLayout3
            comp.GridLayout3 = uigridlayout(comp.GridLayout2);
            comp.GridLayout3.ColumnWidth = {40, '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            comp.GridLayout3.RowHeight = {'1x'};
            comp.GridLayout3.Padding = [0 0 0 0];
            comp.GridLayout3.Interruptible = 'off';
            comp.GridLayout3.Layout.Row = 2;
            comp.GridLayout3.Layout.Column = 1;

            % Create CalculateButton
            comp.CalculateButton = uibutton(comp.GridLayout3, 'push');
            comp.CalculateButton.ButtonPushedFcn = matlab.apps.createCallbackFcn(comp, @CalculateButtonPushed, true);
            comp.CalculateButton.Layout.Row = 1;
            comp.CalculateButton.Layout.Column = 3;
            comp.CalculateButton.Text = 'Calculate';

            % Create AutoUpdateButton
            comp.AutoUpdateButton = uibutton(comp.GridLayout3, 'state');
            comp.AutoUpdateButton.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @AutoUpdateButtonValueChanged, true);
            comp.AutoUpdateButton.Interruptible = 'off';
            comp.AutoUpdateButton.Text = 'AutoUpdate';
            comp.AutoUpdateButton.Layout.Row = 1;
            comp.AutoUpdateButton.Layout.Column = 4;

            % Create Label
            comp.Label = uilabel(comp.GridLayout3);
            comp.Label.HorizontalAlignment = 'right';
            comp.Label.Layout.Row = 1;
            comp.Label.Layout.Column = 2;
            comp.Label.Text = '';

            % Create statusLamp
            comp.statusLamp = uilamp(comp.GridLayout3);
            comp.statusLamp.BusyAction = 'cancel';
            comp.statusLamp.Interruptible = 'off';
            comp.statusLamp.Layout.Row = 1;
            comp.statusLamp.Layout.Column = 1;

            % Create statusLabel
            comp.statusLabel = uilabel(comp.GridLayout3);
            comp.statusLabel.Interruptible = 'off';
            comp.statusLabel.Layout.Row = 1;
            comp.statusLabel.Layout.Column = 2;
            comp.statusLabel.Text = '';

            % Create ToWorkspaceButton
            comp.ToWorkspaceButton = uibutton(comp.GridLayout3, 'push');
            comp.ToWorkspaceButton.ButtonPushedFcn = matlab.apps.createCallbackFcn(comp, @ToWorkspaceButtonPushed, true);
            comp.ToWorkspaceButton.Tooltip = {''};
            comp.ToWorkspaceButton.Layout.Row = 1;
            comp.ToWorkspaceButton.Layout.Column = 5;
            comp.ToWorkspaceButton.Text = 'To Workspace';

            % Create SaveButton
            comp.SaveButton = uibutton(comp.GridLayout3, 'push');
            comp.SaveButton.Enable = 'off';
            comp.SaveButton.Layout.Row = 1;
            comp.SaveButton.Layout.Column = 6;
            comp.SaveButton.Text = 'Save';

            % Create LoadButton
            comp.LoadButton = uibutton(comp.GridLayout3, 'push');
            comp.LoadButton.Enable = 'off';
            comp.LoadButton.Layout.Row = 1;
            comp.LoadButton.Layout.Column = 7;
            comp.LoadButton.Text = 'Load';

            % Create Button3
            comp.Button3 = uibutton(comp.GridLayout3, 'push');
            comp.Button3.Enable = 'off';
            comp.Button3.Visible = 'off';
            comp.Button3.Layout.Row = 1;
            comp.Button3.Layout.Column = 8;
            comp.Button3.Text = 'Button3';
            
            % Execute the startup function
            initialiseApp(comp)
        end
    end
end