%% Example - Transition Dipole Moments
clear all
addpath("Molecules\");addpath("Atoms\")
C = Constants;
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot,'defaultAxesBox','on')
set(groot,'defaultAxesLineWidth',0.5)
%% Define Fields
maxN = 2;
B.value = 542; %Gauss
B.dir = [0,0,1];
B.scaling = 1e-4; %to convert to Tesla
E.value = 420; %V/cm
E.dir = [0,1,0];
E.scaling = 1e2;%to convert to V/m
I.value = 0; %mW/cm^2
I.dir = [0,0,1];
I.pol = [1,0,0];
I.scaling = 1e7;%to convert to W/m^2
Fields.B = B;
Fields.E = E;
Fields.I = I;
%% Define Hamiltonian
Mol = KRb(40,87);
Ham = Hamiltonian(Molecule = Mol, Fields=Fields, maxN=maxN);
H = Ham.hyperfine.total + B.value*B.scaling*Ham.zeeman + I.value*I.scaling*Ham.ac_stark + E.value*E.scaling*Ham.dc_stark;
Base = Ham.Basis;
nStates = Base.NStates;
%% Diagonalise 
[V,D] = eig(full(H/C.h),'vector'); %no need to use sparse matrices for matrices smaller than 1000x1000. Divide by Plancks constant to get energy in Hz
[energies,sortIdx] = sort(real(D));
states = V(:,sortIdx);
clear V D
%%
UCBasis = Base.getStates('all');
stateIdx = 7;
stateComp = round(squeeze(states(:,stateIdx)).^2,6); %state 7 usually |-4,0.5>
[statesCompTable, ~] = sortrows([array2table(stateComp,"VariableNames","Comp") , UCBasis], "Comp","descend");
statesCompTable(1:5,:)
%%
initState = states(:,stateIdx)'; %row vec
Ediff = (energies - energies(stateIdx))*1e-6;%in MHz
TDM = nan(nStates,3);
for p = 1:3
    TDM(:,p) = Mol.d0/Constants.D*real(initState * Ham.dipoleOperator{p} * states); %in Debye
end

NSelect = [0,1];
selectStates = any(UCBasis.N==NSelect,2);
% TDM = TDM(selectStates,:);
% Ediff = Ediff(selectStates);
maxTDMs = max(TDM,[],2);
[~,sortIdx] = sort(maxTDMs,"descend");

sigma = char(931);
TDMTable = round(array2table([Ediff(sortIdx),TDM(sortIdx,:),sortIdx],"VariableNames",["E","s-","p","s+","StateIdx"]),4);
TDMTable(1:10,:)
%%
stateComp = abs(squeeze(states(:,55))).^2; %state 7 usually |-4,0.5>
[statesCompTable, ~] = sortrows([array2table(stateComp,"VariableNames","Comp") , UCBasis(:,["N","mN","mi1","mi2"])], "Comp","descend");
round(statesCompTable(1:5,:),5)

%% Something to label states by (N,mF)_index like in Cornish?
cs = cumsum(statesCompTable.Comp);
%check number of states necessary to capture 99% of the wavefunction
reqStates = find(diff(cs > 0.99));
(statesCompTable(1:reqStates,["mN","mi1","mi2"]))

(statesCompTable(1:reqStates,:))
%% Plot Results
N = Base.getStates("N");
figure(12); clf;
t = tiledlayout(2,1, "TileSpacing","tight");
nexttile(2); %N=0 manifold
plot(x,energyMap(:,N==0)*1e-6, color=[1,1,1]*0.3); hold on;
xlabel("B (G)")
ylabel(sprintf("E (MHz)"))
ylim([-1.5,1])
xlim([min(x),max(x)])
axis square

nexttile(1); %N=1 manifold
plot(x,energyMap(:,N==1)*1e-6, color=[1,1,1]*0.3);
ylabel(sprintf("E (MHz)"))
xticklabels([])
ylim([978.5,981.5])
xlim([min(x),max(x)])
axis square

title(t,"Zeeman Map")
subtitle(t, sprintf("I=%.2f W/m, E=%.2f V/m", Fields.I.value,Fields.E.value))

%% 
Ham.hyperfine.total
