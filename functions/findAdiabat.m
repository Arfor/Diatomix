function [diabat, diabatStates, adiabat, adiabatStates, diabatSorting] = findAdiabat(x,y,states, xIdx,yIdx)
%[diabat, diabatStates, adiabat, adiabatStates] = findAdiabat(x,y,states, xIdx,yIdx)
% x = array of x values
% y = spectrum for certain number of states
% states = array containing the state-composition of every state in y
% xIdx,yIdx = starting indices to look for adiabat 
% Doesnt work when E-field is not in z-direction?
s = size(states);
nStates = s(end);
adiabat = y(:,yIdx); %follow xth energy state, stores energies
adiabatStates = abs(squeeze(states(:,:,yIdx))).^2;
%Checks adjacent states to determine least change in state
diabatSorting = nan(1,length(x)); %stores yIdx, not energy
diabat = nan(1,length(x)); %stores yIdx, not energy

diabatSorting(xIdx) = yIdx;
diabat(xIdx) = y(xIdx,yIdx);
diabatStates = nan(size(adiabatStates));
diabatStates(xIdx,:) = (squeeze(states(xIdx,:,yIdx)));
% [~, sIdx] = sortrows(diabatStates(xIdx,:)'.^2,"descend");
adjC = 3; %set maximum number of adjacent states to look at
for k=xIdx:(length(x)-1) %forward loop
    stateNow = (squeeze(states(k,:,diabatSorting(k))));
    stateSearch = diabatSorting(k)-adjC:diabatSorting(k)+adjC; 
    stateSearch(stateSearch<=0)=[];stateSearch(stateSearch>nStates)=[];
    statesNext = conj(squeeze(states(k+1,:,stateSearch)));
    [~,S]=sort(abs(sum((statesNext.*stateNow'),1))); %computes overlap integral
    idx = stateSearch(S(end)); %restrict search to adjacent states
    if abs(idx-diabatSorting(k))>3
        idx = stateSearch(S(end-1));
    end
    diabatSorting(k+1) = idx;
    diabat(k+1) = y(k+1,idx);
    diabatStates(k+1,:) = (squeeze(states(k+1,:,idx)));
end
for k=(xIdx-1):-1:1 %backward loop
    stateIdx = diabatSorting(k+1);
    stateNow = (squeeze(states(k+1,:,stateIdx)));
    stateSearch =stateIdx-adjC:stateIdx+adjC; 
    stateSearch(stateSearch<=0)=[];stateSearch(stateSearch>nStates)=[];
    statesPrevious = conj(squeeze(states(k,:,stateSearch))); 
    [~,S]=sort(abs(sum((statesPrevious.*stateNow'),1))); %maximum overlap
    idx = stateSearch(S(end));
    if abs(idx-stateIdx)>3
        idx = stateSearch(S(end-1));
    end
    diabatSorting(k) = idx;
    diabat(k) = y(k,idx);
    diabatStates(k,:) = (squeeze(states(k,:,idx)));
end    
end

