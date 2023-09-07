function [adiabat, adiabatStates, diabat, diabatStates] = findAdiabat(x,y,states, xIdx,yIdx)
%[adiabat, adiabatStates, diabat, diabatStates] = findAdiabat(x,y,states, xIdx,yIdx)
% x = array of x values
% y = spectrum for certain number of states
% states = array containing the state-composition of every state in y
% xIdx,yIdx = starting indices to look for adiabat  
s = size(states);
nStates = s(end);
%Checks adjacent states to determine least change in state
diabat = y(:,yIdx); %follow xth energy state, stores energies
diabatStates = abs(squeeze(states(:,:,yIdx))).^2;

adiabatYidx = nan(1,length(x)); %stores yIdx, not energy
adiabat = nan(1,length(x)); %stores yIdx, not energy

adiabatYidx(xIdx) = yIdx;
adiabat(xIdx) = y(xIdx,yIdx);
adiabatStates = nan(size(diabatStates));
adiabatStates(xIdx,:) = (squeeze(states(xIdx,:,yIdx)));
% [~, sIdx] = sortrows(diabatStates(xIdx,:)'.^2,"descend");
adjC = 3; %set maximum number of adjacent states to look at
for k=xIdx:(length(x)-1) %forward loop
    stateNow = abs(squeeze(states(k,:,adiabatYidx(k))));
    stateSearch =adiabatYidx(k)-adjC:adiabatYidx(k)+adjC; stateSearch(stateSearch<=0)=[];
    statesNext = abs(squeeze(states(k+1,:,stateSearch)));
    [~,S]=sort(sum((statesNext.*stateNow'),1)); %computes overlap integral
    idx = stateSearch(S(end)); %restrict search to adjacent states
    if abs(idx-adiabatYidx(k))>3
        idx = stateSearch(S(end-1));
    end
    adiabatYidx(k+1) = idx;
    adiabat(k+1) = y(k+1,idx);
    adiabatStates(k+1,:) = (squeeze(states(k+1,:,idx)));
end
for k=(xIdx-1):-1:1 %backward loop
    stateIdx = adiabatYidx(k+1);
    stateNow = abs(squeeze(states(k+1,:,stateIdx)));
    stateSearch =stateIdx-adjC:stateIdx+adjC; 
    stateSearch(stateSearch<=0)=[];stateSearch(stateSearch>nStates)=[];
    statesPrevious = abs(squeeze(states(k,:,stateSearch))); 
    [~,S]=sort(sum((statesPrevious.*stateNow'),1)); %maximum overlap
    idx = stateSearch(S(end));
    if abs(idx-stateIdx)>3
        idx = stateSearch(S(end-1));
    end
    adiabatYidx(k) = idx;
    adiabat(k) = y(k,idx);
    adiabatStates(k,:) = (squeeze(states(k,:,idx)));
end    

end

