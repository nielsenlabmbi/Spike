function Y = alignWvs(X)

spkMin = min(X,[],2);
Y = nan(size(X));

for m = 1:length(spkMin)

    minOff = find(X(m,:)==spkMin(m))-round(size(X,2)/2);
    if isempty(minOff)
       continue 
    end
    Y(m,:) = circshift(X(m,:),minOff*-1);

end

end