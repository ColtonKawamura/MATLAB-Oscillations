function colts_find_peaks(y)

n = length(y);
    colts_peaks = [];
    for i = 2:n-1
        if y(i) > y(i-1) && y(i) > y(i+1)
            colts_peaks = [colts_peaks, i];
        end
    end
end