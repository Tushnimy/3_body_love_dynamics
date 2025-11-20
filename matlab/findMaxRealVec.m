function maxReal = findMaxRealJac(complexVector, bool)
    % Extract the real parts of the complex numbers
    realParts = real(complexVector);
    if bool ==1
        maxReal = min(realParts);
    end
    if bool ==0
    % Find the maximum real part
        maxReal = max(realParts);
    end
end
