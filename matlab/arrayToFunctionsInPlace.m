function arr = arrayToFunctionsInPlace(arr,symVar)    
    % Loop through the elements of the input array
    for i = 1:numel(arr)
        % Create a function handle for the symbolic expression
        func_handle = matlabFunction(arr(i), 'Vars', symVar);
        
        % Replace the symbolic expression with the function handle
        arr(i) = func_handle;
    end
end