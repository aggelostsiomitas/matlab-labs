function [local_max, local_min] = find_local(x, y)

    df = diff(y) ./ diff(x); 

    local_max_logical = (df(1:end-1) > 0) & (df(2:end) < 0);  
    local_min_logical = (df(1:end-1) < 0) & (df(2:end) > 0); 

    local_max_indexes = find(local_max_logical); 
    local_min_indexes = find(local_min_logical); 

  
    local_max = x(local_max_indexes + 1);  
    local_min = x(local_min_indexes + 1);  
end
