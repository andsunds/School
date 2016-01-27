%% This is a tets of the case-of-end_case structure in matlab
clc
dim=3;

switch logical(true)
    case dim>=3 
        disp(3)
        dim=dim-1;
    case dim>=2 
        disp(2)
    case dim>=1 
        disp(1)
end
    