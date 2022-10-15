function [result] = simprule(x,y)
%% Simpson's Rule Calculator (simprule): Numeric Integration with Simpson's Rule
% Takes inputs for x (independent variable) and y (dependent variable)
% Outputs one value for numeric integration
%
% Author: Caleb Bristol
% Collaborators: N/A
% Date: 08/31/21

%% Define Variable to Store Summation

sum = 0;

%% Run Numeric Integration

for i = 1:2:length(x)-2
    
    %Define delta variable to reduce clutter
    
    delta = (x(i+2) - x(i)) / 6;
    
    %Individual (non-composite) value for Simpson's Rule
    simpInd = delta * (y(i) + 4*y(i+1) + y(i+2));
    
    %Sum individual components to make Composite Simpson's Rule
    sum = sum + simpInd;
end


%% Return Integration Result

result = sum;

end

