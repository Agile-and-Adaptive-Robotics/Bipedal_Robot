function [output1,output2] = generatedFunction_minimizeFlx1_withReuse(input1, input2, input3)
% GENERATEDFUNCTION_MINIMIZEFLX1_WITHREUSE Compute function values and
% reuse previous values whenever possible.
% 
% [OUTPUT1,OUTPUT2] = generatedFunction_minimizeFlx1_withReuse(INPUT1,
% INPUT2, INPUT3) returns the function values computed at the inputs. The
% function returns the values stored in persistent variables if the inputs
% are equal to their previous values.
% 
% Auto-generated by prob2struct on 08-May-2024 02:45:29
% 
% 

% Previous input values.
persistent lastInput1;
persistent lastInput2;
persistent lastInput3;
% Previous output values.
persistent lastOutput1;
persistent lastOutput2;

if isequal(input1,lastInput1) && isequal(input2,lastInput2) && isequal(input3,lastInput3)
    % The function is being evaluated with the previous input values.
    % Return the stored outputs.
    output1 = lastOutput1;
    output2 = lastOutput2;
else
    % The function is being evaluated at new inputs.
    % Evaluate the function and store the current input and output values.
    [output1,output2] = minimizeFlx(input1, input2, input3);

    % Save the input and output values.
    lastInput1 = input1;
    lastInput2 = input2;
    lastInput3 = input3;
    lastOutput1 = output1;
    lastOutput2 = output2;
end

end