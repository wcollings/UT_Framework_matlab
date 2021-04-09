function update_lib( libfile, parameter, value )
%UPDATE_LIB takes in a specific libfile and a specified parameter within a
%subcircuit to have the indicated value. 
%   libfile contains the .lib file name. parameter is the parameter to 
%   update and value is the value to which parameter should be updated. The 
%   value 'value' must be a number with no characters or a string, i.e. the 
%   value 400 micro must be entered as 0.0004, 400e-6, or '400u'. Note that
%   '400u' is entered as a string as it contains a non-numerical character. 

fid=fopen(libfile,'r');                                                     % Open libfile file as readable
i = 1;                                                                      % Initialize index to count lines
cur_line = fgetl(fid);                                                      % Read in the first line of simfile
LIBtxt = {};
while ischar(cur_line)                                                      % Until the end of the .lib file is reached
    cur_line = strtrim(cur_line);                                           % Trim whitespace
    if isempty(cur_line)                                                    % If the current line is empty
        LIBtxt{i} = '';                                                     % Insert empty string into LIBtxt cell array
        
    else                                                                    % If there are some characters in the current line
        target = ['.param ', parameter];                                    % The target is '.param parameter'
        if regexpi(cur_line, target) == 1                                   % If the line begins with '.param parameter'

            cur_val = extractBetween(cur_line, '=', ';');                   % The current value is between the equals sign and the comment
            
            if isempty(cur_val)                                             % If there is no comment
                lib_line = strsplit(cur_line, '=');                         % Split at the equals sign
                new_line = [strtrim(lib_line(1)), '=', num2str(value)];     % Replace the old value with the new value
                
            else                                                            % If there is a comment
                lib_line = strsplit(cur_line, cur_val);                     % Split the string using the current value as the delimeter
                new_line = [lib_line(1), num2str(value), lib_line(2)];      % Insert the new value in place of the old value

            end
            new_line = join(new_line);                                      % Join the new line
            new_line = new_line{:};                                         % Extract new_line string from cell
            
            LIBtxt{i} = new_line;                                           % Insert the new string into the Libtxt cell array

        else                                                                % Otherwise, it is not the target line
            LIBtxt{i} = cur_line;                                           % Insert unmodified string into LIBtxt cell array
        end
    
    end
   
    cur_line = fgetl(fid);                                                  % Read in the next line of the .lib file
    i = i + 1;                                                              % Increment the line counter
end
fclose(fid);                                                                % Close readable file


fid = fopen(libfile, 'w');                                                  % Reopen model_asc as writable
for i = 1:numel(LIBtxt)                                                     % For each recorded row
    if i == numel(LIBtxt)                                                   % If the row is the last recorded row (second to the appended -1)
        fprintf(fid,'%s', LIBtxt{i});                                       % Write the row
        break                                                               % break out of the loop
    else                                                                    % For all other rows
        fprintf(fid,'%s\n', LIBtxt{i});                                     % Write the row and return
    end                                                                     %
end                                                                         %
fclose(fid);                                                                % Close the writable file
end
