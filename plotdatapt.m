function output_txt = plotdatapt(obj,event_obj)
% Display data cursor position in a data tip
% obj          Currently not used
% event_obj    Handle to event object
% output_txt   Data tip text, returned as a character vector or a cell array of character vectors

pos = event_obj.Position;

x = get(get(event_obj,'Target'),'XData');
y = get(get(event_obj,'Target'),'YData');

% Find index
index_x = find(x == pos(1));
index_y = find(y == pos(2));
index = intersect(index_x,index_y);

% Set output text
output_txt = {['X: ',num2str(pos(1),4)], ...
    ['Y: ',num2str(pos(2),4)], ...
    ['Index: ', num2str(index)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end

%***********************************************************%

function formattedValue = formatValue(value,event_obj)
% If you do not want TeX formatting in the data tip, uncomment the line below.
% event_obj.Interpreter = 'none';
if strcmpi(event_obj.Interpreter,'tex')
    valueFormat = ' \color[rgb]{0 0.6 1}\bf';
    removeValueFormat = '\color[rgb]{.25 .25 .25}\rm';
else
    valueFormat = ': ';
    removeValueFormat = '';
end
formattedValue = [valueFormat num2str(value,4) removeValueFormat];
