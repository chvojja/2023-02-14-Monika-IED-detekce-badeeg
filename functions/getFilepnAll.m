function y = getFilepnAll(startpath, prompt)
% getFilepnAll  Get file path names
%   y ... cell array of path names to files
%   

[filen, filep] = uigetfile(startpath, prompt, 'MultiSelect', 'on');

if isa(filen, 'double') % when nothing selected, filen=filep=0
    disp('No files selected');
    y = [];
    return
end
y = fullfile(filep, filen);
if ~iscell(y), y = {y}; end
y=y';
end