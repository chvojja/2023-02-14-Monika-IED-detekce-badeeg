function [label_folder_name] = parserootpath_monika(str)


idx = strfind(str,'*'); 
if ~isempty(idx)
    str = str(1:idx(1)-1);
end

str = regexpi(str, '(.*?)(\\*)$', 'tokens', 'once');
str = str{1};

[~,label_folder_name,~] = fileparts(str);





