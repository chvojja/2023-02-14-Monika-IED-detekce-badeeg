function out = parsefilepath_monika(path,folder_name_mice)
% path = '\\neurodata2\Large data\Monika 2p\VIP_tdT\\**\*.h5'; % example input
% path = '\neurodata2\Large data\Monika 2p\VIP_tdT';

 %fpn_allC=fullfile({dirContents_h5(:).folder},{dirContents_h5(:).name})';
%path = '\\neurodata2\Large data\Monika 2p\VIP_tdT WT mTOR\585F\2p 20210905 585 F1\IED_jancaspike_monika_7_andHFOs\file_0001-lbl3.mat';

out.mouse_number = NaN;
out.label_folder_name = '';
out.Fnumber = NaN;


tokens = regexpi(path, ['\\.*?(?<=' folder_name_mice ').*?\\(\d{3})\D'], 'tokens'); 
if ~isempty(tokens)
    out.mouse_number = str2num(tokens{1}{1});
end




[~, out.label_folder_name, ~] = fileparts(fileparts(path));

tokens = regexpi(path, ['[fF](\d+)\D*\s*\d*\\' out.label_folder_name], 'tokens');
if ~isempty(tokens)
    out.Fnumber = str2num(tokens{1}{1});
end



