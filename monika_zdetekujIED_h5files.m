%% do not edit
disp('Welcome')
disp('This is Jan Chvojka speaking...')
disp('This is the end of the speach. Thank you.')
addpath('functions');

%% select .h5 files    = programmers can edit
fpn_allC = getFilepnAll('*.h5','Select .h5 files'); % for manual selection of a few .h5 files
% % Use this for a broad selection of all .h files in all subfolders in a directory instead
% dirContents_h5=dir(['\\neurodata2\Large data\Monika 2p\VIP_tdT' '\**\*.h5']);
% fpn_allC=fullfile({dirContents_h5(:).folder},{dirContents_h5(:).name})';

%% settings  = Monika can edit
zesileni_AMsystemsZesilovac = 1000;
prisnost_detektoru = 10;  % 3 az 10
det_name = ['IED_jancaspike_monika_' num2str(prisnost_detektoru)];
% R G B color of the labels
color = [1 0.5 0.8]; % barbie pink
%color = [1 0 0]; % rosy red


%% do not edit
Nfiles = numel(fpn_allC);
for i=1:Nfiles
    fp = fpn_allC{i};
    disp(['Processing file: ' fp])
    [signal,fs,chNames,start_dt] = loadh5ondrej(fp,zesileni_AMsystemsZesilovac);
    Nsamples = size(signal,1);

    det_settings = ['-k1 ' num2str(prisnost_detektoru) ' -k2 ' num2str(prisnost_detektoru) ' -k3 0 ']; 
    [lbl3] = jancaspike_detect_signal2lbl3struct(Signal = signal, Fs = fs, ChNames = chNames, FilePath = fp,  StartDT = start_dt, DetName = det_name , DetSettings =  det_settings ,  Color = color );
    
    % center the detections
    lbl3 = lbl3_center2min(signal,fs,lbl3);
    
    % save the bitchmens
    if size( lbl3.lblSet , 1 ) > 0  % if detected something
    % save the lbl3 structure
    path_new = [char(lbl3.sigInfo.FilePath(1)) '\' det_name];
    mkdir(path_new);
    fpname = [path_new '\' filename(char(lbl3.sigInfo.FileName(1))) '-lbl3.mat'];
    save(fpname, '-struct', 'lbl3');
    end
end
rmpath('functions');
disp(['Finished, the lbl3 labels are in a folder: ' det_name] );