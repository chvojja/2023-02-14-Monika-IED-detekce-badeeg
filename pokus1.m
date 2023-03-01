 



%%
fp = '\\neurodata2\Large data\Monika 2p\VIP_tdT\534M\2p 20210831 534 F7\file_0006.h5';
fp = '\\neurodata3\Lab Neuro Optophys\Ca imaging\Monika 2p\spontánní záchvaty 2p\490\file_0044.h5';

fp = 'D:\temp\temp_monika sandbox\file_0044.h5';


%%
disp('Welcome')
disp('This is Jan Chvojka speaking...')
disp('This is the end of the speach. Thank you.')
fpn_allC = getFilepnAll('*.h5','Select .h5 files');

[signal,fs,chNames,start_dt] = loadh5ondrej(fp,1000);
Nsamples = size(signal,1);

%% detect the hobbits

det_settings = '-k1 10 -k2 10 -k3 0 -n 10000' ; %'-k1 5.50 -k2 5.50 -k3 0 -w 5000 -n 4000'
det_settings = '-k1 10 -k2 10 -k3 0 ' ; 
det_name = 'IED_jancaspike_monika1';

[lbl3] = jancaspike_detect_signal2lbl3struct(Signal = signal, Fs = fs, ChNames = chNames, FilePath = fp,  StartDT = start_dt, DetName = det_name , DetSettings =  det_settings ,  Color = [1 0.5 0.8] );

% center the detections
lbl3 = lbl3_center2min(signal,fs,lbl3);

%% save the bitchmens
if size( lbl3.lblSet , 1 ) > 0  % if detected something
% save the lbl3 structure
path_new = [char(lbl3.sigInfo.FilePath(1)) '\' det_name];
mkdir(path_new);
fpname = [path_new '\' filename(char(lbl3.sigInfo.FileName(1))) '-lbl3.mat'];
save(fpname, '-struct', 'lbl3');
end

 

