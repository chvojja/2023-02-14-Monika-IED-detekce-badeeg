

% clear all
%% do not edit
disp('Welcome')
disp('This is Jan Chvojka speaking...')
disp('This is the end of the speach. Thank you.')
addpath('functions');

%% select .h5 files    = programmers can edit
%%fpn_allC = getFilepnAll('*.h5','Select .h5 files'); % for manual selection of a few .h5 files
% % Use this for a broad selection of all .h files in all subfolders in a directory instead


 dirContents_h5=dir(['\\neurodata2\Large data\Monika 2p\VIP_tdT WT mTOR\' '\**\*.h5']);
 fpn_allC=fullfile({dirContents_h5(:).folder},{dirContents_h5(:).name})';
% i = 90 viptdt 

 %dirContents_h5=dir(['\\neurodata2\Large data\Monika 2p\VIP_tdT\' '\**\*.h5']);
 %fpn_allC=fullfile({dirContents_h5(:).folder},{dirContents_h5(:).name})';
%  dirContents_h5=dir(['\\neurodata2\Large data\Monika 2p\VIP_tdT\397 F\2p 20220209 397 F1' '\**\*.h5']);
%  fpn_allC=fullfile({dirContents_h5(:).folder},{dirContents_h5(:).name})';
%  dirContents_h5=dir(['\\neurodata2\Large data\Monika 2p\VIP_tdT\397 F\2p 20220216 397 F2\' '\**\*.h5']);
%  fpn_allC=fullfile({dirContents_h5(:).folder},{dirContents_h5(:).name})';%1063

 %dirContents_h5=dir(['D:\temp\temp_monika sandbox\' '\**\*.h5']);

 %fpn_allC=fullfile({dirContents_h5(:).folder},{dirContents_h5(:).name})';

 %fpn_allC = {'D:\temp\temp_monika sandbox\file_0044.h5'}
%% settings  = Monika can edit
zesileni_AMsystemsZesilovac = 1000;
prisnost_detektoru = 7;  % 3 az 10
det_folder_name = ['IED_jancaspike_monika_' num2str(prisnost_detektoru) '_andHFOs'];
% R G B color of the labels
color = [1 0.5 0.8]; % barbie pink
%color = [1 0 0]; % rosy red

color_r = [0.1 0.8 0.2];
color_fr = [1 0.1 0.3];
hd = HFORMSdetector_staba_chvojka_simplified_v2;



%% do not edit
Nfiles = numel(fpn_allC);
Nfiles_detectedSomething = 0;
Nfiles_skipped =  0;
for i=53:Nfiles
    fp = fpn_allC{i};
    disp(['Processing file: ' fp])
    try
        [signal,fs,chNames,start_dt] = loadh5ondrej(fp,zesileni_AMsystemsZesilovac);
        Nsamples = size(signal,1);
        file_len_min = (Nsamples/fs)/60;
        minFileLength_min = 4;
        failedtoload = false;
        if file_len_min<minFileLength_min
            badl_length = true;

        else
            badl_length = false;
        end
    catch
        failedtoload = true;
        badl_length = false;
    end

    if badl_length 
        disp(['File skipped because it was shorter: ' num2str(file_len_min) ' than minimum of: ' num2str(minFileLength_min) ' seconds. Thank you for understanding.'])
           Nfiles_skipped = Nfiles_skipped  + 1;
    else
        if ~failedtoload
            disp(['File duration ' num2str(file_len_min) ])
    
    
    
    
    
            % hfos
            % select channels to detect
            chSelectedInd = find(contains(chNames,{'FCD','C','contra','lesion','Lesion','L'}));
    
            % detect
            hd.params = hd.paramverse.ripples.HFOpaper2023;
            hd.run(signal, fs, chSelectedInd);
            %hd.run(signal(:,chSelectedInd), fs, chSelectedInd);
            % replace the channel numbering to match original signal variable
            %channel_mapfun = @(x) chSelectedInd(x,1);
            %hd.outT.chan = channel_mapfun(hd.outT.chan);
    
            [lbl3_hfos_r] = LBL3.fill_by_jancadetections( Signal = signal,...
                                    Fs = fs,...
                                    Pos = hd.outT.pos,...
                                    Dur = hd.outT.dur,...
                                    Chan = hd.outT.chan,...
                                    ChNames = chNames,...
                                    FilePath = fp, ...
                                    StartDT = start_dt,...
                                    LabelClassName = 'ripples.HFOpaper2023',...
                                    ChannelMode = 'one',...
                                    LabelType = 'roi',...
                                    Color = color_r );
    
            % detect
            hd.params = hd.paramverse.fripples.HFOpaper2023;
            hd.run(signal, fs, chSelectedInd);
    
            %hd.run(signal(:,chSelectedInd), fs, chSelectedInd);
            % replace the channel numbering to match original signal variable
            %channel_mapfun = @(x) chSelectedInd(x,1);
            %hd.outT.chan = channel_mapfun(hd.outT.chan);
    
            [lbl3_hfos_fr] = LBL3.fill_by_jancadetections( Signal = signal,...
                                    Fs = fs,...
                                    Pos = hd.outT.pos,...
                                    Dur = hd.outT.dur,...
                                    Chan = hd.outT.chan,...
                                    ChNames = chNames,...
                                    FilePath = fp, ...
                                    StartDT = start_dt,...
                                    LabelClassName = 'fripples.HFOpaper2023',...
                                    ChannelMode = 'one',...
                                    LabelType = 'roi',...
                                    Color = color_fr );
             % ieds
              det_settings = ['-k1 ' num2str(prisnost_detektoru) ' -k2 ' num2str(prisnost_detektoru) ' -k3 0 ']; 
             [lbl3_ieds] = jancaspike_detect_signal2lbl3struct(Signal = signal, Fs = fs, ChNames = chNames, FilePath = fp,  StartDT = start_dt, DetName = det_folder_name , DetSettings =  det_settings ,  Color = color );
             % center the detections
             lbl3_ieds = lbl3_center2min(signal,fs,lbl3_ieds);
    
    
    
            % save the bitchmens
            %if size( lbl3_ieds.lblSet , 1 ) > 0  |  size( lbl3_hfos.lblSet , 1 ) > 0  % if detected something in at least one
            
            lbl3 = LBL3.merge([lbl3_ieds; lbl3_hfos_r; lbl3_hfos_fr; ]); %lbl3_hfos_r]);
        
            if size(lbl3.lblSet,1)>0
                status = LBL3.save2file_nexttosignal(lbl3,det_folder_name);
                Nfiles_detectedSomething = Nfiles_detectedSomething  + 1;
                disp(['Detected and saved']);
            else
                disp('nothing detected not saved');
            end
        end
    end
end
rmpath('functions');
disp(['Finished, the lbl3 labels are in a folder: ' det_folder_name] );
disp(['Skipped: ' num2str(Nfiles_skipped) ', found detections in: ' num2str(Nfiles_detectedSomething) ', analyzed total: ' num2str(Nfiles) ' files.']);