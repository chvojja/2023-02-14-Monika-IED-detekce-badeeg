

% clear all
%% do not edit
disp('Welcome')
disp('This is Jan Chvojka speaking again...')
disp('This is the end of the speach. Thank you.')
addpath('functions');

%% select .h5 files    = programmers can edit
%%fpn_allC = getFilepnAll('*.h5','Select .h5 files'); % for manual selection of a few .h5 files
% % Use this for a broad selection of all .h files in all subfolders in a directory instead

path_monika = ['\\neurodata2\Large data\Monika 2p\VIP_tdT WT mTOR\' '\**\*lbl3.mat'];
[folder_name_mice] = parserootpath_monika(str);

dirContents_h5=dir(path_monika);

fpn_allC=fullfile({dirContents_h5(:).folder},{dirContents_h5(:).name})';

% fp = fpn_allC{10};


%% settings  = Monika can edit
settings.lblnames.ied = 'IED_jancaspike_monika_7';
settings.lblFolderName = 'IED_jancaspike_monika_7_andHFOs';
settings.lblnames.hfor = 'ripples.HFOpaper2023';
settings.lblnames.hfofr = 'fripples.HFOpaper2023';

settings.ch_name_treatment = 'FCD';
settings.ch_name_ctrl = 'CTRL';
settings.ch_names_treatment_possible = {'FCD','lesion','Lesion','L'};
settings.ch_names_ctrl_possible= {'contra','C','lat'};

%% do not edit
Nfiles = numel(fpn_allC);

%  Gather table with label files
FilePathName = fpn_allC;
LabelFolderName = cell(Nfiles,1);
ID = zeros(Nfiles,1);
Fnumber = NaN(Nfiles,1);
MouseNumber = NaN(Nfiles,1);
RootFolderNameMice = cell(Nfiles,1);

for i=1:Nfiles

    file_info_monika = parsefilepath_monika( fpn_allC{i} ,folder_name_mice);
    
    LabelFolderName{i} = file_info_monika.label_folder_name;
    Fnumber(i) = file_info_monika.Fnumber;
    MouseNumber(i) = file_info_monika.mouse_number;
    RootFolderNameMice{i} = folder_name_mice;
    ID(i) = i;

end

T_label_files = table(ID,MouseNumber,Fnumber,LabelFolderName,FilePathName,RootFolderNameMice);
T_label_files = categorify(T_label_files);

%% Load info from label files
% NumIED_TREAT = NaN(Nfiles,1);
% NumFR = NaN(Nfiles,1);
% NumR = NaN(Nfiles,1);
for i=1:Nfiles

    fpname = char( T_label_files.FilePathName(i) );
    lbl3 = load(fpname);
    
    % go through channels in lblSet
    chnums  = double(unique(lbl3.lblSet.Channel));
    chnums=chnums(:)';
    for chnum=chnums
        % get lesion or contra
        chName =  lbl3.sigInfo.ChName( chnum );
        if contains(chName,settings.ch_names_treatment_possible)
            chName = settings.ch_name_treatment;
        elseif contains(chName,settings.ch_names_ctrl_possible)
            chName = settings.ch_name_ctrl;
        else

            continue % skip channel
        end
        
        % count detections for each type
        t1 = lbl3.lblSet(lbl3.lblSet.ClassName == settings.lblnames.ied & lbl3.lblSet.Channel==chnum , : );
        N = size(t1,1);
        T_label_files.(['ied_' chName])(i) = N;

        % count detections for each type
        t1 = lbl3.lblSet(lbl3.lblSet.ClassName == settings.lblnames.hfor & lbl3.lblSet.Channel==chnum , : );
        N = size(t1,1);
        T_label_files.(['hfor_' chName])(i) = N;

        % count detections for each type
        t1 = lbl3.lblSet(lbl3.lblSet.ClassName == settings.lblnames.hfofr & lbl3.lblSet.Channel==chnum , : );
        N = size(t1,1);
        T_label_files.(['hfofr_' chName])(i) = N;
    end


end


%%
%  Parse and plot label files
mice_numbers = unique(T_label_files.MouseNumber);
Nmice = length(mice_numbers);

T_mice = table;
for i=1:Nmice
    m_number=mice_numbers(i);
    T_label_mouse = T_label_files(T_label_files.MouseNumber == m_number & T_label_files.LabelFolderName == settings.lblFolderName ,  :);
    f_numbers = unique(T_label_files.Fnumber);
    %T_label_mouse = sortrows(T_label_mouse,'Fnumber');
    f_numbers = sort(f_numbers);
    f_numbers = f_numbers(:)';
    %T_f = table;
    T_f = table('Size',[0,10],'VariableTypes',repmat({'double'},1,10),'VariableNames',{'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10'});
    for fnum = f_numbers
        

        %selected_label  = settings.lblFolderName;
        selected_label  ='IED_jancaspike_monika_7';
        group_name = 'ied_FCD';
        T_f{group_name,fnum} = sum( T_label_files{ T_label_files.MouseNumber == m_number & T_label_files.LabelFolderName == selected_label & T_label_files.Fnumber == fnum, group_name} );
        group_name = 'ied_CTRL';
        T_f{group_name,fnum} = sum( T_label_files{ T_label_files.MouseNumber == m_number & T_label_files.LabelFolderName == selected_label & T_label_files.Fnumber == fnum, group_name} );
    
        selected_label  = settings.lblFolderName;
        %selected_label  ='IED_jancaspike_monika_7';
        group_name = 'hfor_FCD';
        T_f{group_name,fnum} = sum( T_label_files{ T_label_files.MouseNumber == m_number & T_label_files.LabelFolderName == selected_label & T_label_files.Fnumber == fnum, group_name} );
        group_name = 'hfor_CTRL';
        T_f{group_name,fnum} = sum( T_label_files{ T_label_files.MouseNumber == m_number & T_label_files.LabelFolderName == selected_label & T_label_files.Fnumber == fnum, group_name} );


        selected_label  = settings.lblFolderName;
        %selected_label  ='IED_jancaspike_monika_7';
        group_name = 'hfofr_FCD';
        T_f{group_name,fnum} = sum( T_label_files{ T_label_files.MouseNumber == m_number & T_label_files.LabelFolderName == selected_label & T_label_files.Fnumber == fnum, group_name} );
        group_name = 'hfofr_CTRL';
        T_f{group_name,fnum} = sum( T_label_files{ T_label_files.MouseNumber == m_number & T_label_files.LabelFolderName == selected_label & T_label_files.Fnumber == fnum, group_name} );
            
    end

    T_mice{i,'MouseNumber'} = m_number;
    T_mice(i,'Table_F') = {T_f};

    figure; 
  
    feature_name = 'ied';
    feature_name_TREAT = [feature_name '_FCD'];
    feature_name_CTRL = [feature_name '_CTRL'];
    max_y = max(   max(T_f{feature_name_TREAT,:}) , max(T_f{feature_name_CTRL,:})     ) +eps;
    subplot(6,1,1);
    hold on;
    title(['Stats for mouse number: ' num2str(m_number) ]);

    plot_area(1:10,T_f{feature_name_TREAT,:},{'FaceColor', 'r', 'EdgeColor', 'none'})
    set(gca,'YLim',[0 max_y]);
    ylabel(feature_name_TREAT, 'Interpreter', 'latex');
    subplot(6,1,2);
    plot_area(1:10,T_f{feature_name_CTRL,:},{'FaceColor', 'b', 'EdgeColor', 'none'})
    set(gca,'YLim',[0 max_y]);
   % ylabel(feature_name_CTRL);
    ylabel(feature_name_CTRL, 'Interpreter', 'latex');



    feature_name = 'hfor';
    feature_name_TREAT = [feature_name '_FCD'];
    feature_name_CTRL = [feature_name '_CTRL'];
    max_y = max(   max(T_f{feature_name_TREAT,:}) , max(T_f{feature_name_CTRL,:})     );
    subplot(6,1,3);
    plot_area(1:10,T_f{feature_name_TREAT,:},{'FaceColor', 'r', 'EdgeColor', 'none'})
    set(gca,'YLim',[0 max_y]);
    ylabel(feature_name_TREAT, 'Interpreter', 'latex');
    subplot(6,1,4);
    plot_area(1:10,T_f{feature_name_CTRL,:},{'FaceColor', 'b', 'EdgeColor', 'none'})
    set(gca,'YLim',[0 max_y]);
   % ylabel(feature_name_CTRL);
    ylabel(feature_name_CTRL, 'Interpreter', 'latex');
    

    feature_name = 'hfofr';
    feature_name_TREAT = [feature_name '_FCD'];
    feature_name_CTRL = [feature_name '_CTRL'];
    max_y = max(   max(T_f{feature_name_TREAT,:}) , max(T_f{feature_name_CTRL,:})     );
    subplot(6,1,5);
    plot_area(1:10,T_f{feature_name_TREAT,:},{'FaceColor', 'r', 'EdgeColor', 'none'})
    set(gca,'YLim',[0 max_y]);
    ylabel(feature_name_TREAT, 'Interpreter', 'latex');
    subplot(6,1,6);
    plot_area(1:10,T_f{feature_name_CTRL,:},{'FaceColor', 'b', 'EdgeColor', 'none'})
    set(gca,'YLim',[0 max_y]);
   % ylabel(feature_name_CTRL);
    ylabel(feature_name_CTRL, 'Interpreter', 'latex');

    pause
    saveas(gcf, [ num2str(m_number) '.png'])
    close all;
    disp('done')


end

