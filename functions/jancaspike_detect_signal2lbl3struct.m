function [lbl3]=jancaspike_detect_signal2lbl3struct(nv)
arguments
    nv.Signal   % required
    nv.Fs;      % required
    % optionals
    nv.FilePath = "";
    nv.ChNames = [];
    nv.StartDT = datetime(0,'ConvertFrom','datenum');
    nv.Subject = "UnknownAlien";
    nv.DetName = "JancaSpikeDetectorDefaultSettings";
    nv.DetSettings = '-k1 3.65 -k2 3.65 -k3 0';
    nv.Color = [1 0 0];
    nv.Value = 5;
end

[Nsamples,Nchannels] = size(nv.Signal);

% create artificial ChNames if not available
if isempty(nv.ChNames)
    nv.ChNames =  arrayfun(@num2str, 1:Nchannels, 'UniformOutput', 0);
end
nv.ChNames = string(nv.ChNames);


%%
% Create lbl3 tables
label_class_name = string( nv.DetName ); 
nv.FilePath = string( nv.FilePath );

% sigInfo
nv.StartDT.Format = 'dd-MM-yyyy HH:mm:ss' ;

[p,fn,ext] = fileparts(char(nv.FilePath));
FileName = string([fn  ext]);
FilePath = string([p]);
Subject = string( nv.Subject );
SigStart = nv.StartDT;
SigEnd = nv.StartDT + (Nsamples/nv.Fs)/(3600*24);
Fs = nv.Fs;

lbl3.sigInfo = [];
Nchan = numel(nv.ChNames);
for ich = 1:Nchan
    ChName = string( nv.ChNames{ich} );
    lbl3.sigInfo = [lbl3.sigInfo; table(FileName, FilePath, Subject, ChName, SigStart, SigEnd, Fs)    ];
end

% lblDef
ClassName = string( label_class_name );
ChannelMode = "one";
LabelType = "roi";
LabelType = "point";
Color = string( num2str(nv.Color) );
lbl3.lblDef = table(ClassName, ChannelMode, LabelType, Color);
lbl3.lblDef.ChannelMode = categorical(lbl3.lblDef.ChannelMode);
lbl3.lblDef.LabelType = categorical(lbl3.lblDef.LabelType);


% Detect IEDs

out = spike_detector_hilbert_v23(nv.Signal,nv.Fs,nv.DetSettings); 
Ndets = size( out.pos , 1);
%
% add to lblSet
ClassName = categorical(repmat(label_class_name,Ndets,1)); % categorical( repmat({ label_class_name},Ndets,1) );
Channel = int16(out.chan);
pos_dt = nv.StartDT + out.pos/(3600*24);
pos_dt.Format = 'dd-MM-yyyy HH:mm:ss' ;

Start = pos_dt-0.012/(3600*24);
End = pos_dt-0.002/(3600*24);
Value = nv.Value*ones(Ndets,1);
Comment = categorical(NaN(Ndets,1));
Selected = true( Ndets,1 );
ID = int64([1:Ndets]');
SignalFile = categorical( repmat(nv.FilePath,Ndets,1) );

lbl3.lblSet = table(ClassName, Channel, Start, End, Value, Comment, Selected, ID, SignalFile );


end
