function [lbl3]=filllbl3(nv)
arguments
    % Requied
    nv.Pos
    nv.Dur
    nv.Chan

    % optionals

    % lblDef
    nv.LabelClassName = 'UnknownLabelType',...
    nv.ChannelMode = "one";
    nv.LabelType = "roi";
    nv.Color = [1 0 0];

    % sigInfo
    nv.Signal = []; 
    nv.Fs = 0; 
    nv.FilePath = "";
    nv.ChNames = [];
    nv.StartDT = datetime(0,'ConvertFrom','datenum');
    nv.Subject = "UnknownAlien";
    
    % lblSet
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

% lblDef
ClassName = string( nv.LabelClassName );
Color = string( num2str(nv.Color) );
ChannelMode = string( nv.ChannelMode );
LabelType = string( nv.LabelType );
lbl3.lblDef = table(ClassName, ChannelMode, LabelType, Color);
lbl3.lblDef.ChannelMode = categorical(lbl3.lblDef.ChannelMode);
lbl3.lblDef.LabelType = categorical(lbl3.lblDef.LabelType);

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





% lblSet
Ndets = length(nv.Pos);
ClassName = categorical(repmat(ClassName,Ndets,1)); % categorical( repmat({ label_class_name},Ndets,1) );
Channel = int16(nv.Chan);
pos_dt = nv.StartDT + nv.Pos/(3600*24);
pos_dt.Format = 'dd-MM-yyyy HH:mm:ss' ;

switch nv.LabelType
    case "roi"
        Start = pos_dt;
        End = pos_dt + nv.Dur/(3600*24);
    case "point"
        Start = pos_dt;
        End = pos_dt;

end


Value = nv.Value*ones(Ndets,1);
Comment = categorical(NaN(Ndets,1));
Selected = true( Ndets,1 );
ID = int64([1:Ndets]');
SignalFile = categorical( repmat(nv.FilePath,Ndets,1) );

lbl3.lblSet = table(ClassName, Channel, Start, End, Value, Comment, Selected, ID, SignalFile );



end
