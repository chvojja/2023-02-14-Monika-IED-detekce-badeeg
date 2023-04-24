function lbl3 = lbl3_center2min(signal,fs,lbl3)
% this only works with LabelType = "point"


dt_sec = 0.06; % margin around detection for centering
%dt_sec = 0;
for i = 1:size(lbl3.lblSet,1)
    ch_num = lbl3.lblSet.Channel(i);
    start_rel_sec = 24*3600*datenum(  lbl3.lblSet.Start(i)-lbl3.sigInfo.SigStart(ch_num)    );
    %end_rel_sec = 24*3600*datenum(  lbl3.lblSet.End(i)-lbl3.sigInfo.SigStart(ch_num)    );
    sI = round( fs*(start_rel_sec-dt_sec) );
    eI = round( fs*(start_rel_sec+dt_sec) );
    
    sig_ied = signal(sI:eI,ch_num);

%     plot(sig_ied)
%     pause

    % process the detection and center it around minimum
    sig_ied = filtfilt([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]/15,1,sig_ied);

    Nsignal = length(sig_ied);
    Nfade = round(Nsignal/5);
    sig_ied_w = -sig_ied.*fadeinoutwin(Nsignal,Nfade,@blackman)';

    [pks,locs,w,p] = findpeaks(sig_ied_w,'SortStr','descend');
    if numel(locs)>1
        offsetI = Nsignal/2-locs(1);
    else
        offsetI =0;
    end

%     % center the original labels for roi
%     lbl3.lblSet.Start(i)=lbl3.lblSet.Start(i)-offsetI/fs/24/3600;
%     lbl3.lblSet.End(i)=lbl3.lblSet.End(i)-offsetI/fs/24/3600;

    % center the original labels for point
    lbl3.lblSet.Start(i)=lbl3.lblSet.Start(i)-offsetI/fs/24/3600;
    lbl3.lblSet.End(i)=lbl3.lblSet.Start(i);

end

function y = fadeinoutwin(L,Lfades,fun)
%FADEINOUT Summary of this function goes here
% L ... length
% Lfades ... length of fade in - outs

wb = fun(2*Lfades)';

y = ones(1,L);
y(1:Lfades) = wb(1:Lfades);
y(end-Lfades+1:end) = wb(Lfades+1:end);
end


end