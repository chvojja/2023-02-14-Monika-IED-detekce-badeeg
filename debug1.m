clear all;
load('temp')
%%


s= downsample_byfs(s, fs,5000);
% hold on;
% plot(stretch_vector(s2',length(s))');
% plot(s)

%%





s_emg = stretch_vector(sum(    stft_mag(s, fs, 0.1,[400 800])'   ), length(s))';

q_thr = quantile( s_emg , 0.85);


emg_cand=s_emg>q_thr; % boolean of samples , candidtes

%%
emg_cand = emg_cand(:)';
emg_cand(1)=0;
emg_cand(end)=0;

enlarge_sec = 1.2;
min_emg_width_sec = 0.3;
short_protection_sec = 0.5;

emg_dil = imdilate(  emg_cand , ones(1,round(enlarge_sec*fs))  ); % enlarge

emg_close=imclose(emg_dil,  ones(1,round(short_protection_sec*fs))  ); % remove shorter than minAcceptLen


%
% hold on;
% 
% % plot current
% det_inds = out.det_inds;
% t_dets = zeros(1,length(s));
% t_dets(det_inds) = true;
% stem(s2t(s,fs),t_dets);

%
hold on;
plot(s2t(s,fs),[100*s'; s_emg'     ]' );

plot(  s2t(s,fs), q_thr * ones(size(s))  )

plot(  s2t(s,fs), 100*hfo_close )

plot(  s2t(s,fs), -100*emg_cand -100 )
plot(  s2t(s,fs), -100*emg_dil -200 )
plot(  s2t(s,fs), -100*emg_close -300 )




%%



s_emg_sorted = sort(s_emg);
range_percent = [0.1 0.5];

range_inds = round(range_percent * numel(s_emg));
range_vals = s_emg_sorted(range_inds);

s_emg_caped = s_emg;
s_emg_caped(s_emg>range_vals(2) | s_emg<range_vals(1) ) = NaN;
s_emg_caped = fillMissingValues(s_emg_caped);

medWin_sec = 10;
medWin_N = round(medWin_sec*fs);
s_emg_caped_med = medfilt1(s_emg_caped, medWin_N);

plot(s_emg_caped_med)

