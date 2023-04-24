classdef HFORMSdetector_staba_chvojka_simplified_v2 < handle
    %PARAMSHFOSTABACHVOJKA 
    % only works with one dim signal 

    % Input parameters
    % (nv.Signal,nv.Fs) and current_params 
    %

    % Output parameters
    % This should output solely: the position of the detections in seconds , duration in seconds and channel number (index)
    % out.pos
    % out.dur
    % out.chan % multichannel not implemented
    
    
    properties
        paramverse
        params
        %out % this is structure similar to janca detector
        outT % this is table

        s
        fs
    end

    
    
    methods 
        function self = HFORMSdetector_staba_chvojka_simplified_v2()
          
%         arguments
%             nv.SignalCols = [];
%             nv.SignalRows = [];
%             nv.Fs;
%         end
% 
%         self.setSignal(nv)


        

        % Ripples - default        
        self.paramverse.ripples.default.n_std_rms=2.6;  %2.8
        self.paramverse.ripples.default.freq_bounds = [50 200];
        self.paramverse.ripples.default.rmsLen_ms=6;
        self.paramverse.ripples.default.join_gap_ms=15; %2
        self.paramverse.ripples.default.minAcceptLen_ms=20;  %ms %4.5
        self.paramverse.ripples.default.minPeaks = 9;
        self.paramverse.ripples.default.peakPromRatio = 0.25;  
        self.paramverse.ripples.default.fstops = [35 45 250 300];

        % Ripples - hfo paper
        self.paramverse.ripples.HFOpaper2023.n_std_rms=2.2;  %2.8 %D
        self.paramverse.ripples.HFOpaper2023.freq_bounds = [40 250];
        self.paramverse.ripples.HFOpaper2023.rmsLen_ms=5; %5
        self.paramverse.ripples.HFOpaper2023.join_gap_ms=15; %2 %15
        self.paramverse.ripples.HFOpaper2023.minAcceptLen_ms=18;  %ms %4.5
        self.paramverse.ripples.HFOpaper2023.minPeaks = 9;
        self.paramverse.ripples.HFOpaper2023.peakPromRatio = 0.25;  
        self.paramverse.ripples.HFOpaper2023.fstops = [35 45 250 300];
        
        
        % FastRipples - default        
        self.paramverse.fripples.default.n_std_rms=4.4;  %2.8
        self.paramverse.fripples.default.freq_bounds = [250 900];
        self.paramverse.fripples.default.rmsLen_ms=3;
        self.paramverse.fripples.default.join_gap_ms=1.5; %2
        self.paramverse.fripples.default.minAcceptLen_ms=5;  %ms %4.5
        self.paramverse.fripples.default.minPeaks = 13;
        self.paramverse.fripples.default.peakPromRatio = 0.25;
        self.paramverse.fripples.default.fstops = [250 300 800 900];
        
        % FastRipples - HFOpaper 
        self.paramverse.fripples.HFOpaper2023.n_std_rms=3;  %2.8  53.8
        self.paramverse.fripples.HFOpaper2023.freq_bounds = [250 900];
        self.paramverse.fripples.HFOpaper2023.rmsLen_ms=3;
        self.paramverse.fripples.HFOpaper2023.join_gap_ms=3; %2 5 2.2
        self.paramverse.fripples.HFOpaper2023.minAcceptLen_ms=4;  %ms %4.5
        self.paramverse.fripples.HFOpaper2023.minPeaks = 13;
        self.paramverse.fripples.HFOpaper2023.peakPromRatio = 0.25;
        self.paramverse.fripples.HFOpaper2023.fstops = [250 300 800 900];

  
        end
        
%         function self = setParams(self,params)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             self.params = params;
%         end

        function self = run(self,signal,fs,chnums) % Signal = signal, Fs = fs 
            self.outT = [];
            if nargin<4
                chnums=1:size(signal,2);
            end
            chnums = chnums(:);
            
            self.fs = fs;

            % assure even number of samples
            if mod(numel(signal), 2)
                self.s = signal(1:end-1,:);
            else
                self.s = signal; % each column is a signal
            end



           %b_disp = true;
          % sOne = self.s(100000:104999,2);
      
        %  Nsamples = size(self.s,1);
%            sOne = self.s(Nsamples/2 -1:end,2);
          %Nchans = size(self.signal,2);
          for chnum = chnums'

              new_s = self.s(:,chnum);
              new_fs = fs;
              fsMax = 5000;
              if self.fs > fsMax
                  new_s = downsample_byfs(self.s(:,chnum), self.fs,fsMax);
                  new_fs = fsMax;
              end

    
               out = self.RMSdetector_staba_chvojka_simplified_v2(new_s,new_fs,self.params.fstops, self.params.n_std_rms, self.params.freq_bounds , self.params.rmsLen_ms, self.params.join_gap_ms, self.params.minAcceptLen_ms, self.params.minPeaks, self.params.peakPromRatio);
              % if ~isempty(out.pos)
                   out.chan = chnum*ones(size(out.pos));
                   self.outT = [self.outT; struct2table(out)];
             %  end
          end



           disp('finished')
          % lbl3 = filllbl3(ClassName = label_class_name, Channel = int16(out.chan), Start, End)

        end


 
        function self = setSignal(self,nv)
            if ~isempty(nv.SignalCols)
                self.s = nv.SignalCols;
            else
                self.s=nv.SignalRows(:); % make sure s is column
            end
        end


        function y = filterfft2_1D(self,x,magnitude_vector)
        % From signal x of N points generates spectrum by fft2()
        % Then, it applies magnitude vector of N points on the spectrum.
        
        
       
        % num_rows = size(x,1);
        % num_cols = size(x,2);
        %[X,Y] = meshgrid(1:num_cols,1:num_rows);
        freq_domain = fft2(x);
        freq_domain_shifted=fftshift(freq_domain);
        freq_pass_window = ones(size(x));
        % freq_pass_window_center_x = floor(size(freq_pass_window,2)/2)+1;
        % freq_pass_window_center_y = floor(size(freq_pass_window,1)/2)+1;
        
        
        freq_pass_window = freq_pass_window.*magnitude_vector;
        
        % plotbode(Magnitude = mag_filter(2500:end), Fs =5000);
        % pause
        
        windowed_freq_domain_shifted = freq_domain_shifted.*freq_pass_window;
        adjusted_freq_domain = ifftshift(windowed_freq_domain_shifted);
        im_2 = ifft2(adjusted_freq_domain);
        y = im_2;
        
        end



        function out = RMSdetector_staba_chvojka_simplified_v2(self,s,fs,fstops, n_std_rms, freq_bounds , rmsLen_ms, join_gap_ms, minAcceptLen_ms, minPeaks, peakPromRatio, Nover, b_disp)
        arguments 
            self
            s;  % input signal each column is channel
            fs; % sampling freq

            % optional
            fstops = [35 45 250 300];
            n_std_rms = 3; % kind of sensitivity
            freq_bounds = [40 250];
            rmsLen_ms = 3;
            join_gap_ms = 15;
            minAcceptLen_ms = 18;
            minPeaks = 8; % number of oscillations
            peakPromRatio = 0.25; %0.25; % oscillatorines

            % optional default
            Nover = 20; % default xcorr oversampling
            b_disp = false; %true; %false;
        end
        out = [];

    
        if ~isempty(freq_bounds)
            f_max = freq_bounds(2);
            f_min=freq_bounds(1);
        end
        
        %% filtrace 
        %fs = 5000;
        magFilter = gaussmagbp(fstops, fs,numel(s));
        s = s(1:length(magFilter));
    
        fd = self.filterfft2_1D(s, magFilter ) ;
        %%
%         plot(s);
%         hold on;
%         plot(fd)
        
        %%
        rms_length=round(rmsLen_ms*10^-3*fs); % delka segmentu = 3ms ze ktere se pocita rms
        rms_fd=zeros(size(fd));
        
        
        seg_down=floor((rms_length-1)/2);
        seg_up=ceil((rms_length-1)/2);
        for i=1:size(s,1)
            seg_start=i-seg_down; %nastaveni zacatku segmentu
            seg_end=i+seg_up; %nastaveni konce segmentu
            if seg_start<1 %korekce u zacatku signalu
                seg_start=1;
            end
            if seg_end>size(s,1) %korekce u konce signalu
                seg_end=size(s,1);
            end
            rms_fd(i,:)=sqrt(  mean(fd(seg_start:seg_end,:).^2));   % rms value of filtred signal fd of the size rms_length 
        end

        %% correction for emg artefacts

        s_emg = stretch_vector(sum(    stft_mag(s, fs, 0.1,[400 800])'   ), length(s))';
        
 
s_emg_sorted = sort(s_emg);
range_percent = [0.05 0.3];  % [0.1 0.5]; 

range_inds = round(range_percent * numel(s_emg));
range_vals = s_emg_sorted(range_inds);

s_emg_caped = s_emg;
s_emg_caped(s_emg>range_vals(2) | s_emg<range_vals(1) ) = NaN;
s_emg_caped = fillMissingValues(s_emg_caped');

medWin_sec = 5;
medWin_N = round(medWin_sec*fs);
emg_thr = 1.5*medfilt1(s_emg_caped, medWin_N)';


%         q_thr = quantile( s_emg , 0.8);
        
        
        emg_cand=s_emg>emg_thr; % boolean of samples , candidtes
        
        emg_cand = emg_cand(:)';
        emg_cand(1)=0;
        emg_cand(end)=0;
        
        enlarge_sec = 1.2;
        short_protection_sec = 0.5;
        
        emg_dil = imdilate(  emg_cand , ones(1,round(enlarge_sec*fs))  ); % enlarge
        
        emg_close=imclose(emg_dil,  ones(1,round(short_protection_sec*fs))  ); % remove shorter than minAcceptLen
        emg_close = logical(emg_close);
        rms_fd(emg_close)=NaN;
                   
        %% prahovani
        mean_rms_fd=nanmean(rms_fd); % mean rms of the whole filtered signal
        std_rms_fd=nanstd(rms_fd); % mean std of the whole filtered signal
       % rms_fd(isnan(rms_fd))=0;
       rms_fd_b = rms_fd ;

       if numel(find(~isnan(rms_fd)))<medWin_N
           rms_fd = zeros(size(rms_fd));
       end
      
        rms_fd = fillMissingValues(rms_fd)';
%plot(rms_fd)

       % threshold_rms_fd=repmat((mean_rms_fd+n_std_rms*std_rms_fd),size(fd,1),1);
        medWin_sec = 5;
        medWin_N = round(medWin_sec*fs);
        
       mean_rms_fd = 0.8*medfilt1(rms_fd, medWin_N)  + 0.2*mean_rms_fd;
     %   mean_rms_fd = movmean(rms_fd, medWin_N) ;
        std_rms_fd = movstd(rms_fd,medWin_N);
         threshold_rms_fd=(mean_rms_fd+n_std_rms*std_rms_fd);

       
        hfo_cand=rms_fd_b>threshold_rms_fd; % boolean of samples , HFO candidtes
        
        %% kontrola delky, kratke se vyhodi, blizke spoji
        join_gap_samples=round(join_gap_ms*10^-3*fs); % gap in indexes
        minAcceptLen_samples=round(minAcceptLen_ms*10^-3*fs); % min len in indexes
        
        hfo_cand(1)=0;
        hfo_cand(end)=0;
        
        hfo_open=imopen(hfo_cand,strel('line',round(fs/f_max),90)); % remove shorter than 1 period of fmax
        hfo_close=imclose(hfo_open,strel('line',join_gap_samples,90)); % join closer than join_gap
        hfo_close=imopen(hfo_close,strel('line',minAcceptLen_samples,90)); % remove shorter than minAcceptLen
%         
%          hold on;
% %         scc=0.2;
% %         plot(scc*fd);
%          plot(90*emg_close)
%          plot(emg_thr);
%          plot(s_emg);
%          plot(200*s+800)
%          plot(rms_fd);
%          plot( threshold_rms_fd)
% %         plot(scc*s);
% %         plot(hfo_cand)
%          plot(168*hfo_close)
%          disp('')
% %         
%         medWin_sec = 1;
%         medWin_N = round(medWin_sec*fs);
%         plot( medfilt1(rms_fd, medWin_N) );
%         disp('')
        % close(gcf); 
           
        % v hfo_close jsou všechny nad trhresholdem, s minlength a spojeny vedlejsi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         subplot(2,1,1);
%         hold on;
%         plot(0.4*hfo_close,'k');
%         plot(rms_fd,'r')
%         plot(fd,'b')
%         plot(threshold_rms_fd,'k')
%         if sum(hfo_cand)~=0
%         nonzerohfoc=find(hfo_cand);
%         tt=zeros(size(fd)); tt(nonzerohfoc(1):nonzerohfoc(1)+minAcceptLen)=max(fd);
%         plot(tt,'g');
%         end
%         subplot(2,1,2);
%         plot(s);
%         close(gcf);
        
        
        if ~isempty(hfo_close)
        
        %% 
        % abs_fd=abs(fd);
        % m_abs_fd=mean( abs_fd ); % rectified filtered
        % s_abs_fd=std( abs_fd );
        %  % pocet vrcholku ktere maji byt vyssi jak ntimes std
        % vector_threshold2=repmat((m_abs_fd+n_std_peak*s_abs_fd),size(fd,1),1);
        
        %% Measuring frequency
        
        % onset_offsets=get_OOI(double(hfo_close));
        [ detections_len , onsetsInd, offsetsInd ] = s2intervals( double(hfo_close) );
        onset_offsets = [ onsetsInd offsetsInd  ];
        
        count_detections=numel( onsetsInd );
        hfo_freqs=zeros(count_detections,1);
        hfo_pwr=zeros(count_detections,1);
        passedInds = true(count_detections,1);
        peaksIndsC=cell(count_detections,1); %false(count_OOI,length(fd));
        figdataC=cell(count_detections,1);
        
        failedOnC=cell(count_detections,1); % debugging info
        diagC=cell(count_detections,1); % debugging info
        

        for i=1:count_detections
            start=onset_offsets(i,1); stop=onset_offsets(i,2);

            hfo_pwr(i) = mean( rms_fd( start : stop ) ); 


            % start autocorelation freq detection
            xOver=fd(start:stop);
            xOver=xOver-mean(xOver); 
        
            % we oversample fs
            xOver = interpft(xOver,Nover*numel(xOver));
            fsOver = Nover*fs;
            
            if ~isempty(freq_bounds)
                [fp,fi,~,p]=findpeaks(xcorr(xOver,xOver),'SortStr','descend','MinPeakDistance', round(fsOver/f_max) ); %,'MinPeakWidth',widthSec);
            else
                [fp,fi,~,p]=findpeaks(xcorr(xOver,xOver),'SortStr','descend'); %
            end
            [fi_linear,I] = sort(fi,'ascend');
            fp_linear = fp(I);
            p_linear = p(I);
            fi_linear_middleInd = floor(numel(fi)/2)+1; % index of 
        
            % ultimate condition
            if numel(fp)>3
                % conditions
                if numel(fp)<minPeaks  % if test number of peaks
                   passedInds(i)=false; % this HFO failed
                   failedOnC{i} = [ failedOnC{i} '   '  'minPeaks ' num2str(numel(fp)) ' expected: ' num2str(max(minPeaks,2))];
                end
                diagC{i} = [ diagC{i} '   '  'minPeaks ' num2str(numel(fp)) ' expected: ' num2str(max(minPeaks,2))];
            
                if ~isempty( peakPromRatio ) % if test prominance
                   prominanceSecondVsFirst=p_linear( fi_linear_middleInd+1 ) /p(1);
                   if prominanceSecondVsFirst<peakPromRatio
                      passedInds(i)=false; % this HFO failed
                      failedOnC{i} = [ failedOnC{i} '   '  'peakPromRatio ' num2str(prominanceSecondVsFirst) ' expected: ' num2str(peakPromRatio) ];
                   else
                       diagC{i} = [ diagC{i} '   '  'peakPromRatio ' num2str(prominanceSecondVsFirst) ' expected: ' num2str(peakPromRatio) ];
                   end 
                end
                
                % measure freq
                hfo_freqs(i) = fsOver/min(abs(fi(2:end)-fi(1)));
                if ~isempty(freq_bounds)
                    if hfo_freqs(i)<f_min || hfo_freqs(i)>f_max
                       passedInds(i)=false; % this HFO failed
                       failedOnC{i} = [ failedOnC{i} '   '  'f bounds ' num2str( hfo_freqs(i) ) ];
                    else
                        diagC{i} = [ diagC{i} '   '  'f bounds ' num2str( hfo_freqs(i) ) ];
                    end
                end

               
            
               
               % record oscillation peaks
               peaksIndsC{i} = [  start + round(fi_linear(1:fi_linear_middleInd)/Nover)   ];
            else
               passedInds(i)=false; % this HFO failed
            end
        
        

         % Visualization
         if b_disp %%&& selectedIDx(i)
        
            % plot all or only hfos that passed
            if passedInds(i)
                b_plotit=true;
            else
                b_plotit=false;
            end
            %b_plotit=true;
        
            if b_plotit
        
%                 figurefull;
%                 subplot(3,1,1);
%                 hold on;
%                 plot(0.4*hfo_close,'g');
%                 plot(0.2*hfo_cand,'k');
%                 plot(rms_fd,'r')
%                 plot(fd,'b')
%                 plot(threshold_rms_fd,'k')
%                 dd=zeros(size(threshold_rms_fd,1),1);
%                 dd(start:stop)=0.4;
%                 %plot(dd,'y')
%                 
%                 subplot(3,1,2);
%                 plot(s);
%                 
%                 subplot(3,1,3)
%                 findpeaks(xcorr(xOver,xOver),'SortStr','descend','MinPeakDistance', round(fsOver/f_max) )
%                 %plotstockwell(s,fs,1.8,[300 800],10,'linear')
%                
%             
            
                %findpeaks(xcorr(xOver,xOver),'SortStr','descend','MinPeakDistance', round(fsOver/f_max) )
            
                contextN = 5000;
                startCtx = start - contextN;
                stopCtx = stop + contextN;

                figurefull;
                subplot(3,1,1);
                hold on;
                plot(0.4*hfo_close(startCtx:stopCtx),'g');
                plot(0.2*hfo_cand(startCtx:stopCtx),'k');
                plot(rms_fd(startCtx:stopCtx),'r');
                plot(fd(startCtx:stopCtx),'b');
                plot(threshold_rms_fd(startCtx:stopCtx),'k');

                dd=zeros(size(threshold_rms_fd,1),1);
               % dd(start:stop)=0.4;
                %plot(dd,'y')
                
                subplot(3,1,2);
                plot(s(startCtx:stopCtx));
                
                subplot(3,1,3)
                findpeaks(xcorr(xOver,xOver),'SortStr','descend','MinPeakDistance', round(fsOver/f_max) )
                %plotstockwell(s,fs,1.8,[300 800],10,'linear')
                
                title( { [  'accepted: ' num2str( passedInds(i) )  '   freq '  num2str(hfo_freqs(i))  '   ' failedOnC{i}]  ,  diagC{i} });
            
    
                pause
            
                % save figure;
                fig_frame=getframe(gcf);
                figdataC{i}=fig_frame.cdata;
                close(gcf); 
        
            end
           
        
         else % no fig ata
             figdataC{i}=[];
         end
        
        end
        
        out=struct;
        
        % related to individual detections
        out.pos = onset_offsets(passedInds,1)/fs; % save detections in onset offset format
        out.dur = detections_len(passedInds,:)/fs;
   
        out.peaksIndsC = peaksIndsC(passedInds,:);
        out.freq=hfo_freqs(passedInds); % save frequency
        out.powerMeanRms=hfo_pwr(passedInds); % save frequency

        out.figdataC = figdataC(passedInds);


        %plot(s2t(s,fs),[s'; 0.01*s_emg     ] );
        %subplot(2,1,1); plot(s2t(s,fs),s); subplot(2,1,2); stem(s2t(s,fs),  stretch_vector(sum(    stft_mag(s, fs, 0.1,[400 800])'   ), length(s))  )
        %plot(rms_fd)
        
        else % no output
            out = [];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [ intervals_len , onsetsIdx, offsetsIdx ] = s2intervals(s)
        %S2INTERVALS 
        % This ignores unfinished intervals e.g non zeros at the beggining or end of the s signal
        
        % n = lengths of consecutive intervals
        % is = start indexes of intervals
        % ie = end indexes if intervals
        s = logical(s); % anything non zero is considered to be an interval
        dx = diff(s);
        
        %     if x(end)<=0
           onsetsIdx = find(dx>0)+1; 
           offsetsIdx = find(dx<0); 
        
           if ~isempty(onsetsIdx) && ~isempty(offsetsIdx) 
        
               if offsetsIdx(1) < onsetsIdx(1) % if first interval does not start from zero
                   offsetsIdx = offsetsIdx(2:end);
               end
            
               if offsetsIdx(end) < onsetsIdx(end) % if last interval  is not finished
                   onsetsIdx = onsetsIdx(1:end-1);
               end
           end
        
           intervals_len = offsetsIdx-onsetsIdx+1;
         
        end

        end %%%%%%%%%%%%%%%%%% end of rms detector fcn
end


end

