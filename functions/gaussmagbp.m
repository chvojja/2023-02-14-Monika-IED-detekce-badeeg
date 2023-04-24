function y = gaussmagbp(f_3dB,fs, Ns)

% works only fs is even 
gwa = gaussup( round(Ns*f_3dB(2)/fs) , round(Ns*(f_3dB(2)-f_3dB(1))/fs)  );
Ngwup = numel(gwa);
Nfsh = Ns/2+1;
Ntop = round((Ns*(f_3dB(3)-f_3dB(2))/fs) );
gwd = gaussdown( Nfsh-Ngwup-Ntop ,  round(Ns*(f_3dB(4)-f_3dB(3))/fs) );

mag = [  gwa  ones(1,Ntop)  gwd  ];

mag_filter = [fliplr(mag)  mag(2:end-1)]';
y = mag_filter;

    function y = gaussdown(Nmax,samples_3db)
        alpha = 2.5;
        N=samples_3db*128/30;
        %N =128; 
        n = 0:Nmax-1;
        stdev = (N-1)/(2*alpha);
        y = exp(-1/2*(n/stdev).^2);
        %plot(y);

    end

    function y = gaussup(Nmax,samples_3db)
        y = flip( gaussdown(Nmax,samples_3db) );
    end

end