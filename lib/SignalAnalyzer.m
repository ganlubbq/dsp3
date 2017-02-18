classdef SignalAnalyzer < WDeviceHandle
    %SIGNALANALYZER Summary of this class goes here
    %   Detailed explanation goes here
    
    %   Copyright2011 default
    properties
        colormap        = 'hot'
        datamode        = 'real'
        plot_type       = '2D Line'
        signal_type     = 'ElectricalSignal'
		symrate         = 10e9
		sps             = 1
    end
    
    methods
        
        function this = SignalAnalyzer(varargin)
            SetVariousProp(this, varargin{:})
        end
        
        function reset(~)
        end

        function Output(this, x)
            % prepare the signal
			am = x;
            switch this.signal_type
                case 'ElectricalSignal'
                    % square-law and normalize for evelope detection
                    am_sqr = abs(am).^2;
%                     am_nlz = am_sqr / max(am_sqr);
                    am_nlz = am.' / max(abs(am));
                    % normalized field for scatter-plot
                    am_cmp = am.' / max(abs(am));
                case 'OpticalSignal'
                    % square-law, sum and normalize for evelope detection
                    am_sqr = sum(abs(am).^2);
                    am_nlz = am_sqr / max(am_sqr);
                    % normalized field for scatter-plot
                    am_cmp = am.' / max(max(abs(am)));
            end
            
            % creat an eye_scope and plot
            Rs = this.symrate;
            Ns = this.sps;
            if strcmp(this.datamode, 'real') && strcmp(this.signal_type, 'ElectricalSignal')
                eye = this.eyediag(Rs, Ns, 1.2, -1.2, this.datamode, this.plot_type);
                update(eye, am_nlz);
            elseif strcmp(this.datamode, 'real') && strcmp(this.signal_type, 'OpticalSignal')
                eye = this.eyediag(Rs, Ns, 1.2, -0.1, this.datamode, this.plot_type);
                update(eye, am_nlz);
            elseif strcmp(this.datamode, 'complex')
                eye = this.eyediag(Rs, Ns, 1.2, -1.1, this.datamode, this.plot_type);
                data = am_cmp(1:end);
                update(eye, data./max([max(real(data)) max(imag(data))]));
            end
            
            if strcmpi(this.plot_type,'2D Color')
                % specify the colormap of the eye_scope
                cmap = colormap(this.colormap);
                cmap(1,:) = [0 0 0];    % black backgroup
                plot(eye, cmap)
            end
            
%             if strcmpi(this.plot_type,'2D Line')
%                 plot(eye, 'k-');
%             end
            
            % plot the real signal in vector, note that both polarizations
            % of the signal have been treated independently
            h1 = scatterplot(am_cmp(1:end),1,0,'c-');
            hold on
            if Ns>1
                h2 = scatterplot(am_cmp(1:end),Ns,Ns/2-1,'o',h1);
            else
                h2 = scatterplot(am_cmp(1:end),Ns,Ns,'o',h1);
            end
            h22 = get(h2,'Children');
            h222 = get(h22,'Children');
            set(h222,'MarkerSize', 9, 'MarkerFaceColor','r', 'MarkerEdgeColor','g')
            title('Constellation')
            grid on
            hold off
                        
            % deal with spectrum
            Nfft = 2^nextpow2(length(am_cmp));
            SamplingRate = Rs * Ns;
            switch this.signal_type
                case 'OpticalSignal'
                    Amplitude = fftshift(fft(sum(am_cmp,2), Nfft)./Nfft);
                    frequency = SamplingRate *linspace(-0.5,0.5,Nfft);
                    strTitle = 'OSA output';
                case 'ElectricalSignal'
                    Amplitude = fft(sum(am_cmp,2), Nfft)./Nfft;
                    frequency = SamplingRate *linspace(0,0.5,Nfft/2);
                    Amplitude = Amplitude(1:end/2);
                    strTitle = 'ESA output';
            end
            h3 = figure;
            plot( frequency, 10*log10(abs(Amplitude).^2) ); 
            xlabel('Frequency (Hz)');
            ylabel('Power (dB)');
            title(strTitle)
            grid on
            
            % manage  windows
            managescattereyefig(h1,eye,'right')
            h2_p = get(h2,'Position');
            h3_p = get(h3,'Position');
            h3_p(1) = h2_p(1) + h2_p(3) + 10;
            h3_p(2) = h2_p(2);
            set(h3,'Position',h3_p)
        end
		
    end
    
    methods (Static)
	
        function eyeDiag = eyediag( symrate, nsamp, maxa, mina, mode, type)
            eyeDiag 						= commscope.eyediagram;
            eyeDiag.SamplesPerSymbol  		= nsamp;
            eyeDiag.SamplingFrequency 		= nsamp*symrate;
            eyeDiag.MaximumAmplitude  		= maxa;
            eyeDiag.MinimumAmplitude  		= mina;
            eyeDiag.PlotType          		= type;
            eyeDiag.AmplitudeResolution 	= 0.01;
            eyeDiag.NumberOfStoredTraces 	= 400;
            eyeDiag.ColorScale        		= 'Log';
            eyeDiag.OperationMode     		= mode;
            eyeDiag.RefreshPlot       		= 'on';
            eyeDiag.SymbolsPerTrace   		= 2;
        end
        
    end
    
end

