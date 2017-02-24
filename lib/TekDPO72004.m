function [varargout] = TekDPO72004(h,cmdStr,parmCell)
% TekDPO72004 Scope control programing
%
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 

switch cmdStr
    
	case 'start'
		set(h,'EOSMode','read&write');
		set(h,'InputBufferSize',10000000);
		set(h,'OutputBufferSize',10000000);
		fopen(h);
        varargout{1} = query(h,'*IDN?');
        fprintf(h,'ACQuire:SAMPlingmode RT');
        fprintf(h,'HORizontal:MODE MANual');
        fprintf(h,['HORizontal:MODE:RECOrdlength ' num2str(parmCell{1})]);
        fprintf(h,['HORizontal:MODE:SAMPLERate ' num2str(parmCell{2})]);
        fprintf(h,'ACQuire:STOPAfter RUNSTop');
        fprintf(h,'ACQuire:STATE RUN');  pause(1.0);
		fprintf(h,'ACQuire:STATE STOP');
        fprintf(h,'ACQuire:STOPAfter sequence');
        fprintf(h,'ACQuire:STATE 1');
        
    case 'stop'
        fclose(h);
		
	case 'readWaveform'
        fprintf(h,'ACQuire:STATE 1');
        nch = length(parmCell)-1;
        for ii = 1:nch
            fprintf(h,['SELect:' parmCell{ii} ' ON']);
            fprintf(h,[parmCell{ii} ':position 0']);
        end
        yscale = str2double(query(h,'WFMOutpre:YMULt?'));
        for ii = 1:nch
            fprintf(h,['DATa:SOUrce ' parmCell{ii}]);
            fprintf(h,':DATa:ENCdg RIBinary');
            fprintf(h,'WFMOutpre:BYT_Nr 1'); % 1 or 2
            fprintf(h,[':DATa:STARt ' num2str(1) ';STOP ' num2str(parmCell{end})]);
            fprintf(h,'WFMOutpre?');
            fprintf(h,'CURVe?');
            tmp_data = binblockread(h,'int8');
            varargout{ii} = tmp_data * yscale;
        end
        clrdevice(h)
        flushinput(h)
        
    case 'setdeskew'
		fprintf(h,['CH',num2str(parmCell{1}),':DESKEW ',num2str(parmCell{2})]);
        pause(0.1)
        
	case 'debug'
		disp('TekDPO72004 debug'); keyboard;
		
	otherwise
		disp('Un-supported command in TekDPO72004.m');
end

return