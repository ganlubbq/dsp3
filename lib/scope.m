% DESCRIPTION
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
% 
% Copyright 2015 default

function data = scope(Scope)

% Scope IP address such as 'TCPIP::10.71.103.173::INSTR'
rsrcName        = Scope.rsrcName;
vendor          = Scope.vendor;
InBuffersize    = Scope.InBuffersize;
OutBuffersize   = Scope.OutBuffersize;
SampleRate      = Scope.SampleRate;
Timeout         = Scope.Timeout;
Channel         = Scope.ChannelNo;
Points          = Scope.Points;
VerticalScale   = Scope.Vertical;

%Trigger Params.
Trigger.Type    = Scope.Trigger.Type;
Trigger.Level   = Scope.Trigger.Level;
Trigger.Source  = Scope.Trigger.Source;
Trigger.Time    = Scope.Trigger.Time;

switch vendor
    case 'NI'
    case 'Tektronix'
    case 'Keysight'
    otherwise
        error('unsupported instrument vendor in scope.m!!!');
end

data = zeros(length(Channel),Points); % output data from Scope

%% Interface configuration and instrument connection

% Find a VISA-TCPIP object.
Obj = instrfind('Type', 'visa-tcpip', 'RsrcName', rsrcName, 'Tag', '');

% Create the VISA-TCPIP object if it does not exist
% otherwise use the object that was found.
if isempty(Obj)
    Obj = visa(vendor, rsrcName);
else
    fclose(Obj);
    Obj = Obj(1);
end

% Set the buffer size
Obj.InputBufferSize = InBuffersize;
Obj.OutputBufferSize = OutBuffersize;

% Set the timeout value
Obj.Timeout = Timeout;

% Set the Byte order
Obj.ByteOrder = 'littleEndian';

fopen(Obj);
flushinput(Obj);
flushoutput(Obj);

instinfo = query(h,'*IDN?')

%% Instrument control and data retreival
% Now control the instrument using SCPI commands. refer to the instrument
% programming manual for your instrument for the correct SCPI commands for
% your instrument.

% Scope SETUP
fprintf(Obj,'*RST'); % Reset the instrument
fprintf(Obj,'AUTOSet EXECute'); % Autoscale
ChannelNum = length(Channel); 
for index = 1:ChannelNum
    fprintf(Obj,['SELect:CH' num2str(Channel(index)) ' ON']); % Specify Displayed Channel
end

% Horizontal Settings
%MANUAL mode can change the sample rate and the record length. the scale is read only.
fprintf(Obj,'HORizontal:MODE MANUAL'); 
fprintf(Obj,['HORIZONTAL:MODE:SCALE ' num2str(Points/10/SampleRate)]);
fprintf(Obj,['HORIZONTAL:MODE:SAMPLERATE ' num2str(SampleRate)]);
fprintf(Obj,['HORIZONTAL:MODE:RECORDLENGTH ' num2str(Points)]);

% Acquiring Settings
fprintf(h,'ACQuire:SAMPlingmode RT');
fprintf(h,'ACQuire:STOPAfter RUNSTop');
fprintf(h,'ACQuire:STATE RUN');  pause(1.0);
fprintf(h,'ACQuire:STATE STOP');
fprintf(h,'ACQuire:STOPAfter sequence');
fprintf(h,'ACQuire:STATE 1');

% Vertical Settings
for index = 1:ChannelNum
    fprintf(Obj,['CH' num2str(Channel(index)) ':SCALE ' num2str(VerticalScale(index))]);
end
fprintf(Obj,'CH1:DESKEW 0e-12');%-2210

if Scope.BandWidth < 33e9
    fprintf(Obj,['CH1:BANDWIDTH ' num2str(Scope.BandWidth)]);
    fprintf(Obj,['CH3:BANDWIDTH ' num2str(Scope.BandWidth)]);
    fprintf(Obj,'CH1:BANDWIDTH:ENHANCED AUTO');
    fprintf(Obj,'CH3:BANDWIDTH:ENHANCED AUTO');
    fprintf(Obj,'CH1:BANDWIDTH:ENHANCED:FORCE ON');
    fprintf(Obj,'CH3:BANDWIDTH:ENHANCED:FORCE ON');
end
    
% TRIGGER
% Set up trigger Type.
fprintf(Obj,'TRIGGER:A:TYPE PULSE');
fprintf(Obj, ['TRIGGER:A:PULSE:CLASS ' Trigger.Type]); 

% Set up trigger A level. 
fprintf(Obj, ['TRIGGER:A:LEVEL ' num2str(Trigger.Level)]);

% Set up trigger source (CH1/CH2/CH3/CH4).
fprintf(Obj, ['TRIGGER:A:PULSE:SOURCE ' Trigger.Source]); 

% fprintf(Obj,'TRIGGER:A:PULSE:WINDOW:WHEN OCCURS');
% fprintf(Obj,'TRIGGER:A:PULSE:WINDOW:TYPE EXITSWindow');
% fprintf(Obj,'TRIGGER:A:PULSE:WINDOW:THRESHOLD:HIGH 25e-3');
% fprintf(Obj,'TRIGGER:A:PULSE:WINDOW:THRESHOLD:LOW -25e-3');

fprintf(Obj, 'TRIGGER:A:PULSE:TIMEOUT:QUALIFY OCCURS');% Set up trigger qualifier (occurs / logic). 
fprintf(Obj, ['TRIGGER:A:PULSE:TIMEOUT:TIME ' num2str(Trigger.Time)]);% Set up trigger timeout time. 
fprintf(Obj, 'TRIGGER:A:PULSE:TIMEOUT:POLARITY STAYSLOW'); % Set up trigger polarity, high / low.

% Set up the horizontal trigger position.
fprintf(Obj,'HORIZONTAL:POSITION 2');
% tmp = query(Obj, 'TRIGGER:A:PULSE?')

% OUTPUT FORMAT
pause(0.28);
fprintf(Obj,'ACQuire:STATE OFF'); % STOP
for index = 1:ChannelNum
    fprintf(Obj,['DATa:SOUrce CH' num2str(Channel(index))]);
    fprintf(Obj,':DATa:ENCdg RIBinary');% ASCIi|FAStest|RIBinary|RPBinary|FPBinary|SRIbinary|SRPbinary|SFPbinary
    fprintf(Obj,'WFMOutpre:BYT_Nr 1');% 1 or 2
%   fprintf(Obj,'ACQuire:STATE OFF'); % STOP
%     horizLen = str2num(query(Obj,'HORIZONTAL:MODE:RECORDLENGTH?'));
    fprintf(Obj,[':DATA:START ' num2str(1) ';STOP ' num2str(Points)]);
    fprintf(Obj,'CURVE?');
    data1 = binblockread(Obj,'int8');
    data(index,:) = data1;
end

X_Zero = str2num(query(Obj,'WFMOUTPRE:XZERO?'));% get x zero
Sampling_Interval = str2num(query(Obj,'WFMOUTPRE:XINCR?'));% get the sampling interval 
Trigger_Point = str2num(query(Obj,'WFMOUTPRE:PT_OFF?'));% get the trigger point within the record
Y_Multiple = str2num(query(Obj,'WFMOUTPRE:YMULT?'));
Y_Offset = str2num(query(Obj,'WFMOUTPRE:YOFF?'));
Y_Zero = str2num(query(Obj,'WFMOUTPRE:YZERO?'));
data = Y_Multiple*(data - Y_Offset) - Y_Zero;
% fprintf(Obj,'HEADer On');% Enable header for query
% query(Obj,'WFMOutpre?')% Query the output setting
% fprintf(Obj,'HEADer Off');% disable header for query
fclose(Obj);
flushinput(Obj);
flushoutput(Obj);
delete(Obj);

return

