%the excution starts
clear;
close all;
main();

% Description:
%Create the data stream for each Mudulation type 
%the Data Stream length depands on the number of bits in the symbol
function DataStreamArr = DataCreation(ModulationType)
    NumOfSymbols = 1000;
    DataStreamArr = randi([0 1], 1, ModulationType*NumOfSymbols);
end


% Description:
%map the bits to the corresponding symbols 
function MappedArr = Mapper(DataStream, Positions,NumOfSymbols)
    if length(DataStream) == NumOfSymbols 
        MappedArr = DataStream * 2 - 1;
    
    else
        MappedArr = zeros(1,NumOfSymbols);
        ModulationType = length(DataStream)/NumOfSymbols;
        for i = 1:ModulationType:length(DataStream)
            symbol = DataStream(i:i+ModulationType-1);  % every 2 bits is a symbol
            decimal_value = bin2dec(num2str(symbol)); % getting the decimal value of the symbol
            MappedArr(1,fix(i/ModulationType)+1) = Positions(1,decimal_value+1); % map to the corresponding symbol
        end
     end
end




% Description:
% Define the demapper function (get the data Stream from the received symbols)
% NOTE: After creating the channel, replace every mapped data by the mapped data after the channel
% Define the demapper function (get the data Stream from the received symbols)
% NOTE: After creating the channel, replace every mapped data by the mapped data after the channel
function DeMappedArr = DeMapper(ModulationType, RecievedArr ,Positions,NumOfSymbols)

    % Perform demapping based on the modulation type
    if ModulationType == length(RecievedArr)/NumOfSymbols
        DeMappedArr = fix(RecievedArr / 2 + 0.5); % -1->0 & 1->1

    else 
        DeMappedArr = zeros(1, NumOfSymbols * ModulationType);
        for i = 1:NumOfSymbols
            symbol_index = find(Positions == RecievedArr(1,i));
            TempList = dec2bin(symbol_index - 1, ModulationType); % Represent the symbol decimal value in n binary bits
            DeMappedArr((i - 1) * ModulationType + 1 : i * ModulationType) = TempList - '0'; % Append the bits to the demapped version as integers
        end
    end
end


% Description:
% Entry point for the code's main logic
function main()
    NumOfSymbols = 1000;
    %number of bit represents every modulation type
    BPSK = 1;
    PSK8 = 3;
    QPSK = 2;
    QAM16 = 4;
    
    % Define symbols positions for QPSK modulation
    
    BPSK_Positions=[1 ,-1];
    QPSK_Positions = [-1-1i , -1+1i , 1-1i , 1+1i];
    
    % Define symbols positions for 16-QAM modulation
    QAM16_Positions = [-3-3i , -3-1i , -3+3i , -3+1i , -1-3i , -1-1i , -1+3 , -1+1i , ...
                       3-3i , 3-1i , 3+3i , 3+1i , 1-3i , 1-1i , 1+3i , 1+1i];
    
    
    % Define symbols positions for 8-PSK modulation
    PSK8_Positions = zeros(1, 8);
    tolerance = 1e-15;
    for i = 1:8
        angle = (i - 1) * pi / 4;
        cos_value = cos(angle);
        sin_value = sin(angle);
        % Check if the cosine value is close to zero due to the inaccuracy
        if abs(cos_value) < tolerance
            cos_value = 0;
        end
        % Check if the sine value is close to zero due to the inaccuracy
        if abs(sin_value) < tolerance
            sin_value = 0;
        end
        PSK8_Positions(i) = cos_value + 1i * sin_value; % Store complex value
    end
    
    
    
    
    BPSK_DataStream = DataCreation(BPSK);
    BPSK_Mapped = Mapper(BPSK_DataStream,BPSK_Positions,NumOfSymbols);
    BPSK_DeMapped = DeMapper(BPSK,BPSK_Mapped, BPSK_Positions,NumOfSymbols);
    disp(isequal(BPSK_DataStream, BPSK_DeMapped));
    
    QPSK_DataStream = DataCreation(QPSK);
    QPSK_Mapped = Mapper(QPSK_DataStream,QPSK_Positions,NumOfSymbols);
    QPSK_DeMapped = DeMapper(QPSK,QPSK_Mapped, QPSK_Positions,NumOfSymbols);
    disp(isequal(QPSK_DataStream, QPSK_DeMapped));
    
    PSK8_DataStream = DataCreation(PSK8);
    PSK8_Mapped = Mapper(PSK8_DataStream,PSK8_Positions,NumOfSymbols);
    PSK8_DeMapped = DeMapper(PSK8,PSK8_Mapped, PSK8_Positions,NumOfSymbols);
    disp(isequal(PSK8_DataStream, PSK8_DeMapped));
    
    
    QAM16_DataStream = DataCreation(QAM16);
    QAM16_Mapped = Mapper(QAM16_DataStream,QAM16_Positions,NumOfSymbols);
    QAM16_DeMapped = DeMapper(QAM16,QAM16_Mapped, QAM16_Positions,NumOfSymbols);
    disp(isequal(QAM16_DataStream, QAM16_DeMapped));

end


