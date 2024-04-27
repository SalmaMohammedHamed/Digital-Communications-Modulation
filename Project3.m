%the excution starts
clear;
close all;


% Description:
%Create the data stream for each Mudulation type 
%the Data Stream length depands on the number of bits in the symbol
function DataStreamArr = DataCreation(ModulationType)
    NumOfSymbols = 1000;
    DataStreamArr = randi([0 1], 1, ModulationType*NumOfSymbols);
end


% Description:
%map the bits to the corresponding symbols 
function MappedArr = Mapper(ModulationType,DataStream, Positions,NumOfSymbols)
    if ModulationType == 1
        MappedArr = DataStream * 2 - 1;
    else
        MappedArr = zeros(1,NumOfSymbols);
        
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
function DeMappedArr = DeMapper(ModulationType, RecievedArr,NumOfSymbols)
    %number of bit represents every modulation type
    BPSK = 1;
    PSK8 = 3;
    QPSK = 2;
    QAM16 = 4;
    DeMappedArr = zeros(1,NumOfSymbols*ModulationType);
    % Perform demapping based on the modulation type
    if ModulationType == BPSK
        for j = 1 : length(RecievedArr)
            if (RecievedArr(j) > 0)
            DeMappedArr(j) = 1;
            else
            DeMappedArr(j) = 0;
            end
        end

    elseif ModulationType == QPSK
        TempList= [0 0];
        for j = 1 : NumOfSymbols
            if(real(RecievedArr(j))>=0 &&imag(RecievedArr(j))>0)
                TempList=[1, 1];
            elseif(real(RecievedArr(j))<=0 &&imag(RecievedArr(j))>0)
                TempList=[0, 1];
            elseif(real(RecievedArr(j))>=0 &&imag(RecievedArr(j))<0)
                TempList=[1, 0];
            elseif(real(RecievedArr(j))<=0 &&imag(RecievedArr(j))<=0)
                TempList=[0, 0];
            end
            DeMappedArr(ModulationType*j - 1 : j * ModulationType) = TempList ; % Append the bits to the demapped version as integers
        end


    elseif ModulationType == PSK8
        m=1;
        angles=zeros(1,8);
        for k=1:2:15
        angles(m)=k*pi/8;
        m=m+1;
        end
        for j = 1 : NumOfSymbols
            TempList = [0 0 0];
            if (angle(RecievedArr(j)) >= angles(8) && angle(RecievedArr(j)) <= angles(1) )
                TempList = [0 0 0];
            elseif ( angle(RecievedArr(j)) >= angles(1) && angle(RecievedArr(j)) <=angles(2) )
                TempList = [0 0 1];
            elseif ( angle(RecievedArr(j)) >= angles(2) && angle(RecievedArr(j)) <= angles(3))
                TempList = [0 1 1];
            elseif ( angle(RecievedArr(j)) >=angles(3)&& angle(RecievedArr(j)) <= angles(4) )
                TempList = [0 1 0];
            elseif ( angle(RecievedArr(j)) >= angles(4) && angle(RecievedArr(j)) <= angles(5) )
                TempList = [1 1 0];
            elseif ( angle(RecievedArr(j)) >= angles(5) && angle(RecievedArr(j)) <= angles(6) )
                TempList = [1 1 1];
            elseif ( angle(RecievedArr(j)) >= angles(6) && angle(RecievedArr(j)) <= angles(7) )
                TempList = [1 0 1];
            elseif ( angle(RecievedArr(j)) >= angles(7) && angle(RecievedArr(j)) <= angles(8) )
                TempList = [1 0 0];
            end
            DeMappedArr((j - 1) * ModulationType + 1 : j * ModulationType) = TempList; % Append the bits to the demapped version as integers
        end

    elseif ModulationType == QAM16
        TempList = [0 0 0 0];
        for j = 1:NumOfSymbols
            if imag(RecievedArr(j))>=2
                TempList(3:4)=[1 0];
            elseif imag(RecievedArr(j))>=0 && imag(RecievedArr(j))<2
                TempList(3:4)=[1 1];
            elseif imag(RecievedArr(j))<0 && imag(RecievedArr(j))>=-2
                TempList(3:4)=[0 1];
            elseif imag(RecievedArr(j))<-2
                TempList(3:4)=[0 0];
            end
            if real(RecievedArr(j))>=2
                TempList(1:2)=[1 0];
            elseif real(RecievedArr(j))>=0 && real(RecievedArr(j))<2
                TempList(1:2)=[1 1];
            elseif real(RecievedArr(j))<0 && real(RecievedArr(j))>=-2
                TempList(1:2)=[0 1];
            elseif real(RecievedArr(j))<-2
                TempList(1:2)=[0 0];
            end
            DeMappedArr((j - 1) * ModulationType + 1 : j * ModulationType) = TempList;
        end
    end
    
    
end


% Description:
% Entry point for the code's main logic

    NumOfSymbols = 1000;
    %number of bit represents every modulation type
    BPSK = 1;
    PSK8 = 3;
    QPSK = 2;
    QAM16 = 4;
    %just a number to represent BFSK
    BFSK = 5;

    
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
    BPSK_Mapped = Mapper(BPSK,BPSK_DataStream,BPSK_Positions,NumOfSymbols);
    BPSK_DeMapped = DeMapper(BPSK,BPSK_Mapped,NumOfSymbols);
    disp(isequal(BPSK_DataStream, BPSK_DeMapped));
    
    QPSK_DataStream = DataCreation(QPSK);
    QPSK_Mapped = Mapper(QPSK,QPSK_DataStream,QPSK_Positions,NumOfSymbols);
    QPSK_DeMapped = DeMapper(QPSK,QPSK_Mapped,NumOfSymbols);
    result=isequal(QPSK_DataStream,QPSK_DeMapped);
    disp(result);
    
    PSK8_DataStream = DataCreation(PSK8);
    PSK8_Mapped = Mapper(PSK8,PSK8_DataStream,PSK8_Positions,NumOfSymbols);
    PSK8_DeMapped = DeMapper(PSK8,PSK8_Mapped,NumOfSymbols);
    disp(isequal(PSK8_DataStream, PSK8_DeMapped));

    


    QAM16_DataStream = DataCreation(QAM16);
    QAM16_Mapped = Mapper(QAM16,QAM16_DataStream,QAM16_Positions,NumOfSymbols);
    QAM16_DeMapped = DeMapper(QAM16,QAM16_Mapped,NumOfSymbols);
    disp(isequal(QAM16_DataStream, QAM16_DeMapped));
    


