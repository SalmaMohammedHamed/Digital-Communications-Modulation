%the excution starts
clear;
close all;
main();

% Description:
%Create the data stream for each Mudulation type 
%the Data Stream length depands on the number of bits in the symbol
function DataStreamArr = DataCreation(NumOfBits)
    DataStreamArr = randi([0 1], 1, NumOfBits);
end


% Description:
%map the bits to the corresponding symbols 
function MappedArr = Mapper(ModulationType,DataStream, Positions,NumOfSymbols)
    BPSK = 1;
    PSK8 = 3;
    QPSK = 2;
    QAM16 = 4;
    BFSK=5;
    if ModulationType == BPSK
        MappedArr = DataStream * 2 - 1;
    elseif ModulationType==BFSK
        MappedArr = zeros(1,NumOfSymbols);
        for j = 1:NumOfSymbols
            if DataStream(j)==1
                MappedArr(j)=Positions(1);
            elseif DataStream(j)==0
                MappedArr(j)=Positions(2);
            end
        end
    elseif ModulationType== PSK8 || ModulationType==QPSK || ModulationType==QAM16
        MappedArr = zeros(1,NumOfSymbols);
        
        for j = 1:ModulationType:length(DataStream)
            symbol = DataStream(j:j+ModulationType-1);  % every n bits is a symbol
            decimal_value = bin2dec(num2str(symbol)); % getting the decimal value of the symbol
            MappedArr(1,fix(j/ModulationType)+1) = Positions(1,decimal_value+1); % map to the corresponding symbol
        end
    else 
        error("Wronge Moduation Type");
     end
end




% Description:
% Define the demapper function (get the data Stream from the received symbols)
% NOTE: After creating the channel, replace every mapped data by the mapped data after the channel
function DeMappedArr = DeMapper(ModulationType, RecievedArr,SymbolsArr,NumOfBits)
    %number of bit represents every modulation type
    BPSK = 1;
    PSK8 = 3;
    QPSK = 2;
    QAM16 = 4;
    BFSK=5;
    DeMappedArr = zeros(1,NumOfBits);
    NumOfSymbols = NumOfBits/ModulationType;
    % Perform demapping based on the modulation type
    if ModulationType == BPSK
        for j = 1 : length(RecievedArr)
            if (RecievedArr(j) > 0)
            DeMappedArr(j) = SymbolsArr(1);
            else
            DeMappedArr(j) = SymbolsArr(2);
            end
        end

    elseif ModulationType == QPSK
        TempList= [0 0];
        for j = 1 : NumOfSymbols
            if(real(RecievedArr(j))>=0 &&imag(RecievedArr(j))>0)
                TempList=SymbolsArr(1,:);
            elseif(real(RecievedArr(j))<=0 &&imag(RecievedArr(j))>0)
                TempList=SymbolsArr(2,:);
            elseif(real(RecievedArr(j))>=0 &&imag(RecievedArr(j))<0)
                TempList=SymbolsArr(3,:);
            elseif(real(RecievedArr(j))<=0 &&imag(RecievedArr(j))<=0)
                TempList=SymbolsArr(4,:);
            end
            DeMappedArr(ModulationType*j - 1 : j * ModulationType) = TempList ; % Append the bits to the demapped version as integers
        end


    elseif ModulationType == PSK8

        PSK8Angles=zeros(1,8);
        i=1;
        for j=1:4
            PSK8Angles(j)=pi*(i)/8;
            i=i+2;
        end
        for j = 1 : NumOfSymbols
            Angle=angle(RecievedArr(j));
           TempList = [0 0 0];
            % Map the angle to binary bits
           if (Angle > -PSK8Angles(1)) && (Angle < PSK8Angles(1))
                TempList=SymbolsArr(1,:);
           elseif (Angle > PSK8Angles(1)) && (Angle < PSK8Angles(2))
                TempList=SymbolsArr(2,:);
        
            elseif (Angle > PSK8Angles(2)) && (Angle < PSK8Angles(3))
                TempList=SymbolsArr(3,:);
            elseif (Angle > PSK8Angles(3)) && (Angle < PSK8Angles(4))
                TempList=SymbolsArr(4,:);
            elseif (Angle > PSK8Angles(4)) || (Angle <  -PSK8Angles(4))
                TempList=SymbolsArr(5,:);
            elseif (Angle > -PSK8Angles(4)) && (Angle < -PSK8Angles(3))
                TempList=SymbolsArr(6,:);
            elseif (Angle > -PSK8Angles(3)) && (Angle < -PSK8Angles(2))
                TempList=SymbolsArr(7,:);
           elseif (Angle > -PSK8Angles(2)) && (Angle < -PSK8Angles(1))
                TempList=SymbolsArr(8,:);
            end
            DeMappedArr((j - 1) * ModulationType + 1 : j * ModulationType) = TempList; % Append the bits to the demapped version as integers
        end

    elseif ModulationType == QAM16
        TempList = [0 0 0 0];
        for j = 1:NumOfSymbols
            if imag(RecievedArr(j))>=2
                TempList(3:4)=SymbolsArr(1,:);
            elseif imag(RecievedArr(j))>=0 && imag(RecievedArr(j))<2
                TempList(3:4)=SymbolsArr(2,:);
            elseif imag(RecievedArr(j))<0 && imag(RecievedArr(j))>=-2
                TempList(3:4)=SymbolsArr(3,:);
            elseif imag(RecievedArr(j))<-2
                TempList(3:4)=SymbolsArr(4,:);
            end
            if real(RecievedArr(j))>=2
                TempList(1:2)=SymbolsArr(1,:);
            elseif real(RecievedArr(j))>=0 && real(RecievedArr(j))<2
                TempList(1:2)=SymbolsArr(2,:);
            elseif real(RecievedArr(j))<0 && real(RecievedArr(j))>=-2
                TempList(1:2)=SymbolsArr(3,:);
            elseif real(RecievedArr(j))<-2
                TempList(1:2)=SymbolsArr(4,:);
            end
            DeMappedArr((j - 1) * ModulationType + 1 : j * ModulationType) = TempList;
        end
        
    elseif ModulationType==BFSK
        Angle1=pi/4;
        Angle2=3*pi/4;
        for j = 1 : length(RecievedArr)
            if (angle(RecievedArr(j)) < Angle1 && angle(RecievedArr(j))> -Angle2)
            DeMappedArr(j) = SymbolsArr(1);
            else
            DeMappedArr(j) = SymbolsArr(2);
            end
        end
    end
     
end




% Description:
% Entry point for the code's main logic
function main
    NumOfBits=240000;
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
    QPSK_NotGreyCodedPositions = [-1-1i , -1+1i , 1+1i , 1-1i];
    % Define symbols positions for 16-QAM modulation
    QAM16_Positions = [-3-3i , -3-1i , -3+3i , -3+1i , -1-3i , -1-1i , -1+3i , -1+1i , ...
                       3-3i , 3-1i , 3+3i , 3+1i , 1-3i , 1-1i , 1+3i , 1+1i];
    
    % Define symbols positions for 8PSK modulation
    PSK8_Positions = zeros(1, 8);
    PSK8_Positions(1)=1.0 + 0.0i;
    PSK8_Positions(2)=0.707106781186548 + 0.707106781186547i;
    PSK8_Positions(4)=0.0 + 1.0i;
    PSK8_Positions(3)=-0.707106781186547 + 0.707106781186548i;
    PSK8_Positions(5)=0.707106781186547 -0.707106781186548i;
    PSK8_Positions(6)=0.0 - 1.0i;
    PSK8_Positions(7)=-1.0 + 0.0i;
    PSK8_Positions(8)=-0.707106781186547 -0.707106781186548i;
    %symbol positions for BFSK 
    BFSK_Positions=[ 1.0+0.0i,0.0+1.0i];


    %demapping arrays (every region equivalent symbol)
    BPSK_SympolsArr=[1 0];
    QPSK_SympolsArr=[[1 1]; [0 1]; [1 0]; [0 0]];
    QPSK_NotGreyCodedSympolsArr=[[1 0]; [0 1]; [1 1]; [0 0]];
    PSK8_SympolsArr=[[0 0 0];[0 0 1];[0 1 1];[0 1 0];[1 1 0];[1 1 1];[1 0 1];[1 0 0];];
    QAM16_SympolsArr=[[1 0];[1 1];[0 1];[0 0]];
    BFSK_SympolsArr=[1 0];
    
    
    BPSK_DataStream = DataCreation(NumOfBits);
    BPSK_Mapped = Mapper(BPSK,BPSK_DataStream,BPSK_Positions,NumOfBits/BPSK);
    BPSK_DeMapped = DeMapper(BPSK,BPSK_Mapped,BPSK_SympolsArr,NumOfBits);%after creating the channel use the Mapped symbols after the channel effect
    disp(isequal(BPSK_DataStream, BPSK_DeMapped));
    
    QPSK_DataStream = DataCreation(NumOfBits);
    QPSK_Mapped = Mapper(QPSK,QPSK_DataStream,QPSK_Positions,NumOfBits/QPSK);
    QPSK_DeMapped = DeMapper(QPSK,QPSK_Mapped,QPSK_SympolsArr,NumOfBits);
    result=isequal(QPSK_DataStream,QPSK_DeMapped);
    disp(result);
    
    QPSK_NotGreyCodedDataStream = DataCreation(NumOfBits);
    QPSK_NotGreyCodedMapped = Mapper(QPSK,QPSK_NotGreyCodedDataStream,QPSK_NotGreyCodedPositions,NumOfBits/QPSK);
    QPSK_NotGreyCodedDeMapped = DeMapper(QPSK,QPSK_NotGreyCodedMapped,QPSK_NotGreyCodedSympolsArr,NumOfBits);
    result=isequal(QPSK_NotGreyCodedDataStream,QPSK_NotGreyCodedDeMapped);
    disp(result);
    
    PSK8_DataStream = DataCreation(NumOfBits);
    PSK8_Mapped = Mapper(PSK8,PSK8_DataStream,PSK8_Positions,NumOfBits/PSK8);
    PSK8_DeMapped = DeMapper(PSK8,PSK8_Mapped,PSK8_SympolsArr,NumOfBits);
    disp(isequal(PSK8_DataStream, PSK8_DeMapped));
    


    QAM16_DataStream = DataCreation(NumOfBits);
    QAM16_Mapped = Mapper(QAM16,QAM16_DataStream,QAM16_Positions,NumOfBits/QAM16);
    QAM16_DeMapped = DeMapper(QAM16,QAM16_Mapped,QAM16_SympolsArr,NumOfBits);
    disp(isequal(QAM16_DataStream, QAM16_DeMapped));
    

    BFSK_DataStream = DataCreation(NumOfBits);
    BFSK_Mapped = Mapper(BFSK,BFSK_DataStream,BFSK_Positions,NumOfBits);
    BFSK_DeMapped = DeMapper(BFSK,BFSK_Mapped,BFSK_SympolsArr,NumOfBits);
    disp(isequal(BFSK_DataStream, BFSK_DeMapped));
end