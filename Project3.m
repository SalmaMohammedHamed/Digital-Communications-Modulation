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
            if (RecievedArr(j) >= 0)
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
%   calculateBER calculates the Bit Error Rate (BER) and theoretical BER for a given set of Signal-to-Noise Ratio (SNR) values.
%   It takes the SNR values, energy per bit (Eb), noise factor, mapped signal, symbol array, and data stream as inputs, and
%   returns the simulated BER and theoretical BER.
function [BER_sim, BER_ther] = calculateBER(SNR, Eb, noise_bits, ModulationType, mapped_signal, symbol_array, data_stream, NumOfBits)
    BER_sim = zeros(size(SNR));
    BER_ther = zeros(size(SNR));
    N = length(SNR);
    PSK8 = 3;
    BFSK = 5;
    QAM16 = 4;
    for i = 1:N
        N0 = Eb / (10^(SNR(i) / 10));
        noise = noise_bits * sqrt(N0/2);
        received_signal = mapped_signal + noise;
        demapped_signal = DeMapper(ModulationType, received_signal, symbol_array, NumOfBits);
        
        calc_error = sum(demapped_signal ~= data_stream);
        BER_sim(i) = calc_error / NumOfBits;
        N1 = Eb / N0;
        if ModulationType == PSK8
            BER_ther(i)=(1/3)*erfc(sqrt(1/N0)*sin(pi/8));
        elseif ModulationType == BFSK
            BER_ther(i) = 0.5 * erfc(sqrt(N1/2));
        elseif ModulationType == QAM16
            BER_ther(i) = (3/8) * erfc(sqrt(N1/2.5));
        else
            BER_ther(i) = 0.5 * erfc(sqrt(N1));
        end
    end
end
% Description:
%   Plotting BER theortical and simulated to compare between them
function plotBER(SNR, BER_sim, BER_ther, title_text)
    figure();
    semilogy(SNR, BER_sim, 'b');
    hold on;
    semilogy(SNR, BER_ther, 'r--');
    title(title_text);
    xlabel('Eb/No (dB)');
    ylabel('BER');
    legend('practical', 'theoretical');
end
% Description:
%   Plotting BER theortical and simulated for all modulations on one plot
function plotAll(SNR, BER_sim_bpsk, BER_ther_bpsk, BER_sim_qpsk, BER_ther_qpsk, BER_sim_qpsk_ng, ...
 BER_sim_psk8, BER_ther_psk8, BER_sim_QAM, BER_ther_QAM, title_text)
    
    figure();
    semilogy(SNR, BER_sim_bpsk, 'LineWidth', 2); % Thickened line for BPSK simulated
    hold on;
    semilogy(SNR, BER_ther_bpsk, '--', 'LineWidth', 2); % Thickened dashed line for BPSK theoretical
    hold on;
    semilogy(SNR, BER_sim_qpsk, 'LineWidth', 2); % Thickened line for QPSK simulated
    hold on;
    semilogy(SNR, BER_ther_qpsk, '--', 'LineWidth', 2); % Thickened dashed line for QPSK theoretical
    hold on;
    semilogy(SNR, BER_sim_qpsk_ng, 'LineWidth', 2); % Thickened line for QPSK Not Grey simulated
    hold on;
    semilogy(SNR, BER_sim_psk8, 'LineWidth', 2); % Thickened line for PSK8 simulated
    hold on;
    semilogy(SNR, BER_ther_psk8, '--', 'LineWidth', 2); % Thickened dashed line for PSK8 theoretical
    hold on;
    semilogy(SNR, BER_sim_QAM, 'LineWidth', 2); % Thickened line for QAM16 simulated
    hold on;
    semilogy(SNR, BER_ther_QAM, '--', 'LineWidth', 2); % Thickened dashed line for QAM16 theoretical
    
    hold off;
    title(title_text);
    xlabel('Eb/No (dB)');
    ylabel('BER');
    legend('BPSK Simulated', 'BPSK Theoretical', 'QPSK Simulated', 'QPSK Theoretical', ...
        'QPSK Not Grey Simulated', 'PSK8 Simulated', 'PSK8 Theoretical', ...
        'QAM16 Simulated', 'QAM16 Theoretical');
    ylim([1e-6, 1]); % Set y-axis limits from 10^-4 to 1
end


%Descroption:
%creating a whole ensample of BFSK
function BFSK_EnsampleArrMapped = BFSK_Ensample(NumRealizations,NumRVs)
    Tb=10;
    Eb=1;
    t=1:1:Tb;
    delf=1/Tb;
    BFSK_EnsampleArr = zeros(NumRealizations,NumRVs);
    BFSK_EnsampleArrMapped=zeros(NumRealizations,NumRVs);
    for j = 1:NumRealizations
        BFSK_EnsampleArr(j,:) = randi([0 1], 1, NumRVs);
        
    end
    for i=1:1:NumRealizations
        for j=1:1:NumRVs
          BFSK_EnsampleArrMapped(i,(j-1)*Tb+1:(j)*Tb)=BFSK_EnsampleArr(i,j)*sqrt(2*Eb/Tb)+(~BFSK_EnsampleArr(i,j))*...
              sqrt(2*Eb/Tb)*(cos(2*pi*delf*t)+1i*sin(2*pi*delf*t));
        end
    end
end

%Description:
%calculating the Satatistical autocorrelation
    function StatAutoArr = StatAuto(ensample,NumRealizations,NumRVs)
    StatAutoArr = zeros(NumRVs, NumRVs);
    for t = 1:NumRVs  % This loop sweeps over t.
        for taw = 0:NumRVs-t % This loop sweeps over taw for a certain t.
            sum = 0;
            for i = 1:NumRealizations
                sum_mat = ensample(i,t) .* ensample(i,t+taw); % Autocorrelation
                sum = sum + sum_mat; % Sum of autocorrelated values of certain t and taw in each realization.
            end
            StatAutoArr(taw+1,t) = sum / NumRealizations;
        end
    end
end

%Description:
%BFSK PSD
function BFSK_PSD (StatAutoArr)
    % Get the length of the autocorrelation vector
    N = length(StatAutoArr(:, 1));
    
    % Perform FFT on the autocorrelation vector
    Sx = fft(StatAutoArr(:, 1), N);
    
    % Define the frequency axis
    n = -N/2:N/2-1;
    fs = 10; % Sampling frequency
    
    % Plot the PSD
    figure;
    plot(fs * n / N, fftshift(abs(Sx)), 'LineWidth', 1);
    title('PSD for BFSK');
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density');
    ylim([0, 5]);
    xlim([-3 ,3]);
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
    
    %%----------------------- BPSK -------------- %%
    BPSK_DataStream = DataCreation(NumOfBits);
    BPSK_Mapped = Mapper(BPSK,BPSK_DataStream,BPSK_Positions,NumOfBits/BPSK);
    
    noise_bits=randn(size(BPSK_Mapped));
    SNR=-4:14;
    Eb=1;
    [BER_sim, BER_ther] = calculateBER(SNR, Eb, noise_bits, BPSK, BPSK_Mapped, BPSK_SympolsArr, BPSK_DataStream, NumOfBits);
    plotBER(SNR, BER_sim, BER_ther, 'BER Performance for BPSK');


    
    %%---------------------QPSK--------------------%%
    QPSK_DataStream = DataCreation(NumOfBits);
    QPSK_Mapped = Mapper(QPSK,QPSK_DataStream,QPSK_Positions,NumOfBits/QPSK);
    noise_qpsk=randn(size(QPSK_Mapped))+randn(size(QPSK_Mapped))*1i;
    Eb=1;
    SNR=-4:14;
    [BER_sim2, BER_ther2] = calculateBER(SNR, Eb, noise_qpsk, QPSK, QPSK_Mapped, QPSK_SympolsArr, QPSK_DataStream, NumOfBits);
    plotBER(SNR, BER_sim2, BER_ther2, 'BER Performance for QPSK');

    %%---------------------QPSK not grey ----------%%
    QPSK_NotGreyCodedDataStream = DataCreation(NumOfBits);
    QPSK_NotGreyCodedMapped = Mapper(QPSK,QPSK_NotGreyCodedDataStream,QPSK_NotGreyCodedPositions,NumOfBits/QPSK);
    noise_not_grey_qbsk=randn(size(QPSK_NotGreyCodedMapped))+randn(size(QPSK_NotGreyCodedMapped))*1i;
    Eb=1;
    SNR=-4:14;
    [BER_sim3, BER_ther3] = calculateBER(SNR, Eb, noise_not_grey_qbsk, QPSK ,QPSK_NotGreyCodedMapped, QPSK_NotGreyCodedSympolsArr, QPSK_NotGreyCodedDataStream, NumOfBits);
    plotBER(SNR, BER_sim3, BER_ther3, 'BER Performance for QPSK Not Grey');

    %%-----------------PSK8-------------------%%    
    PSK8_DataStream = DataCreation(NumOfBits);
    PSK8_Mapped = Mapper(PSK8,PSK8_DataStream,PSK8_Positions,NumOfBits/PSK8);
    noise_PSK8=randn(size(PSK8_Mapped))+randn(size(PSK8_Mapped))*1i;
    Eb=1/3;
    SNR=-4:14;
    [BER_sim4, BER_ther4] = calculateBER(SNR, Eb, noise_PSK8, PSK8, PSK8_Mapped, PSK8_SympolsArr, PSK8_DataStream, NumOfBits);
    plotBER(SNR, BER_sim4, BER_ther4, 'BER Performance for PSK8');


    %%-------------------QAM16--------------%%
    QAM16_DataStream = DataCreation(NumOfBits);
    QAM16_Mapped = Mapper(QAM16,QAM16_DataStream,QAM16_Positions,NumOfBits/QAM16);
    noise_QAM=randn(size(QAM16_Mapped))+randn(size(QAM16_Mapped))*1i;
    Eb=2.5;
    SNR=-4:14;
    [BER_sim5, BER_ther5] = calculateBER(SNR, Eb, noise_QAM, QAM16, QAM16_Mapped, QAM16_SympolsArr, QAM16_DataStream, NumOfBits);
    plotBER(SNR, BER_sim5, BER_ther5, 'BER Performance for QAM16');
    
    %%------------------plot All on one plot ---%
    plotAll(SNR, BER_sim, BER_ther,BER_sim2, BER_ther2, BER_sim3, BER_sim4, BER_ther4, BER_sim5, BER_ther5, "A plot for all modulations BER");
    %%--------------------BFSK-------------%%
    BFSK_DataStream = DataCreation(NumOfBits);
    BFSK_Mapped = Mapper(BFSK,BFSK_DataStream,BFSK_Positions,NumOfBits);
    noise_BFSK=randn(size(BFSK_Mapped))+randn(size(BFSK_Mapped))*1i;
    Eb=1;
    SNR=-4:14;
    [BER_sim6, BER_ther6] = calculateBER(SNR, Eb, noise_BFSK, BFSK, BFSK_Mapped, BFSK_SympolsArr, BFSK_DataStream, NumOfBits);
    plotBER(SNR, BER_sim6, BER_ther6, 'BER Performance for BFSK');

    %%------------BFSK PSD---------------%%
    NumRealizations=500;
    NumRVs=700;
    BFSK_EnsampleArr=BFSK_Ensample(NumRealizations,NumRVs);
    StatAutoArr=StatAuto(BFSK_EnsampleArr,NumRealizations,NumRVs);
    BFSK_PSD(StatAutoArr);

end
