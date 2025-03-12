%% EX 4 -------------------------------------------------------------------

% Taxa de amostragem
fs = 48000;

% Diretório com os sinais
dirData = 'data';
userId = '40';

%{
% Loop que importa e reproduz um exemplo de cada dígito
for digito = 0:9

    % Nome do arquivo
    file = [num2str(digito), '_', userId, '_0.wav'];

    % Caminho completo para o arquivo
    f = fullfile(dirData, userId, file);
    
    % Testar se o caminho está correto
    disp(f);

    % Importar o sinal de áudio
    [x, ~] = audioread(f);
    
    % Reproduzir o sinal de áudio
    sound(x, fs);

    % Aguardar 2 segundos entre áudios
    pause(2);

    % Fazer plot do sinal de áudio
    figure;
    
    % Vetor de tempo
    t = (0:length(x)-1) / fs;

    plot(t, x);
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    title(['Dígito: ', num2str(digito)]);
end
%}

%% EX 4.1 -----------------------------------------------------------------

% Pares de dígitos
ParesDigitos = [0, 1; 2, 3; 4, 5; 6, 7; 8, 9];

%{
% Loop para comparar os pares de dígitos
for i = 1:size(ParesDigitos, 1)

    % Obter dígitos
    digito1 = ParesDigitos(i, 1);
    digito2 = ParesDigitos(i, 2);

    % Fazer os nomes dos arquivos
    file1 = [num2str(digito1), '_', userId, '_0.wav'];
    file2 = [num2str(digito2), '_', userId, '_0.wav'];

    % Fazer os caminhos para os arquivos
    f1 = fullfile(dirData, userId, file1);
    f2 = fullfile(dirData, userId, file2);

    % Importar os sinais de áudio
    [x1, ~] = audioread(f1);
    [x2, ~] = audioread(f2);

    % Fazer plot dos sinais de áudio
    figure;

    % Vetor de tempo para o primeiro sinal
    t1 = (0:length(x1)-1) / fs;
    % Vetor de tempo para o segundo sinal
    t2 = (0:length(x2)-1) / fs;
    
    %[M,I] = max(x1);
    %maximo1= M
    %tempo1= t1(I)   

    subplot(2, 1, 1);
    plot(t1, x1);
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    title(['Dígito ', num2str(digito1)]);

    %[M,I] = max(x2);
    %maximo2= M
    %tempo2= t2(I) 

    subplot(2, 1, 2);
    plot(t2, x2);
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    title(['Dígito ', num2str(digito2)]);

    % Aguardar 2 segundos entre plots
    pause(2);
end
%}


%% EX 4.2 -----------------------------------------------------------------

% Definir cores para cada dígito
colors = parula(10);

%{
% Vetores para armazenar as características dos sinais
energia = zeros(10, 50);
amplitudeMaxima = zeros(10, 50);

% Loop para percorrer os dígitos de 0 a 9
for digito = 0:9
    % Loop para percorrer os 50 arquivos de áudio para cada dígito
    for i = 0:49
        % Construir o nome do arquivo
        file = [num2str(digito), '_', userId, '_', num2str(i), '.wav'];

        % Caminho para o arquivo
        f = fullfile(dirData, userId, file);

        % Ler o sinal de áudio
        [x, ~] = audioread(f);

        % Calcular as características do sinal
        energia(digito+1, i+1) = sum(x.^2);  % Energia
        amplitudeMaxima(digito+1, i+1) = max(abs(x));  % Amplitude máxima
    end
    
    % Plotar as características do dígito atual
    scatter(energia(digito+1, :), amplitudeMaxima(digito+1, :), 50, colors(digito+1, :), 'filled');
    hold on;
end

% Configurar o gráfico
xlabel('Energia');
ylabel('Amplitude Máxima');
title('Características Temporais');
legend('0', '1', '2', '3', '4', '5', '6', '7', '8', '9');
%}

%{
% Criar vetores para armazenar as características dos sinais
digitos = zeros(1, 500);
amplitudeMaxima = zeros(1, 500);

% Loop para percorrer os dígitos de 0 a 9
for digito = 0:9
    % Loop para percorrer os 50 arquivos de áudio para cada dígito
    for repeticao = 0:49
        % Construir o nome do arquivo
        file = [num2str(digito), '_', userId, '_', num2str(repeticao), '.wav'];

        % Caminho completo para o arquivo
        f = fullfile(dirData, userId, file);

        % Ler o sinal de áudio
        [x, ~] = audioread(f);

        % Calcular a amplitude máxima do sinal
        maxAmplitude = max(abs(x));  % Amplitude máxima

        % Armazenar os valores no vetor correspondente
        index = digito * 50 + repeticao + 1;
        digitos(index) = digito;
        amplitudeMaxima(index) = maxAmplitude;
    end
end

% Plotar as características dos dígitos - Max Amplitude
scatter(digitos, amplitudeMaxima, 50, colors(digitos+1, :), 'Marker', 'o', 'MarkerFaceColor', 'none');

xlabel('Dígitos');
ylabel('Amplitude Máxima');
title('Caracteristicas Temporais - Amplitude Máxima');
%}

%{
% Criar vetores para armazenar as características dos sinais
digitos = zeros(1, 500);
energia = zeros(1, 500);

% Loop para percorrer os dígitos de 0 a 9
for digito = 0:9
    % Loop para percorrer os 50 arquivos de áudio para cada dígito
    for repeticao = 0:49
        % Construir o nome do arquivo
        file = [num2str(digito), '_', userId, '_', num2str(repeticao), '.wav'];

        % Caminho completo para o arquivo
        f = fullfile(dirData, userId, file);

        % Ler o sinal de áudio
        [x, ~] = audioread(f);

        % Calcular a energia do sinal
        energy = sum(x.^2);  % Energia

        % Armazenar os valores no vetor correspondente
        index = digito * 50 + repeticao + 1;
        digitos(index) = digito;
        energia(index) = energy;
    end
end

% Plotar as características dos dígitos - Energy
scatter(digitos, energia, 50, colors(digitos+1, :), 'Marker', 'o', 'MarkerFaceColor', 'none');

xlabel('Dígitos');
ylabel('Energia');
title('Características Temporais - Energia');
%}


%% EX 4.3 -----------------------------------------------------------------

%% EX 5 -------------------------------------------------------------------

%{
dirData = 'C:\Users\claud\Ambiente de Trabalho\UC\2º ano\2 sem\ATD\PROJETO\data\40';
numAmostras = 50;
fs = 8000;

espectros = zeros(numAmostras, fs/2+1, 10);
mediana = zeros(fs/2+1, 10);
primeiroQuartil = zeros(fs/2+1, 10);
terceiroQuartil = zeros(fs/2+1, 10);

% Loop through the digits 0 to 9
for digito = 0:9
    % Initialize matrix to store the amplitude spectra for the current digit signals
    espectrosDigito = zeros(numAmostras, fs/2+1);
    
    % Loop through the 50 audio files for each digit
    for repeticao = 0:numAmostras-1
        % Build the file name
        file = [num2str(digito), '_40_', num2str(repeticao), '.wav'];
        
        % Complete file path
        f = fullfile(dirData, file);

        % Read the audio signal
        audioData = audioread(f);
        window = rectwin(length(audioData));
        
        % Compute the Fourier transform using FFT
        fft_signal = fft(audioData .* window);
        
        % Compute the single-sided amplitude spectrum
        amplitude_spectrum = abs(fft_signal / length(audioData));
        amplitude_spectrum = amplitude_spectrum(1:fs/2+1);
        
        % Store the amplitude spectrum of the current signal
        espectrosDigito(repeticao+1, :) = amplitude_spectrum';
    end
    
    % Compute the median, first quartile, and third quartile of the amplitude spectra
    mediana(:, digito+1) = median(espectrosDigito);
    primeiroQuartil(:, digito+1) = quantile(espectrosDigito, 0.25);
    terceiroQuartil(:, digito+1) = quantile(espectrosDigito, 0.75);

    % Store the amplitude spectra of the current digit
    espectros(:, :, digito+1) = espectrosDigito;
    
    % Plot the median, first quartile, and third quartile for the current digit
    subplot(4, 3, digito+1);
    plot(0:(fs/2), mediana(:, digito+1), 'b', 'LineWidth', 2);
    hold on;
    plot(0:(fs/2), primeiroQuartil(:, digito+1), 'r--', 'LineWidth', 1);
    plot(0:(fs/2), terceiroQuartil(:, digito+1), 'k--', 'LineWidth', 1);
    xlabel('Frequência (Hz)');
    ylabel('Amplitude');
    title(['Dígito ', num2str(digito)]);
    legend('Mediana', '1º Quartil', '3º Quartil');
    hold off;
end

% Adjust the subplot layout
sgtitle('Espectros de Amplitude - Dígitos 0 a 9');
%}

%% 5.1

%{
folder = 'C:\Users\claud\Ambiente de Trabalho\UC\2º ano\2 sem\ATD\PROJETO\data\40';
n = 50;
mf = 8000;

% Define the window types
tiposJanela = {'Retangular', 'Hann', 'Hamming', 'Blackman'};

% Generate the file names
for i = 1:10
    figure; % Create a new figure for each digit
    
    for janela = 1:numel(tiposJanela)
        spectrum = zeros(50, mf);
        
        for j = 1:n
            fileName = [num2str(i-1), '_40_', num2str(j-1), '.wav'];
            fileName = fullfile(folder, fileName); %fullfile creates the full path of a file from 2 strings

            [audioData, sampleRate] = audioread(fileName);
            
            % Select the window type
            switch tiposJanela{janela}
                case 'Retangular'
                    window = rectwin(length(audioData));
                case 'Hann'
                    window = hann(length(audioData));
                case 'Hamming'
                    window = hamming(length(audioData));
                case 'Blackman'
                    window = blackman(length(audioData));
                otherwise
                    error('Invalid window type');
            end

            % Compute the Fourier transform using FFT
            fft_signal = fft(audioData .* window);
            
            % Compute the single-sided amplitude spectrum
            amplitude_spectrum = abs(fft_signal / length(audioData));
            amplitude_spectrum = amplitude_spectrum(1: floor(length(amplitude_spectrum)/2)+1);
            amplitude_spectrum = amplitude_spectrum(1: mf);
            spectrum(j, 1:length(amplitude_spectrum))  = amplitude_spectrum';
        end

        q2 = median(spectrum);
        q1 = quantile(spectrum, 0.25);
        q3 = quantile(spectrum, 0.75);

        % Create a subplot for each window type
        subplot(2, 2, janela);
        
        % Plot the amplitude spectrum
        plot(q2, 'b', 'LineWidth', 2);
        hold on;
        plot(q1, 'r--', 'LineWidth', 1);
        plot(q3, 'k--', 'LineWidth', 1);
        
        xlabel('Frequência (Hz)');
        ylabel('Amplitude');
        title(tiposJanela{janela}, 'Interpreter', 'none'); % Display the window type as the title
        
        grid on;
        hold off;
    end
    
    sgtitle(['Espectro de Amplitude - Dígito ', num2str(i-1)]);
end
%}

%% 5.2

%{
folder = 'C:\Users\claud\Ambiente de Trabalho\UC\2º ano\2 sem\ATD\PROJETO\data\40';
n = 50;
mf = 8000;

% Choose two pairs of digits for comparison (e.g., 0 and 1, and 2 and 3)
pair1 = [0, 1];
pair2 = [2, 3];

% Plotting the spectral characteristics for pair 1
figure;
subplot(2, 1, 1);
for digit = pair1
    spectrum = zeros(n, mf);
    for repetition = 1:n
        fileName = [num2str(digit), '_40_', num2str(repetition-1), '.wav'];
        fileName = fullfile(folder, fileName);

        [audioData, ~] = audioread(fileName);
        window = rectwin(length(audioData));

        % Compute the Fourier transform using FFT
        fft_signal = fft(audioData .* window);

        % Compute the single-sided amplitude spectrum
        amplitude_spectrum = abs(fft_signal / length(audioData));
        amplitude_spectrum = amplitude_spectrum(1:floor(length(amplitude_spectrum)/2)+1);
        amplitude_spectrum = amplitude_spectrum(1:mf);
        spectrum(repetition, 1:length(amplitude_spectrum)) = amplitude_spectrum';
    end

    plot(mean(spectrum), 'LineWidth', 2);
    hold on;
end
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(['Spectral Characteristics - Pair ', num2str(pair1(1)), ' and ', num2str(pair1(2))]);
legend(['Digit ', num2str(pair1(1))], ['Digit ', num2str(pair1(2))]);
grid on;
hold off;

% Plotting the spectral characteristics for pair 2
subplot(2, 1, 2);
for digit = pair2
    spectrum = zeros(n, mf);
    for repetition = 1:n
        fileName = [num2str(digit), '_40_', num2str(repetition-1), '.wav'];
        fileName = fullfile(folder, fileName);

        [audioData, ~] = audioread(fileName);
        window = rectwin(length(audioData));

        % Compute the Fourier transform using FFT
        fft_signal = fft(audioData .* window);

        % Compute the single-sided amplitude spectrum
        amplitude_spectrum = abs(fft_signal / length(audioData));
        amplitude_spectrum = amplitude_spectrum(1:floor(length(amplitude_spectrum)/2)+1);
        amplitude_spectrum = amplitude_spectrum(1:mf);
        spectrum(repetition, 1:length(amplitude_spectrum)) = amplitude_spectrum';
    end

    plot(mean(spectrum), 'LineWidth', 2);
    hold on;
end
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(['Spectral Characteristics - Pair ', num2str(pair2(1)), ' and ', num2str(pair2(2))]);
legend(['Digit ', num2str(pair2(1))], ['Digit ', num2str(pair2(2))]);
grid on;
hold off;
%}

%% 5.3

%{
folder = 'C:\Users\claud\Ambiente de Trabalho\UC\2º ano\2 sem\ATD\PROJETO\data\40';

% Definir o número total de dígitos
numDigitos = 10;

% Definir o número de amostras para cada dígito
numAmostras = 50;

% Definir a taxa de amostragem
fs = 8000;

% Inicializar matriz para armazenar os espectros de amplitude de cada dígito
espectros = zeros(numDigitos, fs/2+1);

% Loop para percorrer os dígitos de 0 a 9
for digito = 0:numDigitos-1
    % Inicializar matriz para armazenar os espectros de amplitude dos sinais do dígito atual
    espectrosDigito = zeros(numAmostras, fs/2+1);
    
    % Loop para percorrer as amostras de áudio para cada dígito
    for amostra = 1:numAmostras
        % Construir o nome do arquivo
        fileName = [num2str(digito), '_40_', num2str(amostra-1), '.wav'];
        filePath = fullfile(folder, fileName);
        
        % Ler o sinal de áudio
        [audioData, ~] = audioread(filePath);
        
        % Calcular o espectro de amplitude
        espectro = abs(fft(audioData))./length(audioData);
        
        % Armazenar o espectro de amplitude da amostra atual
        espectrosDigito(amostra, :) = espectro(1:fs/2+1);
    end
    
    % Calcular a média dos espectros de amplitude do dígito atual
    mediaEspectrosDigito = mean(espectrosDigito, 1);
    
    % Armazenar a média dos espectros de amplitude do dígito atual
    espectros(digito+1, :) = mediaEspectrosDigito;
end

% Plotar os espectros de amplitude para cada dígito
figure;
frequencies = 0:(fs/2);
for digito = 0:numDigitos-1
    plot(frequencies, espectros(digito+1, :));
    hold on;
end
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Características Espectrais dos Dígitos');
legend('0', '1', '2', '3', '4', '5', '6', '7', '8', '9');
grid on;
hold off;
%}

%% 5.4

%{
% Dados das características obtidas
amplitudeMinima = [4.168e-06, 7.645e-07, 3.469e-06, 2.688e-06, 1.126e-06, 1.422e-06, 2.278e-06, 3.835e-06, 2.523e-06, 8.889e-07];
amplitudeMaxima = [0.0004348, 0.0004628, 0.0003758, 0.0002597, 0.0003621, 0.0002169, 0.0001569, 0.0002716, 0.0002768, 0.0004020];
amplitudeMedia = (amplitudeMinima + amplitudeMaxima) / 2;

% Plotar superfícies de decisão em 3D
figure;
scatter3(amplitudeMinima, amplitudeMaxima, amplitudeMedia, 50, [0:9], 'filled');
xlabel('Amplitude Mínima');
ylabel('Amplitude Máxima');
zlabel('Amplitude Média');
title('Superfícies de Decisão');
colorbar;

% Plotar superfícies de decisão em 2D
figure;
subplot(1, 2, 1);
scatter(amplitudeMinima, amplitudeMaxima, 50, [0:9], 'filled');
xlabel('Amplitude Mínima');
ylabel('Amplitude Máxima');
title('Superfícies de Decisão (2D)');
colorbar;

subplot(1, 2, 2);
scatter(amplitudeMinima, amplitudeMedia, 50, [0:9], 'filled');
xlabel('Amplitude Mínima');
ylabel('Amplitude Média');
title('Superfícies de Decisão (2D)');
colorbar;
%}

%% 6