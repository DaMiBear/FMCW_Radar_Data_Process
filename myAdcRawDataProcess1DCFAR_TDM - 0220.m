clear;
close all;

global Params;
global dataSet;
Params.NChirp = 128;                %一帧的chirp个数
Params.NChan = 4;                   %RxAn数,ADC通道数
Params.NSample = 256;               %每个chirp ADC采样数
Params.Fs = 10e6;                    %采样频率
Params.c = 3.0e8;                   %光速
Params.startFreq = 77.18e9;            %起始频率 ????
Params.freqSlope = 29.9817e12;      %chirp的斜率
Params.bandwidth = (1.5351e9);        %带宽!!!!---真实带宽
Params.lambda=Params.c/Params.startFreq;   %雷达信号波长
Params.Tc = 160e-6;                 %chirp周期 ????
Params.NChan_V = 2*4;               %MIMO中包括虚拟天线总数 Ntx*Nrx
Params.NChirp_V = Params.NChirp/2;  %分离出虚拟天线的数据后chirp会减半
Params.Tc_V = 2*Params.Tc;          %有虚拟天线，chirp的周期要加倍
Params.dopplerPFA = 0.03;           %多普勒维cfarPFA
Params.rangePFA = 0.05;             %距离维cfarPFA

%打开文件 fid用于储存文件句柄值,r表示只读 b表示二进制形式打开
Params.fid_rawData = fopen('adc_data_TDM_10M.bin','rb');
Params.dataSizeOneFrame = Params.NSample*4*Params.NChirp*Params.NChan;

%处理数据
dp_updateFrameData(1);  %第i帧(1-8)

% 分离出虚拟天线的数据
dp_separateVirtualRx();

%2D FFT
%对dataSet.rawFrameData执行距离FFT 
dataSet.radarCubeData = processingChain_rangeFFT(1);

% [X,Y] = meshgrid(Params.c*(0:Params.NSample-1)*Params.Fs/2/Params.freqSlope/Params.NSample, ...
%     (-Params.NChirp/2:Params.NChirp/2 - 1)*Params.lambda/Params.Tc/Params.NChirp/2);    
% mesh(X,Y,abs(reshape(dataSet.radarCubeData(:,1,:),Params.NChirp,Params.NSample)));
% title('1d fft 结果');
% 执行多普勒FFT
dataSet.radarCubeData = processingChain_dopplerFFT(1);


%对NChirp_V个接收天线Chan求和
dataSet.sumFFT2D_radarCubeData = single(zeros(Params.NChirp_V,Params.NSample));   %创建一个NChirp_V行 NSample列的空二维数组
for chanNum = 1:Params.NChan_V
    %把三维数组变为二维数组
    FFT2D_radarCubeData = reshape(dataSet.radarCubeData(:,chanNum,:),Params.NChirp_V,Params.NSample);   
    %求模->取log->相加
    dataSet.sumFFT2D_radarCubeData = (abs(FFT2D_radarCubeData)) + dataSet.sumFFT2D_radarCubeData;
%     mesh(X,Y,10*log10(abs(reshape(dataSet.radarCubeData(:,chanNum,:),Params.NChirp,Params.NSample))));
%     figure;
end

dopplerDimCfarThresholdMap = zeros(size(dataSet.sumFFT2D_radarCubeData));  %创建一个二维矩阵存放doppler维cfar后的结果
dopplerDimCfarResultMap = zeros(size(dataSet.sumFFT2D_radarCubeData));


[X,Y] = meshgrid(Params.c*(0:Params.NSample-1)*Params.Fs/2/Params.freqSlope/Params.NSample, ...
                (-Params.NChirp_V/2:Params.NChirp_V/2 - 1)*Params.lambda/Params.Tc_V/Params.NChirp_V/2);    

% 多普勒维度进行CFAR
for i = 1:Params.NSample
    dopplerDim = reshape(dataSet.sumFFT2D_radarCubeData(:,i),1,Params.NChirp_V);  %变成一行数据
    [cfar1D_Arr,threshold] = ac_cfar1D(12,2,Params.dopplerPFA,dopplerDim);  %进行1D cfar
    dopplerDimCfarResultMap(:,i) = cfar1D_Arr.'; 
    dopplerDimCfarThresholdMap(:,i) = threshold.';
end

mesh(X,Y,(20*log10(dataSet.sumFFT2D_radarCubeData)));
xlabel('距离(m)');ylabel('速度(m/s)');zlabel('信号幅值dB');
title('四天线求和后的2D_FFT结果');
figure;

mesh(X,Y,20*log10(dopplerDimCfarThresholdMap));
xlabel('距离(m)');ylabel('速度(m/s)');zlabel('信号幅值dB');
title('doppler维度CFAR门限图');
figure;

mesh(X,Y,(dopplerDimCfarResultMap));
xlabel('距离(m)');ylabel('速度(m/s)');zlabel('信号幅值dB');
title('doppler维度CFAR判决结果');
figure;

%沿着doppler维度方向寻找在doppler维cfar判决后为1的结果
saveMat = zeros(size(dataSet.sumFFT2D_radarCubeData));
for range = 1:Params.NSample
    indexArr = find(dopplerDimCfarResultMap(:,range)==1);
    objDopplerArr = [indexArr;zeros(Params.NChirp_V - length(indexArr),1)];   %补充长度
    saveMat(:,range) = objDopplerArr; %保存doppler下标
end
% 保存有物体的doppler坐标
objDopplerIndex = unique(saveMat);  % unqiue是不重复的返回数组中的数

% 根据之前doppler维的cfar结果对应的下标saveMat，对相应的速度进行range维度的CFAR
rangeDimCfarThresholdMap = zeros(size(dataSet.sumFFT2D_radarCubeData));  %创建一个二维矩阵存放range维cfar后的结果
rangeDimCfarResultMap = zeros(size(dataSet.sumFFT2D_radarCubeData));
i = 1;
while(i<=length(objDopplerIndex))
    if(objDopplerIndex(i)==0)   % 因为数组里面有0,防止下面j取到0
        i = i + 1;
        continue;
    else    %根据速度下标进行range CFAR
        j = objDopplerIndex(i);     % 获得物体所在的行
        rangeDim = reshape(dataSet.sumFFT2D_radarCubeData(j, :),1,Params.NSample);  %变成一行数据
        % tip 这个PFA很迷啊,如果设置的低一些,在进行分支聚集的时候,可能的结果是没有检测到物体
        % 因为在进行rangeCFAR的时候,把附近的最大值给滤掉了,那这样在进行峰值聚集的时候,判决结果对应的最大值一直是小的
        % ╮(╯▽╰)╭
        [cfar1D_Arr,threshold] = ac_cfar1D(3,2,Params.rangePFA,rangeDim);  %进行1D cfar
        rangeDimCfarResultMap(j,:) = cfar1D_Arr; 
        rangeDimCfarThresholdMap(j,:) = threshold;
        i = i + 1;
        
        plot(20*log10(rangeDim));hold on;
        plot(20*log10(threshold));
        title(['距离维度门限dopplerIndex=',num2str(j)]);
        xlabel('距离(m)');ylabel('信号幅值dB');
        figure;
        
    end
end
% plot((rangeDim));hold on;
% plot(threshold);
mesh(X,Y,(rangeDimCfarResultMap));
xlabel('距离(m)');ylabel('速度(m/s)');zlabel('信号幅值');
title('rangeCFAR之后判决结果(峰值聚集前)');
xlim([0     Params.c*(Params.NSample-1)*Params.Fs/2/Params.freqSlope/Params.NSample]);
ylim([(-Params.NChirp_V/2)*Params.lambda/Params.Tc_V/Params.NChirp_V/2    (Params.NChirp_V/2 - 1)*Params.lambda/Params.Tc_V/Params.NChirp_V/2]);
figure;

% 进行峰值聚焦
[objDprIdx,objRagIdx] = peakFocus(rangeDimCfarResultMap);
objDprIdx(objDprIdx==0)=[]; %去掉后面的0
objRagIdx(objRagIdx==0)=[];
% 根据物体的点,计算速度和距离
objSpeed = ( objDprIdx - Params.NChirp_V/2 - 1)*Params.lambda/Params.Tc_V/Params.NChirp_V/2;
objRange = single(Params.c*(objRagIdx-1)*Params.Fs/2/Params.freqSlope/Params.NSample);
plot(objRange,objSpeed,'*r');
xlabel('距离(m)');ylabel('速度(m/s)');
title('峰值聚集后');
xlim([0     Params.c*(Params.NSample-1)*Params.Fs/2/Params.freqSlope/Params.NSample]);
ylim([(-Params.NChirp_V/2)*Params.lambda/Params.Tc_V/Params.NChirp_V/2    (Params.NChirp_V/2 - 1)*Params.lambda/Params.Tc_V/Params.NChirp_V/2]);



% 如果有物体则进行角度FFT
if(~isempty(objDprIdx))
    % 进行多普勒补偿
     processingChain_dopplerCompensation(objDprIdx,objRagIdx,objSpeed)

    dataSet.angleFFTOut = processingChain_angleFFT(1,objDprIdx,objRagIdx);
end







%% 2D CFAR部分，暂时不用
% [cfar_radarCubeData,tMap] = ac_cfar(12,4,10,4,0.0001,dataSet.sumFFT2D_radarCubeData);


% mesh(X,Y,abs(reshape(dataSet.radarCubeData(:,1,:),Params.NChirp,Params.NSample)));
%三维视图
%调整坐标系,标出距离和速度
% [X,Y] = meshgrid(Params.c*(0:Params.NSample-1)*Params.Fs/2/Params.freqSlope/Params.NSample, ...
%                 (-Params.NChirp/2:Params.NChirp/2 - 1)*Params.lambda/Params.Tc/Params.NChirp/2);    
% mesh(X,Y,(dataSet.sumFFT2D_radarCubeData));
%           title("log前");  figure;
% mesh(X,Y,20*log10(dataSet.sumFFT2D_radarCubeData));
% xlabel('距离(m)');ylabel('速度(m/s)');zlabel('信号幅值');
% title('2D FFT(4个天线相加后)');
% xlim([0     Params.c*(Params.NSample-1)*Params.Fs/2/Params.freqSlope/Params.NSample]);
% ylim([(-Params.NChirp/2)*Params.lambda/Params.Tc/Params.NChirp/2    (Params.NChirp/2 - 1)*Params.lambda/Params.Tc/Params.NChirp/2]);
% hold on;
% 
% mesh(X,Y,20*log10(tMap));
% xlabel('距离(m)');ylabel('速度(m/s)');zlabel('信号幅值');
% title('门限视图 & 2dfft');
% xlim([0     Params.c*(Params.NSample-1)*Params.Fs/2/Params.freqSlope/Params.NSample]);
% ylim([(-Params.NChirp/2)*Params.lambda/Params.Tc/Params.NChirp/2    (Params.NChirp/2 - 1)*Params.lambda/Params.Tc/Params.NChirp/2]);
% figure;
% 
% 
% mesh(X,Y,(cfar_radarCubeData));
% xlabel('距离(m)');ylabel('速度(m/s)');zlabel('判决结果1或0');
% title('2D FFT cfar后');
% xlim([0     Params.c*(Params.NSample-1)*Params.Fs/2/Params.freqSlope/Params.NSample]);
% ylim([(-Params.NChirp/2)*Params.lambda/Params.Tc/Params.NChirp/2    (Params.NChirp/2 - 1)*Params.lambda/Params.Tc/Params.NChirp/2]);





% rangeProfileData = dataSet.radarCubeData(1, 1, :);
% chFFT = rangeProfileData(:);    %把三维数据变成一列数据
% channelData = 20 * log10 (abs(chFFT));
% plot(channelData);

%% 函数

% 加载一帧的数据并且执行距离FFT
% 输入:frameIdx - 帧索引(从1开始)
% 输出:dataSet.rawFrameData - 要进行FFT的三维矩阵
%      dataSet.radarCubeData - 执行距离FFT后的三维矩阵
function dp_updateFrameData(frameIdx)
    global Params;
    global dataSet;
    
    rawDataComplex = dp_loadOneFrameData(Params.fid_rawData,Params.dataSizeOneFrame,frameIdx);
    dataSet.rawDataUint16 = uint16(rawDataComplex);
    % time domain data y value adjustments
    %把那种65534的数表示为负数-2
    %如果二进制的数首位为1(即大于等于2的15次方)那么这个数是负数，再减去65536
    timeDomainData = rawDataComplex - ( rawDataComplex >=2.^15).* 2.^16;
    frameData = dp_reshape2LaneLVDS(timeDomainData);

    %咱们目前iqSwapSel = 0
    %frameData的第一列I加1i*第二列Q变为一个复数(其中1i表示一个虚数单位,为了以防定义了变量i混淆)
    frameCplx = frameData(:,1) + 1i*frameData(:,2); 

    %初始化存放多个帧的三位数组 (Nchirp,NChan,NSample)
    frameComplex = single(zeros(Params.NChirp, Params.NChan, Params.NSample));

    % 咱们的chInterleave = 1 
    % non-interleave data
    % 把之前存放复数的矩阵frameCplx reshape为NSample*NChan行,NChirp列的二维数组,再取转置
    % 因为之前的frameCplx是以每个chirp为'周期'的
    % 这样的temp矩阵就是 行:某个chirp的数据(Rx按顺序排列的ADC采样数据) 列:0_0-255 1_0-255 ```3_0-255 
    temp = reshape(frameCplx, [Params.NSample * Params.NChan, Params.NChirp]).';
    % 把temp中的每一行reshape为NSample*NChan的二维矩阵,转置后再放入frameComplex的NChirp维度中
    for chirp=1:Params.NChirp                            
    frameComplex(chirp,:,:) = reshape(temp(chirp,:), [Params.NSample, Params.NChan]).';
    end
    % 保存到全局变量
    dataSet.rawFrameData = frameComplex;
end

%从bin文件中加载一帧的数据
% 输入:fid_rawData - bin文件句柄
%      dataSizeOneFrame - 一帧的数据大小(字节)
%      frameIdx - 帧索引
% 输出:rawData - 一帧的数据
function [rawData] = dp_loadOneFrameData(fid_rawData, dataSizeOneFrame, frameIdx)
    % 重新定位文件位置
    % 从头开始往后偏移(frameIdx-1)*dataSizeOneFrame个字节
    fseek(fid_rawData,(frameIdx-1)*dataSizeOneFrame,'bof');
    
    try   %如果文件打开成功
    %data:读取的数据 
    %count:读取数据的数量 
    %[M,N]读取数据到M*N的数组里，数据按列存放 262144是一帧的数个数 
    %262144 = 每个chirpADC采样数*2(I&Q)*Rx数*chirp个数 
    %uint16 数据格式为十六位无符号整型
    %这里有个疑问:uint16输出的结果是和例程输出的一样，但是并没有像手算的那样把连续2个aabb颠倒为bbaa
    %我也不知道为什么要uint16=>single,可能是为了加快计算速度吧,single是比double计算速度快的
    rawData = fread(fid_rawData,dataSizeOneFrame/2,'uint16=>single');
    catch
        error("error reading binary file");
    end
    fclose(fid_rawData);
end

% 把数据变为复数形式
function [frameData] = dp_reshape2LaneLVDS(rawData)
    % Convert 2 lane LVDS data to one matrix 
    rawData4 = reshape(rawData, [4, length(rawData)/4]);    %把rawData矩阵变为[M,N]的矩阵,按列排列(4个为1列)牛X牛X！！

    rawDataI = reshape(rawData4(1:2,:), [], 1);     %rawData4(1:2,:)表示取矩阵rawData4的1,2行(1,2行放的都是I数据)
                                                    %取了1,2行之后再reshape为一列I数据
                                                    
    rawDataQ = reshape(rawData4(3:4,:), [], 1);     %rawData4(3:4,:)表示取3,4行,再把3,4行数据reshape为一列Q数据
    
    frameData = [rawDataI, rawDataQ];   %把上面的一列I数据和Q数据合并为一个n行2列的矩阵(一帧的数据，128个chirp)
end


% 分离出虚拟天线的数据
function  dp_separateVirtualRx()
    global Params;
    global dataSet;
    
    NChirp = Params.NChirp;
    NChan = Params.NChan;
%     NRangeBin = Params.NSample;
    % 把chirp坐标为2 4 6....偶数的找出来放到对应的虚拟天线5678的位置
    for i = 2:2:NChirp      
        for j = 1:NChan
            dataSet.rawFrameData(i-1,j+4,:) = dataSet.rawFrameData(i,j,:);
        end
    end
    
    % 提取完毕后删除chirp坐标为2 4 6...偶数的行   这样会导致chirp数减半变为每一帧有原来一般的个数比如128/2
    dataSet.rawFrameData(2:2:NChirp,:,:) = [];
end
% 对数据进行距离FFT
function [radarCubeData] = processingChain_rangeFFT(rangeWinType)
    global Params;
    global dataSet;
    
    NChirp = Params.NChirp_V;
    NChan = Params.NChan_V;
    NRangeBin = Params.NSample;
    
    % 窗函数选择
    switch rangeWinType
        case 1 %hann
            win = hann(NRangeBin);
        case 2 %blackman
            win = blackman(NRangeBin);
        otherwise
            win = rectwin(NRangeBin);
    end
    %初始化一个三维数组存放FFT结果
    radarCubeData = single(zeros(NChirp,NChan,NRangeBin));
    for chirpIdx=1:NChirp
        for chIdx = 1: NChan
            frameData(1,:) = dataSet.rawFrameData(chirpIdx,chIdx,:);    
            frameData = fft(frameData .* win', NRangeBin);      %进行NRangeBin点的range-FFT
            radarCubeData(chirpIdx,chIdx,:) = frameData(1,:);   %把FFT的结果放入radarCubeData
        end
    end
end

% 对数据进行多普勒FFT
function [radarCubeData] = processingChain_dopplerFFT(rangeWinType)
    global Params;
    global dataSet;
    
    NChirp = Params.NChirp_V;
    NChan = Params.NChan_V;
    NRangeBin = Params.NSample;
    
    % 窗函数选择
    switch rangeWinType
        case 1 %hann
            win = hann(NChirp);
        case 2 %blackman
            win = blackman(NChirp);
        otherwise
            win = rectwin(NChirp);
    end
    %初始化一个三维数组存放FFT结果
    radarCubeData = single(zeros(NChirp,NChan,NRangeBin));
    for chIdx=1:NChan
        for rangeIdx = 1: NRangeBin
            frameData(1,:) = dataSet.radarCubeData(:,chIdx,rangeIdx);    
            frameData = fftshift(fft(frameData .* win', NChirp));      %进行NChirp点的doppler-FFT
                                                                       %fftshift移动零频点到频谱中间---不是很懂
            radarCubeData(:,chIdx,rangeIdx) = frameData(1,:);   %把FFT的结果放入radarCubeData
        end
    end
end

% 输入:NTrainRange - 距离维度训练单元数量
% NGuardRange - 距离维度保护单元数量
% NTrainDoppler - 多普勒维度
% NGuardDoppler
% PFA - 期望虚晶率
% inputRDM - range doppler map
% 输出: cfar_radarCubeData - 判决结果
% tMap - 门限图
function [cfar_radarCubeData,tMap] = ac_cfar(NTrainRange,NGuardRange,NTrainDoppler,NGuardDoppler,PFA,inputRDM)
    global Params; 
    cfar_radarCubeData = zeros(Params.NChirp_V,Params.NSample);   %存放cfar后的结果,
    tMap = NaN*zeros(Params.NChirp_V,Params.NSample);                 %存放门限图
    %限制CUT的范围
    for i = NTrainRange+NGuardRange+1 : Params.NSample-NTrainRange-NGuardRange
        for j = NTrainDoppler+NGuardDoppler+1 : Params.NChirp_V-NTrainDoppler-NGuardDoppler
            %存放每次迭代的训练单元的噪音值
            noiseLevel = zeros(1,1);
            %在指定范围内(guardCell之外)遍历当前CUT训练单元
            for p = i-NTrainRange-NGuardRange : i+NTrainRange+NGuardRange
                for q = j-NTrainDoppler-NGuardDoppler : j+NTrainDoppler+NGuardDoppler
                    if(abs(i-p)>NGuardRange || abs(j-q)>NGuardDoppler)
                        %因为咱们的RDM横坐标(列)是range,纵坐标(行)是doppler
                        noiseLevel = noiseLevel + (inputRDM(q,p));  %第q行,第p列,
                    end
                end
            end
            
            %根据噪音的平均值计算门限
            totalNTrain = (2*(NTrainRange+NGuardRange)+1)*(2*(NTrainDoppler+NGuardDoppler)+1)-(2*NGuardRange+1)*(2*NGuardDoppler+1);
            a = totalNTrain*((PFA^(-1/totalNTrain))-1);
            threshold = (noiseLevel / totalNTrain);
            threshold =  a*threshold;
            
            tMap(j,i) = threshold;
            CUT = inputRDM(j,i);
            
           %判决
            if(CUT < threshold)
                cfar_radarCubeData(j,i) = 0;
            else
                cfar_radarCubeData(j,i) = 1; 
            end
            
        end
    end
end

function [cfar1D_Arr,threshold] = ac_cfar1D(NTrain,NGuard,PFA,inputArr)
    cfar1D_Arr = zeros(size(inputArr));
    threshold = zeros(size(inputArr));

    totalNTrain = 2*(NTrain);
    a = totalNTrain*((PFA^(-1/totalNTrain))-1);
    %求平均值
    for i = NTrain+NGuard+1:length(inputArr)-NTrain-NGuard
        avg = mean([inputArr((i-NTrain-NGuard):(i-NGuard-1))   inputArr((i+NGuard+1):(i+NTrain+NGuard))]);
        threshold(1,i) = a.*avg;
        %根据threshold比较
        if(inputArr(i) < threshold(i))
            cfar1D_Arr(i) = 0;
        else
            cfar1D_Arr(i) = 1;
        end
    end
    
end
 
% 输入: inputCfarResMat - 进行峰值聚焦的二维矩阵,即进行range维CFAR判决后得到的结果矩阵
% 输出: row - 物体的行坐标(对应速度)
% column - 物体的列坐标(对应距离)
function [row,column] = peakFocus(inputCfarResMat)
    global dataSet;
    global Params;
    j = 1;
    row = zeros([1 Params.NSample]);
    column = zeros([1 Params.NSample]);
    [d,r] = find(inputCfarResMat==1);   %寻找进行range维cfar后的判决为1的坐标
    for i = 1 : length(d)
        peakRow = d(i);
        peakColumn = r(i);
        peak = dataSet.sumFFT2D_radarCubeData(peakRow,peakColumn);  %待验证的峰值
        % 在附近的3*3矩阵中的数进行比较,如果中间的数是最大值,就判定为1  tip:现在知道变量名太长的后果了 ┭┮﹏┭┮
        % 把附近的8+1个数取出来
        % 根据之前进行的2次cfar,因为有TrainCell和GuardCell，所以不会碰到边缘
        tempArr =[dataSet.sumFFT2D_radarCubeData(peakRow-1,peakColumn-1) , dataSet.sumFFT2D_radarCubeData(peakRow-1,peakColumn) ,  dataSet.sumFFT2D_radarCubeData(peakRow-1,peakColumn+1), ...
                  dataSet.sumFFT2D_radarCubeData(peakRow,peakColumn-1)   ,                     peak                             ,  dataSet.sumFFT2D_radarCubeData(peakRow,peakColumn+1), ...
                  dataSet.sumFFT2D_radarCubeData(peakRow+1,peakColumn-1) , dataSet.sumFFT2D_radarCubeData(peakRow+1,peakColumn) ,  dataSet.sumFFT2D_radarCubeData(peakRow+1,peakColumn+1)] ;    
        truePeak = max(tempArr);     % 寻找最大值
        if(truePeak == peak)         %如果中间的是最大值就保存当前的坐标
            row(j) = peakRow;
            column(j) = peakColumn;
            j = j+1;
        end
    end
end


% 多普勒补偿
% 输入：
% objDprIndex - 物体对应的多普勒维坐标
% objRagIndex - 物体对应的距离维坐标
% Vest - 雷达测出的物体速度
function processingChain_dopplerCompensation(objDprIndex,objRagIndex,Vest)
    global Params;
    global dataSet;
    for n = 1:length(objDprIndex)
        deltePhi = 4*pi*Vest(n)*Params.Tc/Params.lambda;  %多普勒相移
        %进行补偿
        dataSet.radarCubeData(objDprIndex(n),5:8,objRagIndex(n)) = dataSet.radarCubeData(objDprIndex(n),5:8,objRagIndex(n)).*exp(-1i*1*deltePhi);
    end
    
end

% 输入：rangeWinType - 窗函数选择
% objDprIndex - 物体对应的多普勒维坐标
% objRagIndex - 物体对应的距离维坐标
% 输出：angleFFTOut - 进行角度FFT输出的结果（二维，行 - 第几个物体，列 - FFT结果）
function [angleFFTOut] = processingChain_angleFFT(rangeWinType,objDprIndex,objRagIndex)

    global Params;
    global dataSet;
    angleFFTNum = 180;
    %NChirp = Params.NChirp_V;
    NChan = Params.NChan_V;
    %NRangeBin = Params.NSample;
    NObject = length(objDprIndex);  %检测的物体数目
    % 窗函数选择
    switch rangeWinType
        case 1 %hann
            win = hann(NChan);
        case 2 %blackman
            win = blackman(NChan);
        otherwise
            win = rectwin(NChan);
    end
    %初始化一个二维数组存放FFT结果
    angleFFTOut = single(zeros(NObject,angleFFTNum));
    for i=1:NObject
        frameData(1,:) = dataSet.radarCubeData(objDprIndex(i),:,objRagIndex(i));   
        % 有个问题，这里到底是几个点的FFT
        % 几个点，结果差别不大，另外不用加窗！！！！
        frameFFTData = fftshift(fft(frameData, angleFFTNum));      %进行NChan点的angler-FFT
                                                                   %fftshift移动零频点到频谱中间---不是很懂
                                                                   %有+Π和-Π
        angleFFTOut(i,:) = frameFFTData(1,:);   %把FFT的结果放入angleFFTOut
        
%         figure;plot(abs(angleFFTOut(i,:)));
        maxIndex= find(abs(angleFFTOut(i,:))==max(abs(angleFFTOut(i,:))));
        angle = asin((maxIndex - angleFFTNum/2 - 1)*2/angleFFTNum) * 180/pi;
        fprintf("object%d angle:%.2f°\n",i,angle);
    end

end


