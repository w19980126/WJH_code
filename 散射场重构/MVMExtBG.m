function BG = MVMExtBG(TiffInformation,varargin)
% median-value method (MVM) to extract stationary background
%
% syntax
%   BG = MVMExtBG(TiffInformation)
%       从TiffInformation提示的信息中提取图片并进行背景提取
%   BG = MVMExtBG(TiffInformation,startFrame,endFrame)
%       从TiffInformation提示的从startFrame到endFrame的图片中提取背景
%
% input:
%       TiffInformation         图片信息，可以是图片保存文件夹，也可以是图片序列
%       varargin                输入两个变量，第一个是起始帧，第二个是结束帧
% output:
%       BG                      提取到的背景

    if ischar(TiffInformation)
        tiffpath = TiffInformation;
        tiffs = dir(fullfile(tiffpath,'*.tif*'));
        if nargin == 1
            startFrame = 1;
            endFrame = length(tiffs);
        elseif nargin == 3
            startFrame = varargin{1};
            endFrame = varargin{2};
        end
        Num = endFrame - startFrame + 1;    % 用于求BG的图片张数
        I0 = double(imread(fullfile(tiffpath,tiffs(1).name)));
        [m,n] = size(I0);
        RawFigs = zeros([m n Num]);
        for ii = startFrame:endFrame
            RawFigs(:,:,ii) = double(imread(fullfile(tiffpath,tiffs(ii).name)));
        end
    else
        RawFigs = TiffInformation;
        if nargin == 1
            startFrame = 1;
            endFrame = size(RawFigs,3);
        elseif nargin == 3
            startFrame = varargin{1};
            endFrame = varargin{2};
        end
        RawFigs = RawFigs(:,:,startFrame:endFrame);
    end
    BG = median(RawFigs,3);
end

