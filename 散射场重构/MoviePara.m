classdef MoviePara
    % 这个类预设Movie S1视频中会用到的参数
    %   此处显示详细说明
    
    properties(Hidden)
        color                               % 大图背景颜色
        pstn                                % 整个图片的位置
        pstn1                               % 强度曲线图片的位置
        pstn2                               % SPRM视频的位置
        pstn3                               % 重构视频的位置
        caxis_SPRM                          % SPRM图片的colorbar动态范围
        caxis_rec                           % 重构图片的colorbar动态范围
        xtick                               % 曲线图的x轴刻度位置
        xticklabel                          % 曲线图的x轴刻度标志
        ytick                               % 曲线图的y轴刻度位置
        yticklabel                          % 曲线图的y轴刻度标志
        xlim                                % 曲线图x轴范围
        ylim                                % 曲线图y轴范围
        xlabel                              % x轴label
        ylabel                  
        
        int_spr                             % 颗粒SPR强度
        int_rec                             % 颗粒重构强度
        recpath                             % 重构图片路径
        sprpath                             % SPRM图片路径 
        rectiffs                            % 重构图片文件
        sprtiffs                            % SPRM图片文件
        SPR_BG
        
        x                                   % 第二个颗粒的x坐标
        y                                   % 第二个颗粒的y坐标
        frame                               % 第二个颗粒出现的帧数
    end
    
    methods
        function obj = MoviePara(recpath,sprpath,datapath)
            obj.color = [0.94,0.94,0.94];
            obj.pstn = [285.0000  149.0000  887.2000  560.0000];
            obj.pstn1 = [0.1200    0.6500    0.8000    0.3300];
            obj.pstn2 = [0.0500    0.0200    0.4500    0.4500];
            obj.pstn3 = [0.5000    0.0200    0.4500    0.4500];
            obj.caxis_SPRM = [-8000 10000];
            obj.caxis_rec = [0 2500];
            obj.xtick = linspace(0,10000,6);
            obj.xticklabel = {'0','20','40','60','80','100'};
            obj.ytick = [0,0.5,1];
            obj.yticklabel = {'0','0.5','1'};
            obj.xlim = [0 10000];
            obj.ylim = [0 1.35];
            obj.xlabel = 'Time (s)';
            obj.ylabel = 'Intensity (a.u.)';
            
            S = load(datapath);                                                     % 从datapath里面读取通过imageJ统计的数据，包括第二个颗粒的SPRM强度和重构强度
                                                                                    % 以及第二个颗粒的位置信息，用来在视频中标注其位置
            obj.x = S.pstn{2}(:,2);                                                 % 第二个颗粒的x坐标
            obj.y = S.pstn{2}(:,3);                                                 % 第二个颗粒的y坐标
            obj.frame = S.pstn{2}(:,4);                                             % 第二个颗粒出现的帧数
            int_rec = S.rec(:,2);
            int_spr = S.spr(:,2);
            obj.int_rec = int_rec/max(int_rec);
            [~,ind] = max(int_rec);
            obj.int_spr = int_spr/int_spr(ind);
            
            obj.recpath = recpath;
            obj.sprpath = sprpath;
            rectiffs = dir(fullfile(recpath,'*.*tif*'));
            obj.rectiffs = struct_sort(rectiffs);
            sprtiffs = dir(fullfile(sprpath,'*.*tif*'));
            obj.SPR_BG = double(imread(fullfile(sprpath,sprtiffs(1250).name)));
            sprtiffs(1:1250) = [];
            obj.sprtiffs = sprtiffs;
        end
        function plot(obj,ii)
            recImage = double(imread(fullfile(obj.recpath,obj.rectiffs(ii).name)));
            sprImage = double(imread(fullfile(obj.sprpath,obj.sprtiffs(ii).name))) - obj.SPR_BG;
            recData = obj.int_rec(1:ii);
            sprData = obj.int_spr(1:ii);
            FigureCreate(obj,ii,recData,sprData,recImage,sprImage);
        end
    end
    
end
%% 局部函数
function FigureCreate(obj,ii,recData,sprData,recImage,sprImage)
    figure('Position',obj.pstn,'Color',obj.color,'Visible','on')
    ax1 = axes('Position',obj.pstn1);
    hold on
    plot(ax1,recData,'linewidth',2,'color',[0.64,0.38,0.77])
    plot(ax1,sprData,'LineWidth',2,'color',[1.00,0.89,0.07])
    legend({'Reconstruction','Raw Data'});
    legend('boxoff')
    ax1.XLim = obj.xlim;
    ax1.YLim = obj.ylim;
    ax1.XTick = obj.xtick;
    ax1.XTickLabel = obj.xticklabel;
    ax1.YTick = obj.ytick;
    ax1.YTickLabel = obj.yticklabel;
    xlabel(ax1,obj.xlabel,'FontSize',15,'FontWeight','bold');
    ylabel(ax1,obj.ylabel,'FontSize',15,'FontWeight','bold');
    ax1.FontSize = 15;
    ax1.FontWeight = 'bold';
    ax1.Box = 'on';

    ax2 = axes('Position',obj.pstn2,'yDir','reverse');
    hold on
    axis(ax2,'off','equal')
    imagesc(ax2,sprImage)
    colormap(ax2,'gray')
    caxis(ax2,obj.caxis_SPRM)
    title(ax2,'Raw images','FontSize',15,'FontWeight','bold');
    line(ax2,[100 100+5000/74],[380 380],'LineWidth',6,'Color','w');  
    text(72,420,'5 \mum','Color','w','FontSize',20,'FontWeight',"bold")
    text(480,420,[num2str((ii-1)/100) ' s'],'Color','w','FontSize',20,'FontWeight',"bold")

    ax3 = axes('Position',obj.pstn3,'yDir','reverse');
    hold on
    axis(ax3,'off','equal')
    imagesc(ax3,recImage)
    colormap(ax3,'gray')
    caxis(ax3,obj.caxis_rec)
    title(ax3,'Reconstructed images','FontSize',15,'FontWeight','bold');
    line(ax3,[100 100+5000/74],[380 380],'LineWidth',6,'Color','w');
    text(72,420,'5 \mum','Color','w','FontSize',20,'FontWeight',"bold")
    text(480,420,[num2str((ii-1)/100) ' s'],'Color','w','FontSize',20,'FontWeight',"bold")

    if ismember(ii,obj.frame)     
        k = find(ii == obj.frame);
        plot(ax2,obj.x(k),obj.y(k),'o','MarkerSize',13,'MarkerEdgeColor','r','LineWidth',3);      
        plot(ax3,obj.x(k),obj.y(k),'o','MarkerSize',13,'MarkerEdgeColor','r','LineWidth',3);      
    end
    
end

function sorted_struct = struct_sort(raw_struct)
% 本函数用来对结构体进行重新排序
% 通过dir读入的文件结构体，排列顺序完全是按字符排序的，也不怪他，毕竟是机器，但
% 没按原始文件顺序读入就很难受了，这里重新按照数字大小排序
    L = length(raw_struct);
    names = zeros(L,1);
    for ii = 1:length(raw_struct)
        temp = raw_struct(ii).name;
        temp = split(temp,'.');
        names(ii) = str2double(temp{1});
    end
    [~,ind] = sort(names);
    sorted_struct = raw_struct(ind);
end