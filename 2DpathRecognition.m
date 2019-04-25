%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright: Copyright (c) 2018
%Created on 2018-6-16 
%Author:MengDa (github:pilidili)
%Version 1.0 
%Title: 2DpathRecognition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;
disp('Initializing...');
t1=clock;%初始化部分计时开始
%%%%%%%%导入和压缩图片%%%%%%%%
I_origin=imread('f.jpg');
I_Singelpic_Size=1024;
I_resize_Size=512;
I_origin=imresize(I_origin,I_resize_Size/I_Singelpic_Size);
I_origin_Size=size(I_origin);
spop_start=[I_origin_Size(2)-I_resize_Size,(I_origin_Size(1)-I_resize_Size)/2];%第一个摄像区原点坐标（Singel Picture Original Point）
I=I_origin(spop_start(2)+1:spop_start(2)+I_resize_Size,spop_start(1)+1:I_origin_Size(2),:);%读取第一个摄像区
%%%%%%%%开辟识别线端点存储空间%%%%%%%%
lpm=zeros(ceil(I_resize_Size/10),2,ceil(I_origin_Size(2)/I_resize_Size)+1);
%%%%%%%%开辟各个摄像区原点坐标存储空间%%%%%%%%
opm=zeros(ceil(I_origin_Size(2)/I_resize_Size)+1,2);
opm(1,:)=spop_start;
%%%%%%%%开辟各个摄像区处理时间的存储空间%%%%%%%%
time_m=zeros(1,ceil(I_origin_Size(2)/I_resize_Size)+1);
%%%%%%%%生成路径判断条件矩阵%%%%
unit_scale=35;%采样尺寸
x=[(1:(unit_scale-1))-0.5,(unit_scale-1)*ones(1,unit_scale-1),unit_scale-0.5-(1:(unit_scale-1))];
y=[floor(unit_scale/2)*ones(1,unit_scale-1),unit_scale/2-(1:(unit_scale-1)),-floor(unit_scale/2)*ones(1,unit_scale-1)];
k_switch=[y./x,-inf];%生成判断斜率区间
x_delta=[-(0:unit_scale-2),-(unit_scale-1)*ones(1,unit_scale),(0:unit_scale-2)-unit_scale+2];
y_delta=[-floor(unit_scale/2)*ones(1,unit_scale),(1:unit_scale-2)-floor(unit_scale/2),floor(unit_scale/2)*ones(1,unit_scale)];
P_delta=[x_delta',y_delta'];%生成区间对应坐标变化量
d_d=dcm(unit_scale);%生成斜率求解卷积矩阵
t2=clock;
time_int=etime(t2,t1);%初始化部分计时结束


clc
disp('calculating...');
t1=clock;%计算部分计时开始
%%%%%%%%起点定位%%%%%%%%
%%%% x坐标定位 %%%%
qtdc_size_x=32;%区域分裂采样尺寸
qtdc_storage_x=zeros(qtdc_size_x,I_resize_Size);%开辟区域分裂结果存储空间
for ii=1:I_resize_Size/qtdc_size_x %区域分裂算法沿x采样计算
    pp=qtdc_size_x*(ii-1)+1:qtdc_size_x*(ii);
    I_q=qtdecomp(I(1:qtdc_size_x,pp),0.05);
    qtdc_storage_x(:,pp)=I_q;
end
sum_x=sum(qtdc_storage_x); %区域分裂结果对行求和降维
sum_x=edge(sum_x,'roberts',1); %利用roberts函数（梯度算子）检测边界
sum_x=bwareaopen(sum_x,50); %利用连通域函数去除小面积的边界
x_start=find(sum_x,1,'first'); %搜索出起点x坐标


%%%% y坐标定位 %%%%
qtdc_size_y=32;
if x_start>qtdc_size_y  %%防错判断
    range=x_start-qtdc_size_y+1:x_start;
    qtdc_storage_y=zeros(I_resize_Size,qtdc_size_y);  %开辟区域分裂结果存储空间
    for ii=1:I_resize_Size/qtdc_size_y  %区域分裂算法沿y采样计算
        pp=qtdc_size_y*(ii-1)+1:qtdc_size_y*(ii);
        I_q_y=qtdecomp(I(pp,range),0.05);
        qtdc_storage_y(pp,:)=I_q_y;
    end
    sum_y=sum(qtdc_storage_y,2);%区域分裂结果对行求和降维

    y_sum_fft=fft(sum_y);  %低通滤波
    y_sum_fft(10:end)=0;
    sum_y=real(ifft(y_sum_fft));
else %%错误丢出
    clc;
    disp('---ERROR---');
    disp('Start point not found!');
    return;
end
[maxsum,y_start]=max(sum_y); %求极大值获得起点x坐标


%%%% 起点坐标采集与坐标系转换 %%%%
point_mini=[x_start,y_start];%设置起点(单摄像区域坐标系)
origin_big=spop_start;%设置当前摄像区域的原点（工作台坐标系）



jj=1;%摄像区域序号
while origin_big(1)> 0 && origin_big(2)> 0  && ...
        origin_big(1)<=(I_origin_Size(2)-I_resize_Size) && origin_big(2)<=(I_origin_Size(1)-I_resize_Size)
    
    I=I_origin(origin_big(2)+1:origin_big(2)+I_resize_Size,origin_big(1)+1:origin_big(1)+I_resize_Size,:);    
    lpm(1,:,jj)=point_mini+origin_big;
    ii=1;

    %%%%%%%%识别裂缝%%%%%%%%
    while point_mini(1)>spop_start(1)+x_start-2802 && point_mini(2)>unit_scale  && ...
            point_mini(1)<=I_resize_Size && point_mini(2)<=(I_resize_Size-floor(unit_scale/2))
        ii=ii+1;
        x_range=point_mini(1)-unit_scale+1:point_mini(1);
        y_range=point_mini(2)-floor((unit_scale-1)/2):point_mini(2)+floor(unit_scale/2);
        I_current=I(y_range,x_range,:);
        if point_mini(1)+origin_big(1)<1055 %分区域切换阈值
            I_current=1-im2bw(I_current,0.25);
            I_current(1,1)=1;
        else
            I_current=1-im2bw(I_current,0.4);
        end
        if sum(sum(I_current))==0 || (sum(sum(I_current.*[ones(unit_scale,ceil(unit_scale*2/3)),zeros(unit_scale,floor(unit_scale*1/3))]))==0)%%若错失曲线 向下（y方向）找曲线
            x_range=point_mini(1)-floor((unit_scale-1)/2):point_mini(1)+floor(unit_scale/2);
            y_range=point_mini(2):point_mini(2)+unit_scale-1;
            I_current=I_origin(y_range+origin_big(2),x_range+origin_big(1),:);
            I_current=imrotate(I_current, -90);
            I_current=1-im2bw(I_current,0.4);
            if sum(sum(I_current))==0 %%若仍错失曲线 向上（-y方向）找曲线
                x_range=point_mini(1)-floor((unit_scale-1)/2):point_mini(1)+floor(unit_scale/2);
                y_range=point_mini(2)-unit_scale+1:point_mini(2);
                I_current=I_origin(y_range+origin_big(2),x_range+origin_big(1),:);
                I_current=imrotate(I_current,90);
                I_current=1-im2bw(I_current,0.4);
                if sum(sum(I_current))==0  %%若依然错失曲线 报错退出
                    clc;
                    disp('---ERROR---');
                    disp('Missing the curve!');
                    return;
                else %向上找的情况
                    p=sum(sum(I_current.*d_d))/sum(sum(I_current));
                    Pixal_out=find(k_switch<p(1),1,'first');
                    pos_delta=fliplr(P_delta(Pixal_out,:)).*[-1,1];
                end
            else  %向下找的情况
                p=sum(sum(I_current.*d_d))/sum(sum(I_current));
                Pixal_out=find(k_switch<p(1),1,'first');
                pos_delta=fliplr(P_delta(Pixal_out,:)).*[1,-1];
            end
        else %正常情况
            p=sum(sum(I_current.*d_d))/sum(sum(I_current));
            Pixal_out=find(k_switch<p(1),1,'first');
            pos_delta=P_delta(Pixal_out,:);
        end
        point_mini=point_mini+pos_delta;
        lpm(ii,:,jj)=point_mini+origin_big;
    end
    t2=clock;
    time_m(jj)=etime(t2,t1);%计时一次
    t1=clock;
    
    origin_big=point_mini+origin_big-[I_resize_Size,I_resize_Size/2];%更新下一个循环区域的原点位置(工作台坐标系)
    point_mini=[I_resize_Size,I_resize_Size/2]; %更新下一个循环的初始点位置(单摄像区域坐标系)
    jj=jj+1;
    opm(jj,:)=origin_big;
end
opm=opm(1:end-1,:);

%%%%%%%%作图输出结果%%%%%%%%
clc
disp('plotting...');
t1=clock;%绘图部分计时开始
imshow(I_origin);
hold on;
for ii=1:size(lpm,3); %画出识别出的裂缝线段
    eee=find(lpm(:,1,ii)==0,1,'first')-1;
    line(lpm(1:eee,1,ii),lpm(1:eee,2,ii),'linewidth',3,'color',[0 0 1]);
    text(opm(ii,1)+I_resize_Size/2,opm(ii,2)-30,[num2str(time_m(ii)*1000),'ms'],'color',[1 1 0],'fontsize',13);
end
for ii=1:size(opm,1) %框出摄像头区域
    line([opm(ii,1)+1,opm(ii,1)+I_resize_Size,opm(ii,1)+I_resize_Size,opm(ii,1),opm(ii,1)+1],[opm(ii,2)+1,opm(ii,2)+1,opm(ii,2)+I_resize_Size,opm(ii,2)+I_resize_Size,opm(ii,2)+1],'color',[1 1 0])
    for jj=1:size(lpm,1) %框出每个采样的小区域
    line([lpm(jj,1,ii)-unit_scale+1,lpm(jj,1,ii)-unit_scale+1,lpm(jj,1,ii),lpm(jj,1,ii),lpm(jj,1,ii)-unit_scale+1],[lpm(jj,2,ii)-floor(unit_scale/2),lpm(jj,2,ii)+floor(unit_scale/2),lpm(jj,2,ii)+floor(unit_scale/2),lpm(jj,2,ii)-floor(unit_scale/2),lpm(jj,2,ii)-floor(unit_scale/2)],'color',[0 1 0])
    end
end
%%%%给出起点、终点坐标%%%%
sp=[spop_start(1)+x_start,spop_start(2)+y_start];
ep=point_mini+origin_big;
text(sp(1),sp(2),['起点坐标（',num2str(sp),'）'],'fontsize',15,'color',[1 0 1]);
text(ep(1),ep(2),['终点坐标(',num2str(ep),')'],'fontsize',15,'color',[0 0 1]);
%%%%粉红框框出起点位置%%%%
line([lpm(1,1,1)-floor(unit_scale/2),lpm(1,1,1)-floor(unit_scale/2),lpm(1,1,1)+floor(unit_scale/2),lpm(1,1,1)+floor(unit_scale/2),lpm(1,1,1)-floor(unit_scale/2)],[lpm(1,2,1)-floor(unit_scale/2),lpm(1,2,1)+floor(unit_scale/2),lpm(1,2,1)+floor(unit_scale/2),lpm(1,2,1)-floor(unit_scale/2),lpm(1,2,1)-floor(unit_scale/2)],'color',[1 0 1],'linewidth',5)
%%%%蓝色框框出终点位置%%%%
line([ep(1)-floor(unit_scale/2),ep(1)-floor(unit_scale/2),ep(1)+floor(unit_scale/2),ep(1)+floor(unit_scale/2),ep(1)-floor(unit_scale/2)],[ep(2)-floor(unit_scale/2),ep(2)+floor(unit_scale/2),ep(2)+floor(unit_scale/2),ep(2)-floor(unit_scale/2),ep(2)-floor(unit_scale/2)],'color',[0 0 1],'linewidth',5)
axis equal
t2=clock;
time_plot=etime(t2,t1);%绘图部分计时结束
title(['初始化部分总耗时:',num2str(sum(time_int)*1000),'ms, ','计算部分总耗时:',num2str(sum(time_m)*1000),'ms, ','绘图部分总耗时:',num2str(sum(time_plot)*1000),'ms'],'fontsize',15)


clc
disp('finished!');