function [ d_d ] = dcm( scale )
%%%%本函数用于此次建模中生成估算图像线条斜率的卷积矩阵
%%%%scale为需求矩阵尺寸，矩阵为方阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright: Copyright (c) 2018
%Created on 2018-6-16 
%Author:MengDa (github:pilidili)
%Version 1.0 
%Title: function dcm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear;clc;
format long
% scale=35;
d_d=zeros(scale);
count=zeros(scale);
x=[-(1:scale-2),-(scale-1)*ones(1,scale),(0:scale-3)-scale+2];
y=[-floor(scale/2)*ones(1,scale-1),(1:scale-2)-floor(scale/2),floor(scale/2)*ones(1,scale-1)];
d=y./x;
for ii=1:floor(size(d,2)/2)
    for jj=1:(scale-1)
        ppp=ceil(scale/2-jj*d(ii));
        if ppp<scale && ppp>0
            count(ppp,scale-jj)=count(ppp,scale-jj)+1;
            d_d(ppp,scale-jj)=d_d(ppp,scale-jj)+d(ii);
        end
    end
end


count(1:floor(scale/2),scale)=1;
d_d(1:floor(scale/2),scale)=1e10;

count(floor(scale/2)+2:end,:)=flipud(count(1:floor(scale/2),:));
count(floor(scale/2)+1,:)=1;
d_d(floor(scale/2)+2:end,:)=-flipud(d_d(1:floor(scale/2),:));
d_d(floor(scale/2)+1,:)=d(ceil(size(d,2)/2));


count(count==0)=inf;
d_d=d_d./count;
d_d=atan(d_d);
% d_d(2:end-1,2:end-1)=0;
% imagesc(count);

% x=[-(0:scale-2),-(scale-1)*ones(1,scale),(0:scale-2)-scale+2]+scale-0.5;
% y=[-floor(scale/2)*ones(1,scale),(1:scale-2)-floor(scale/2),floor(scale/2)*ones(1,scale)]+scale/2;
% for ii=0:scale
%     line([0,scale],[ii ii])
% end
% for ii=0:scale
%     line([ii ii],[0,scale])
% end
% for ii=1:scale*3-2
%     line([scale-0.5,x(ii)],[scale/2,y(ii)])
% end


% axis equal