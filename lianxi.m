function varargout = lianxi(varargin)
% LIANXI MATLAB code for lianxi.fig
%      LIANXI, by itself, creates a new LIANXI or raises the existing
%      singleton*.
%
%      H = LIANXI returns the handle to a new LIANXI or the handle to
%      the existing singleton*.
%
%      LIANXI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LIANXI.M with the given input arguments.
%
%      LIANXI('Property','Value',...) creates a new LIANXI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lianxi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lianxi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lianxi

% Last Modified by GUIDE v2.5 20-Dec-2019 17:59:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lianxi_OpeningFcn, ...
                   'gui_OutputFcn',  @lianxi_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before lianxi is made visible.
function lianxi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lianxi (see VARARGIN)

% Choose default command line output for lianxi
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lianxi wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lianxi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function caiyangpinlv_Callback(hObject, eventdata, handles)
% hObject    handle to caiyangpinlv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of caiyangpinlv as text
%        str2double(get(hObject,'String')) returns contents of caiyangpinlv as a double


% --- Executes during object creation, after setting all properties.
function caiyangpinlv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to caiyangpinlv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function caiyangdianshu_Callback(hObject, eventdata, handles)
% hObject    handle to caiyangdianshu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of caiyangdianshu as text
%        str2double(get(hObject,'String')) returns contents of caiyangdianshu as a double


% --- Executes during object creation, after setting all properties.
function caiyangdianshu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to caiyangdianshu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in shengka.
function shengka_Callback(hObject, eventdata, handles)
% hObject    handle to shengka (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shengka
%声卡部分
set(handles.shengka,'value',1);  %声卡按钮打开
set(handles.luyinshijian,'enable','on');%录音按钮打开
set(handles.kaishiluyin,'enable','on');%开始录音按钮打开
%函数发生器部分
set(handles.hanshufashengqi,'value',0);%函数发生器按钮关闭，以下一系列与声卡无关的均关闭
set(handles.hundie,'enable','off');%混叠关闭
set(handles.xuanzeboxingcaidan,'enable','off');
set(handles.pinlv1,'enable','off');
set(handles.fudu1,'enable','off');
set(handles.xiangwei1,'enable','off');
set(handles.shengchengboxing,'enable','off');
%WAV文件部分
set(handles.wavwenjian,'value',0);
set(handles.edit66,'enable','off');
set(handles.shengdao,'enable','off');
set(handles.dakaiwenjian,'enable','off');




function luyinshijian_Callback(hObject, eventdata, handles)
% hObject    handle to luyinshijian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of luyinshijian as text
%        str2double(get(hObject,'String')) returns contents of luyinshijian as a double


% --- Executes during object creation, after setting all properties.
function luyinshijian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to luyinshijian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in kaishiluyin.
function kaishiluyin_Callback(hObject, eventdata, handles)
% hObject    handle to kaishiluyin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Fs = str2double(get(handles.caiyangpinlv,'String'));%取采样频率并取值
y = audiorecorder(Fs, 16, 1);%创建录音机对象
recordblocking(y,str2double(get(handles.luyinshijian,'string')));%与记录相同，但是直到记录完成才返回控制
 handles.boxing =  getaudiodata(y);%将记录的音频数据返回到MATLAB工作区
 handles.boxing =  handles.boxing*100;
handles.Nn  = size(handles.boxing);%取采样点数
set(handles.caiyangdianshu,'string',num2str(handles.Nn));%将生成的采样点数值放到采样点数里
guidata(hObject,handles);
axes(handles.axes1);%选第一个框
plot(handles.boxing);%画出音频波形
title('初始波形图');





% --- Executes during object creation, after setting all properties.
function kaishiluyin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kaishiluyin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pinyufenxi.
function pinyufenxi_Callback(hObject, eventdata, handles)
% hObject    handle to pinyufenxi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fs1 = str2double(get(handles.caiyangpinlv,'string'));
N1 = str2double(get(handles.caiyangdianshu,'string'));

if get(handles.suoyoudian,'value')==1 %分析所有点
    kaishidian1 = 1;  %开始的采样点数
    jieshudian1 = N1;  %结束的采样点数
else
    kaishidian1 = str2double(get(handles.kaishidian,'string'));% 所需的开始点
    jieshudian1 = str2double(get(handles.jieshudian,'string'));%所需的结束点
end

if str2double(get(handles.jieshudian,'string'))>N1   %确定声卡分析结束点范围
 msgbox('结束点输入数值应小于生成的采样点数，请输入正确分析范围！');
  return; 
 end 

sample = handles.boxing(kaishidian1:jieshudian1);%提取出待分析样本，将其存入sample中
f = linspace(0,Fs1/2,(jieshudian1-kaishidian1+1)/2);%以采样频率做离散化的间隔
Y = fft(handles.boxing,jieshudian1-kaishidian1+1);%
[C,I] = max(abs(Y));%获得幅值最大的点并取其对应的下标I
set(handles.zhouqi2,'string',1/f(I));%将值放入周期里
set(handles.pinlv2,'string',f(I));%将值放入频率里
Y = Y(1:(jieshudian1-kaishidian1+1)/2);%只取Y的前半部分

axes(handles.axes2);%选图2
plot(f,2*sqrt(Y.*conj(Y)));%画幅值谱
axes(handles.axes3);%选图3
plot(f,angle(Y));%相位谱
axes(handles.axes4);%选图4
plot(f,real(Y));%实频谱
axes(handles.axes5);%选图5
plot(f,imag(Y));%虚频谱
axes(handles.axes6);%选图6
plot(f,abs(Y).^2);%功率谱
chazhi=jieshudian1-kaishidian1+1;%布莱克曼窗区间选取
 juxing1=blackman(chazhi);%布莱克曼窗函数
 y=juxing1.*sample ;   %窗函数相乘
 [H,W]=dtft(y,jieshudian1);%dtft
  axes(handles.axes7);%选图7
  plot(W/pi,20*log(abs(H)));%加窗后的db特性曲线。

 %获取基波
Xabs=abs(Y);
Xabs(1) = 0; %直流分量置0
for i= 1:(jieshudian1-kaishidian1+1)%基波变换
[Amax,index]=max(Xabs);
if(Xabs(index-1) > Xabs(index+1))
a1 = Xabs(index-1) / Xabs(index);
r1 = 1/(1+a1);
k01 = index -1;
else
a1 = Xabs(index) / Xabs(index+1);
r1 = 1/(1+a1);
k01 = index;
end
Fn = (k01+r1-1)*Fs1/(jieshudian1-kaishidian1+1); %基波频率
An = 2*pi*r1*Xabs(k01)/((jieshudian1-kaishidian1+1)*sin(r1*pi)); %基波幅值
Pn = phase(Y(k01))-pi*r1; %基波相角 单位弧度
Pn = mod(Pn(1),pi);
end
set(handles.jibofuzhi,'string',An);%基波幅值输出
set(handles.jibopinlv,'string',Fn);%基波频率输出
set(handles.jiboxiangwei,'string',Pn);%基波相位输出
t=[1:0.001:(jieshudian1-kaishidian1+1)];%基波时间范围
jibohanshu=An*sin((2*pi/Fn)*t+Pn);%基波函数
axes(handles.axes8);%选图8
plot(jibohanshu);%基波特性曲线。
 

% --- Executes on button press in shiyufenxi.
function shiyufenxi_Callback(hObject, eventdata, handles)
% hObject    handle to shiyufenxi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Fs0 = str2double(get(handles.caiyangpinlv,'String'));%取采样频率并取值
N0 = str2double(get(handles.caiyangdianshu,'string'));%获取采样点数

T = zeros(size(handles.boxing));%过零监测的周期向量
amp = zeros(size(handles.boxing));%过零监测的幅值向量
n = zeros(size(handles.boxing));%过零点数

if  str2double(get(handles.kaishidian,'String'))>Fs0  %确定开始点范围
    msgbox('开始点输入数值应小于采样频率，请输入正确分析范围！');
    return; 
end

if str2double(get(handles.jieshudian,'string'))>N0   %确定声卡分析结束点范围
 msgbox('结束点输入数值应小于生成的采样点数，请输入正确分析范围！');
  return; 
 end 


if get(handles.hanshufashengqi,'value')==1 && get(handles.suoyoudian,'value')==1
     kaishidian = 1;  %开始的采样点数
     jieshudian = N0;  %结束的采样点数
end
 if get(handles.suoyoudian,'value')==1 && get(handles.shengka,'value')==1
     kaishidian = 1;  %开始的采样点数
     jieshudian = handles.Nn;  %结束的采样点数
 end
 if get(handles.suoyoudian,'value')==1 && get(handles.wavwenjian,'value')==1
     kaishidian = 1;  %开始的采样点数
     jieshudian = handles.Nnn;  %结束的采样点数
 end
 if get(handles.suoyoudian,'value')==0
     kaishidian = str2double(get(handles.kaishidian,'String'));% 所需的开始点
     jieshudian = str2double(get(handles.jieshudian,'string'));%所需的结束点
 end

% if  jieshudian>(2*Fs0+1)
%     msgbox('输入结束点的数值应小于二倍的采样频率，请输入正确分析范围！');
%     return;
% end
 
axes(handles.axes1);%选第一个框
 plot(handles.boxing);%画出波形
% plot([kaishidian:jieshudian],handles.boxing);%画出波形

xianglinagweizhi = 1;  %过零点数
for  i=kaishidian+2:jieshudian-2  %过零点数计算
      if handles.boxing(i-2)<0 && handles.boxing(i)>=0 && handles.boxing(i-1)<0 && handles.boxing(i+2)>0 && handles.boxing(i+1)>0 %判断是否是过零点
         if handles.boxing(i)==0
            n(xianglinagweizhi) = i; %过零点计数
         else
               n(xianglinagweizhi) = (i-0.001);
         end
         xianglinagweizhi = xianglinagweizhi+1; %过零点数加1
      end
end

xianglinagweizhi = xianglinagweizhi-1;
%计算频率和周期
for i = 1:xianglinagweizhi-1
    T(i) = n(i+1)-n(i);%相邻两个过零点相减
end

pinlv1=(Fs0/mean(T(1:xianglinagweizhi-1)));%计算频率
set(handles.zhouqi,'string',1/pinlv1);%周期值
set(handles.pinlv,'string',num2str(pinlv1));%频率估计值

max = 0;%预设最大值
min = 0;%预设最小值
for  i=kaishidian+2:jieshudian-2  %最大最小值计算
        if max<handles.boxing(i) %如果下一个的值比max大
            max = handles.boxing(i);%值给max
            amp(xianglinagweizhi) = max;%幅值等于最大加最小除以2
        end
        if min>handles.boxing(i)%如果下一个的值比min小
            min = handles.boxing(i);%值给min
        end
end

% amp(xianglinagweizhi) = (max-min)/2;%幅值等于最大加最小除以2
set(handles.fuzhi,'string',num2str(mean((max-min)/2)));%幅值输出
phase = 2*pi*(1-(n(1:xianglinagweizhi-1)-1)./T(1:xianglinagweizhi-1)+floor((n(1:xianglinagweizhi-1)-1)./T(1:xianglinagweizhi-1)));%计算相位
set(handles.xiangwei,'string',num2str(mean(phase)));%相位输出
set(handles.fengzhi,'string',max);%峰值输出
set(handles.junzhi,'string',mean(handles.boxing(kaishidian:jieshudian)));%均值输出
set(handles.junfangzhi,'string',mean(handles.boxing(kaishidian:jieshudian).^2));%均方值输出
set(handles.fangcha,'string',std(handles.boxing(kaishidian:jieshudian))^2);%方差输出




function zhouqi_Callback(hObject, eventdata, handles)
% hObject    handle to zhouqi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zhouqi as text
%        str2double(get(hObject,'String')) returns contents of zhouqi as a double


% --- Executes during object creation, after setting all properties.
function zhouqi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zhouqi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pinlv_Callback(hObject, eventdata, handles)
% hObject    handle to pinlv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pinlv as text
%        str2double(get(hObject,'String')) returns contents of pinlv as a double


% --- Executes during object creation, after setting all properties.
function pinlv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pinlv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fuzhi_Callback(hObject, eventdata, handles)
% hObject    handle to fuzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fuzhi as text
%        str2double(get(hObject,'String')) returns contents of fuzhi as a double


% --- Executes during object creation, after setting all properties.
function fuzhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fuzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xiangwei_Callback(hObject, eventdata, handles)
% hObject    handle to xiangwei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xiangwei as text
%        str2double(get(hObject,'String')) returns contents of xiangwei as a double


% --- Executes during object creation, after setting all properties.
function xiangwei_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xiangwei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fengzhi_Callback(hObject, eventdata, handles)
% hObject    handle to fengzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fengzhi as text
%        str2double(get(hObject,'String')) returns contents of fengzhi as a double


% --- Executes during object creation, after setting all properties.
function fengzhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fengzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function junzhi_Callback(hObject, eventdata, handles)
% hObject    handle to junzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of junzhi as text
%        str2double(get(hObject,'String')) returns contents of junzhi as a double


% --- Executes during object creation, after setting all properties.
function junzhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to junzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function junfangzhi_Callback(hObject, eventdata, handles)
% hObject    handle to junfangzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of junfangzhi as text
%        str2double(get(hObject,'String')) returns contents of junfangzhi as a double


% --- Executes during object creation, after setting all properties.
function junfangzhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to junfangzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fangcha_Callback(hObject, eventdata, handles)
% hObject    handle to fangcha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fangcha as text
%        str2double(get(hObject,'String')) returns contents of fangcha as a double


% --- Executes during object creation, after setting all properties.
function fangcha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fangcha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kaishidian_Callback(hObject, eventdata, handles)
% hObject    handle to kaishidian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kaishidian as text
%        str2double(get(hObject,'String')) returns contents of kaishidian as a double


% --- Executes during object creation, after setting all properties.
function kaishidian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kaishidian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function jieshudian_Callback(hObject, eventdata, handles)
% hObject    handle to jieshudian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jieshudian as text
%        str2double(get(hObject,'String')) returns contents of jieshudian as a double


% --- Executes during object creation, after setting all properties.
function jieshudian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jieshudian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in suoyoudian.
function suoyoudian_Callback(hObject, eventdata, handles)
% hObject    handle to suoyoudian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')==0.0     
    set(handles.kaishidian,'Enable','on');    %开始点按钮开启
    set(handles.jieshudian,'Enable','on');  %结束点按钮开启
else
   set(handles.kaishidian,'String','1','Enable','off'); %开始点 1
    set(handles.jieshudian,'String',get(handles.caiyangdianshu,'String'),'Enable','off');%结束点 采样点数
end
% Hint: get(hObject,'Value') returns toggle state of suoyoudian



%function kaishidian_Callback(hObject, eventdata, handles)
% hObject    handle to kaishidian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kaishidian as text
%        str2double(get(hObject,'String')) returns contents of kaishidian as a double



function zhouqi2_Callback(hObject, eventdata, handles)
% hObject    handle to zhouqi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zhouqi2 as text
%        str2double(get(hObject,'String')) returns contents of zhouqi2 as a double


% --- Executes during object creation, after setting all properties.
function zhouqi2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zhouqi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pinlv2_Callback(hObject, eventdata, handles)
% hObject    handle to pinlv2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pinlv2 as text
%        str2double(get(hObject,'String')) returns contents of pinlv2 as a double


% --- Executes during object creation, after setting all properties.
function pinlv2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pinlv2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
 function shiyufenxi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shiyufenxi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in shengchengboxing.
function shengchengboxing_Callback(hObject, eventdata, handles)
% hObject    handle to shengchengboxing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shengchengboxing

Fs2 = str2double(get(handles.caiyangpinlv,'String'));
N2 = str2double(get(handles.caiyangdianshu,'String'));

xuanzeboxing = get(handles.xuanzeboxingcaidan,'value');%根据选择菜单中的选项选择那个波形
pinlv2 = str2double(get(handles.pinlv1,'String'));%获取幅度值
fudu2 = str2double(get(handles.fudu1,'String'));%获取频率值
xiangwei2= str2double(get(handles.xiangwei1,'String'));%获取相位值
x = linspace(0,N2/Fs2,N2)';
switch xuanzeboxing
    case 1
        y = fudu2*sin(2*pi*pinlv2*x+xiangwei2);
    case 2
        y = fudu2*sign(sin(2*pi*pinlv2*x+xiangwei2));
    case 3
        y = fudu2*sawtooth(2*pi*pinlv2*x+xiangwei2,0.5);
    case 4
        y = fudu2*sawtooth(2*pi*pinlv2*x+xiangwei2);
    case 5
        y = fudu2*(2*rand(size(x))-1);
    otherwise
        errordlg('请选择正确波形');
end

if get(handles.hundie,'value')==0.0
    handles.boxing = y;
else
    handles.boxing = handles.boxing+y;
end

guidata(hObject,handles);
axes(handles.axes1);%选图1做图
plot(handles.boxing);

% --- Executes on button press in hanshufashengqi.
function hanshufashengqi_Callback(hObject, eventdata, handles)
% hObject    handle to hanshufashengqi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hanshufashengqi

%声卡部分
set(handles.shengka,'value',0);  %声卡按钮打开
set(handles.luyinshijian,'enable','off');%录音按钮打开
set(handles.kaishiluyin,'enable','off');%开始录音按钮打开
%函数发生器部分
set(handles.hanshufashengqi,'value',1);%函数发生器按钮关闭，以下一系列与声卡无关的均关闭
set(handles.hundie,'enable','on');%混叠关闭
set(handles.xuanzeboxingcaidan,'enable','on');
set(handles.pinlv1,'enable','on');
set(handles.fudu1,'enable','on');
set(handles.xiangwei1,'enable','on');
set(handles.shengchengboxing,'enable','on');
%WAV文件部分
set(handles.wavwenjian,'value',0);
set(handles.edit66,'enable','off');
set(handles.shengdao,'enable','off');
set(handles.dakaiwenjian,'enable','off');


function fudu1_Callback(hObject, eventdata, handles)
% hObject    handle to fudu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fudu1 as text
%        str2double(get(hObject,'String')) returns contents of fudu1 as a double


% --- Executes during object creation, after setting all properties.
function fudu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fudu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in xuanzeboxingcaidan.
function xuanzeboxingcaidan_Callback(hObject, eventdata, handles)
% hObject    handle to xuanzeboxingcaidan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xuanzeboxingcaidan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xuanzeboxingcaidan



% --- Executes during object creation, after setting all properties.
function xuanzeboxingcaidan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xuanzeboxingcaidan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pinlv1_Callback(hObject, eventdata, handles)
% hObject    handle to pinlv1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pinlv1 as text
%        str2double(get(hObject,'String')) returns contents of pinlv1 as a double


% --- Executes during object creation, after setting all properties.
function pinlv1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pinlv1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xiangwei1_Callback(hObject, eventdata, handles)
% hObject    handle to xiangwei1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xiangwei1 as text
%        str2double(get(hObject,'String')) returns contents of xiangwei1 as a double


% --- Executes during object creation, after setting all properties.
function xiangwei1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xiangwei1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hundie.
function hundie_Callback(hObject, eventdata, handles)
% hObject    handle to hundie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hundie


% --- Executes on button press in dakaiwenjian.
function dakaiwenjian_Callback(hObject, eventdata, handles)
% hObject    handle to dakaiwenjian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[wenjianming00,wenjianlujing]=uigetfile('*.wav','wavfile');%通过函数将输入的wav文件找到
set(handles.wenjianming,'string',wenjianming00);
wenjianming0 = audioread(get(handles.wenjianming,'String'));%选择文件
shengdao = str2double(get(handles.shengdao,'String'));%选择声道
handles.boxing = wenjianming0(:,shengdao);%将波形数据存入
 [handles.Nnn zz] = size(handles.boxing);
 set(handles.caiyangdianshu,'string',num2str(handles.Nnn));%将数据放到采样点数里
guidata(hObject,handles);%获取数据
axes(handles.axes1);%选图1做图
plot(handles.boxing);%画图
% set(handles.caiyangpinlv,'string',Fs);

% --- Executes on button press in wavwenjian.
function wavwenjian_Callback(hObject, eventdata, handles)
% hObject    handle to wavwenjian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wavwenjian
%声卡部分
set(handles.shengka,'value',0);  %声卡按钮打开
set(handles.luyinshijian,'enable','off');%录音按钮打开
set(handles.kaishiluyin,'enable','off');%开始录音按钮打开
%函数发生器部分
set(handles.hanshufashengqi,'value',0);%函数发生器按钮关闭，以下一系列与声卡无关的均关闭
set(handles.hundie,'enable','off');%混叠关闭
set(handles.xuanzeboxingcaidan,'enable','off');
set(handles.pinlv1,'enable','off');
set(handles.fudu1,'enable','off');
set(handles.xiangwei1,'enable','off');
set(handles.shengchengboxing,'enable','off');
%WAV文件部分
set(handles.wavwenjian,'value',1);
set(handles.edit66,'enable','on');
set(handles.shengdao,'enable','on');
set(handles.dakaiwenjian,'enable','on');




function wenjianming_Callback(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit66 as text
%        str2double(get(hObject,'String')) returns contents of edit66 as a double


% --- Executes during object creation, after setting all properties.
function wenjianming_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function shengdao_Callback(hObject, eventdata, handles)
% hObject    handle to shengdao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shengdao as text
%        str2double(get(hObject,'String')) returns contents of shengdao as a double


% --- Executes during object creation, after setting all properties.
function shengdao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shengdao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function jibopinlv_Callback(hObject, eventdata, handles)
% hObject    handle to jibopinlv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jibopinlv as text
%        str2double(get(hObject,'String')) returns contents of jibopinlv as a double


% --- Executes during object creation, after setting all properties.
function jibopinlv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jibopinlv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function jibofuzhi_Callback(hObject, eventdata, handles)
% hObject    handle to jibofuzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jibofuzhi as text
%        str2double(get(hObject,'String')) returns contents of jibofuzhi as a double


% --- Executes during object creation, after setting all properties.
function jibofuzhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jibofuzhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function jiboxiangwei_Callback(hObject, eventdata, handles)
% hObject    handle to jiboxiangwei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jiboxiangwei as text
%        str2double(get(hObject,'String')) returns contents of jiboxiangwei as a double


% --- Executes during object creation, after setting all properties.
function jiboxiangwei_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jiboxiangwei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
