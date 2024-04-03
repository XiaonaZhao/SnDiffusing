% Detailed process for cases

%% input image sequence
%selected sample 1
tifPath = 'H:\☆SnSe2\☆selected\2017-06-12 A5_SnSe_-1.1V_30s_50 fps_0.5M Na2SO4☆\';
[tifDir, ~] = ReadPicFileNames(tifPath);
frame.begin = 192;
frame.end = 323; % 323
frame.fps = 50;
validDir = tifDir(frame.begin:frame.end);  %changed

%selected sample 2
tifPath = 'H:\☆SnSe2\☆selected\2017-09-12 B6_SnSe_-1.3V_10s_0.5M Na2SO4_106fps\';
[tifDir, ~] = ReadPicFileNames(tifPath);
frame.begin = 248;
frame.end = 762; % 762
frame.fps = 106;
validDir = tifDir(frame.begin:frame.end);  %changed

%selected sample 3
tifPath = 'H:\☆SnSe2\☆selected\2018-12-05 A3_SnSe2_200fps_0.5M Na2SO4_-1.0V 10s☆\';
[tifDir, ~] = ReadPicFileNames(tifPath);
frame.begin = 548; % 543, 548
frame.end = 1755; % 1755
frame.fps = 200;
validDir = tifDir(frame.begin:frame.end);  %changed

%% image load
L = frame.end - frame.begin + 1;
tif0 = double(imread(fullfile(tifPath, tifDir(1).name)));
tif = cell(L, 1); % why not 'tif = cell(L, 1);'

for ii = 1:L
    temp = double(imread(fullfile(tifPath, validDir(ii).name))) - tif0;
    % tif{ii, 1} = imboxfilt(temp, 7);
    tif{ii, 1} = temp;
end
%%
%diffusion coefficient
%=================fit for Sample1\2\3===================
% figure('color', 'w')
saveRoute = 'H:\MATLAB\Result_SnSe2\Sample2';
radii = zeros(L, 1);
speed = 1/frame.fps;
D = zeros(L-1, 1);
dD = zeros(L-1, 1);
for ii = 2:L
    % ii = 400;
    % I = tif{ii, 1} - tif{ii-1, 1};
    I = tif{ii, 1};
    % imshow(I, 'DisplayRange', [], 'InitialMagnification', 'fit');
    I = imboxfilt(I, 3);
    % I1 = I.*(I>0);
    % imshow(I1, 'DisplayRange', [], 'InitialMagnification', 'fit');
    I2 = imbinarize(I,'adaptive');
    % imshow(I2, 'DisplayRange', [], 'InitialMagnification', 'fit');

    imLabel = bwlabel(I2);
    stats = regionprops(imLabel, 'Area', "MajorAxisLength", "MinorAxisLength");
    area = cat(1, stats.Area);
    index = find(area == max(area));% Find the index of the smallest connected domain
    I3 = ismember(imLabel, index);
    h = imshow(I3, 'DisplayRange', [], 'InitialMagnification', 'fit');

    drawnow; % =======draw image=======
    refreshdata(h);

    diameters = mean([stats(index).MajorAxisLength stats(index).MinorAxisLength]);
    radii(ii, 1) = diameters/2;
    T = (ii-1)*speed;
    D(ii-1, 1) = (((0.00005)/200*radii(ii, 1))^2)/(4*T);

    dR = radii(ii, 1) - radii(ii-1, 1);
    dD(ii-1, 1) = (((50e-6)/200*dR).^2)/(4*speed);

    % title({['\rm2-D Diffusion with \itD\rm = ',num2str(dD),' m^2s^-^1'];...
    %     ['\rmtime (\itt\rm) = ',num2str(roundn(T, -3)),' s']})
    % title({['\rm\itD\rm = ',num2str(dD),' m^2s^-^1, \itt\rm = ',num2str(roundn(T, -3)),' s']})
    title({['\itt\rm = ',num2str(T, '%.3f'),' s, \rm\itD\rm = ',...
        num2str(D(ii-1, 1), 3),' m^2s^-^1, \rm\itdD\rm = ' num2str(dD(ii-1, 1), 3), ' m^2s^-^1']})
    scalebar; % sample1, 0.14; sample2, 0.11; sample3, 0.07

    figPath = [saveRoute '\num_' num2str(ii, '%04d')];
    saveas(h, figPath, 'jpg')

end

%Fick's Law: R^2 = 4*D*t;
% dR = diff(radii);
% dt = speed;
% D = (((50e-6)/200*dR).^2)/(4*dt); % hmmt, pike with 1.5x convertor on 60x Object

% %% This is not important.
% load('H:\MATLAB\Result_SnSe2\Figures\DiffusionCoefficient.mat')
% 
% % 2017-06-12 A5_SnSe_-1.1V_30s_50 fps_0.5M Na2SO4☆
% figure('color', 'w');
% Log10D1 = log10(D1);
% plot(t1(1:end-1), Log10D1(1:end-1))
% 
% % 2017-09-12 B6_SnSe_-1.3V_10s_0.5M Na2SO4_106fps
% figure('color', 'w');
% Log10D2 = log10(D2);
% plot(t2(1:end-1), Log10D2(1:end-1))
% 
% 
% % 2018-12-05 A3_SnSe2_200fps_0.5M Na2SO4_-1.0V 10s☆ ====Sample3====
% figure('color', 'w');
% Log10D3 = log10(D3);
% plot(t3(1:end-1), Log10D3(1:end-1))
% 
% 
% % -1.0 V, -1.1 V, -1.3 V
% figure('color', 'w');
% plot(t3(1:end-1), Log10D3(1:end-1)) % -1.0 V
% hold on
% plot(t1(1:end-1), Log10D1(1:end-1))
% plot(t2(1:end-1), Log10D2(1:end-1))
% legend '-1.0 V' '-1.1 V' '-1.3 V'
% xlabel('Time (s)')
% ylabel('the log10 of real-time Diffusion coefficient (log10(m2/s))')
% hold off
% 
% %%
% %merge images, pde calculation
% load('center.mat')
% load('around.mat')
% 
% ucell = ucell1;
% nt = size(ucell1, 1);
% for ii = 1:nt
%     temp = ucell{ii, 1};
%     temp(201:300, 201:300) = ucell2{ii, 1};
%     ucell{ii, 1} = temp;
% end
% 
% 
% saveRoute = 'H:\MATLAB\Result_SnSe2\All';
% vis = 8;
% for it = 1:nt
%     z = ucell{it, 1};
%     h = imshow(z, 'DisplayRange',[], 'InitialMagnification', 'fit');
%     colormap parula
%     axis on
%     xlabel('Spatial co-ordinate (x) pixel \rightarrow')
%     ylabel('{\leftarrow} Spatial co-ordinate (y) pixel')
%     title({['\rm2-D Diffusion with \itD\rm  = ',num2str(vis/1e10),' m^2s^-^1'];...
%         ['\rmtime (\itt\rm) = ',num2str((it-1)*0.02),' s']})
%     scalebar
%     pause(0.1)
%     figPath = [saveRoute '\' num2str(it) '_all'];
%     saveas(h, figPath, 'jpg')
% end

%%
%output videos
picFolder = 'H:\MATLAB\Result_SnSe2\Sample1\';
[~, picNames] = ReadPicFileNames(picFolder); %====== pay attention to the format of pictures=======
video = VideoWriter('simulation1_lockedclim', 'MPEG-4'); % Prepare a .mp4 file
video.FrameRate = 4; % The same as the fps
open(video);
for ii = 1:size(picNames,1)  % The frame number for video, 1:10:size(picNames,1)
    frame = imread(fullfile(picFolder, picNames{ii}));
    writeVideo(video, frame);
end
close(video);