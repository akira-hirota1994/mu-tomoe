
%%% Read table("~.csv")
%data_table = readtable("test6.csv");
data_table = readtable("simultaneous_event_candidates_table_MU+Tomo-e6.csv");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION SETTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FILE OUTPUT(.pngの出力
% !sh mv_imagefile.sh %shスクリプト実行 
option_setting_image = 0; %yes:1, no:0


%ANALIZE SINGLE IMAGE FILE
option_setting_single = 1; %yes:1, no:0



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get file name and angle info %%
ID = data_table(:,"ID");ID = table2array(ID);
image_filename = "sub" + ID + ".png";
rotimage_filename = "sub" + ID + "_rot.png";
photometry_filename = "sub" + ID + "_phot.png";
ID ="d_fits/sub" + ID + ".fits";
angle = data_table(:,"Angle");
angle = table2array(angle);
angle_rad = deg2rad(angle);

if option_setting_single == 1
    disp "AAA";
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% angular velocity [deg/sec] %%%
azimuth_1 = data_table(:,"azimuth_start_at_kiso"); %az of start point
azimuth_1 = table2array(azimuth_1);
azimuth_2 = data_table(:,"azimuth_end_at_kiso"); %az of end point
azimuth_2 = table2array(azimuth_2);
elevation_1 = data_table(:,"elevation_start_at_kiso"); %el of start point
elevation_1 = table2array(elevation_1);
elevation_2 = data_table(:,"elevation_end_at_kiso"); %el of end point
elevation_2 = table2array(elevation_2);
duration = data_table(:,"duration"); %duration
duration = table2array(duration);
flight_distance = sqrt(abs(azimuth_1 - azimuth_2).^2 + ...
    abs(elevation_1 - elevation_2).^2);
angular_velocity = flight_distance./duration;
sensor_fov = 1.19/3600; %tomoeの1pixelの視野(1.19"×1.19")
through_pixel = angular_velocity./sensor_fov;
through_pixel = 0.5*through_pixel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data setting %%%
center_x = data_table(:,"X");
center_y = data_table(:,"Y");
center_x = table2array(center_x);
center_y = table2array(center_y);
slope_angle = 90 - angle;
slope = tan(deg2rad(slope_angle));
y_intercept = center_y - (slope.*center_x);
x1 = data_table(:,"x_1");
x2 = data_table(:,"x_2");
x1 = table2array(x1);
x2 = table2array(x2);
y1 = slope.*x1 + y_intercept;
y2 = slope.*x2 + y_intercept; 
x_r = [x1,x2] ;y_r = [y1,y2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ループ処理（回転、rangeを求めて加算、ガウシアン）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% とりあえずy軸に±20pixel範囲をとる %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% imageの範囲を出る場合はずらす  %%%%%%%%%%%%
a = 20*cos(deg2rad(90- slope_angle));
b = 20*sin(deg2rad(90- slope_angle));
w1y1 = y1 - b; w1x1 = x1 + a;
w2y1 = y1 + b; w2x1 = x1 - a;
w1y2 = y2 - b; w1x2 = x2 + a;
w2y2 = y2 + b; w2x2 = x2 - a;
for i = 1:size(ID)
    if angle(i,1) < 90
        while (w2x1(i,1) < 0) || (w1y1(i,1) < 0)
            x1(i,1) = x1(i,1) + 1;
            y1(i,1) = slope(i,1)*x1(i,1) + y_intercept(i,1);
            w1y1(i,1) = y1(i,1) - b(i,1); w1x1(i,1) = x1(i,1) + a(i,1);
            w2y1(i,1) = y1(i,1) + b(i,1); w2x1(i,1) = x1(i,1) - a(i,1);
        end
        while (w2y2(i,1) > 282) || (w1x2(i,1) > 500)
            x2(i,1) = x2(i,1) - 1;
            y2(i,1) = slope(i,1)*x2(i,1) + y_intercept(i,1);
            w1y2(i,1) = y2(i,1) - b(i,1); w1x2(i,1) = x2(i,1) + a(i,1);
            w2y2(i,1) = y2(i,1) + b(i,1); w2x2(i,1) = x2(i,1) - a(i,1);
        end
    elseif angle(i,1) > 90
        while (w2y1(i,1) > 282) || (w1x1(i,1) < 0)
            x1(i,1) = x1(i,1) + 1;
            y1(i,1) = slope(i,1)*x1(i,1) + y_intercept(i,1);
            w1y1(i,1) = y1(i,1) - b(i,1); w1x1(i,1) = x1(i,1) + a(i,1);
            w2y1(i,1) = y1(i,1) + b(i,1); w2x1(i,1) = x1(i,1) - a(i,1);
        end
        while (w2x2(i,1) > 500) || (w1y2(i,1) < 0)
            x2(i,1) = x2(i,1) - 1;
            y2(i,1) = slope(i,1)*x2(i,1) + y_intercept(i,1);
            w1y2(i,1) = y2(i,1) - b(i,1); w1x2(i,1) = x2(i,1) + a(i,1);
            w2y2(i,1) = y2(i,1) + b(i,1); w2x2(i,1) = x2(i,1) - a(i,1);
        end
    end
end    
l1 = [w1x1,w2x1] ;l2 = [w1y1,w2y1];
l3 = [w1x1,w1x2] ;l4 = [w1y1,w1y2];
l5 = [w2x1,w2x2] ;l6 = [w2y1,w2y2];
l7 = [w1x2,w2x2] ;l8 = [w1y2,w2y2];
sz = size(ID);
LoopNum = 1871; %waitbar
t = 0; %waitbar
flux = zeros(sz); %途中から再開する場合はこの二行消す
line_intensity = zeros(sz);
h = waitbar(0,'running');
for i = 1:sz
    t = t + 1; %waitbar
    h = waitbar(i/LoopNum,h,['Number of Loop = ',num2str(t)]); %waitbar
    fitsfile_name = ID(i,1);
    data = fitsread(fitsfile_name); % reading fits file
    data = data(:,:,1); % 差分後の画像
    rotate_data = imrotate(data, -1*angle(i,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% midpoint %%%
    position_x = center_x(i,1) - size(data,2)/2;
    position_y = center_y(i,1) - size(data,1)/2;
    rotate_x = position_x*cos(angle_rad(i,1)) - position_y*sin(angle_rad(i,1));
    rotate_y = position_x*sin(angle_rad(i,1)) + position_y*cos(angle_rad(i,1));
    position_rotate_x = size(rotate_data,2)/2 + rotate_x;
    position_rotate_y = size(rotate_data,1)/2 + rotate_y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 両端の点  %%%
    range_x1 = x1(i,1) - size(data,2)/2;
    range_x2 = x2(i,1) - size(data,2)/2;
    range_y1 = y1(i,1) - size(data,1)/2;
    range_y2 = y2(i,1) - size(data,1)/2;
    rotate_x1 = range_x1*cos(angle_rad(i,1)) - range_y1*sin(angle_rad(i,1));
    rotate_y1 = range_x1*sin(angle_rad(i,1)) + range_y1*cos(angle_rad(i,1));
    rotate_x2 = range_x2*cos(angle_rad(i,1)) - range_y2*sin(angle_rad(i,1));
    rotate_y2 = range_x2*sin(angle_rad(i,1)) + range_y2*cos(angle_rad(i,1));
    range_rotate_x1 = size(rotate_data,2)/2 + rotate_x1;
    range_rotate_y1 = size(rotate_data,1)/2 + rotate_y1;
    range_rotate_x2 = size(rotate_data,2)/2 + rotate_x2;
    range_rotate_y2 = size(rotate_data,1)/2 + rotate_y2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 矩形の頂点  %
    range_w1x1 = w1x1(i,1) - size(data,2)/2;
    range_w2x1 = w2x1(i,1) - size(data,2)/2;
    range_w1x2 = w1x2(i,1) - size(data,2)/2;
    range_w2x2 = w2x2(i,1) - size(data,2)/2;
    range_w1y1 = w1y1(i,1) - size(data,1)/2;
    range_w2y1 = w2y1(i,1) - size(data,1)/2;
    range_w1y2 = w1y2(i,1) - size(data,1)/2;
    range_w2y2 = w2y2(i,1) - size(data,1)/2;
    rotate_w1x1 = range_w1x1*cos(angle_rad(i,1)) - range_w1y1*sin(angle_rad(i,1));
    rotate_w2x1 = range_w2x1*cos(angle_rad(i,1)) - range_w2y1*sin(angle_rad(i,1));
    rotate_w1x2 = range_w1x2*cos(angle_rad(i,1)) - range_w1y2*sin(angle_rad(i,1));
    rotate_w2x2 = range_w2x2*cos(angle_rad(i,1)) - range_w2y2*sin(angle_rad(i,1));
    rotate_w1y1 = range_w1x1*sin(angle_rad(i,1)) + range_w1y1*cos(angle_rad(i,1));
    rotate_w2y1 = range_w2x1*sin(angle_rad(i,1)) + range_w2y1*cos(angle_rad(i,1));
    rotate_w1y2 = range_w1x2*sin(angle_rad(i,1)) + range_w1y2*cos(angle_rad(i,1));
    rotate_w2y2 = range_w2x2*sin(angle_rad(i,1)) + range_w2y2*cos(angle_rad(i,1));
    range_rotate_w1x1 = size(rotate_data,2)/2 + rotate_w1x1;
    range_rotate_w2x1 = size(rotate_data,2)/2 + rotate_w2x1;
    range_rotate_w1x2 = size(rotate_data,2)/2 + rotate_w1x2;
    range_rotate_w2x2 = size(rotate_data,2)/2 + rotate_w2x2;
    range_rotate_w1y1 = size(rotate_data,1)/2 + rotate_w1y1;
    range_rotate_w2y1 = size(rotate_data,1)/2 + rotate_w2y1;
    range_rotate_w1y2 = size(rotate_data,1)/2 + rotate_w1y2;
    range_rotate_w2y2 = size(rotate_data,1)/2 + rotate_w2y2;
    
    rotl1 = [range_rotate_w1x1,range_rotate_w2x1];
    rotl2 = [range_rotate_w1y1,range_rotate_w2y1];
    rotl3 = [range_rotate_w1x1,range_rotate_w1x2];
    rotl4 = [range_rotate_w1y1,range_rotate_w1y2];
    rotl5 = [range_rotate_w2x1,range_rotate_w2x2];
    rotl6 = [range_rotate_w2y1,range_rotate_w2y2];
    rotl7 = [range_rotate_w1x2,range_rotate_w2x2];
    rotl8 = [range_rotate_w1y2,range_rotate_w2y2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% フラックス求める %%%
    if w1x1(i,1) < 0 || w2x1(i,1) < 0 || w1x2(i,1) > 500 || w2x2(i,1) > 500 ||...
            w1y1(i,1) < 0 || w2y1(i,1) > 282 || w1y2(i,1) < 0 || w2y2(i,1) > 282 ...
            || w1x1(i,1) == w1x2(i,1) || w1y1(i,1) == w1y2(i,1) ||...
            fitsfile_name == "d_fits/sub661714_118.fits" ||...
            fitsfile_name == "d_fits/sub705614_052.fits"
        
        flux(i,1) = 0;
        disp("This frame is out of range.");
    else
        data3 = rotate_data(range_rotate_w1y1:range_rotate_w1y2 ...
            ,range_rotate_w2x1:range_rotate_w1x1);
        y = mean(data3);
        x = range_rotate_w2x1:range_rotate_w1x1;  
        f = fit(x.',y.','gauss1');
        coeff = coeffvalues(f);
        fun = @(x) coeff(1)*exp(-((x-coeff(2))/coeff(3)).^2);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        option_setting_image = 1;
        plot(x,y); %ガウシアン画像
        title(fitsfile_name);
        hold on;
        plot(f,x,y); 
        hold off;
        grid on
        %saveas(gcf,photometry_filename(i,1));
        clf;
        %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% photometry %%
%2パターン(関数の積分,普通にフラックスで足し合わせる)
        fwhm = coeff(3)/sqrt(2);
        aperture_radius = 1.5*fwhm;
        aperture_x1 = position_rotate_x - aperture_radius;
        aperture_x2 = position_rotate_x + aperture_radius;
        integer = integral(fun,aperture_x1,aperture_x2);
        line_intensity(i,1) = integer/aperture_radius;%[ADU/pix]
        flux(i,1) = line_intensity(i,1)*through_pixel(i,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        imshow(rotate_data); %回転後の画像
        hold on
        title(fitsfile_name);
        scatter(range_rotate_x1,range_rotate_y1,50,"red","filled")
        scatter(range_rotate_x2,range_rotate_y2,50,"red","filled")
        line(rotl1,rotl2,"Color","green","Linewidth",2,"Linestyle","--")
        line(rotl3,rotl4,"Color","green","Linewidth",2,"Linestyle","--")
        line(rotl5,rotl6,"Color","green","Linewidth",2,"Linestyle","--")
        line(rotl7,rotl8,"Color","green","Linewidth",2,"Linestyle","--")
        hold off
        axis equal
        saveas(gcf,rotimage_filename(i,1));
        %}
    end
    %{
    clf;
    imshow(data);
    title(fitsfile_name); %回転前の画像
    hold on
    scatter(x1(i,1),y1(i,1),50,"red","filled")
    scatter(x2(i,1),y2(i,1),50,"red","filled")
    line(l1(i,:),l2(i,:),"Color","green","Linewidth",2,"Linestyle","--")
    line(l3(i,:),l4(i,:),"Color","green","Linewidth",2,"Linestyle","--")
    line(l5(i,:),l6(i,:),"Color","green","Linewidth",2,"Linestyle","--")
    line(l7(i,:),l8(i,:),"Color","green","Linewidth",2,"Linestyle","--")
    hold off
    axis equal
    saveas(gcf,image_filename(i,1));
    %}
end
close(h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 元のデータテーブルに加える %%%
photometry_data = horzcat(line_intensity,flux);
photometry_data = array2table(photometry_data,'VariableNames',...
   {'line_intensity','flux_meteor'});
data_table = horzcat(data_table,photometry_data);
%writetable(data_table,"test7.csv");
writetable(data_table,"simultaneous_event_candidates_table_MU+Tomo-e7.csv");
disp("Complete!")