%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section : SETTING  

%%%  FILE OUTPUT(.pdfの出力) 
option_show_image = 1; %yes:1, no:0
option_save_image = 1; %yes:1, no:0

%%%  ANALIZE SINGLE IMAGE FILE
option_single_image = 1; %yes:1, no:0
ID_number = "572534_007"; %571226_095 , 571342_058(recheck),572534_007
pixel_width = 10; %default is 20pix

%%% ANALYZIE SETTING
re_analyze = 0; %yes:1, no:0
version = "main"; % 'test' or 'main'


%%%%% TO DO LIST %%%%%%
%・2回目以降で全部やるやつはしばらく放置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section : 色々値引っ張ってくる

if option_single_image == 1 && re_analyze == 1 %単体
    if version == "test"
        data_table = readtable("test6.csv");
        recheck_list =readtable("test_recheck_photometry.csv");
        data_table_overwrite = readtable("test8.csv");
    elseif version == "main"
        data_table = readtable("simultaneous_event_candidates_table_MU+Tomo-e6.csv");
        recheck_list =readtable("recheck_photometry.csv");
        data_table_overwrite = readtable("simultaneous_event_candidates_table_MU+Tomo-e8.csv");
    end
    ID = "d_fits/sub" + ID_number + ".fits";
    ID_list = data_table(:,"ID");ID_list = table2array(ID_list);
    ID_recheck_list = recheck_list(:,"ID");
    ID_recheck_list = table2array(ID_recheck_list);
    image_filename = "sub" + ID_number + ".pdf";
    %rotimage_filename = "sub" + ID_number + "_rot.png";
    %photometry_filename = "sub" + ID_number + "_phot.png";
    
    index = [];index_recheck = [];%IDの該当するindexナンバーを取得
    for i = 1:size(ID_list)
        if ID_number == ID_list(i,1)
            index = vertcat(index,i);
        end
    end
    for i = 1:size(ID_recheck_list)
        if ID_number == ID_recheck_list(i,1)
            index_recheck = vertcat(index_recheck,i);
        end
    end
    angle = data_table(index,"Angle");angle = table2array(angle);
    angle_rad = deg2rad(angle);
    center_x = data_table(index,"X");center_y = data_table(index,"Y");
    center_x = table2array(center_x);center_y = table2array(center_y);
    slope_angle = 90 - angle;
    slope = tan(deg2rad(slope_angle));
    y_intercept = center_y - (slope.*center_x);
    x1 = data_table(index,"x_1");x2 = data_table(index,"x_2");
    x1 = table2array(x1);x2 = table2array(x2);
    y1 = slope.*x1 + y_intercept;y2 = slope.*x2 + y_intercept; 
    x_r = [x1,x2] ;y_r = [y1,y2];
    azimuth_1 = data_table(index,"azimuth_start_at_kiso");
    azimuth_1 = table2array(azimuth_1);
    azimuth_2 = data_table(index,"azimuth_end_at_kiso");
    azimuth_2 = table2array(azimuth_2);
    elevation_1 = data_table(index,"elevation_start_at_kiso");
    elevation_1 = table2array(elevation_1);
    elevation_2 = data_table(index,"elevation_end_at_kiso");
    elevation_2 = table2array(elevation_2);
    duration = data_table(index,"duration");
    duration = table2array(duration);
    loop_number = size(index);
%{    
elseif option_single_image == 0 && re_analyze == 1; %二回目-全部
    data_table = readtable("test6.csv");
    recheck_list =readtable("recheck_photometry.csv");
    data_table_overwrite = readtable("test8.csv");   
    
    ID_list = data_table(:,"ID");ID_list = table2array(ID_list);
    ID_recheck_list = recheck_list(:,"ID");ID_recheck_list = table2array(ID_recheck_list);
    image_filename = "sub" + ID_recheck_list + ".pdf";
    ID ="d_fits/sub" + ID_recheck_list + ".fits";
    pixel_width = recheck_list(:,"pixel_width_modified");
    pixel_width = table2array(pixel_width);
    pixel_width = pixel_width(1,1);
    
    %IDの該当するindexナンバーを取得
    index = [];
    for i = 1:size(ID_recheck_list)
        for j = 1:size(ID_list)
            if table2array(ID_recheck_list(i,1)) == table2array(ID_list(j,1))
                index = vertcat(index,j);
            end
        end
    end
    
    angle = data_table(index,"Angle");angle = table2array(angle);
    angle_rad = deg2rad(angle);
    center_x = data_table(index,"X");center_y = data_table(index,"Y");
    center_x = table2array(center_x);center_y = table2array(center_y);
    slope_angle = 90 - angle;
    slope = tan(deg2rad(slope_angle));
    y_intercept = center_y - (slope.*center_x);
    x1 = data_table(index,"x_1");x2 = data_table(index,"x_2");
    x1 = table2array(x1);x2 = table2array(x2);
    y1 = slope.*x1 + y_intercept;y2 = slope.*x2 + y_intercept; 
    x_r = [x1,x2] ;y_r = [y1,y2];

    azimuth_1 = data_table(index,"azimuth_start_at_kiso");
    azimuth_1 = table2array(azimuth_1);
    azimuth_2 = data_table(index,"azimuth_end_at_kiso");
    azimuth_2 = table2array(azimuth_2);
    elevation_1 = data_table(index,"elevation_start_at_kiso");
    elevation_1 = table2array(elevation_1);
    elevation_2 = data_table(index,"elevation_end_at_kiso");
    elevation_2 = table2array(elevation_2);
    duration = data_table(index,"duration");
    duration = table2array(duration);
    loop_number = size(index);
%}    
elseif option_single_image == 0 && re_analyze == 0 %初回
    if version == "test"
        data_table = readtable("test6.csv");
    elseif version == "main"
        data_table = readtable("simultaneous_event_candidates_table_MU+Tomo-e6.csv");
    end
    ID = data_table(:,"ID");ID = table2array(ID);
    ID_list = data_table(:,"ID");ID_list = table2array(ID_list);
    image_filename = "sub" + ID + ".pdf";
    %rotimage_filename = "sub" + ID + "_rot.png";
    %photometry_filename = "sub" + ID + "_phot.png";
    ID ="d_fits/sub" + ID + ".fits";
    
    angle = data_table(:,"Angle");angle = table2array(angle);
    angle_rad = deg2rad(angle);
    center_x = data_table(:,"X");center_y = data_table(:,"Y");
    center_x = table2array(center_x);center_y = table2array(center_y);
    slope_angle = 90 - angle;
    slope = tan(deg2rad(slope_angle));
    y_intercept = center_y - (slope.*center_x);
    x1 = data_table(:,"x_1");x2 = data_table(:,"x_2");
    x1 = table2array(x1);x2 = table2array(x2);
    y1 = slope.*x1 + y_intercept;y2 = slope.*x2 + y_intercept; 
    x_r = [x1,x2] ;y_r = [y1,y2];

    azimuth_1 = data_table(:,"azimuth_start_at_kiso");
    azimuth_1 = table2array(azimuth_1);
    azimuth_2 = data_table(:,"azimuth_end_at_kiso");
    azimuth_2 = table2array(azimuth_2);
    elevation_1 = data_table(:,"elevation_start_at_kiso");
    elevation_1 = table2array(elevation_1);
    elevation_2 = data_table(:,"elevation_end_at_kiso");
    elevation_2 = table2array(elevation_2);
    duration = data_table(:,"duration");
    duration = table2array(duration);
    pixel_width = 20;
    loop_number = size(ID);
end

recheck_ID = [];angle_original = [];pixel_width_original = [];
flux_meteor_original = [];line_intensity_original = [];
pixel_width_modified = [];flux_meteor_modified = [];
line_intensity_modified = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section : angular velocity [deg/sec]

flight_distance = sqrt(abs(azimuth_1 - azimuth_2).^2 + ...
    abs(elevation_1 - elevation_2).^2);
angular_velocity = flight_distance./duration;
sensor_fov = 1.19/3600; %tomoeの1pixelの視野(1.19"×1.19")
through_pixel = angular_velocity./sensor_fov;
through_pixel = 0.5*through_pixel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section : ループ処理（回転、rangeを求めて加算、ガウシアン）
% imageの範囲を出る場合はずらす

a = pixel_width*cos(deg2rad(90- slope_angle));
b = pixel_width*sin(deg2rad(90- slope_angle));
w1y1 = y1 - b; w1x1 = x1 + a;
w2y1 = y1 + b; w2x1 = x1 - a;
w1y2 = y2 - b; w1x2 = x2 + a;
w2y2 = y2 + b; w2x2 = x2 - a;
for i = 1:loop_number
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
LoopNum = loop_number(1,1);t = 0; %waitbar
flux = zeros(loop_number); 
line_intensity = zeros(loop_number);
h = waitbar(0,'running');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section : midpoint
for i = 1:loop_number %%%途中からループはここをいじる
    t = t + 1; 
    h = waitbar(i/LoopNum,h,['Number of Loop = ',num2str(t)]); %waitbar
    if option_single_image == 1
        fitsfile_name = ID(1,1);
    else    
        fitsfile_name = ID(i,1);
    end
    data = fitsread(fitsfile_name);
    data = data(:,:,1); % 差分後の画像
    rotate_data = imrotate(data, -1*angle(i,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    position_x = center_x(i,1) - size(data,2)/2;
    position_y = center_y(i,1) - size(data,1)/2;
    rotate_x = position_x*cos(angle_rad(i,1)) - position_y*sin(angle_rad(i,1));
    rotate_y = position_x*sin(angle_rad(i,1)) + position_y*cos(angle_rad(i,1));
    position_rotate_x = size(rotate_data,2)/2 + rotate_x;
    position_rotate_y = size(rotate_data,1)/2 + rotate_y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section : 両端の点
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
% section : 矩形の頂点
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
% section : フラックス求める
    if w1x1(i,1) < 0 || w2x1(i,1) < 0 || w1x2(i,1) > 500 || w2x2(i,1) > 500 ||...
            w1y1(i,1) < 0 || w2y1(i,1) > 282 || w1y2(i,1) < 0 || w2y2(i,1) > 282 ...
            || w1x1(i,1) == w1x2(i,1) || w1y1(i,1) == w1y2(i,1) ||...
            fitsfile_name == "d_fits/sub661714_118.fits" ||...
            fitsfile_name == "d_fits/sub705614_052.fits"
        line_intensity(i,1) = 0;
        flux(i,1) = 0;
        disp("This frame is out of range.");
        
        if option_single_image == 1 && re_analyze == 1
            %{
            ここいらない
            pixel_width = array2table(pixel_width);
            line_intensity = array2table(line_intensity(i,1));
            flux = array2table(flux(i,1));
            recheck_list(index,"pixelwidth_modified") = pixel_width;
            recheck_list(index,"line_intensity_modified") = line_intensity;
            recheck_list(index,"flux_modified") = flux;
            %}
        %{
        elseif option_single_image == 0 && re_analyze == 1;
            %いらない
            %line_intensity = array2table(line_intensity(i,1));
            %flux = array2table(flux(i,1));
            %recheck_list(index,"pixelwidth") = pixel_width;
            %recheck_list(index,"line_intensity_modified") = line_intensity;
            %recheck_list(index,"flux_modified") = flux;
        %} 
        elseif option_single_image == 0 && re_analyze == 0 %recheck_list
            recheck_ID = vertcat(recheck_ID,ID_list(i,1));
            angle_original = vertcat(angle_original,angle(i,1));
            pixel_width_original = vertcat(pixel_width_original,pixel_width);
            line_intensity_original = vertcat(line_intensity_original,line_intensity(i,1));
            flux_meteor_original = vertcat(flux_meteor_original,flux(i,1));
            pixel_width_modified = vertcat(pixel_width_modified,0);
            line_intensity_modified = vertcat(line_intensity_modified,0);
            flux_meteor_modified = vertcat(flux_meteor_modified,0);
        
        end
        
    else
        data3 = rotate_data(range_rotate_w1y1:range_rotate_w1y2 ...
            ,range_rotate_w2x1:range_rotate_w1x1);
        y = mean(data3);
        x = range_rotate_w2x1:range_rotate_w1x1;  
        f = fit(x.',y.','gauss1');
        coeff = coeffvalues(f);
        fun = @(x) coeff(1)*exp(-((x-coeff(2))/coeff(3)).^2);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section :  photometry
%2パターン(関数の積分,普通にフラックスで足し合わせる)
        fwhm = coeff(3)/sqrt(2);
        aperture_radius = 1.5*fwhm;
        aperture_x1 = position_rotate_x - aperture_radius;
        aperture_x2 = position_rotate_x + aperture_radius;
        integer = integral(fun,aperture_x1,aperture_x2);
        line_intensity(i,1) = integer/aperture_radius;%[ADU/pix]
        flux(i,1) = line_intensity(i,1)*through_pixel(i,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section : photometryイメージ
        if option_show_image == 1
            figure('visible', 'on');
            subplot(2,2,1)
            option_show_image = 1;
            plot(x,y); %ガウシアン画像
            title(fitsfile_name);
            hold on;
            plot(f,x,y); 
            hold off;
            grid on
            %if option_save_image == 1
                %saveas(gcf,photometry_filename(i,1));
            %end
            subplot(2,2,2)
            imshow(rotate_data); %回転後の画像
            hold on
            title(fitsfile_name);
            scatter(range_rotate_x1,range_rotate_y1,5,"red","filled")
            scatter(range_rotate_x2,range_rotate_y2,5,"red","filled")
            line(rotl1,rotl2,"Color","green","Linewidth",0.5,"Linestyle","-.")
            line(rotl3,rotl4,"Color","green","Linewidth",0.5,"Linestyle","-.")
            line(rotl5,rotl6,"Color","green","Linewidth",0.5,"Linestyle","-.")
            line(rotl7,rotl8,"Color","green","Linewidth",0.5,"Linestyle","-.")
            hold off
            axis equal
            
            subplot(2,2,4),imshow(data);
            title(fitsfile_name); %回転前の画像
            hold on
            scatter(x1(i,1),y1(i,1),5,"red","filled")
            scatter(x2(i,1),y2(i,1),5,"red","filled")
            line(l1(i,:),l2(i,:),"Color","green","Linewidth",0.5,"Linestyle","-.")
            line(l3(i,:),l4(i,:),"Color","green","Linewidth",0.5,"Linestyle","-.")
            line(l5(i,:),l6(i,:),"Color","green","Linewidth",0.5,"Linestyle","-.")
            line(l7(i,:),l8(i,:),"Color","green","Linewidth",0.5,"Linestyle","-.")
            hold off
            axis equal
            if option_save_image == 1
                saveas(gcf,image_filename(i,1));
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %画像一枚一枚別にする場合はここ編集
        %{
        if option_show_image == 1
            subplot(2,2,3),imshow(rotate_data); %回転後の画像
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
            if option_save_image == 1
                saveas(gcf,rotimage_filename(i,1));
            end
            clf;
        end
        %}
        %{
        if option_show_image == 1
            subplot(2,2,4),imshow(data);
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
            if option_save_image == 1
                saveas(gcf,image_filename(i,1));
            end
        end
        %}
    end
end
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section : 元のデータテーブルに加える
if option_single_image == 1 && re_analyze == 1 %単体&二回目
    
    line_intensity_modified = array2table(line_intensity);
    flux_meteor_modified = array2table(flux);
    pixel_width_modified = array2table(pixel_width);
    data_table_overwrite(index,"line_intensity") = line_intensity_modified;
    data_table_overwrite(index,"flux_meteor") = flux_meteor_modified;
    
    if version == "test"
        writetable(data_table_overwrite,"test8.csv");
    elseif version == "main"
        writetable(data_table_overwrite,"simultaneous_event_candidates_table_MU+Tomo-e8.csv");
    end
    %recheck_listに入っている or 入ってない
    if isempty(index_recheck) %入ってない
        data_table_past = readtable("test7.csv");
        ID_number = array2table(ID_number,'VariableNames',{'ID'});
        angle_original = array2table(angle,'VariableNames',{'angle_original'});
        line_intensity_original = data_table_past(index,"line_intensity");
        line_intensity_original.Properties.VariableNames{'line_intensity'}...
            = 'line_intensity_original';
        flux_meteor_original = data_table_past(index,"flux_meteor");
        flux_meteor_original.Properties.VariableNames{'flux_meteor'}...
            = 'flux_meteor_original';
        pixel_width_original =array2table(20);
        pixel_width_original.Properties.VariableNames{'Var1'} = 'pixel_width_original';
        
        pixel_width_modified = array2table(pixel_width,'VariableNames',...
            {'pixel_width_modified'});
        line_intensity_modified = array2table(line_intensity(i,1),...
            'VariableNames',{'line_intensity_modified'});
        flux_meteor_modified = array2table(flux(i,1),'VariableNames',...
            {'flux_meteor_modified'});
        T = horzcat(ID_number,angle_original,line_intensity_original,...
            flux_meteor_original,pixel_width_original,pixel_width_modified,...
            line_intensity_modified,flux_meteor_modified);
        recheck_list = vertcat(recheck_list,T);
        if version == "test"
            writetable(recheck_list,"test_recheck_photometry.csv");
        elseif version == "main"
            writetable(recheck_list,"recheck_photometry.csv");
        end
    else %入ってる
        recheck_list(index_recheck,"pixel_width_modified") = pixel_width_modified;
        recheck_list(index_recheck,"line_intensity_modified") = line_intensity_modified;
        recheck_list(index_recheck,"flux_meteor_modified") = flux_meteor_modified;
        if version == "test"
            writetable(recheck_list,"test_recheck_photometry.csv");
        elseif version == "main"
            writetable(recheck_list,"recheck_photometry.csv");
        end
    end
%{
elseif option_single_image == 0 && re_analyze == 1 %2回目以降&全部
    photometry_data = horzcat(line_intensity,flux);
    photometry_data = array2table(photometry_data,'VariableNames',...
       {'line_intensity_modified','flux_meteor_modified'});
    recheck_table = horzcat(recheck_table,photometry_data)
    writetable(recheck_table,"recheck_photometry.csv");
%}    
elseif option_single_image == 0 && re_analyze == 0 %初回
    photometry_data = horzcat(line_intensity,flux);
    photometry_data = array2table(photometry_data,'VariableNames',...
       {'line_intensity','flux_meteor'});
    data_table = horzcat(data_table,photometry_data);
    if version == "test"
        writetable(data_table,"test7.csv");
        !cp test7.csv test8.csv
    elseif version == "main"
        writetable(data_table,"simultaneous_event_candidates_table_MU+Tomo-e7.csv");
        !cp simultaneous_event_candidates_table_MU+Tomo-e7.csv simultaneous_event_candidates_table_MU+Tomo-e8.csv
    end
    angle_original = array2table(angle_original);
    pixel_width_original = array2table(pixel_width_original);
    line_intensity_original = array2table(line_intensity_original);
    flux_meteor_original = array2table(flux_meteor_original);
    
    pixel_width_modified = array2table(pixel_width_modified);
    line_intensity_modified = array2table(line_intensity_modified);
    flux_meteor_modified = array2table(flux_meteor_modified);
    
    recheck_list = horzcat(recheck_ID,angle_original,line_intensity_original,...
        flux_meteor_original,pixel_width_original,pixel_width_modified,...
        line_intensity_modified,flux_meteor_modified);
    recheck_list.Properties.VariableNames{'Var1'} = 'ID';
    if version == "test"
        writetable(recheck_list,"test_recheck_photometry.csv");
    elseif version == "main"
        writetable(recheck_list,"recheck_photometry.csv");
    end
end
close;
disp("Complete!")
