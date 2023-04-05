% analyze_svf.m 
%
% This MATLAB code is designed to perform the Soot Volume Fraction (SVF)
% measurements for droplet combustion experiments. The program is designed
% to automatically find SVF for multiple axes within a given region. Public
% functions within the image process toolbox: ?imfindcircles.m? and
% ?improfile.m?, as well as the financial toolbox: ?tsmovavg.m? are used in
% the program. 
%
% The original image and attenuated image will be converted into grayscale
% intensity matrix with 8-bit unsigned integer. The user inputs a direction
% of analysis angle by doing two mouse input, start and end accordingly.
% The program separates the input range equally with spokes and obtain
% the intensity of both original image and attenuated image. The
% intensity profiles are sent to ?calAbel.m? to calculate the SVF using
% three-point Abel inversion. The results are stored into a matrix and
% could be exported for further analysis. 
%
% The program output the SVF profile image, spokes image, and a .CSV 
% document with SVF value through each spokes.
%
%
%                       2023 Cornell University
% (c) Yuhao Xu, Yiren Shen, C. Thomas Avedisian, Michael C. Hicks, Mun Y. Choi
% If using analyze_svf in your research, please cite the following paper: 
% Y. Xu, Y. Shen, C.T. Avedisian, M.C. Hicks, M.Y. Choi, Quantitative Investigation of
% Sooting Dynamics in Droplet Combustion Using Automated Image Analysis Algorithms, Fuel, 2023
%
% For user support and report of bugs, contact ys672 at cornell.edu 
%*******************************************************************
 
close all
clear all

%*************************************************************************
%****************************INPUT SECTION********************************
% Input image names, ex: E066C01A_00149, E066C01A_00153 to be analyzed,
% enter 'E066C01A_00' for foldername and [149 153] for imgfiles. The
% original image (t0) of the droplet must be the first file in the list. 
%
% Please keep this program in the same directory with the images. 
foldername = 'E073C09A_00';
imgfiles = [134 195];
% the number of image in the file list to be analyzed
n=2;
% estimate the size of droplet in pixel
dropletRadMin = 10;     %usually pick the size of true droplet -20 pixel
dropletRadMax = 50;     %usually pick the size of true droplet +20 pixel
% Abel inversion attenuation calibration factor
%  based on wavelength, resolution and soot properity
calibration=30.265/1024/1000; 
% enter the number of single-sided ROI (ie, number of rays on one side)
SpoksNumber = 4;


%*************************************************************************
%*********************************Code************************************

imgnum = size(imgfiles);
filename = strings(imgnum) + foldername;
 
for i = 1:imgnum(1,2)
    
    filename(i) = filename(i) + imgfiles(i) + '.tif';
 
end
 

% read the original image and the soot image
original=uint8((imread(char(filename(1)))./256));
sootimg=uint8((imread(char(filename(n)))./256));
 
% identification of the center C and radius R of droplet in original image
figure 
imshow(original);
title('Reference Image');
[CenterOrig,rOrig]=imfindcircles(original,[dropletRadMin,dropletRadMax],'ObjectPolarity','dark');
viscircles(CenterOrig,rOrig);
CenterOrig=round(CenterOrig);
fprintf('Original image center (%d, %d), radius %d \n', CenterOrig(1), CenterOrig(2), rOrig)
fprintf('Please check the validity of droplet boundary marked in red \n')
pause(3) % give user time to check identified droplet boundary
 
% identification of the center C and radius R of droplet in the soot image
figure
imshow(sootimg);
title('Attenuated Image');
[CenterSoot,rS]=imfindcircles(sootimg,[dropletRadMin,dropletRadMax],'ObjectPolarity','dark');
viscircles(CenterSoot,rS);
CenterSoot=round(CenterSoot);
fprintf('Sooty image center (%d, %d), radius %d \n', CenterSoot(1), CenterSoot(2), rS)
fprintf('Please check the validity of droplet boundary marked in red \n')
pause(3) % give user time to check identified droplet boundary
 
i2=length(CenterSoot);
 

 
%% single spoke analysis
if SpoksNumber==1
 
    close all
    
    % obtain the spoke of interest in the soot image by inputting one point
    % through mouse click. A(x2,y2)
    fprintf('************ ROI SELECTION ************************ \n')
    fprintf('Please use mouse click to input the endpoint of ROI \n')
    figure
    hold on
    imshow(sootimg)
    [x2,y2]=ginput(1);
    x2=round(x2);
    y2=round(y2);
    linearrayx2=[CenterSoot(1), x2];
    linearrayy2=[CenterSoot(2), y2];
    line(linearrayx2,linearrayy2,'Color','r','LineWidth',1)
    hold off
    
    dx=x2-CenterSoot(1);
    dy=y2-CenterSoot(2);
    x1=CenterOrig(1)+dx;
    y1=CenterOrig(2)+dy;
    linearrayx1=[CenterOrig(1), x1];
    linearrayy1=[CenterOrig(2), y1];
 
    centerdist=sqrt((CenterSoot(1)-CenterOrig(1))^2+(CenterSoot(2)-CenterOrig(2))^2);
 
    if centerdist >= (rOrig/5)
        originalprofile=improfile(original,linearrayx2,linearrayy2);
        figure
        hold on
        imshow(original)
        line(linearrayx2,linearrayy2,'Color','r','LineWidth',1)
        hold off
    else
       originalprofile=improfile(original,linearrayx1,linearrayy1);
       figure
       hold on
       imshow(original)
       line(linearrayx1,linearrayy1,'Color','r','LineWidth',1)
       hold off
    end
 
    % calculate the SVF profile through the spoke.
    % using function ?improfile? to obtain the intensity of the image.
    %
    % the calculation is based on Lambert Law, which offers the
    % relationship between attenuation of light and the transmittance of
    % material in the path. 
    % 
    % the intensity of original image and soot image is compared to get the
    % transmittance (named as 'sootprofile'):
    %                   T = I/I0 = exp[-u*L]
    % the field distribution of Soot Volume Fraction (P_fv) is described as:
    %                   P_fv = -(lamda/Ke)*ln(I/I0)
    % with three-point Abel inversion (funciton ?calAbel?), final soot
    % volume fraction (sootintensity) is obtained. 
    sootimgprofile=improfile(sootimg,linearrayx2,linearrayy2);
    sootimgprofile=tsmovavg(sootimgprofile,'s',4,1);
    originalprofile=tsmovavg(originalprofile,'s',4,1);
    
    sootprofile=sootimgprofile./originalprofile;    % I/I0
    P_fv=-(653*10^(-9)/8.6)*log(sootprofile);       % P_fv 
    P_fv(1:5)=0;
    sootintensity=calAbel(P_fv,calibration);    % Abel inversion
   
    dist=0:1:(length(sootintensity)-1);
    dist=dist*calibration;
 
    rSmm=rS*calibration;
    rOmm=rOrig*calibration;
    dist1=dist/rOmm;
    dist2=dist/rSmm;
 
    tfdist=(dist1>=1);
    sootintensity=smooth(smooth(sootintensity));
    sootintensity=sootintensity.*double(tfdist)';
 
    % output file
    savingtable= {7};
    savingtable{1,1}={'Original Intensity'};
    savingtable{1,2}='SootImg Intensity';
    savingtable{1,3}='Soot Volume Fraction';
    savingtable{1,4}='Original droplet Center (x,y)';
    savingtable{1,5}='SootImg droplet center (x,y)';
    savingtable{1,6}='dx (pix)';
    savingtable{1,7}='dy (pix)';
    savingtable{2,1}=originalprofile;
    savingtable{2,2}=sootimgprofile;
    savingtable{2,3}=sootintensity;
    savingtable{2,4}=CenterOrig';
    savingtable{2,5}=CenterSoot';
    savingtable{2,6}=dx;
    savingtable{2,7}=dy;
 
    dist2=dist2';
    date=char(datetime('now','TimeZone','local','Format','d-MMM-y-HH_mm_ss'));
    filename=char([date 'SVFOutputFile.txt']);
    save(filename)
 
    figure
    plot(dist2,sootintensity)
    xlabel('r/r_i')
    ylabel('f_v [ppm]')
    title('Soot Volume Fraction on given line segment')
 
%% multiple spoke analysis
else  
    close all
    figure
    hold on
 
    % obtain the range of interest in the soot image by inputting two points
    % through mouse. A(x2,y2), B(X3,y3)
    % the range must be inputted in anti-clock direction. 
    fprintf('************ ROI SELECTION ************************ \n')
    fprintf('Please use two mouse clicks to input the endpoint of ROI \n')
    
    imshow(sootimg)
    [x2,y2]=ginput(1);
    [x3,y3]=ginput(1);
    x2=round(x2); x3=round(x3);
    y2=round(y2); y3=round(y3);
 
    % show the range of interest
    linearrayx2=[CenterSoot(1), x2];
    linearrayy2=[CenterSoot(2), y2];
    line(linearrayx2,linearrayy2,'Color','r','LineWidth',1)
    
    linearrayx3=[CenterSoot(1), x3];
    linearrayy3=[CenterSoot(2), y3];
    line(linearrayx3,linearrayy3,'Color','b','LineWidth',1)
    hold off
    pause
    fprintf('Paused, click any key to continue\n')
    
    % separate the range of interest equally with spokes. 
    % the crossing angle in the range of interest is calculated through Law
    % of cosines. 
    distancebetorigin=sqrt((CenterSoot(1)-CenterOrig(1))^2+(CenterSoot(2)-CenterOrig(2))^2);
    
    if distancebetorigin >= (0.2*rOrig)
        CenterAnalysis=CenterSoot;
    else
        CenterAnalysis=CenterOrig;
    end
    
    Radius=((x2-CenterAnalysis(1))^2+(y2-CenterAnalysis(2))^2)^0.5;
    R2=Radius;                            % length of OA
    R3=((x3-CenterAnalysis(1))^2+(y3-CenterAnalysis(2))^2)^0.5; 
                                          % length of OB
    d=((x2-x3)^2+(y2-y3)^2)^0.5;          % length of AB
    CosAngValue=(R2^2+R3^2-d^2)/(2*R2*R3);
    CrossAngle=acos(CosAngValue);
    DeltAngle=CrossAngle/(SpoksNumber-1); % angle between each two spokes
    AngleIni=atan((CenterAnalysis(2)-y2)/(x2-CenterAnalysis(1)));
                                          % initial angle between OA and x axis
   
    % the profile of SVF vs. radius
        figure
        hold on
        xlabel('r/r_i')
        ylabel('f_v [ppm]')
        title('Soot Volume Fraction on given line segment')
        STORE={};
        
        currAngle=AngleIni;
        currAnalysis1=CenterAnalysis(1);
        currAnalysis2=CenterAnalysis(2);
        
    % calculate the start point and end point of each spokes. each spoke
    % has a another one located in the opposite of the range of interest. 
    for alpha=1:(2*SpoksNumber)
        if alpha <= SpoksNumber
            EndpointsX(alpha)=sqrt(Radius^2/((tan(currAngle))^2+1))+CenterAnalysis(1);
            EndpointsY(alpha)=CenterAnalysis(2)-sqrt(Radius^2-(EndpointsX(alpha)...
                -CenterAnalysis(1))^2);
            currAngle=currAngle+DeltAngle;
            
        else if alpha > SpoksNumber
                EndpointsX(alpha)= CenterAnalysis(1)-abs(CenterAnalysis(1)-EndpointsX(alpha-SpoksNumber));
                EndpointsY(alpha)= CenterAnalysis(2)+abs(CenterAnalysis(2)-EndpointsY(alpha-SpoksNumber));
            end
        end
    end
  
    % calculate the SVF profile through each spoke.
    % using function ?improfile? to obtain the intensity in image.
    %
    % the calculation is based on Lambert Law, which offers the
    % relationship between attenuation of light and the transmittance of
    % material in the path. 
    % 
    % the intensity of original image and soot image is compared to get the
    % transmittance (named as 'sootprofile'):
    %                   T = I/I0 = exp[-u*L]
    % the integration of the field distribution of Soot Volume Fraction (P_fv) 
    % is described as:
    %                   P_fv = -(lamda/Ke)*ln(I/I0)
    % with three-point Abel inversion (function ?calAbel?), final soot
    % volume fraction (sootintensity) is obtained. 
    for alpha=1:(2*SpoksNumber)
        
            linearrayX=[CenterAnalysis(1),EndpointsX(alpha)];
            linearrayY=[CenterAnalysis(2),EndpointsY(alpha)];
            samplerate=round(Radius);
            
            originalprofile=improfile(original,linearrayX,linearrayY,samplerate);
            sootimgprofile=improfile(sootimg,linearrayX,linearrayY,samplerate);
            sootimgprofile=smooth(tsmovavg(sootimgprofile,'s',4,1));
            originalprofile=smooth(tsmovavg(originalprofile,'s',4,1));
            sootprofile=sootimgprofile./originalprofile;    % I/I0
            P_fv=-(653*10^(-9)/8.6)*log(sootprofile);       % P_fv 
            P_fv(1:5)=0;
            sootintensity=calAbel(P_fv,calibration);         % Abel inversion
 
            dist=0:1:(length(sootintensity)-1);
            dist1=dist/rOrig;
            dist2=dist/rS;                                  % r/rs
            tfdist=(dist1>=1.2);
            sootintensity=sootintensity.*double(tfdist)';
            STORE{alpha}=sootintensity;                     % SVF results
            plot(dist2,sootintensity)                       % SVF profile
      
    end
 
    hold off
    
    figure
    hold on
    title('ROI')
    line(linearrayx2,linearrayy2,'Color','r','LineWidth',1)
    line(linearrayx3,linearrayy3,'Color','b','LineWidth',1)
  
    image(sootimg)
    for i=1:alpha
        plot(EndpointsX(i),EndpointsY(i),'k*','MarkerSize',10)
    end
 
    hold off
    
  
end
 
% data saving
% this MATLAB code saves the image of spokes location, the profile of
% different spokes, and a .CSV document with SVF value through each spoke. 
close all 
 
dist2=dist2';   
 
% draw the endpoint of each spoke in the soot image and save the image. 
figure
imshow(sootimg)
hold on
title('ROI') 
for i=1:(2*SpoksNumber)
    if i==1
        plot(EndpointsX(i),EndpointsY(i),'r.','MarkerSize',8)
        text(EndpointsX(i)+7,EndpointsY(i)-7,num2str(i),'HorizontalAlignment','center','FontSize',6)
    else if i==2*SpoksNumber
        plot(EndpointsX(i),EndpointsY(i),'b.','MarkerSize',8)
        text(EndpointsX(i)+7,EndpointsY(i)-7,num2str(i),'HorizontalAlignment','center','FontSize',6)
        else    
        plot(EndpointsX(i),EndpointsY(i),'k.','MarkerSize',8)
        text(EndpointsX(i)+7,EndpointsY(i)-7,num2str(i),'HorizontalAlignment','center','FontSize',6)
        end
    end
    
end
 
savefig(strcat(num2str(imgfiles(n)),'(spokes)'))    % save spokes image
SVFdata = cell2mat(STORE);
[rows,cols] = size(SVFdata);
 
for b = 1:(2*SpoksNumber)
   
    for a = 1:rows
        % remove negative svf (unphysical data pt)
        if SVFdata(a,b) < 0
            SVFdata(a,b) = 0;
            
        end
        
    end
    
end
 
% plot the profile of each spoke, and save the profile image. 
figure
hold on
 
 for b = 1:(2*SpoksNumber)
 
     plot(dist2,SVFdata(:,b));
     
 end
 
 legend('show')
 
 xlabel('r/r_i')
 ylabel('SVF [ppm]')
 title('SFV Profiles Along Spokes')
 
 savefig(strcat(num2str(imgfiles(n)),'(profiles)')) % save profile image
 
 
 % saving the SVF value into .CSV document. The program will ask the user
 % to choose the profiles they want to average, which normally are the data
 % of spokes go through perfect sphere range. 
 % 
 % a additional .CSV document named as 'image#(averaged profiles)' will be
 % created to record the number of spokes that are averaged.
 fprintf('*************** SELECT PROFILE AND EXPORT*************** \n');
 fprintf('Input the index of profiles according to ROI image index, \n');
 fprintf('use vectorized format, eg [1,2,4] representing selecting index 1,2, and 4.  \n');
 fprintf('All ROIs will be exported to $test_number$.csv  \n');
 selection = input('Choose index of profiles to average, in vector format \n');
 
 [rows, cols] = size(selection);
 
 clear SVFmax; clear locSVFmax;
 
 for b = 1:cols
     
     [SVFmax(b), indx] = max(SVFdata(:,selection(1,b)));
     
     locSVFmax(b) = dist2(indx);
     
 end
 
 meanSVFmax(n-1) = mean(SVFmax);
 meanSVFloc(n-1) = mean(locSVFmax);
 scaledTime(n-1) = (n-1)*0.05;
 
 clear profiles; clear aveProfile;
 
 for b = 1:cols
     
     profiles(:,b) = SVFdata(:,selection(1,b));
     
 end
 
aveProfile = mean(profiles,2);
exporttdata = [dist2 SVFdata aveProfile];
 
% write the SVF data into image#.csv
dlmwrite(strcat(num2str(imgfiles(n)),'.csv'),exporttdata,'precision',8);
 
masterDATA = [meanSVFmax' meanSVFloc' scaledTime'];
dlmwrite('Mean Max SVF Data.csv',masterDATA,'precision',8);
 
%write the number of spokes averaged into image#(averaged profiles).csv
dlmwrite(strcat(num2str(imgfiles(n)),'(averaged profiles).csv'),selection);
 
 
%end

%******************************************************************
%*************************SUB FUNCTIONS****************************
%
%


function SVF=calAbel(P,calibration)
    % The MATLAB code of ?calAbel? is provided below.
    % calAbel.m
    %
    % The code is designed to obtain the Soot Volume Fraction(SVF) through
    % three-point Abel inversion. The mathematical technique is recommended by
    % Dasch (1992) for one-dimensional tomography. 
    % The field distribution fv(ri) can be found at a finite domain with
    % operator coefficients in the inversion. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main: Start
 
    [row,col] = size(P);
    D = zeros(col,col);
 
    for i = 1:row
        for j = 1:row       % D operator index start from 0
            D(i,j) = OP_D(i-1,j-1);
        end
    end
 
    inv_Mat = zeros(size(P));
     F = [];
     
    for i=1:row
        Fj=0;
        for j=1:row
           Fj = Fj+D(i,j)* P(j);
        end
        F(i)=Fj;
    end
 
    F=F';
    F=F*10^6;
    F=F/calibration/2;
    
    SVF=F;
 
    function D = OP_D(i,j)
        
    % Calculate three-point abel inversion operator Di,j.
    % The index i,j start from 0.
    % The formula followed Dasch 1992 (Applied Optics) which contains several typos.
    % One correction is done in function OP1 follow Martin's PhD thesis.
    global PI
    PI = 3.14159265;
 
    if j<i-1
        D = 0;
    elseif j==i-1
        D = OP0(i,j+1)-OP1(i,j+1);
    elseif j==i
        D = OP0(i,j+1)-OP1(i,j+1)+2*OP1(i,j);
    elseif i==0&&j==1
        D = OP0(i,j+1)-OP1(i,j+1)+2*OP1(i,j)-2*OP1(i,j-1);
    elseif j>=i+1
        D = OP0(i,j+1)-OP1(i,j+1)+2*OP1(i,j)-OP0(i,j-1)-OP1(i,j-1);
    end
 
    function I0 = OP0(i,j)        % Define operator OP0
        if j<i || (j==i&&i==0)
            I0 = 0;
        elseif (j==i&&i~=0)
            I0 = log((((2*j+1)^2-4*i^2)^0.5+2*j+1)/(2*j))/(2*PI);
        elseif j>i
            I0 = log((((2*j+1)^2-4*i^2)^0.5+2*j+1)/(((2*j-1)^2-4*i^2)^0.5+2*j-1))/(2*PI);
        end
    end
 
    function I1 = OP1(i,j)       % Define operator OP1
        if j<i
            I1 = 0;
        elseif j==i
            I1 = ((2*j+1)^2-4*i^2)^0.5/(2*PI)-2*j*OP0(i,j);
        elseif j>i
            I1 = (((2*j+1)^2-4*i^2)^0.5-((2*j-1)^2-4*i^2)^0.5)/(2*PI)-2*j*OP0(i,j); 
        end
 
    end
 
    end
end

