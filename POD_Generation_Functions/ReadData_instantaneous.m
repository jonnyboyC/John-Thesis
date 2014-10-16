function[xi,yi,ui,vi,Nx,Ny] = ReadData_instantaneous(ni)
%%%%%% Airfoil velocity plotter %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% For phaselocked movies only for this experiment %%%%%

%%%%%% see date.xls for details
% clear all; clc; clf; close all;

% Creating directories for processed data
% warning off

%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = .952;  % Tunnel constant
c = 8; % inches
Lscale = c*.0254;  % meters
patm = 98600; % Pa
phases = 8;
countph = 0; % counter for phases
rect = [-70 -60 515 335]; % rectangular window to create movie includes title, labels and colorbar
%[left bottom width height]
% ni=1930;
%%%%%%%%%%%  Read data from excel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [data, act] = xlsread([dscribe, date], 'DataSheet'); % read data from excel sheet, "act" not used here
% run = data(:,1); %chronological run number
% re = data(:,2); % Re/1000 tested
% ia = data(:,3); % Angle of attack (deg)
% q = 249*data(:,4)/k;  % Dynamic Pressure (Pa)
% dp = 249*data(:,5);  % Static Pressure (Pa)
% temp = (data(:,6) + 459.67)/1.8; % Temperature (K)
% % data(:,7) is tunnel setting
% % data(:,8) is actuator seen as text ("act")
% xx = 1000*data(:,9); % Actuator Location
% % data(:,10) is DC input voltage to pulser
% vp = data(:,11); % pulser voltage
% fp = data(:,12); % pulser frequency
% ph = data(:,15);
% clear data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear
fprintf(1, 'Please choose a directory for the data files\n');
d = uigetdir('', 'Choose Source Image Directory');
files1 = [ dir([d '\*.vc7'])];% dir([d '\*.'])];%
nf=length(files1);
if nf==0
    files1 = [ dir([d '\*.im7'])];% dir([d '\*.'])];%
    nf=length(files1);
end
if ni < nf
    nf=ni;
end

% tt=['Test_01'; 'Test_02'; 'Test_03'; 'Test_04'; 'Test_05'; 'Test_06';...
%     'Test_07'; 'Test_08'; 'Test_09'; 'Test_10'; 'Test_11'; 'Test_12';...
%     'Test_13'; 'Test_14'; 'Test_15'; 'Test_16'; 'Test_17'; 'Test_18';...
%     'Test_19'; 'Test_20'; 'Test_21'; 'Test_22'; 'Test_23'; 'Test_24';...
%     'Test_25'; 'Test_26'; 'Test_27'; 'Test_28'; 'Test_29'; 'Test_30';...
%     'Test_31'; 'Test_32'; 'Test_33'; 'Test_34'; 'Test_35'; 'Test_36';...
%     'Test_37'; 'Test_38'; 'Test_39'; 'Test_40'; 'Test_41'; 'Test_42';...
%     'Test_46'; 'Test_47'; 'Test_45'; 'Test_46'; 'Test_47'; 'Test_48';...
%     'Test_49'; 'Test_50'; 'Test_51'; 'Test_52'; 'Test_53'];

fn=1;
% for t = 1:53
%     if strfind(d, tt(t,:)) > 0
%         fn=t+1;
%         break
%     end
% end

% Processing data
for j = fn%:4%length(re)
    for i = 1:nf
        if ~files1(i).isdir
            n1      = files1(i).name;
            c(i) = str2num(n1(3:6));
            fprintf(1, 'Processing file %s\n', n1);
            files2 = dir([d '\' n1 '*']);
            
            fn1       = [d '\' n1];
            %         j
            %         length(re)
            
            %%% Preparing airfoil to plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Coordinates of airfoil tail from raw data
            % Fine tune in Matlab to look best
            
%             if ia(j) == 10
%                 xactual = -111; yactual = -46; %mm tail coordinates from raw data (-108,-65)
%                 fudge = .5; % angle of attack adjustment to fit the data
%             elseif ia(j) == 12
%                 xactual = -110; yactual = -49; %mm tail coordinates from raw data (-108,-65)
%                 fudge = .5; % angle of attack adjustment to fit the data
%             elseif ia(j) == 14
%                 xactual = -112; yactual = -53.5; %mm tail coordinates from raw data (-108,-65)
%                 fudge = 1; % angle of attack adjustment to fit the data
%             elseif ia(j) == 16
%                 xactual = -110; yactual = -56; %mm tail coordinates from raw data (-108,-65)
%                 fudge = 1.5; % angle of attack adjustment to fit the data
%             elseif ia(j) == 18
%                 xactual = -111; yactual = -60; %mm tail coordinates from raw data (-108,-65)
%                 fudge = 2; % angle of attack adjustment to fit the data
%             elseif ia(j) == 19
%                 xactual = -112; yactual = -62; %mm tail coordinates from raw data (-108,-65)
%                 fudge = 2.5; % angle of attack adjustment to fit the data
%             elseif ia(j) == 20
%                 xactual = -112; yactual = -64; %mm tail coordinates from raw data (-108,-65)
%                 fudge = 2.5; % angle of attack adjustment to fit the data
%             end
%             
%             %%% Loading and finding airfoil coordinates
%             load(['NACA0015.mat']); % c:\airfoil\load NACA 0015 coordinates
%             coordinates = coordinates/100; % removing percentagte
%             
%             % Rotating airfoil
%             coordinates(:,1) = coordinates(:,1)*cos(-(ia(j)+fudge)*pi/180) - coordinates(:,2)*sin(-(ia(j)+fudge)*pi/180);
%             coordinates(:,2) = coordinates(:,1)*sin(-(ia(j)+fudge)*pi/180) + coordinates(:,2)*cos(-(ia(j)+fudge)*pi/180);
%             
%             % Finding the airfoil tail after rotation
%             r = sqrt(coordinates(:,1).^2 + coordinates(:,2).^2);
%             rmaxi = find(r == max(r));
%             xtip = coordinates(rmaxi,1);
%             ytip = coordinates(rmaxi,2);
%             
%             % Correction factor
%             xcorrect = xtip + xactual/(Lscale*1000);
%             ycorrect = ytip - yactual/(Lscale*1000);
%             
            % Run nomenclature and calcs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             runcode = ['r',num2str(re(j)),'i',num2str(ia(j)),'a',act{j+1,8},'x',num2str(xx(j)),'v',num2str(vp(j)*10),'f',num2str(fp(j)),'p',num2str(ph(j))]; % parameter description for each test;
%             uset = sqrt(498*287*temp(j)*(k*q(j)/249)/(k*(patm-dp(j)))); % Normalization velocity
%             rho = ((patm - dp(j)))/(287*temp(j));
            
            % Loading PIV data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Loop to read and reshape PIV data
%             if j == 1
                lavdata = readimx([d '\' n1]);
%             elseif j < 11
%                 lavdata = readimx(['PIVdata\RawData\Test_0',num2str(j-1),'\PostProc\TimeMeanQF_Vector\B00001_Avg V.VC7']);
%             else
%                 lavdata = readimx(['PIVdata\RawData\Test_',num2str(j-1),'\PostProc\TimeMeanQF_Vector\B00001_Avg V.VC7']);
%             end
            
            [x,y,u,v] = showimx(lavdata); % spatial x,y in mm
            Nx=lavdata.Nx; Ny=lavdata.Ny;
            clear lavdata
            
%             x = -1*x/Lscale/1000 + xcorrect; %x/c
%             y = y/Lscale/1000 + ycorrect; %y/c
%             u = -1*u/uset; % dimens
%             v = v/uset;
            
            if (x(1,1)>x(end,1) && y(1,1)>y(1,end))
                for ix =1:size(x,1)
                    for iy = 1:size(y,2)
                        xn(ix,iy)=x(end-ix+1,end-iy+1);
                        yn(ix,iy)=y(end-ix+1,end-iy+1);
                        un(ix,iy)=u(end-ix+1,end-iy+1);
                        vn(ix,iy)=v(end-ix+1,end-iy+1);
                    end
                end
            elseif (x(1,1)>x(end,1) && y(1,1)<y(1,end))
                for ix =1:size(x,1)
                    xn(ix,:)=x(end-ix+1,:);
                    yn(ix,:)=y(end-ix+1,:);
                    un(ix,:)=u(end-ix+1,:);
                    vn(ix,:)=v(end-ix+1,:);
                end
            elseif (x(1,1)<x(end,1) && y(1,1)>y(1,end))
                for iy = 1:size(y,2)
                    xn(:,iy)=x(:,end-iy+1);
                    yn(:,iy)=y(:,end-iy+1);
                    un(:,iy)=u(:,end-iy+1);
                    vn(:,iy)=v(:,end-iy+1);
                end
            end
            x=xn;
            y=yn;
            u=un;
            v=vn;
            
            
            %loop to set masked area to white
            %             for ii = 1:Nx
            %                 for jj = 1:Ny
            %                     if u(ii,jj) == 0
            %                         u(ii,jj) = NaN;
            %                     end
            %                 end
            %             end
            
            ui{1}(:,:,i)=u;
            vi{1}(:,:,i)=v;
            
            %Plot Vorticity
            set(gcf,'color','w');
            pcolor(x,y,u); shading interp;
%             set(gca,'clim',[0 1.4],'TickLength',[0.025 0.025],'TickDir','out')
            colorbar;
            axis equal;
%             axis([-.05 1.05 -.4 .3]);
%             title(['<T>/U_\infty ',runcode,' ', date,'n',num2str(run(j))],'FontSize', 10,'fontweight','normal')
%             xlabel('x/c','FontSize', 14,'fontweight','bold')
%             ylabel('y/c','FontSize', 14,'fontweight','bold')
%             set(gca,'GridLineStyle',':','FontSize', 14,'fontweight','bold');
%             hold on
%             fill(coordinates(:,1),coordinates(:,2),'k');
            grid on;
            %         pause
            
            %             saveas(gcf, ['PIVdata\ProcessedData\velocity\figures\','TmpN_',runcode,'_',date,'n',num2str(run(j))],'fig');
            %             saveas(gcf, ['PIVdata\ProcessedData\velocity\figures\','TmpN_',runcode,'_',date,'n',num2str(run(j))],'png');
            
            % Loop to create movie of phaselock data, only works if phaselocked
            %     %% images are acquired for a given case (actuator, waveform, etc)
            %     %% sequentially. Order of phases in the case does not matter.
            %     %% Saveas command above is necessary to capture image, could also use
            
%             if ph(j) ~= 0 % only for phaselock data
%                 countph = countph + 1; %counter for phases
%                 h = gca; % get current axis
%                 
%                 % loop to select and order phases
%                 if ph(j) == 1
%                     M(1) = getframe(h,rect);
%                 elseif ph(j) == 2
%                     M(2) = getframe(h,rect);
%                 elseif ph(j) == 3
%                     M(3) = getframe(h,rect);
%                 elseif ph(j) == 4;
%                     M(4) = getframe(h,rect);
%                 elseif ph(j) == 5;
%                     M(5) = getframe(h,rect);
%                 elseif ph(j) == 6;
%                     M(6) = getframe(h,rect);
%                 elseif ph(j) == 7;
%                     M(7) = getframe(h,rect);
%                 elseif ph(j) == 8;
%                     M(8) = getframe(h,rect);
%                 end
%                 
%                 if countph == phases % must have 8 phases to create movie
%                     for l=1:countph
%                         gif_add_frame(M(l),['PIVdata\ProcessedData\velocity\movies\TmpN_',runcode,'_', date,'.gif'])
%                     end
%                     countph = 0;
%                 end
%             end
            
        end
%         pause(.2)
        clf
    end
    close
    xi{1}=x;
    yi{1}=y;
end
