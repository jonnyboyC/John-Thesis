% clear all;

%P   = 100  ;  %(N)
% d   = 0.125;  %(sec)
% tau = 1  ;  %(sec)
%T   = tau  ;  %time window
% fs  = 1024  ;
% nb=8;
if i == 1
    h=sc/mu0*dl;
    fs  = 1/h;  %sampling frequency
    vd0 = floor(size(Y1,1)/(Ns));
    vd1 = 0;% floor(size(Y2,1)/(Ns));
    vd2=  0;% floor(size(Y3,1)/(Ns));
    vd3=  0;% floor(size(Y4,1)/(Ns));
    Nf   = Ns  ;  %number of points per time windowvd*
    T  =  Nf/fs;
    % h   = T/Nf  ;  %time resolution (s)
    df  = 1/T  ;  %frequency resolution
    fh  = fs/2 ;  %nyquist frequency

    t = 0:h:T-h;    %discrete time points
    f = 0:df:fs-df; %discrete frequency points

    %create function x(i) for input signal
    % x=1*sin(2*pi*10*t);%+4*cos(2*pi*400*t);


    % noise=0.1*(rand(1,N)-0.5);
    % x=x+noise;                   %add low amplitude noise to signal

    w=hanning(Nf);                 %uniform window
    stt=1;
    edp=Nf;
    Xt=0;
    Xty_dB=0;
    Xt1=0;
    Xty1_dB=0;
    Xt2=0;
    Xty2_dB=0;
    Xt3=0;
    Xty3_dB=0;
    X_ref = .1;  %reference value of 1 m-s
    Xp_ref =20*10^(-6);

    if vd0 > 0
        for ns=1:vd0
            y_w=Y1(stt:edp,:).*(w*ones(1,n));

            %compute the fast fourier transform of the windowed data
            Xff_y =fft(y_w);
            Xf_y =sqrt(Xff_y.*conj(Xff_y))/N;
            Xt=Xt+Xf_y;

            Xfy_dB = 10*log10((Xf_y/Xp_ref));
            Xty_dB = Xty_dB+Xfy_dB;

            stt=stt+Nf;
            edp=edp+Nf;
        end
        Xt=Xt/vd0;

        np=0;
        for l=1:Nf/2
            if f(l)< 800
                np=np+1;
            end
        end
        %     for l=1:n
        %         maxi=max(Xt(np:end,l));
        %         Xt(:,l)=(Xt(:,l)/maxi);%);
        %     end
        Xty_dB=20*log10((Xt/Xp_ref));
    end

    if vd1 > 0
        stt=1;
        edp=Nf;
        for ns=1:vd1
            y_w3=Y2(stt:edp,:).*(w*ones(1,n));

            %compute the fast fourier transform of the windowed data
            Xff_y3 =fft(y_w3);
            Xf_y3 =sqrt(Xff_y3.*conj(Xff_y3))/N;
            Xt3=Xt3+Xf_y3;

            Xfy3_dB = 10*log10((Xf_y3/Xp_ref));
            Xty3_dB = Xty3_dB+Xfy3_dB;

            stt=stt+Nf;
            edp=edp+Nf;
        end
        Xt3=Xt3/vd1;
        np=0;
        for l=1:Nf/2
            if f(l)< 800
                np=np+1;
            end
        end
        %     for l=1:n
        %         maxi=max(Xt3(np:end,l));
        %         Xt3(:,l)=(Xt3(:,l)/maxi);%);
        %     end
        Xty3_dB=20*log10((Xt3/Xp_ref));
    end

    if vd2 > 0
        stt=1;
        edp=Nf;

        for ns=1:vd2
            y_w1=Y3(stt:edp,:).*(w*ones(1,n));

            %compute the fast fourier transform of the windowed data
            Xff_y1 =fft(y_w1);
            Xf_y1 =sqrt(Xff_y1.*conj(Xff_y1))/N;
            Xt1=Xt1+Xf_y1;

            Xfy1_dB = 10*log10((Xf_y1/Xp_ref));
            Xty1_dB = Xty1_dB+Xfy1_dB;
            stt=stt+Nf;
            edp=edp+Nf;
        end
        Xt1=Xt1/vd2;
        np=0;
        for l=1:Nf/2
            if f(l)< 800
                np=np+1;
            end
        end
        %     for l=1:n
        %         maxi=max(Xt1(np:end,l));
        %         Xt1(:,l)=(Xt1(:,l)/maxi);%);
        %     end
        Xty1_dB=20*log10((Xt1/Xp_ref));
    end

    if vd3 > 0
        stt=1;
        edp=Nf;
        for ns=1:vd3
            y_w2=Y4(stt:edp,:).*(w*ones(1,n));

            %compute the fast fourier transform of the windowed data
            Xff_y2 =fft(y_w2);
            Xf_y2 =sqrt(Xff_y2.*conj(Xff_y2))/N;
            Xt2=Xt2+Xf_y2;

            Xfy2_dB = 10*log10((Xf_y2/Xp_ref));
            Xty2_dB = Xty2_dB+Xfy2_dB;

            stt=stt+Nf;
            edp=edp+Nf;
        end
        Xt2=Xt2/vd3;
        np=0;
        for l=1:Nf/2
            if f(l)< 800
                np=np+1;
            end
        end
        %     for l=1:n
        %         maxi=max(Xt2(np:end,l));
        %         Xt2(:,l)=(Xt2(:,l)/maxi);%);
        %     end
        Xty2_dB=20*log10((Xt2/Xp_ref));
    end
end


% % %compute the analytical fourier transform for the selected frequencies
% % for q=1:N
% %     total = 0;
% %     for n = 1:d/h
% %         total = total + x_w(n)*exp(-j*2*pi*f(q)*n*h);
% %     end
% %     Xa_FT(q)=h*abs(total);
% % end
% % Xa_FT_dB =(Xa_FT/X_ref);
%
% %compute the discrete fourier transform of the windowed data
% for k=1:N
%     Xd(k) = abs(h*sum(x_w.*exp(-j*2*pi*(0:N-1)*(k-1)/N)));
% end
% Xd_dB =(Xd/X_ref);

%plot data
% subplot(223)
% hold on
% plot(t,u_w,'-'); grid; set(gca,'GridLineStyle',':','FontSize', 12);
% plot(t,v_w,'-'); grid; set(gca,'GridLineStyle',':','FontSize', 12);
% xlabel('\bf Time (sec)'); ylabel('\bf Magnitude (m/s)');
%
% subplot(223)
% plot(f,Xf_u*mu0); hold on;
% plot(f,Xf_v*mu0); hold on;
% grid;  set(gca,'GridLineStyle',':','FontSize', 12);
% xlabel('\bf Frequency(Hz)'); ylabel('\bf Magnitude, m/s'); legend('u', 'v');
% Title('FFT of the velocity field')
% xlim([100 4000]);

% subplot(224)
if n > 4
    if mn == 1
        if vd0 > 0
            plot(f,Xty_dB(:,i),'g-'); hold on;
        end
    elseif mn == 2
        if vd1 > 0
            plot(f,Xty3_dB(:,i),'b-'); hold on;
        end
    elseif mn == 3
        if vd2 > 0
            plot(f,Xty1_dB(:,i),'r--');hold on
        end
    else
        if vd3 > 0
            plot(f,Xty2_dB(:,i),'g-.');hold on
        end
    end
else
    if mn == 1
        if vd0 > 0
            plot(f,Xty_dB(:,i),'g-'); hold on;
        end
    elseif mn == 2
        if vd1 > 0
            plot(f,Xty3_dB(:,i),'b-'); hold on;
        end
    elseif mn == 3
        if vd2 > 0
            plot(f,Xty1_dB(:,i),'r--');hold on
        end
    else
        if vd3 > 0
            plot(f,Xty2_dB(:,i),'g-.');hold on
        end
    end
end
% plot(f,Xfp_dB,'-+'); hold on;
grid;  set(gca,'GridLineStyle',':','FontSize', 11);
xlabel('Frequency(Hz)'); ylabel('Magnitude');% legend('FFT');
title(['FFT of a^',num2str(i), '(t)'])
xlim([1000 5000]);

% figure(3)
% plot(f,Xa_FT_dB,'-',f,Xf_dB); hold on;
% grid; set(gca,'GridLineStyle',':','FontSize', 12);
% xlabel('\bf Frequency(Hz)'); ylabel('\bf Magnitude'); legend('FT','FFT');
%
%
% figure(4)
% plot(f,Xf_dB,'-',f,Xd_dB,'-'); hold on;
% grid; set(gca,'GridLineStyle',':','FontSize', 12);
% xlabel('\bf Frequency(Hz)'); ylabel('\bf Magnitude, N-s'); legend('FT','FFT');
