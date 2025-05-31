classdef MMUPE
    properties
        dv
        NFFT_1
        NFFT_2
        NFFT_3
        pulse_1
        pulse_2
        pulse_3
        waveguide
        z_grid
        OutputInterval
        initial_pulse
        linear_factor
        chi2_factor
        chi3_factor
        Normalize_factor
        vg
        v_grid
        PSD_pump_dBm
    end

    methods
        % Define Class
        function obj = MMUPE(pulse_1,pulse_2,pulse_3,waveguide,z_grid,OutputInterval)
            obj.pulse_1 = pulse_1;
            obj.pulse_2 = pulse_2;
            obj.pulse_3 = pulse_3;
            obj.PSD_pump_dBm = 10.*log10(abs(obj.pulse_1.a_v).^2);
            obj.waveguide = waveguide;
            obj.z_grid = z_grid;
            obj.OutputInterval = OutputInterval;
            obj.NFFT_1 = pulse_1.NFFT;
            obj.dv = pulse_1.dv;
            obj.v_grid = obj.pulse_1.v_grid;
            if ~isempty(pulse_2)
                obj.NFFT_2 = pulse_2.NFFT;
                obj.v_grid = [obj.v_grid;pulse_2.v_grid];
                if ~isempty(pulse_3)
                    obj.NFFT_3 = pulse_3.NFFT;
                    obj.v_grid = [obj.v_grid;pulse_3.v_grid];
                else
                    obj.NFFT_3 = 3*pulse_1.NFFT-2*mod(pulse_1.NFFT,2);
                end
            else
                obj.NFFT_2 = 0;
                if ~isempty(pulse_3)
                    obj.NFFT_3 = pulse_3.NFFT;
                    obj.v_grid = [obj.v_grid;pulse_3.v_grid];
                else
                    obj.NFFT_3 = pulse_1.NFFT;
                end
            end
        end
        
        function obj = prepare(obj)
            cspeed = 299792458;
            epsilon0 = 8.854187817e-12;
            [neff_1,beta_1,divvg,Aeff_1] = obj.pulse_1.mode.getmodeProfile('v_grid');
            %[neff_1_v0,~,~,Aeff_1_v0] = obj.pulse_1.mode.getmodeProfile('v0');
            obj.vg = 1/divvg;
            linear_factor_1 = 1i * beta_1;
            obj.linear_factor = linear_factor_1;
            if ~isempty(obj.waveguide.deff)
                chi2_factor_1 = 1i.*(2.*pi.*obj.pulse_1.v_grid)./(2.*neff_1)./cspeed .* obj.waveguide.deff./sqrt(1/2.*epsilon0.*cspeed.*neff_1.*Aeff_1);
                %chi2_factor_1 = 1i.*(2.*pi.*obj.pulse_1.v_grid)./(2.*neff_1_v0)./cspeed .* obj.waveguide.deff./sqrt(1/2.*epsilon0.*cspeed.*neff_1_v0.*Aeff_1_v0);
                chi2_factor_1(obj.pulse_1.mode.discard_v_index) = 0;
                obj.chi2_factor = chi2_factor_1;
            else
                obj.chi2_factor = [];
            end
            if ~isempty(obj.waveguide.chi3)
                chi3_factor_1 = 1i.*1./4.*(2.*pi.*obj.pulse_1.v_grid)./(2.*neff_1)./cspeed.*obj.waveguide.chi3./(1/2.*epsilon0.*cspeed.*neff_1.*Aeff_1);
                chi3_factor_1(obj.pulse_1.mode.discard_v_index) = 0;
                obj.chi3_factor = chi3_factor_1;
            else
                obj.chi3_factor = [];
            end
            pulse_energy = obj.pulse_1.getpulse_energy();
            a_v_1_co = exp(1i.*beta_1.*0).*obj.pulse_1.a_v;
            obj.initial_pulse = a_v_1_co;
            if ~isempty(obj.pulse_2)
                [neff_2,beta_2,~,Aeff_2] = obj.pulse_2.mode.getmodeProfile('v_grid');
                linear_factor_2 = 1i * beta_2;
                obj.linear_factor = [obj.linear_factor;linear_factor_2];
                if ~isempty(obj.waveguide.deff)
                    chi2_factor_2 = 1i.*(2.*pi.*obj.pulse_2.v_grid)./(2.*neff_2)./cspeed .* obj.waveguide.deff./sqrt(1/2.*epsilon0.*cspeed.*neff_2.*Aeff_2);
                    chi2_factor_2(obj.pulse_2.mode.discard_v_index) = 0;
                    obj.chi2_factor = [obj.chi2_factor;chi2_factor_2];
                end
                if ~isempty(obj.waveguide.chi3)
                    chi3_factor_2 = 1i.*1./4.*(2.*pi.*obj.pulse_2.v_grid)./(2.*neff_2)./cspeed.*obj.waveguide.chi3./(1/2.*epsilon0.*cspeed.*neff_2.*Aeff_2);
                    chi3_factor_2(obj.pulse_2.mode.discard_v_index) = 0;
                    obj.chi3_factor = [obj.chi3_factor;chi3_factor_2];
                end
                pulse_energy = pulse_energy + obj.pulse_2.getpulse_energy();
                a_v_2_co = exp(1i.*beta_2.*0).*obj.pulse_2.a_v;
                obj.initial_pulse = [obj.initial_pulse;a_v_2_co];
            end
            if ~isempty(obj.pulse_3)
                [neff_3,beta_3,~,Aeff_3] = obj.pulse_3.mode.getmodeProfile('v_grid');
                linear_factor_3 = 1i * beta_3;
                obj.linear_factor = [obj.linear_factor;linear_factor_3];
                if ~isempty(obj.waveguide.deff)
                    chi2_factor_3 = 1i.*(2.*pi.*obj.pulse_3.v_grid)./(2.*neff_3)./cspeed .* obj.waveguide.deff./sqrt(1/2.*epsilon0.*cspeed.*neff_3.*Aeff_3);
                    chi2_factor_3(obj.pulse_3.mode.discard_v_index) = 0;
                    obj.chi2_factor = [obj.chi2_factor;chi2_factor_3];
                end
                if ~isempty(obj.waveguide.chi3)
                    chi3_factor_3 = 1i.*1./4.*(2.*pi.*obj.pulse_3.v_grid)./(2.*neff_3)./cspeed.*obj.waveguide.chi3./(1/2.*epsilon0.*cspeed.*neff_3.*Aeff_3);
                    chi3_factor_3(obj.pulse_3.mode.discard_v_index) = 0;
                    obj.chi3_factor = [obj.chi3_factor;chi3_factor_3];
                end
                pulse_energy = pulse_energy + obj.pulse_3.getpulse_energy();
                a_v_3_co = exp(1i.*beta_3.*0).*obj.pulse_3.a_v;
                obj.initial_pulse = [obj.initial_pulse;a_v_3_co];
            end
            average_power = pulse_energy/obj.dv;
            obj.Normalize_factor = sqrt(1/average_power);
            obj.initial_pulse = obj.Normalize_factor * obj.initial_pulse;
            obj.chi2_factor = obj.chi2_factor/obj.Normalize_factor*obj.dv;
            obj.chi3_factor = obj.chi3_factor/obj.Normalize_factor.^2*obj.dv.^2;
        end
            
        function [pulse_1_out,pulse_2_out,pulse_3_out] = Simulate(obj)
            tic;
            OutputFunction = @(z, a_v_co, flag) obj.OutputFunction_MMUPE_SHG(z, a_v_co, flag);
            options = odeset('Stats','on','OutputFcn', OutputFunction,'Refine',1,'MaxStep',obj.OutputInterval,'RelTol',1e-4,'AbsTol',1e-8);
            figure('Name','ode');
            results = ode45(@(z,a_v_co) obj.odefunction_MMUPE_SHG(z,a_v_co), obj.z_grid, obj.initial_pulse, options);
            obj.plotresults_MMUPE(results);
            a_v_end_Normalized = results.y(:,end).*exp(obj.linear_factor.*(obj.z_grid(end)-obj.z_grid(1)));
            a_v_end = a_v_end_Normalized/obj.Normalize_factor;
            a_v_1_end = a_v_end(1:obj.NFFT_1);
            pulse_1_out = obj.pulse_1;pulse_1_out = pulse_1_out.FromPulseSpectrum(a_v_1_end,obj.pulse_1.v0,obj.pulse_1.v_grid);
            if ~isempty(obj.pulse_2)
                a_v_2_end = a_v_end(obj.NFFT_1+1:obj.NFFT_1+obj.NFFT_2);
                pulse_2_out = obj.pulse_2;pulse_2_out = pulse_2_out.FromPulseSpectrum(a_v_2_end,obj.pulse_2.v0,obj.pulse_2.v_grid);
            else
                pulse_2_out = [];
            end
            if ~isempty(obj.pulse_3)
                a_v_3_end = a_v_end(obj.NFFT_1+obj.NFFT_2+1:obj.NFFT_1+obj.NFFT_2+obj.NFFT_3);
                pulse_3_out = obj.pulse_3;pulse_3_out = pulse_3_out.FromPulseSpectrum(a_v_3_end,obj.pulse_3.v0,obj.pulse_3.v_grid);
            else
                pulse_3_out = [];
            end
            toc
        end

        function [dA_v_co_dz] = odefunction_MMUPE_SHG(obj,z,a_v)
            dA_v_co_dz = zeros(size(a_v));
            deltaz = z - obj.z_grid(1);
            a_v_Normalized = a_v.*exp(obj.linear_factor.*deltaz);
            a_v_1 = a_v_Normalized(1:obj.NFFT_1);
            if ~isempty(obj.pulse_2)
                a_t_1_padded = myifft(a_v_1,1,obj.NFFT_3);
                a_v_2 = a_v_Normalized(obj.NFFT_1+1:obj.NFFT_1+obj.NFFT_2);
                a_t_2_padded = myifft(a_v_2,1,obj.NFFT_3);
                if ~isempty(obj.pulse_3)
                    a_v_3 = a_v_Normalized(obj.NFFT_1+obj.NFFT_2+1:obj.NFFT_1+obj.NFFT_2+obj.NFFT_3);
                    a_t_3_padded = myifft(a_v_3,1,obj.NFFT_3);
                end
            else
                if ~isempty(obj.pulse_3)
                    a_t_1_padded = myifft(a_v_1,1,obj.NFFT_3);
                    a_v_3 = a_v_Normalized(obj.NFFT_1+obj.NFFT_2+1:obj.NFFT_1+obj.NFFT_2+obj.NFFT_3);
                    a_t_3_padded = myifft(a_v_3,1,obj.NFFT_3);
                else
                    a_t_1_padded = myirfft(a_v_1,obj.pulse_1.idx,1,obj.pulse_1.antialiasing);
                end
            end

            if ~isempty(obj.waveguide.chi3)
                if ~isempty(obj.pulse_2)
                    if ~isempty(obj.pulse_3)
                        %spmstrength = a_t_1_padded.*conj(a_t_1_padded)+a_t_2_padded.*conj(a_t_2_padded)+a_t_3_padded.*conj(a_t_3_padded);
                        spmstrength = a_t_1_padded.*conj(a_t_1_padded);
                        %conv_1_thg = 3.*myfft(conj(a_t_3_padded).*a_t_1_padded.*a_t_1_padded,1,obj.NFFT_1);
                        conv_1_spm = 3.*myfft(spmstrength.*a_t_1_padded,1,obj.NFFT_1);
                        conv_1 = conv_1_spm;
                        conv_2 = zeros(size(a_v_2));
                        conv_3 = zeros(size(a_v_3));
                        %conv_1 = conv_1_spm + conv_1_thg;
                        %conv_2 = 3.*myfft(spmstrength.*a_t_2_padded,1,obj.NFFT_2);
                        %conv_3_thg = myfft(a_t_1_padded.*a_t_1_padded.*a_t_1_padded,1,obj.NFFT_3);
                        %conv_3_spm = 3.*myfft(spmstrength.*a_t_3_padded,1,obj.NFFT_3);
                        %conv_3 = conv_3_thg + conv_3_spm;
                        conv_total = [conv_1;conv_2;conv_3];
                    else
                        spmstrength = a_t_1_padded.*conj(a_t_1_padded)+a_t_2_padded.*conj(a_t_2_padded);
                        conv_1 = 3.*myfft(spmstrength.*a_t_1_padded,1,obj.NFFT_1);
                        conv_2 = 3.*myfft(spmstrength.*a_t_2_padded,1,obj.NFFT_2);
                        conv_total = [conv_1;conv_2];
                    end
                else
                    if ~isempty(obj.pulse_3)
                        a_t_1_padded = myifft(a_v_1,1,obj.NFFT_3);
                        a_v_3 = a_v_Normalized(obj.NFFT_1+obj.NFFT_2+1:obj.NFFT_1+obj.NFFT_2+obj.NFFT_3);
                        a_t_3_padded = myifft(a_v_3,1,obj.NFFT_3);
                        spmstrength = a_t_1_padded.*conj(a_t_1_padded)+a_t_3_padded.*conj(a_t_3_padded);
                        conv_1_thg = 3.*myfft(conj(a_t_3_padded).*a_t_1_padded.*a_t_1_padded,1,obj.NFFT_1);
                        conv_1_spm = 3.*myfft(spmstrength.*a_t_1_padded,1,obj.NFFT_1);
                        conv_1 = conv_1_spm + conv_1_thg;
                        conv_3_thg = myfft(a_t_1_padded.*a_t_1_padded.*a_t_1_padded,1,obj.NFFT_3);
                        conv_3_spm = 3.*myfft(spmstrength.*a_t_3_padded,1,obj.NFFT_3);
                        conv_3 = conv_3_thg + conv_3_spm;
                        conv_total = [conv_1;conv_3];
                    else
                        %SMUPE
                        conv_1 = myrfft(a_t_1_padded.*a_t_1_padded.*a_t_1_padded,obj.pulse_1.idx,1);
                        conv_total = conv_1;
                    end
                end
                dA_v_co_dz = obj.chi3_factor.*conv_total;
            end
            if ~isempty(obj.waveguide.deff)
                if ~isempty(obj.pulse_2)
                    conv_1 = 2.*myfft(a_t_2_padded.*conj(a_t_1_padded),1,obj.NFFT_1);
                    %conv_1 = zeros(size(conv_1));
                    conv_2 = myfft(a_t_1_padded.*a_t_1_padded,1,obj.NFFT_2);
                    %conv_1 = conv(conj(flip(a_v_1,1)),a_v_2,"same");
                    %conv_2 = conv(a_v_1,a_v_1,"full");
                    if ~isempty(obj.pulse_3)
                        conv_1 = conv_1 + 2.*myfft(a_t_3_padded.*conj(a_t_2_padded),1,obj.NFFT_1);
                        conv_2 = conv_2 + 2.*myfft(a_t_3_padded.*conj(a_t_1_padded),1,obj.NFFT_2);
                        conv_3 = 2.*myfft(a_t_1_padded.*a_t_2_padded,1,obj.NFFT_3);
                    else
                        conv_3 = [];
                    end
                    conv_total = [conv_1;conv_2;conv_3];
                    dA_v_co_dz = dA_v_co_dz + obj.poll_sign(z,obj.waveguide.poll_grid) * obj.chi2_factor .* conv_total;
                else
                    if ~isempty(obj.pulse_3)
                        %Do Nothing
                    else
                        %SMUPE
                        conv_1 = myrfft(a_t_1_padded.*a_t_1_padded,obj.pulse_1.idx,1);
                        conv_total = conv_1;
                        dA_v_co_dz = dA_v_co_dz + obj.poll_sign(z,obj.waveguide.poll_grid) * obj.chi2_factor .* conv_total;
                    end
                end
            end
            dA_v_co_dz = dA_v_co_dz.*exp(-obj.linear_factor.*deltaz);
        end

        function [status] = OutputFunction_MMUPE_SHG(obj, z, a_v_co_Normalized, flag)
            % 初始化状态
            status = 0;
            persistent lastOutputTime;
            persistent SimulationLength;
            persistent Energy_1;
            persistent Energy_2;
            persistent Energy_3;
            switch flag
                case 'init'
                    lastOutputTime = -obj.OutputInterval;
                    SimulationLength = z(1);
                    deltaz = z(1) - obj.z_grid(1);
                    a_v_Normalized = a_v_co_Normalized.*exp(obj.linear_factor.*deltaz);
                    a_v = a_v_Normalized/obj.Normalize_factor;
                    a_v_1 = a_v(1:obj.NFFT_1);
                    Energy_1 = sum(abs(a_v_1).^2).*obj.dv;
                    if ~isempty(obj.pulse_2)
                        a_v_2 = a_v(obj.NFFT_1+1:obj.NFFT_1+obj.NFFT_2);
                        Energy_2 = sum(abs(a_v_2).^2).*obj.dv;
                    else
                        Energy_2 = [];
                    end
                    if ~isempty(obj.pulse_3)
                        a_v_3 = a_v(obj.NFFT_1+obj.NFFT_2+1:obj.NFFT_1+obj.NFFT_2+obj.NFFT_3);
                        Energy_3 = sum(abs(a_v_3).^2).*obj.dv;
                    else
                        Energy_3 = [];
                    end
                case []
                    % 检查是否达到输出时间
                    if z - lastOutputTime >= obj.OutputInterval
                        % 输出结果
                        SimulationLength=[SimulationLength,z];
                        deltaz = z - obj.z_grid(1);
                        a_v_Normalized = a_v_co_Normalized.*exp(obj.linear_factor.*deltaz);
                        a_v = a_v_Normalized/obj.Normalize_factor;
                        a_v_1 = a_v(1:obj.NFFT_1);
                        Energy_1 = [Energy_1,sum(abs(a_v_1).^2).*obj.dv];
                        if ~isempty(obj.pulse_2)
                            a_v_2 = a_v(obj.NFFT_1+1:obj.NFFT_1+obj.NFFT_2);
                            Energy_2 = [Energy_2,sum(abs(a_v_2).^2).*obj.dv];
                        end
                        if ~isempty(obj.pulse_3)
                            a_v_3 = a_v(obj.NFFT_1+obj.NFFT_2+1:obj.NFFT_1+obj.NFFT_2+obj.NFFT_3);
                            Energy_3 = [Energy_3,sum(abs(a_v_3).^2).*obj.dv];
                        end
                        obj.MMUPE_ode_plot(SimulationLength,a_v_co_Normalized,Energy_1,Energy_2,Energy_3);
                        % 更新上次输出时间
                        lastOutputTime = z;
                    end
                case 'done'
            end
        end
        
        function MMUPE_ode_plot(obj,SimulationLength,a_v_co_Normalized,Energy_1,Energy_2,Energy_3)
            cspeed = 299792458;
            lambda_grid_1 = cspeed./obj.pulse_1.v_grid;
            deltaz = SimulationLength(end) - obj.z_grid(1);
            a_v_Normalized = a_v_co_Normalized.*exp((obj.linear_factor-1i./obj.vg.*2.*pi.*obj.v_grid).*deltaz);
            a_v = a_v_Normalized/obj.Normalize_factor;
            a_v_1 = a_v(1:obj.NFFT_1);
            a_t_1 = myifft(a_v_1,obj.dv);
            P_t_1 = abs(a_t_1).^2;
            PSD_1 = abs(a_v_1).^2;
            PSD_1_dBm = 10.*log10(PSD_1);
            Energy_Total = Energy_1(1);
            if ~isempty(obj.pulse_2)
                lambda_grid_2 = cspeed./obj.pulse_2.v_grid;
                a_v_2 = a_v(obj.NFFT_1+1:obj.NFFT_1+obj.NFFT_2);
                a_t_2 = myifft(a_v_2,obj.dv);
                P_t_2 = abs(a_t_2).^2;
                PSD_2 = abs(a_v_2).^2;
                PSD_2_dBm = 10.*log10(PSD_2);
                PSD_2_dBm = PSD_2_dBm - max(obj.PSD_pump_dBm,[],'omitnan');
                Energy_Total = Energy_Total+Energy_2(1);
                SH_proportion = 100 * Energy_2./Energy_Total;
            end
            if ~isempty(obj.pulse_3)
                lambda_grid_3 = cspeed./obj.pulse_3.v_grid;
                a_v_3 = a_v(obj.NFFT_1+obj.NFFT_2+1:obj.NFFT_1+obj.NFFT_2+obj.NFFT_3);
                a_t_3 = myifft(a_v_3,obj.dv);
                P_t_3 = abs(a_t_3).^2;
                PSD_3 = abs(a_v_3).^2;
                PSD_3_dBm = 10.*log10(PSD_3);
                PSD_3_dBm = PSD_3_dBm - max(obj.PSD_pump_dBm,[],'omitnan');
                Energy_Total = Energy_Total+Energy_3(1);
                TH_proportion = 100 * Energy_3./Energy_Total;
            end
            PSD_1_dBm = PSD_1_dBm - max(obj.PSD_pump_dBm,[],'omitnan');
            F_proportion = 100 * Energy_1./Energy_Total;
            
            % Plot pulse results
            subplot(2, 2, 1);
            if ~isempty(obj.pulse_2)
                if ~isempty(obj.pulse_3)
                    plot(obj.pulse_1.t_grid.*1e12,P_t_1,obj.pulse_2.t_grid.*1e12,P_t_2,obj.pulse_3.t_grid.*1e12,P_t_3);
                    legend('F','SH','TH');
                else
                    plot(obj.pulse_1.t_grid.*1e12,P_t_1,obj.pulse_2.t_grid.*1e12,P_t_2);
                    legend('F','SH');
                end
            else
                if ~isempty(obj.pulse_3)
                    plot(obj.pulse_1.t_grid.*1e12,P_t_1,obj.pulse_3.t_grid.*1e12,P_t_3);
                    legend('F','SH');
                else
                    plot(obj.pulse_1.t_grid.*1e12,P_t_1);
                    legend('F');
                end
            end
            xlim([-0.5 0.5]);
            xlabel('Time Window (ps)');
            ylabel('Intensity (W)');
            title('Pulse');
            grid on;
            drawnow limitrate;
            
            subplot(2, 2, 2);
            if ~isempty(obj.pulse_2)
                if ~isempty(obj.pulse_3)
                    %plot(lambda_grid_1.*1e9,PSD_1_dBm',lambda_grid_1*1e9,(obj.PSD_pump_dBm - max(obj.PSD_pump_dBm,[],'omitnan'))',lambda_grid_2.*1e9,PSD_2_dBm',lambda_grid_3.*1e9,PSD_3_dBm');
                    plot(obj.pulse_1.v_grid*1e-12,PSD_1_dBm',obj.pulse_1.v_grid*1e-12,(obj.PSD_pump_dBm - max(obj.PSD_pump_dBm,[],'omitnan'))',obj.pulse_2.v_grid*1e-12,PSD_2_dBm',obj.pulse_3.v_grid*1e-12,PSD_3_dBm');
                else
                    plot(lambda_grid_1.*1e9,PSD_1_dBm',lambda_grid_2.*1e9,PSD_2_dBm');
                end
            else
                if ~isempty(obj.pulse_3)
                    plot(lambda_grid_1.*1e9,PSD_1_dBm',lambda_grid_3.*1e9,PSD_3_dBm');
                else
                    plot(lambda_grid_1.*1e9,PSD_1_dBm');
                end
            end
            xlabel('Frequency (THz)');
            ylabel('Density (dB/Hz)');
            xlim([150 850]);
            %ylim([max(obj.PSD_pump_dBm,[],'omitnan')-55,max(obj.PSD_pump_dBm,[],'omitnan')+5]);
            ylim([-55 5]);
            title('Spectrum');
            legend('F','Pump','SH','TH');
            grid on;
            drawnow limitrate;

            subplot(2, 2, [3 4]);
            if ~isempty(obj.pulse_2)
                if ~isempty(obj.pulse_3)
                    plot(SimulationLength * 1E3,F_proportion,SimulationLength * 1E3,SH_proportion,SimulationLength * 1E3,TH_proportion);
                    legend('F%','SH%','TH%');
                else
                    plot(SimulationLength * 1E3,F_proportion,SimulationLength * 1E3,SH_proportion);
                    legend('F%','SH%');
                end
            else
                if ~isempty(obj.pulse_3)
                    plot(SimulationLength * 1E3,F_proportion,SimulationLength * 1E3,TH_proportion);
                    legend('F%','TH%');
                else
                    plot(SimulationLength * 1E3,F_proportion);
                    legend('F%');
                end
            end
            xlabel('Distance (mm)');
            ylabel('η (%)');
            xlim([SimulationLength(1)*1e3,SimulationLength(end)*1e3]);
            title('Energy Evaluation');
            grid on;
            drawnow limitrate;
        end

        function plotresults_MMUPE(obj,results)
            z = results.x;
            cspeed = 299792458;
            deltaz = z - obj.z_grid(1);
            lambda_grid_1 = cspeed./obj.pulse_1.v_grid;
            a_v_Normalized = results.y.*exp((obj.linear_factor-1i./obj.vg.*2.*pi.*obj.v_grid).*deltaz);
            a_v = a_v_Normalized/obj.Normalize_factor;
            a_v_1 = a_v(1:obj.NFFT_1,:);
            a_t_1_padded = zeros(obj.NFFT_3,size(a_v,2));
            for i = 1:size(a_v,2)
                a_t_1_padded(:,i) = myifft(a_v_1(:,i),obj.dv,obj.NFFT_3);
            end
            PSD_1 = abs(a_v_1).^2;
            PSD_1_dBm = 10.*log10(PSD_1);
            I_t_1 = abs(a_t_1_padded).^2;
            I_t = I_t_1;
            Energy_1 = sum(abs(a_v_1).^2.*obj.dv,1);
            Energy_Total = Energy_1(1);
            figure('Name','Simulation Results');
            if ~isempty(obj.pulse_2)
                lambda_grid_2 = cspeed./obj.pulse_2.v_grid;
                a_v_2 = a_v(obj.NFFT_1+1:obj.NFFT_1+obj.NFFT_2,:);
                a_t_2_padded = zeros(obj.NFFT_3,size(a_v,2));
                for i = 1:size(a_v,2)
                    a_t_2_padded(:,i) = myifft(a_v_2(:,i),obj.dv,obj.NFFT_3);
                end
                I_t_2 = abs(a_t_2_padded).^2;
                PSD_2 = abs(a_v_2).^2;
                PSD_2_dBm = 10.*log10(PSD_2);
                PSD_2_dBm = PSD_2_dBm - max(obj.PSD_pump_dBm,[],'omitnan');
                Energy_2 = sum(abs(a_v_2).^2.*obj.dv,1);
                Energy_Total = Energy_Total+Energy_2(1);
                I_t = I_t + I_t_2;
            end
            if ~isempty(obj.pulse_3)
                lambda_grid_3 = cspeed./obj.pulse_3.v_grid;
                a_v_3 = a_v(obj.NFFT_1+obj.NFFT_2+1:obj.NFFT_1+obj.NFFT_2+obj.NFFT_3,:);
                a_t_3_padded = zeros(obj.NFFT_3,size(a_v,2));
                for i = 1:size(a_v,2)
                    a_t_3_padded(:,i) = myifft(a_v_3(:,i),obj.dv,obj.NFFT_3);
                end
                I_t_3 = abs(a_t_3_padded).^2;
                PSD_3 = abs(a_v_3).^2;
                PSD_3_dBm = 10.*log10(PSD_3);
                PSD_3_dBm = PSD_3_dBm - max(obj.PSD_pump_dBm,[],'omitnan');
                Energy_3 = sum(abs(a_v_3).^2.*obj.dv,1);
                Energy_Total = Energy_Total+Energy_3(1);
                I_t = I_t + I_t_3;
            end
            t_grid = linspace(-1/2,1/2,obj.NFFT_3)'./obj.dv;
            PSD_1_dBm = PSD_1_dBm - max(obj.PSD_pump_dBm,[],'omitnan');
            F_proportion = 100 * Energy_1./Energy_Total;
            subplot(2,2,1);
            imagesc(t_grid * 1E12, z * 1E3, I_t');
            axis xy;
            xlim([-0.5 0.5]);
            xlabel('Time τ (ps)');
            ylabel('Distance (mm)');
            colormap('turbo');
            title('Pulse Evaluation');
            colorbar;

            lambda_linear = (320:1:2000)';
            PSD_1_interp = zeros(length(lambda_linear), size(PSD_1, 2));
            for k = 1:size(PSD_1, 2)
                PSD_1_interp(:,k) = interp1(lambda_grid_1*1e9,PSD_1(:,k),lambda_linear,"linear",0);
            end
            PSD_interp = PSD_1_interp;
            if ~isempty(obj.pulse_2)
                PSD_2_interp = zeros(length(lambda_linear), size(PSD_2, 2));
                for k = 1:size(PSD_1, 2)
                    PSD_2_interp(:,k) = interp1(lambda_grid_2*1e9,PSD_2(:,k),lambda_linear,"linear",0);
                end
                PSD_interp = PSD_interp + PSD_2_interp;
            end
            if ~isempty(obj.pulse_3)
                PSD_3_interp = zeros(length(lambda_linear), size(PSD_3, 2));
                for k = 1:size(PSD_3, 2)
                    PSD_3_interp(:,k) = interp1(lambda_grid_3*1e9,PSD_3(:,k),lambda_linear,"linear",0);
                end
                PSD_interp = PSD_interp + PSD_3_interp;
            end
            PSD_interp_dBm = 10.*log10(PSD_interp);
            PSD_interp_dBm = PSD_interp_dBm - max(PSD_interp_dBm,[],'omitnan');
            subplot(2, 2, 2);
            imagesc(lambda_linear, z * 1E3, PSD_interp_dBm');
            axis xy;
            xlabel('λ (nm)');
            ylabel('Simulation Length (mm)');
            colormap('turbo');
            title('Pulse Evaluation');
            clim([max(PSD_1_dBm(:,end),[],'omitnan')-55,max(PSD_1_dBm(:,end),[],'omitnan')+5]);
            h = colorbar;
            h.Label.String = 'Density (dB/Hz)';

            subplot(2,2,3);
            if ~isempty(obj.pulse_2)
                if ~isempty(obj.pulse_3)
                    %plot(lambda_grid_1*1e9,PSD_1_dBm(:,end),lambda_grid_1*1e9,obj.PSD_pump_dBm - max(obj.PSD_pump_dBm,[],'omitnan'),lambda_grid_2*1e9,PSD_2_dBm(:,end),lambda_grid_3*1e9,PSD_3_dBm(:,end));
                    plot(obj.pulse_1.v_grid*1e-12,PSD_1_dBm(:,end),obj.pulse_1.v_grid*1e-12,obj.PSD_pump_dBm - max(obj.PSD_pump_dBm,[],'omitnan'),obj.pulse_2.v_grid*1e-12,PSD_2_dBm(:,end),obj.pulse_3.v_grid*1e-12,PSD_3_dBm(:,end));
                else
                    plot(lambda_grid_1*1e9,PSD_1_dBm(:,end),lambda_grid_2*1e9,PSD_2_dBm(:,end));
                end
            else
                if ~isempty(obj.pulse_3)
                    plot(lambda_grid_1*1e9,PSD_1_dBm(:,end),lambda_grid_3*1e9,PSD_3_dBm(:,end));
                else
                    plot(lambda_grid_1*1e9,PSD_1_dBm(:,end));
                end
            end
            xlabel('Frequency (THz)');
            ylabel('Density (dB/Hz)');
            %xlim([350 1750]);
            xlim([150 850]);
            %ylim([max(PSD_1_dBm(:,end),[],'omitnan')-55,max(PSD_1_dBm(:,end),[],'omitnan')+5]);
            %ylim([max(obj.PSD_pump_dBm,[],'omitnan')-55,max(obj.PSD_pump_dBm,[],'omitnan')+5]);
            ylim([-55 5]);
            title('Spectrum');
            legend('F','Pump','SH','TH');
            grid on;

            subplot(2,2,4);
            if ~isempty(obj.pulse_2)
                SH_proportion = 100 * Energy_2./Energy_Total;
                if ~isempty(obj.pulse_3)
                    TH_proportion = 100 * Energy_3./Energy_Total;
                    plot(z * 1E3,F_proportion,z * 1E3,SH_proportion,z * 1E3,TH_proportion);
                    legend('F%','SH%','TH%')
                else
                    plot(z * 1E3,F_proportion,z * 1E3,SH_proportion);
                    legend('F%','SH%')
                end
            else
                if ~isempty(obj.pulse_3)
                    TH_proportion = 100 * Energy_3./Energy_Total;
                    plot(z * 1E3,F_proportion,z * 1E3,TH_proportion);
                    legend('F%','TH%')
                else
                    plot(z * 1E3,F_proportion);
                    legend('F%')
                end
            end
            xlabel('Distance (mm)');
            ylabel('η (%)');
            title('Energy Evaluation');
            grid on;
            drawnow limitrate;

            figure;
            if ~isempty(obj.pulse_2)
                if ~isempty(obj.pulse_3)
                    %plot(lambda_grid_1*1e9,PSD_1_dBm(:,end),lambda_grid_1*1e9,obj.PSD_pump_dBm - max(obj.PSD_pump_dBm,[],'omitnan'),lambda_grid_2*1e9,PSD_2_dBm(:,end),lambda_grid_3*1e9,PSD_3_dBm(:,end));
                    plot(obj.pulse_1.v_grid*1e-12,PSD_1_dBm(:,end),obj.pulse_1.v_grid*1e-12,obj.PSD_pump_dBm - max(obj.PSD_pump_dBm,[],'omitnan'),obj.pulse_2.v_grid*1e-12,PSD_2_dBm(:,end),obj.pulse_3.v_grid*1e-12,PSD_3_dBm(:,end));
                else
                    plot(lambda_grid_1*1e9,PSD_1_dBm(:,end),lambda_grid_2*1e9,PSD_2_dBm(:,end));
                end
            else
                if ~isempty(obj.pulse_3)
                    plot(lambda_grid_1*1e9,PSD_1_dBm(:,end),lambda_grid_3*1e9,PSD_3_dBm(:,end));
                else
                    plot(lambda_grid_1*1e9,PSD_1_dBm(:,end));
                end
            end
            xlabel('Frequency (THz)');
            ylabel('Density (dB/Hz)');
            %xlim([350 1750]);
            xlim([150 850]);
            %ylim([max(PSD_1_dBm(:,end),[],'omitnan')-55,max(PSD_1_dBm(:,end),[],'omitnan')+5]);
            %ylim([max(obj.PSD_pump_dBm,[],'omitnan')-55,max(obj.PSD_pump_dBm,[],'omitnan')+5]);
            ylim([-55 5]);
            grid on;
            legend('F','Pump','SH','TH');
            % 创建顶部坐标轴
            ax2 = axes('Color', 'none');
            % 设置顶部轴属性
            %ax2.XColor = 'r'; % 顶部轴颜色
            ax2.YColor = 'none'; % 隐藏右侧 Y 轴
            ax2.XAxisLocation = 'top'; % 设置 X 轴在顶部
            %ax2.FontSize = 14;
            xlabel('Wavelength (nm)');
            %ax2.YTick = []; % 隐藏 Y 轴刻度
            % 自定义顶部坐标轴刻度
            ax2.XTick = [cspeed/2000e-9 cspeed/1600e-9 cspeed/1200e-9 cspeed/1000e-9 cspeed/800e-9 cspeed/600e-9 cspeed/500e-9 cspeed/400e-9]; % 设定刻度
            ax2.XTickLabel = {'2000' '1600', '1200', '1000', '800', '600', '500','400'}; % 设定刻度标签
            ax2.XLim = [135e12,780e12];
            %title('Spectrum');
        end
    end

    methods(Static)
        %pollsign
        function [sign] = poll_sign(z,poll_grid)
            index = find(poll_grid >= z,1);
            sign = 1 - 2 * mod(index,2);
        end
    end
end