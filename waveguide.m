classdef waveguide
    properties
        material
        deff
        chi3
        poll_grid
        Dc
        Leff
        z_v_grid
        error
        d
    end

    methods
        % Define Class
        function obj = waveguide(material)
            obj.material = material;
            switch material
                case 'vacuum'
                    obj.deff = [];
                    obj.chi3 = [];
                case 'LiNbO3_eee'
                    obj.deff = 27e-12;
                    obj.chi3 = 5200e-24;
                case 'LiNbO3_eoo'
                    obj.deff = 5e-12;
                    obj.chi3 = 5200e-24;
                case 'YAG'
                    obj.deff = [];
                    obj.chi3 = 800e-24;
                case 'MgO:LiNbO3_eee'
                    obj.deff = 20e-12;
                    obj.chi3 = [];
                otherwise
                    error('Unsupported material');
            end
            obj.poll_grid = [0,Inf];
        end

        % Periodically poll
        function obj = periodicallypoll(obj,poll_range,LAM)
            z = poll_range(1);
            obj.poll_grid(end) = z;
            while z < poll_range(2)
                z = z + LAM / 2;
                obj.poll_grid = [obj.poll_grid, z];
            end
            obj.poll_grid = [obj.poll_grid, Inf];
        end

        function obj = linearchirppoll(obj,poll_range,lambda_poll_range,LAM_order,mode)
            z = poll_range(1);
            obj.poll_grid(end) = z;
            Polling_start = poll_range(1);
            Polling_end = poll_range(2);
            beta_function = mode.beta_function;
            cspeed = 299792458; % 光速，单位为 m/s
            v_min_poll = cspeed / lambda_poll_range(2);
            v_max_poll = cspeed / lambda_poll_range(1);
            z = Polling_start;
            beta_F_min = double(beta_function(v_min_poll));
            beta_SHG_min = double(beta_function(2 * v_min_poll));
            beta_F_max = double(beta_function(v_max_poll));
            beta_SHG_max = double(beta_function(2 * v_max_poll));
            A_start = 2 * pi / (beta_SHG_min - 2 * beta_F_min) * LAM_order;
            A_end = 2 * pi / (beta_SHG_max - 2 * beta_F_max) * LAM_order;
            disp(['A_start: ', num2str(A_start)]);
            disp(['A_end: ', num2str(A_end)]);
            Dg = 2 * pi * (A_start - A_end) / (A_start * A_end * (Polling_end - Polling_start));
            obj.Dc = Dg*ones(size(mode.v_grid));
            obj.Leff = sqrt(2*pi)./sqrt(obj.Dc);
            lambda_grid_SH = cspeed./mode.v_grid;
            for i = 1:length(mode.v_grid)
                if lambda_grid_SH(i)<lambda_poll_range(1)||lambda_grid_SH(i)>lambda_poll_range(2)
                    %obj.Dc(i) = Inf;
                end
            end
            disp(['Dg: ', num2str(Dg)]);
            obj.poll_grid = [];
            while z < Polling_end
                A_z = A_start / (1 + Dg * (z - Polling_start) * A_start / (2 * pi));
                z = z + A_z / 2;
                if z < Polling_end
                    obj.poll_grid = [obj.poll_grid, z]; 
                end
            end
            deltak = double(mode.Gm_function(mode.v_grid));
            deltak_min = beta_SHG_min - 2 * beta_F_min;
            obj.z_v_grid = (deltak-deltak_min)/Dg + Polling_start;
            obj.poll_grid = [obj.poll_grid, Polling_end];
        end

        function obj = linearchirppoll_eoo(obj,poll_range,lambda_poll_range,LAM_order,mode_o,mode_e)
            z = poll_range(1);
            obj.poll_grid(end) = z;
            Polling_start = poll_range(1);
            Polling_end = poll_range(2);
            beta_function_o = mode_o.beta_function;
            beta_function_e = mode_e.beta_function;
            cspeed = 299792458; % 光速，单位为 m/s
            v_min_poll = cspeed / lambda_poll_range(2);
            v_max_poll = cspeed / lambda_poll_range(1);
            z = Polling_start;
            beta_F_min = double(beta_function_o(v_min_poll));
            beta_SHG_min = double(beta_function_e(2 * v_min_poll));
            beta_F_max = double(beta_function_o(v_max_poll));
            beta_SHG_max = double(beta_function_e(2 * v_max_poll));
            A_start = 2 * pi / (beta_SHG_min - 2 * beta_F_min) * LAM_order;
            A_end = 2 * pi / (beta_SHG_max - 2 * beta_F_max) * LAM_order;
            disp(['A_start: ', num2str(A_start)]);
            disp(['A_end: ', num2str(A_end)]);
            Dg = 2 * pi * (A_start - A_end) / (A_start * A_end * (Polling_end - Polling_start));
            obj.Dc = Dg*ones(size(mode_e.v_grid));
            disp(['Dg: ', num2str(Dg)]);
            while z < Polling_end
                A_z = A_start / (1 + Dg * (z - Polling_start) * A_start / (2 * pi));
                z = z + A_z / 2;
                if z < Polling_end
                    obj.poll_grid = [obj.poll_grid, z]; 
                end
            end
            deltak = double(beta_function_e(mode_e.v_grid)-2*beta_function_o(mode_e.v_grid/2));
            figure(2);
            plot(cspeed./mode_e.v_grid.*1e9,2*pi./deltak.*1e6);
            xlim([350 1000]);
            xlabel('λ (nm)');
            ylabel('Period (um)');
            title('Required polling period');
            grid on;
            data = [cspeed./mode_e.v_grid.*1e9,2*pi./deltak.*1e6];
            writematrix(data,'output1.xlsx');
            deltak_min = beta_SHG_min - 2 * beta_F_min;
            obj.z_v_grid = (deltak-deltak_min)/Dg + Polling_start;
            obj.poll_grid = [obj.poll_grid, Inf];
        end

        function obj = custom_function_poll(obj,poll_range,poll_function,duty_ratio)
            z = poll_range(1);
            obj.poll_grid(end) = z;
            Polling_end = poll_range(2);
            while z < Polling_end
                A_z = poll_function(z);
                z = z + A_z* duty_ratio;
                if z < Polling_end
                    obj.poll_grid = [obj.poll_grid, z]; 
                end
                A_z = poll_function(z);
                z = z + A_z* (1-duty_ratio);
                if z < Polling_end
                    obj.poll_grid = [obj.poll_grid, z]; 
                end
            end
            obj.poll_grid = [obj.poll_grid, Inf];
        end

        function obj = linearchirppoll_Arange(obj,poll_range,A_range)
            z = poll_range(1);
            obj.poll_grid(end) = z;
            Polling_start = poll_range(1);
            Polling_end = poll_range(2);
            A_start = A_range(1);
            A_end = A_range(2);
            disp(['A_start: ', num2str(A_start)]);
            disp(['A_end: ', num2str(A_end)]);
            Dg = 2 * pi * (A_start - A_end) / (A_start * A_end * (Polling_end - Polling_start));
            disp(['Dg: ', num2str(Dg)]);
            while z < Polling_end
                A_z = A_start / (1 + Dg * (z - Polling_start) * A_start / (2 * pi));
                z = z + A_z / 2;
                if z < Polling_end
                    obj.poll_grid = [obj.poll_grid, z]; 
                end
            end
            obj.poll_grid = [obj.poll_grid, Inf];
        end

        function obj = adaptivechirppoll(obj,poll_range,lambda_poll_range,LAM_order,v_grid_SH,mode_SH,dz_v,dz_v_threshold)
            cspeed = 299792458;
            lambda_grid_SH = cspeed./v_grid_SH;
            lambda_poll_range = lambda_poll_range./2;
            dz_v(isnan(dz_v)) = Inf;
            dz_v = dz_v./min(dz_v);
            dz_v_raw = dz_v;
            idx = dz_v>dz_v_threshold;
            dz_v(idx) = 0;
            indices = false(size(v_grid_SH));
            for i = 1:size(lambda_poll_range, 1)
                lower_bound = lambda_poll_range(i, 1);
                upper_bound = lambda_poll_range(i, 2);
                indices = indices | (lambda_grid_SH >= lower_bound & lambda_grid_SH <= upper_bound);
            end
            dz_v(~indices) = 0;
            %plot(lambda_grid_SH.*1e9,dz_v);
            z_v = zeros(size(dz_v));
            for i = 1:length(dz_v)
                z_v(i) = sum(dz_v(1:i));
            end
            dz_v = dz_v./max(z_v).*(poll_range(end)-poll_range(1));
            dz_v_raw = dz_v_raw./max(z_v).*(poll_range(end)-poll_range(1));
            z_v = z_v./max(z_v);
            dv = abs(mode_SH.v_grid(2)-mode_SH.v_grid(1));
            obj.Dc = double(mode_SH.GVM_function(mode_SH.v_grid))*2*pi*dv./dz_v_raw;
            obj.Leff = sqrt(2*pi)./sqrt(obj.Dc);
            obj.Dc = double(mode_SH.GVM_function(mode_SH.v_grid))*2*pi*dv./dz_v;
            z_v = (poll_range(end)-poll_range(1)).*z_v + poll_range(1);

            z = poll_range(1);
            obj.poll_grid(end) = z;
            while z < poll_range(end)
                idx = find(z_v > z,1);
                A_z = 2*pi/double(mode_SH.Gm_function(v_grid_SH(idx)))*LAM_order;
                z = z + A_z / 2;
                if z < poll_range(end)
                    obj.poll_grid = [obj.poll_grid, z]; 
                end
            end
            obj.z_v_grid = z_v;
            obj.poll_grid = [obj.poll_grid, Inf];
        end

        function obj = adderror(obj,error_range)
            obj.error = error_range/2+error_range/2.*(rand(size(obj.z_v_grid)));
            obj.error(1:2:end) = 0;
            obj.error(1) = 0;
            obj.error(end) = 0;
            obj.z_v_grid = obj.z_v_grid + obj.error;
        end
    end
end