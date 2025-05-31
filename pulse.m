classdef pulse
    properties
        NFFT
        a_t
        a_v
        v0
        v_c
        dv
        v_grid
        dt
        t_grid
        idx
        mode
        antialiasing
    end
    
    methods(Static)
        %Default v_grid_base
        function [v_grid,dv,v_c] = generate_v_grid(lambda_range,NFFT)
            cspeed = 299792458;
            v_min = cspeed./lambda_range(end);
            v_max = cspeed./lambda_range(1);
            dv = (v_max-v_min)/NFFT;
            idx_min = floor(v_min/dv);
            idx_max = idx_min + NFFT - 1;
            idx_c = idx_min + floor(NFFT/2);
            v_c = idx_c*dv;
            v_grid = (idx_min:idx_max)'.*dv;
        end

        function [v_grid] = generate_v_grid_from_v_c(v_c,dv,NFFT)
            idx_c = v_c/dv;
            idx_min = idx_c - floor(NFFT/2);
            idx_max = idx_min + NFFT - 1;
            v_grid = (idx_min:idx_max)'.*dv;
        end


        function v_grid = v_grid_harmonic(v_grid_F,harmonic_order)
            dv = v_grid_F(2) - v_grid_F(1);
            idx_min = floor(v_grid_F(1)/dv);
            idx_min = harmonic_order*idx_min;
            N = length(v_grid_F);
            NFFT = harmonic_order.*(N-mod(N,2))+mod(N,2);
            idx_max = idx_min + NFFT - 1;
            v_grid = (idx_min:idx_max)'.*dv;
        end
    end

    methods
        % Define Class
        function obj = pulse(v_grid,v0,mode,antialiasing)
            if nargin > 0
                obj.v_grid = v_grid;
                obj.mode = mode;
                obj = obj.complete_pulse;
                obj.v0 = v0;
                obj.a_t = zeros(obj.NFFT,1);
                obj.a_v = myfft(obj.a_t,obj.dv);
            end
            if nargin > 3
                obj.antialiasing = antialiasing;
            else
                obj.antialiasing = 3;
            end
        end

        function obj = complete_pulse(obj)
            obj.NFFT = length(obj.v_grid);
            obj.dv = abs(obj.v_grid(2)-obj.v_grid(1));
            idx_min = floor(obj.v_grid(1)/obj.dv);
            idx_max = idx_min + obj.NFFT - 1;
            obj.idx = [idx_min,idx_max];
            idx_c = idx_min + floor(obj.NFFT/2);
            obj.v_c = obj.dv*idx_c;
            obj.dt = 1 / obj.NFFT / obj.dv;
            obj.t_grid = linspace(-1/2,1/2,obj.NFFT)'.*obj.dt.*obj.NFFT;
        end

        %Generate Pulse From a_t
        function obj = FromPulseShape(obj,a_t,v0,v_grid)
            obj.v0 = v0;
            obj.v_grid = v_grid;
            obj = obj.complete_pulse;
            obj.a_t = a_t;
            obj.a_v = myfft(a_t,obj.dv);
        end

        %Generate Pulse From a_v
        function obj = FromPulseSpectrum(obj,a_v,v0,v_grid)
            obj.v0 = v0;
            obj.v_grid = v_grid;
            obj = obj.complete_pulse;
            obj.a_v = a_v;
            obj.a_t = myifft(a_v,obj.dv);
        end

        %Generate Pulse From Gaussian
        function obj = Gasussian(obj,pulse_energy,t_fwhm,v0,v_grid)
            obj.v0 = v0;
            obj.v_grid = v_grid;
            obj = obj.complete_pulse;
            p_t = 2.^(-(obj.t_grid/(t_fwhm/2)).^2);
            p_t = p_t./(sum(p_t).*obj.dt).*pulse_energy;
            phi_t = 2.*pi.*(v0-obj.v_c).*obj.t_grid;
            obj.a_t = p_t.^0.5 .* exp(1i.*phi_t);
            obj.a_v = myfft(obj.a_t,obj.dv);
            obj.a_v(obj.mode.discard_v_index)=0;
            obj.a_t = myifft(obj.a_v,obj.dv);
        end

        function obj = CW(obj,Power,v0,v_grid)
            obj.v_grid = v_grid;
            obj = obj.complete_pulse;
            diff_grid = abs(v_grid - v0);
            [~,idx_min] = min(diff_grid);
            obj.v0 = v_grid(idx_min);
            %phi_t = 2.*pi.*(v0-obj.v_c).*obj.t_grid;
            %obj.a_t = Power.^0.5 .* exp(1i.*phi_t);
            obj.a_v = zeros(obj.NFFT,1);
            obj.a_v(idx_min) = sqrt(Power)/obj.dv;
            obj.a_v = obj.a_v + obj.vacuum_noise();
            obj.a_v(obj.mode.discard_v_index)=0;
            obj.a_t = myifft(obj.a_v,obj.dv);
        end

        function [E_p] = getpulse_energy(obj)
            E_p = sum(abs(obj.a_v).^2).*obj.dv;
            %E_p = sum(abs(obj.a_t).^2).*obj.dt;
        end

        function [A_v] = getA_v(obj)
            cspeed = 299792458;
            epsilon0 = 8.854187817e-12;
            A_v = obj.a_v./sqrt(epsilon0.*cspeed.*obj.mode.neff.*obj.mode.Aeff./2);
        end

        function obj = readA_v(obj,A_v)
            cspeed = 299792458;
            epsilon0 = 8.854187817e-12;
            obj.a_v = A_v.*sqrt(epsilon0.*cspeed.*obj.mode.neff.*obj.mode.Aeff./2);
            obj.a_t = myifft(obj.a_v,obj.dv);
        end

        function [a_v] = A_v_to_a_v(obj,A_v)
            cspeed = 299792458;
            epsilon0 = 8.854187817e-12;
            a_v = A_v.*sqrt(epsilon0.*cspeed.*obj.mode.neff.*obj.mode.Aeff./2);
        end

        function obj = bandpass_filter_reconstruct(obj,lambda_min,lambda_max)
            cspeed = 299792458;
            v_max = cspeed/lambda_min;
            v_min = cspeed/lambda_max;
            idx_min = find(obj.v_grid>v_min,1);
            idx_max = find(obj.v_grid>v_max,1);
            obj.a_v(1:idx_min) = 0;
            obj.a_v(idx_max:end) = 0;
            idx_center = round((idx_max+idx_min)/2);
            idx_span = idx_max - idx_min;
            NFFT_n = 2^ceil(1+log(idx_span)/log(2));
            if (NFFT_n<obj.NFFT)
                if (idx_center+NFFT_n/2-1>obj.NFFT)
                    obj.v_grid = obj.v_grid(obj.NFFT-NFFT_n:obj.NFFT);
                    obj.a_v = obj.a_v(obj.NFFT-NFFT_n:obj.NFFT);
                else
                    if (idx_center-NFFT_n/2<1)
                        obj.v_grid = obj.v_grid(1:NFFT_n);
                        obj.a_v = obj.a_v(1:NFFT_n);
                    else
                        obj.v_grid = obj.v_grid(idx_center-NFFT_n/2:idx_center+NFFT_n/2-1);
                        obj.a_v = obj.a_v(idx_center-NFFT_n/2:idx_center+NFFT_n/2-1);
                    end
                end
                obj = obj.complete_pulse;

            end
            obj.a_t = myifft(obj.a_v,obj.dv);
        end

        function obj = bandpass_filter(obj,lambda_min,lambda_max)
            cspeed = 299792458;
            v_max = cspeed/lambda_min;
            v_min = cspeed/lambda_max;
            idx_min = find(obj.v_grid>v_min,1);
            idx_max = find(obj.v_grid>v_max,1);
            obj.a_v(1:idx_min) = 0;
            obj.a_v(idx_max:end) = 0;
            obj.a_t = myifft(obj.a_v,obj.dv);
        end

        function obj = add_delay(obj,delay_time)
            delay_idx = round(delay_time./obj.dt);
            obj.a_t = circshift(obj.a_t,delay_idx);
            obj.a_v = myfft(obj.a_t,obj.dv);
        end
        
        function [noise] = vacuum_noise(obj)
            h = 6.62607015e-34;
            S_v = sqrt(h * obj.v_grid / 2);
            noise = randn(size(obj.v_grid));
            noise_freq = fft(noise);
            noise_freq = noise_freq .* S_v;
            noise = ifft(noise_freq);
            noise = obj.A_v_to_a_v(noise);
        end

        function obj = add_dispersion(obj,GDD,wavelength)
            cspeed = 299792458;
            v_center = cspeed/wavelength;
            phase = exp(1i.*GDD./2.*(2.*pi.*(obj.v_grid-v_center)).^2);
            obj.a_v = obj.a_v.*phase;
            obj.a_t = myifft(obj.a_v,obj.dv);
        end

    end
end