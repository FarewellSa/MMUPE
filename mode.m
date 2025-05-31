classdef mode
    properties
        v_grid
        neff
        beta
        Aeff
        neff_v0
        beta_v0
        beta1_v0
        Aeff_v0
        lambda_range
        discard_v_index
        neff_function
        beta_function
        beta1_function
        beta2_function
        Gm_function
        Gmvv1_function
        GVM_function
        dGmvv1_function
    end

    methods
        % Define Class
        function obj = mode(material,neff,muneff,Aeff,muAeff,v0,v_grid,lambda_range)
            cspeed = 299792458;
            flag = [];
            v_range = fliplr(cspeed./lambda_range);
            obj.v_grid = v_grid;
            obj.discard_v_index = find((v_grid < v_range(1)) | (v_grid > v_range(end)));
            switch material
                case 'FromFEM'
                    flag = 1;
                    syms v vnormalized;
                    neff_function(vnormalized) = poly2sym(neff,vnormalized);
                    neff_function(v) = subs(neff_function(vnormalized),vnormalized,(v-muneff(1))/muneff(2));
                    beta_function(v) = 2.*pi.*neff_function(v).*v./cspeed;
                    beta1_function(v) = diff(beta_function,1)/2/pi;
                    beta2_function(v) = diff(beta1_function,1)/2/pi;
                    obj.beta1_function = beta1_function;
                    obj.beta2_function = beta2_function;
                    obj.neff = polyval(neff,v_grid,[],muneff);
                    obj.beta = 2.*pi.*obj.neff.*v_grid./cspeed;
                    obj.neff_v0 = polyval(neff,v0,[],muneff);
                    obj.beta_v0 = 2.*pi.*obj.neff_v0.*v0./cspeed;
                    obj.beta1_v0 = double(beta1_function(v0));
                    obj.Aeff = polyval(Aeff,v_grid,[],muAeff);
                    obj.Aeff_v0 = polyval(Aeff,v0,[],muAeff);
                case 'vacuum'
                    syms v
                    neff_function(v) = 1;
                case 'LiNbO3_bulk_e'
                    syms v
                    lambdaum(v) = 1e6.*cspeed./v;
                    neff_function(v) = sqrt(1+2.9804./(1-0.02047./lambdaum(v).^2)+0.5981./(1-0.0666./lambdaum(v).^2)+8.9543./(1-416.08./lambdaum(v).^2));
                case 'YAG'
                    syms v
                    lambdaum(v) = 1e6.*cspeed./v;
                    neff_function(v) = sqrt(1+2.28200./(1-0.01185./lambdaum(v).^2)+3.27644./(1-282.734./lambdaum(v).^2));
                case 'MgO:LiNbO3_bulk_e'
                    syms v
                    lambdaum(v) = 1e6.*cspeed./v;
                    neff_function(v) = sqrt(5.756+0.0983./(lambdaum(v).^2-0.2020^2)+189.32./(lambdaum(v).^2-12.52^2)-1.32e-2.*lambdaum(v).^2);
                case ohter
                    disp('unsupported material');
            end
            if isempty(flag)
                obj.Aeff = Aeff.*ones(size(v_grid));
                obj.Aeff_v0 = Aeff;
                beta_function(v) = 2.*pi.*neff_function(v).*v./cspeed;
                beta1_function(v) = diff(beta_function,1)/2/pi;
                beta2_function(v) = diff(beta1_function,1)/2/pi;
                obj.beta1_function = beta1_function;
                obj.beta2_function = beta2_function;
                obj.neff = double(neff_function(v_grid));
                obj.beta = 2.*pi.*obj.neff.*v_grid./cspeed;
                obj.neff_v0 = double(neff_function(v0));
                obj.beta_v0 = 2.*pi.*obj.neff_v0.*v0./cspeed;
                obj.beta1_v0 = double(beta1_function(v0));
            end
            if ~isempty(obj.discard_v_index)
                obj.Aeff(obj.discard_v_index) = 0;
                obj.neff(obj.discard_v_index)= 0;
                obj.beta(obj.discard_v_index)= 0;
            end
            if ((v0 < v_range(1)) || (v0 > v_range(end)))
                obj.Aeff_v0 = 0;
                obj.neff_v0 = 0;
                obj.beta_v0 = 0;
                obj.beta1_v0 = 0;
            end
            obj.neff_function = neff_function;
            obj.beta_function = beta_function;
        end

        function [neff,beta,beta1,Aeff] = getmodeProfile(obj,flag)
            switch flag
                case 'v_grid'
                    neff = obj.neff;
                    beta = obj.beta;
                    beta1 = obj.beta1_v0;
                    Aeff = obj.Aeff;
                case 'v0'
                    neff = obj.neff_v0;
                    beta = obj.beta_v0;
                    beta1 = obj.beta1_v0;
                    Aeff = obj.Aeff_v0;
            end
        end

        function obj = calculate_dGmdv_function(obj)
            syms v v1;
            Gm(v) = obj.beta_function(v)-2*obj.beta_function(v/2);
            Gmvv1(v,v1) = obj.beta_function(v)-obj.beta_function(v1)-obj.beta_function(v-v1);
            obj.Gm_function = Gm;
            obj.Gmvv1_function = Gmvv1;
            obj.GVM_function = diff(Gm,1)/2/pi;
            obj.dGmvv1_function = diff(Gmvv1,v1,1);
        end
    end
end