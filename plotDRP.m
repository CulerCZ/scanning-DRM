% function to show DRP in 3d northern hemisphere
function plotDRP(drplist, posInfo, options)
    arguments
        drplist
        posInfo
        options.cMap (1,1) string = "jet"
        options.type (1,1) string = "3d"
    end
    nn = length(posInfo.phi);
    theta_list_original = sort(unique(posInfo.theta),'ascend');
    theta_list = zeros(1,length(theta_list_original)+1);
    for ii = 1:length(theta_list_original)
        if theta_list_original(ii) <= 70
            theta_list(ii) = theta_list_original(ii);
        else
            theta_list(ii) = (theta_list_original(ii) + theta_list_original(ii-1)) / 2;
        end
    end
    theta_list(end) = 90;
    hold on
    for ii = 1:nn
        % datapoints in the lower range
        if posInfo.theta(ii) <= 70
            phi_temp = [posInfo.phi(ii)-5; posInfo.phi(ii)+5; posInfo.phi(ii)+5; posInfo.phi(ii)-5];
            kk = find(theta_list_original == posInfo.theta(ii));
            theta_temp = [theta_list(kk); theta_list(kk); theta_list(kk+1); theta_list(kk+1)];
            [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
            patch(x,y,z,drplist(ii),"EdgeColor","none")
        else
            switch find(theta_list_original == posInfo.theta(ii))
                case 17
                    phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(16))),'ascend');
                    jj = floor(posInfo.phi(ii) / 22.5);
                    phi_lower = phi_list_temp(phi_list_temp > jj*22.5 & phi_list_temp < (jj+1)*22.5);
                    phi_temp = [22.5*jj; phi_lower'; 22.5*(jj+1); 22.5*(jj+1); 22.5*jj];
                    kk = find(theta_list_original == posInfo.theta(ii));
                    theta_temp = [theta_list(kk)*ones(length(phi_lower)+2,1); theta_list(kk+1); theta_list(kk+1)];
                    [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
                    patch(x,y,z,drplist(ii),"EdgeColor","none")
                case 18
                    phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(17))),'ascend');
                    jj = floor(posInfo.phi(ii) / 30);
                    phi_lower = phi_list_temp(phi_list_temp > jj*30 & phi_list_temp < (jj+1)*30);
                    phi_temp = [30*jj; phi_lower'; 30*(jj+1); 30*(jj+1); 30*jj];
                    kk = find(theta_list_original == posInfo.theta(ii));
                    theta_temp = [theta_list(kk)*ones(length(phi_lower)+2,1); theta_list(kk+1); theta_list(kk+1)];
                    [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
                    patch(x,y,z,drplist(ii),"EdgeColor","none")
                case 19
                    phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(18))),'ascend');
                    jj = floor(posInfo.phi(ii) / 90);
                    phi_lower = phi_list_temp(phi_list_temp > jj*90 & phi_list_temp < (jj+1)*90);
                    phi_temp = [90*jj; phi_lower'; 90*(jj+1); 90*(jj+1); 90*jj];
                    kk = find(theta_list_original == posInfo.theta(ii));
                    theta_temp = [theta_list(kk)*ones(length(phi_lower)+2,1); theta_list(kk+1); theta_list(kk+1)];
                    [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
                    patch(x,y,z,drplist(ii),"EdgeColor","none")
            end
        end
    end
    % axis and plot setting
    axis equal
    colormap(options.cMap)
    set(gca,"Visible","off")
end