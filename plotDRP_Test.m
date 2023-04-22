phi = posInfo.phi';
theta = posInfo.theta';
[x,y,z] = sph2cart(phi.*degree, theta.*degree, ones(size(phi)));
figure(Position=[100 100 1000 1000])
scatter3(x,y,z,108,"blue","filled");
axis equal
xlabel("x")
ylabel("y")


%% create patch grid
angle = [phi, theta];
phi_temp = [-5 5 5 -5;5 10 10 5]';
theta_temp = [0 0 5 5;5 5 10 10]';
[x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
figure, patch(x,y,z,'r')

%% craete patch grid (phi, theta) using input pos_info data
nn = length(posInfo.phi);
phi_temp = zeros(4,nn);
theta_temp = zeros(4,nn);
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

for ii = 1:nn
    % datapoints in the lower range
    if posInfo.theta(ii) <= 70
        phi_temp(:,ii) = [posInfo.phi(ii)-5; posInfo.phi(ii)+5; posInfo.phi(ii)+5; posInfo.phi(ii)-5];
        kk = find(theta_list_original == posInfo.theta(ii));
        theta_temp(:,ii) = [theta_list(kk); theta_list(kk); theta_list(kk+1); theta_list(kk+1)];
    else
        switch find(theta_list_original == posInfo.theta(ii))
            case 17
                % phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(16))),'ascend');
                jj = floor(posInfo.phi(ii) / 22.5);
                % phi_lower = phi_list_temp(phi_list_temp > jj*22.5 & phi_list_temp < (jj+1)*22.5);
                phi_temp(:,ii) = 22.5*[jj; jj+1; jj+1; jj];
                kk = find(theta_list_original == posInfo.theta(ii));
                theta_temp(:,ii) = [theta_list(kk); theta_list(kk); theta_list(kk+1); theta_list(kk+1)];
            case 18
                jj = floor(posInfo.phi(ii) / 30);
                phi_temp(:,ii) = 30*[jj; jj+1; jj+1; jj];
                kk = find(theta_list_original == posInfo.theta(ii));
                theta_temp(:,ii) = [theta_list(kk); theta_list(kk); theta_list(kk+1); theta_list(kk+1)];
            case 19
                jj = floor(posInfo.phi(ii) / 90);
                phi_temp(:,ii) = 90*[jj; jj+1; jj+1; jj];
                kk = find(theta_list_original == posInfo.theta(ii));
                theta_temp(:,ii) = [theta_list(kk); theta_list(kk); theta_list(kk+1); theta_list(kk+1)];
        end
    end
    workbar(ii/nn,sprintf("working on %03d / %03d pixels",[ii nn]));
end
%%
[x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree,ones(size(phi_temp)));
figure, patch(x,y,z,drpLib.drpList(1000,:),"EdgeColor","none")
axis equal
