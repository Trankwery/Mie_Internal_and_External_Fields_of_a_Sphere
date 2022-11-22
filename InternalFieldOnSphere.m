function InternalFieldOnSphere
%%  This function calculates evolution of electrical field in and around the sphere.
%      aa    particle radius in nanometers

%      n_out refractive index of medium
% l     number of multipoles in sum
% t     time

%% Forumalas were taken from the Bohren & Huffman "Absorption and Scattering of Light by Small Particles"
% Code written by G. Derkachov & A. Derkachova

% Initial conditions
    n_in = 1.45;                     % droplets refractive index
    n_out = 1;                        % ivironment refractive index
    m = n_in./n_out;              % effective refractive index
    aa = 38780:-0.1:38760  ;  %2815:-1:2e3; %[nm]  % vector of droplet radiuses
    lambda = 805;                 %[nm]   % wavelenght of insident beam
    E_0 =  1;                         % amplitude of insident field
    Step = 1;
    Theta = 0: Step: 180;   % Range of angles changes
    tt = Theta*pi/180;          % Range of angles in radians
    Fi = 0;                            % Range of angle changes (elevation angle)

    S = 'Time estimation...';
    wb = waitbar(0,S);
    set(wb,'position',[447.75,243.75,270.0,56.25]);
    hf = figure('position',[509,448,560,420]);
    ha = axes;
    SEi = zeros(1,length(aa)); % Total energy inside droplet
    for ia = 1:length(aa)
        % Time estimation
            t_tic = tic;
            waitbar(ia./length(aa),wb,S);

            q = 2*pi*aa(ia)*n_out./lambda;       %out of the sphere size parameter

        % Mie scattering coefficients
        [a, b] = MieScatKoeff_Boren_Hufman(q, m);
        [c, d] = MieScatKoeff_Boren_Hufman_cd(q, m);
        nmax = length(a);                                % The number of harmonics

        rr = linspace(0,aa(ia)*1.5,3e2);           % Scanning range linspace(0,aa(ia)*2,3e2)
        EE = zeros(length(rr),length(Theta));  % Initialization of electric field matrix
        wb2 = waitbar(0,'Scanning progress...');
        set(wb2,'position',[449.250,141.0,270.0,56.25]);

        for kkr = length(rr):-1:1 % Scanning loop
            waitbar((1-length(rr)\kkr),wb2)
            % initialisation of matrices
            Er = zeros(length(Theta),length(Fi) );
            Er_inc = Er;
            Ef_inc = Er;
            Et_inc = Er;

            Er_scat = Er;
            Ef_scat = Er;
            Et_scat = Er;

            Er_in = Er;
            Ef_in = Er;
            Et_in = Er;


            k2 = 2*pi*rr(kkr)*n_out./(lambda);
            q_in = 2*pi*rr(kkr)*n_in/lambda;

            ll = 1:nmax;
            El = E_0.*(1i.^ll).*(2*ll+1)./(ll.*(ll+1));

            for ii = 1:length(Theta)
                [Tau, Pi] = MLegandr(tt(ii), nmax); 

                if rr(kkr) >= aa(ia) %------------ The electric field outside the droplet  ----------------------%

                    %------------ Incident field ----------------------%
                    Er_inc(ii,:) = sum(-1i.*El.*(ll.^2+ll).*cosd(Fi).*sind(Theta(ii)).*Pi(ll).'.*fpsi(ll,k2)./(k2.^2));
                    Et_inc(ii,:) = sum( El.*cosd(Fi).*( Pi(ll).'.*fpsi(ll,k2) - 1i.*Tau(ll).'.*dpsi(ll,k2) )./k2 );
                    Ef_inc(ii,:) = sum(-El.*sind(Fi).*( Tau(ll).'.*fpsi(ll,k2) - 1i.*Pi(ll).'.*dpsi(ll,k2) )./k2 );
                    %------------ scattered field ----------------------%
                    Er_scat(ii,:)  = sum(1i.*El.*(ll.^2 + ll).*cosd(Fi).*sind(Theta(ii)).*a(ll).*Pi(ll).'.*fksi(ll,k2)./(k2.^2));
                    Et_scat(ii,:)  = sum( El.*cosd(Fi) .* ( 1i.*a(ll).*Tau(ll).'.*dksi(ll,k2) - b(ll).*Pi(ll).'.*fksi(ll,k2) )./k2 );
                    Ef_scat(ii,:)  = sum( -El.*sind(Fi).*(1i.*a(ll).*Pi(ll).'.*dksi(ll,k2) - b(ll).*Tau(ll).'.*fksi(ll,k2))./k2);
                    
                else  %------------ internal field ----------------------%

                    Er_in(ii,:) = sum(-1i.*El.*(ll.^2+ll).*cosd(Fi).*sind(Theta(ii)).*d(ll).*Pi(ll).'.*fpsi(ll,q_in)./(q_in.^2));
                    Et_in(ii,:) = sum( El.*cosd(Fi).*(c(ll).*Pi(ll).'.*fpsi(ll,q_in) - 1i.*d(ll).*Tau(ll).'.*dpsi(ll,q_in))./q_in);
                    Ef_in(ii,:) = sum( -El.*sind(Fi).*(c(ll).*Tau(ll).'.*fpsi(ll,q_in) - 1i.*d(ll).*Pi(ll).'.*dpsi(ll,q_in))./q_in);
                end
                
                if rr(kkr) >= aa(ia)
                     EE(kkr,ii) =  norm([real( Er_inc(ii) + Er_scat(ii) ),real(Et_inc(ii,:)+Et_scat(ii,:)),(Ef_inc(ii,:)+Ef_scat(ii,:)) ]);
                else
                     EE(kkr,ii) =   norm( [real(Er_in(ii)),real(Et_in(ii)),real(Ef_in(ii))]);
                end

            end;

        end;
        close(wb2);
        SEi(ia) = sum(sqrt(real(Er_in(:)).^2 + real(Et_in(:)).^2 + real(Ef_in(:)).^2));
        
        XX = zeros(length(rr),length(tt));
        YY = XX;
        for ii = 1:length(rr)
            [XX(ii,:),YY(ii,:)] = pol2cart(tt,rr(ii));
        end
        % Drawing and saving data

            Et = zeros(size(EE,1),size(EE,2)*2);
            X = Et;
            Y = Et;

            theta_id = 1:size(EE,2);
            Et(:,theta_id)= real(EE);
            Et(:,theta_id+size(EE,2)) = real(EE(:,end:-1:1));
            X(:,theta_id)= XX;
            X(:,theta_id+size(EE,2)) = XX(:,end:-1:1);
            Y(:,theta_id)= YY;
            Y(:,theta_id+size(EE,2)) = -YY(:,1:end);

            surf(ha,X,Y,Et,'linestyle','none');
           %     colormap(cma);
           set(gca,'fontsize',14);
           xlabel('X[ nm ]'); ylabel('Y[ nm ]');zlabel('E[ arb.u. ]');
           title({['a = ',num2str(aa(ia),'%7.2f'),'[ nm ]; \lambda = 805[ nm ]'], 'm = 1.45; \phi = 0.'})
           axis('equal');view([0,90]);
           drawnow;

    % Saving the result of calculations 
        print(hf,['Images\',num2str(ia,'%04i'),'.png'],'-dpng','-zbuffer','-r600');

        t_toc = toc(t_tic);
        S = sprintf('Calculation time -> %s \n Time left - > %s',...
            datestr(datenum([0,0,0,0,0,t_toc*length(aa)]),31),datestr(datenum([0,0,0,0,0,t_toc*(length(aa)-ia)]),31));
    end
    close(wb)
