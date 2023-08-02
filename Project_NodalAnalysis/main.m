clear
close all
clc

Parameter;

%Temperature along the discretized length of tubing
T_grad = (T2 - T1)/8000;

P_pc = psc_press(Y_g);
T_pc = psc_tempt(Y_g);

% discretizing length
verSection = 0:tub_length:KOP;
MD_c = 1500*pi;                             %Total deviated MD
devSection = KOP:L:(KOP+MD_c);
horSection = (KOP+MD_c):L:(KOP+MD_c+L_h);

entSection = [verSection devSection(:,2:end)]; %Measured depth for entire section
ent_Section = [verSection devSection(:,2:end) horSection];
theta_diff = -0.4775;                             % build rate for deviated section calculated using length of an arc
Rad = 3000;


for i=1:length(ent_Section)
    if i <= 200
        theta(i) = 90;
        vDev(i) =  ent_Section(i) + Rad*sind(theta(i) + 90);
        hDev(i) = ent_Section(i)*cosd(theta(i));

    elseif (i > 200) && (i < 389)
        theta(i) = theta(i-1) + theta_diff;
        vDev(i) = KOP + Rad*(sind(theta(i) + 90));       % actutal TVD depth for entire well section
        hDev(i) = Rad - Rad*(cosd(theta(i) - 90));
    else
        theta(i) = 0;
        vDev(i) = 8000;
        hDev(i) = hDev(i-1) + L;
    end
end



counter = 1;
MAE = 1;
q_g = 1000; %CF/D
z = 1;

    while MAE > 0.00001
        % Project 1
        for k = 2:length(vDev(:,2:388))

            nodeTem(1) = T1;
            nodeP(1) = P1;
            nodeTem(k) = T_grad*(vDev(k) - vDev(k-1)) + nodeTem(k-1);
            nodeTpr(k) = psr_tempt(nodeTem(k), T_pc);
            blocTem(k-1) = (nodeTem(k) - nodeTem(k-1))/log(nodeTem(k)/nodeTem(k-1));
            blocTpr(k-1) = psr_tempt(blocTem(k-1), T_pc);

            nodePpr_(k) = psr_press(nodeP(k-1), P_pc);


            error = 1;
            while error > 0.01
                nrho(k) = (0.27*nodePpr_(k))/(z*nodeTpr(k));
                nZ(k) = 1 + nrho(k)*(A1 + A2/nodeTpr(k) + A3/nodeTpr(k)^3 + A4/nodeTpr(k)^4 + A5/nodeTpr(k)^5) + nrho(k)^2*(A6 + A7/nodeTpr(k) + A8/nodeTpr(k)^2) - A9*nrho(k)^5*(A7/nodeTpr(k) + A8/nodeTpr(k)^2) + A10*(1 + A11*nrho(k)^2)*(nrho(k)^2/nodeTpr(k)^3)*exp(-A11*nrho(k)^2);
                rho_g(k) = (0.27*nodePpr_(k))/(nZ(k)*nodeTpr(k));
                error = (abs((rho_g(k) - nrho(k))/rho_g(k)))*0.5;
                z = (nZ(k) + z)/2;
            end


            nodeBg(k) = gas_fom_vol(nZ(k), nodeTem(k), nodeP(k-1));

            nodeA(k) = m_gA(M_g, nodeTem(k));
            nodeB(k) = m_gB(M_g, nodeTem(k));
            nodeC(k) = m_gC(nodeB(k));
            nodem_g(k) = visz_gas(nodeA(k), nodeB(k), nodeC(k), rho_g(k));

            nodeN_Re(k) = Rey_num_hor(Y_g, q_g, d_tub, nodem_g(k));


            nodef_fact(k) = friction_fact(nodeN_Re(k), rel_rough);  % turbulent condition is assumed to be always true for high pressrue gas flo

            nodes_P(k) = s_P2(Y_g, tub_length, theta(k), nZ(k), nodeTem(k));
            nodeP(k) = P2_cal_ver(nodeP(k-1), nodes_P(k), nodef_fact(k), nZ(k), nodeTem(k), q_g, theta(k), d_tub); %psi

            blocPre(k) = (nodeP(k) + nodeP(k-1))/2;

            blocPpr(k-1) = psr_press(blocPre(k), P_pc);
            error1 = 1;
            while error1 > 0.01
                brho(k-1) = (0.27*blocPpr(k-1))/(z*blocTpr(k-1));
                bZ(k-1) = 1 + brho(k-1)*(A1 + A2/blocTpr(k-1) + A3/blocTpr(k-1)^3 + A4/blocTpr(k-1)^4 + A5/blocTpr(k-1)^5) + brho(k-1)^2*(A6 + A7/blocTpr(k-1) + A8/blocTpr(k-1)^2) - A9*brho(k-1)^5*(A7/blocTpr(k-1) + A8/blocTpr(k-1)^2) + A10*(1 + A11*brho(k-1)^2)*(brho(k-1)^2/blocTpr(k-1)^3)*exp(-A11*brho(k-1)^2);
                brho_g(k-1) = (0.27*blocPpr(k-1))/(bZ(k-1)*blocTpr(k-1));
                error1 = (abs((brho_g(k-1) - brho(k-1))/brho_g(k-1)))*0.5;
                z = (bZ(k-1) + z)/2;
            end

            blocBg(k-1) = gas_fom_vol(bZ(k-1), blocTem(k-1), blocPre(k));
            blocm_g(k-1) = visz_gas(nodeA(k-1), nodeB(k-1), nodeC(k-1), brho_g(k-1));

        end

        % Project 2

        for n = 2:length(horSection)
            nodeP_hor(1) = nodeP(k);
            nodeTem_hor(n) = T2;
            nodeTpr_hor(n) = psr_tempt(nodeTem_hor(n), T_pc);


            nodePpr_hor(n) = psr_press(nodeP_hor(n-1), P_pc);


            nZ_hor(n) = Z_factor(nodeTpr_hor(n), nodePpr_hor(n));

            rho_g_hor(n) = den_g(nodeP_hor(n-1), M_g, nZ_hor(n), R, nodeTem_hor(n));

            nodeBg_hor(n) = gas_fom_vol(nZ_hor(n), nodeTem_hor(n), nodeP_hor(n-1));

            nodeA_hor(n) = m_gA(M_g, nodeTem_hor(n));
            nodeB_hor(n) = m_gB(M_g, nodeTem_hor(n));
            nodeC_hor(n) = m_gC(nodeB_hor(n));

            nodem_g_hor(n) = visz_gas(nodeA_hor(n), nodeB_hor(n), nodeC_hor(n), rho_g_hor(n));

            nodeN_Re_hor(n) = Rey_num_hor(Y_g, q_g, d_tub, nodem_g_hor(n));
            nodef_fact_hor(n) = friction_fact(nodeN_Re_hor(n), rel_rough);

            nodeP_hor(n) = P2_cal_hor(nodeP_hor(n-1), Y_g, nodef_fact(n), nZ_hor(n), nodeTem_hor(n), q_g, d_tub, L);


            P_avg = 3600:-100:2000;

            Q_vlp_pss(1) = q_g;
            Q_vlp(1) = q_g;
%             Q_vlp_tran = q_g;
            Q_ipr_pss(1) = 0;
            c_g(1) = 1/(P_avg(1));
            t_pss = (1200*poro*m_g1*c_g*r_e^2/k_v/24);
            initial_g_in_place = 43560*A*poro*h*(1-Sw)/nodeBg_hor(n);
            time_transient = linspace(1,t_pss,50);
            G_p_pss(1) = 0;
            for j =1:length(P_avg)
%                         for u = 1:length(time_transient)
%                 
%                             Q_ipr_tran(n,u,j) = q_transient(k_v, h, P_avg(j), nodeP_hor(n), poro, nodeTem_hor(n), time_transient(u), c_g, r_w, nodem_g_hor(n));
%                             Q_vlp_tran(n,u,j) = Q_vlp_tran(n-1,u) + Q_ipr_tran(n,u,j);
%                             Q_ipr_hor_tran(u,j) = mean(Q_ipr_tran(n,u,j));
%                             Q_vlp_hor_tran(u,j) = Q_vlp_tran(n,j) - sum(Q_ipr_tran(n,u,j));
%                             G_p_tran(u,j) = initial_g_in_place*(1 - ((P_avg/Z_factor(psr_tempt(T2, T_pc),psr_press(P_avg, P_pc))) / (3600/Z_factor(psr_tempt(T2, T_pc),psr_press(3600, P_pc)))));
%                 
%                         end

                Q_ipr_pss(n,j) = q_cal(nodeP_hor(n), P_avg(j), 30, nZ_hor(n), nodeTem_hor(n), nodem_g_hor(n), r_e, r_w, k_v);
                Q_vlp_pss(n,j) = (Q_vlp_pss(n-1)) + Q_ipr_pss(n,j);

                Q_vlp(n,j) = (Q_vlp(n-1) - Q_ipr_pss(n,j));

                Q_ipr_pss_mean(j) = mean(Q_ipr_pss(n,j));

                G_p_pss(j) = initial_g_in_place*(1 - ((P_avg(j)/Z_factor(psr_tempt(T2, T_pc),psr_press(P_avg(j), P_pc))) / (3600/Z_factor(psr_tempt(T2, T_pc),psr_press(3600, P_pc)))));
                time_pss(j) = (G_p_pss(j)-G_p_pss(1))./Q_ipr_pss_mean(j);

            end


            MAE = (abs((Q_vlp_pss(1) - Q_ipr_pss(n,j))/Q_vlp_pss(1)))*0.5;
            [nodeP_hor(n), nodeP(k), Q_ipr_pss(n,j), q_g]
%             [G_p_pss(j), time_pss(j)]
            q_g = (q_g + Q_ipr_pss(n,j))/2;
            counter = counter + 1;
        end
    end


pressure = [nodeP nodeP_hor];
plot(pressure, 'b-o')
xlabel("Number of discrete elements")
ylabel("Pressure")
title("Well Pressure Profile")

figure
plot(nodeP_hor, 'k-o')
xlabel("Number of discrete elements in Hor section")
ylabel("Pressure")
title("Wellbore Pressure Profile")

figure 
Z= [nZ nZ_hor(1,2:end)];
plot(Z(:,2:end), 'g-o')
xlabel("Number of discrete elements")
ylabel("Z Factor")
title("Z factor profile")

figure %well trajectory
plot (hDev, vDev, 'r-o');
xlabel("Horizontal deviation along Hole ft")
ylabel("Vertical deviation along hole ft")
title("Well Trajectory")


figure %production rate curve versus cummulative production curve
plot(time_pss, Q_ipr_pss_mean, 'b-o')
hold on
plot(time_pss, G_p_pss./10000000, 'r-o')
xlabel("Time (hrs)")
ylabel("Rates")
title("Rate Vs Prod Curve")
hold off

figure %production rate curve
plot(time_pss, Q_ipr_pss_mean, 'b-o')
xlabel("Time (days)")
ylabel("Mean IPR Flowrate (Mcf/d)")
title("IPR rates @ reservoir pressures")


figure %cummulative production curve
plot(time_pss, G_p_pss, 'r-o')
xlabel("Time (hrs)")
ylabel("Cumm Production (Bcf/d)")
title("Cumm Prod for depleting reservoir pressures")


figure %cummulative rate for reservoir pressures
plot(horSection, Q_vlp_pss(:,1), 'r-o')
hold on
plot(horSection, Q_vlp_pss(:,17), 'b-o')
hold off
xlabel("Horiz Deviation")
ylabel("cummulative VLP&IPR rates across section")
title("Cumm VLP&IPR rate curves")


figure % ipr curve at initial reservoir pressure
plot(horSection, Q_ipr_pss, 'r-o')
xlabel("Horizontal Distance")
ylabel("Resservoir pressure drop & IPR rates")
title("IPR rate across horizontal pressures")


figure 
plot(horSection, Q_vlp(:,1), 'r-o')
hold on
plot(horSection, Q_vlp(:,17), 'b-o')
hold off
xlabel("Horiz Deviation")
ylabel("VLP rates across section")
title("VLP curves")

figure %vlp curve at various reservori pressures
plot(Q_vlp(:,1), nodeP_hor, 'r-o')
hold on
plot(Q_vlp(:,17), nodeP_hor, 'b-o')
hold off
xlabel("Rate (SCF/D)")
ylabel("Pressure (psi)")
title("VLP curves")