set(0, 'DefaultLineLineWidth', 2);
figure()
subplot(2,2,1)
tempo = veh_states.Time;
slip_angle = atan2(veh_states.Data(:,1), veh_states.Data(:,2))*180/pi;
plot( tempo, slip_angle);
title('Deslizamento lateral')
xlabel('tempo (s)')
ylabel('\beta (graus)')

subplot(2,2,2)
yaw_rate = veh_states.Data(:,3)*180/pi;
plot( tempo, yaw_rate);
hold on
ref_yaw_rate = ref_states.Data(:,2)*180/pi;
plot( tempo, ref_yaw_rate);
title('Taxa de guinada')
xlabel('tempo (s)')
ylabel('d\psi/dt (graus/s)')

subplot(2,2,3)
roll = veh_states.Data(:,5)*180/pi;
plot( tempo, roll);
hold on
title('Rolagem')
xlabel('tempo (s)')
ylabel('\phi (graus)')

subplot(2,2,4)
plot( tempo, yaw_moment_control_signal.Data(:,1));
hold on
title('Comando')
xlabel('tempo (s)')
ylabel('M_u (Nm)')
