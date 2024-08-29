function OD_Signal(tOut,in,out,sineA,sinef,ind_start,ind_end)


figure()
plot(tOut,in,'k-','LineWidth',2);
hold all;
plot(tOut,out,'b-','LineWidth',2);
hold off
grid minor;
ylabel('[units]')  
xlabel('time [s]')
title('Signal')
legend('input ','output','Location','northwest')
set(gcf,'Position',[151 246 1024 420])

% calculate and plot Modulus and Phase
sig=CalculateModPhase(tOut,in,out,sineA,sinef,ind_start,ind_end);

%calculate and plot WSD

[WSD_mag,WSD_phase,sig]=Calculate_WSD(sinef,sig);

% This is added in GitHib
% Add this 2nd line in RGH
end
