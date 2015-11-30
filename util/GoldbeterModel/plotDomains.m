% Plot vertical lines separating the different domains of stability

clear;
figure(123);
hold on
plot([0.774,0.774],[0,1],'k');
plot([0.289,0.289],[0,1],'k');
plot([0.189,0.189],[0,1],'--k');
plot([0.182,0.182],[0,1],'--k');
plot([0.168,0.168],[0,1],'--k');

hold off
xlim([0 1])
xlabel('$\beta$','Interpreter','LaTex')
set(gca,'ytick',[])
set(gca,'yticklabel',[])