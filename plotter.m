site=3750;
w3=1/2*(U(:,site).^2);

hist
hold off
bar(omega(1:51),height(1:51),1,'FaceColor',[122/255 91/255 117/255]);
hold on
bar(omega(1:51),their*.01/478.573,.5,'FaceColor',[216/255 223/255 230/255]);