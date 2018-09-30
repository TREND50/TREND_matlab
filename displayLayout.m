function [] = DisplayTrigRate(runbeg,runend)
% Displays results from trigrate.mat
% Script used for presentation at 2015 GRAND workshop

SharedGlobals;
DISPLAY = 1;
ants = [101:138 140 148:158];

figure(12)
h1 = scatter(pos_ant(:,1),pos_ant(:,2),'y^','MarkerFaceColor','y');
for i = 1:length(ants)
    text(pos_ant(i,1)+10,pos_ant(i,2)+10,int2str(ants(i)))
end
xlabel('Easting (m)')
ylabel('Northing (m)')
%axis equal
grid on
xlim([-100 2800])
