%% Set Up

% Set Up and remove excluded cells

BT145 = BT145quant;
remove145 = find(BT145.Exclude == 1);
BT145(remove145,:) = [];
BT145(1,:) = [];

BT159 = BT159quant;
remove159 = find(BT159.Exclude == 1);
BT159(remove159,:) = [];
BT159(1,:) = [];

BT179 = BT179quant;
remove179 = find(BT179.Exclude == 1);
BT179(remove179,:) = [];
BT179(1,:) = [];

Kura = Kuramochiquant;
removeKura = find(Kura.Exclude == 1);
Kura(removeKura,:) = [];
Kura(1,:) = [];


%% Calculations

data = {BT145, BT159, BT179, Kura};
CellLines = ["BT145";"BT159";"BT179";"Kuramochi"];
Primary = table(CellLines);
Micro = table(CellLines);


for i = 1:height(Primary)
    Line = data{i};
    Primary.CellCount(i) = height(Line);
    Primary.TotalOldPN(i) = sum(Line.PreexistingPN);
    Primary.TotalNewPN(i) = sum(Line.NewPN);
    Primary.TotalPN(i) = Primary.TotalOldPN(i) + Primary.TotalNewPN(i);
    Primary.PNRate(i) = Primary.TotalPN(i)/height(Line);
    NotRuptured = find(Line.PreexistingPN == 0 & Line.NewPN == 0);
    Primary.NumberRupturedPN(i) = height(Line) - length(NotRuptured);
    Primary.PercentCellsPN(i) = Primary.NumberRupturedPN(i)/height(Line);
    Primary.NewPNRate(i) = Primary.TotalNewPN(i)/height(Line);

    Micro.CellCount(i) = height(Line);
    Micro.TotalOldMN(i) = sum(Line.PreexistingMN);
    Micro.TotalNewMN(i) = sum(Line.NewMN);
    Micro.TotalMN(i) = Micro.TotalOldMN(i) + Micro.TotalNewMN(i);
    Micro.MNRate(i) = Micro.TotalMN(i)/height(Line);
    NotRuptured = find(Line.PreexistingMN == 0 & Line.NewMN == 0);
    Micro.NumberRupturedMN(i) = height(Line) - length(NotRuptured);
    Micro.PercentCellsMN(i) = Micro.NumberRupturedMN(i)/height(Line);
    Micro.NewMNRate(i) = Micro.TotalNewMN(i)/height(Line);

    Primary.PercentPrimary(i) = Primary.TotalPN(i)/(Primary.TotalPN(i) + Micro.TotalMN(i));
    Micro.PercentMicro(i) = Micro.TotalMN(i)/(Primary.TotalPN(i) + Micro.TotalMN(i));
end

%% Graphs
x = categorical(Primary.CellLines);

figure()
hold on
PrimaryRates = [Primary.PNRate Primary.NewPNRate Primary.PercentCellsPN];
bar(x,PrimaryRates)
legend("Frequency of PN Rupture", "Frequency of PN Rupture After t_0","Frequency of Cells with PN Rupture",'FontSize',8,Location="northwest")
xlabel("Cell Line")
ylabel("Frequency")
title("Primary Nuclear Ruptures")

figure()
hold on
MicroRates = [Micro.MNRate Micro.NewMNRate Micro.PercentCellsMN];
bar(x,MicroRates)
legend("Frequency of MN Rupture", "Frequency of MN Rupture After t_0","Frequency of Cells with MN Rupture",'FontSize',8,Location="northwest")
xlabel("Cell Line")
ylabel("Frequency")
title("Micronucleus Ruptures")

figure
hold on
RuptureTypes = [Primary.PercentPrimary Micro.PercentMicro];
bar(x,RuptureTypes,'stacked')
legend('Primary Rupture','Micronucleus Rupture')
title('Rupture Types')
xlabel('Cell Line')
ylabel('Proportion')

