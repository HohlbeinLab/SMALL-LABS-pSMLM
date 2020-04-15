function drawGUIloadbar(perc,textinput,GUIhandles)
cla(GUIhandles.axes1);
hold on
fill(GUIhandles.axes1,[0 1 1 0],[0 0 1 1],'w');
fill(GUIhandles.axes1,[0 perc perc 0],[0 0 1 1],[63 202 142]./255);
set(GUIhandles.axes1,'XLim',[0 1]);
set(GUIhandles.axes1,'YLim',[0 4]);
set(GUIhandles.axes1,'Visible','off');
text(GUIhandles.axes1,0.5,1.4,textinput,'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',14)
drawnow
end