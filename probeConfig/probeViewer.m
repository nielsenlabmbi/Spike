function probeViewer(probeType)


%load probe configuration
if exist(['probeConfig_' probeType],'file')
    probeW=feval(['probeConfig_' probeType]);

    figure
    plot(probeW(:,2),probeW(:,4),'sqr', 'MarkerSize',20);
    hold on
    for i=1:size(probeW,1)
        text(probeW(i,2)-0.5,probeW(i,4),...
            num2str(probeW(i,1)+1),'FontSize',9);
    end

    ylabel('Depth (0=deep)')
    %pbaspect([1 5 1])
    title(probeType)
else
    disp('Incorrect type selected')
end
