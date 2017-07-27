function visDTW( D, Path )

figure;
imagesc(D);
axis xy;
colorbar;
colormap hot;
if nargin > 1
    hold on;
    h1 = line(Path(2,:),Path(1,:),'Color','w','LineWidth',4);
    h2 = line(Path(2,:),Path(1,:),'Color','k','LineWidth',2);
    hold off;
end

end

