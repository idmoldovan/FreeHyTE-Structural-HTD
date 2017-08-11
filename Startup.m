function Startup

%%
% Splash screen sequence
timer = tic;
ShowSplashForSecs = 1;

s = SplashScreen( 'Splashscreen', 'Splash.png', ...
    'ProgressBar', 'on', ...
    'ProgressPosition', 5, ...
    'ProgressRatio', 0.4 );
s.addText( 500, 80, 'Structural HTD', 'FontSize', 35, 'Color', [0 0 0.6] );
s.addText( 500, 120, 'v1.2.1', 'FontSize', 25, 'Color', [0.2 0.2 0.5] );
s.addText( 365, 270, 'Loading...', 'FontSize', 20, 'Color', 'white' );

ElapsedTime = toc(timer);
if ElapsedTime < ShowSplashForSecs
    pause(ShowSplashForSecs - ElapsedTime);
end
delete(s);

fclose('all');

StructDef(2);