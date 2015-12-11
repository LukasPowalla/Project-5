fid=fopen('ergodic100.txt')
    s = textscan(fid,'%f %f')
fclose(fid)

x=s{1};
y=s{2};

y = y/10000;
plot(y,x,'red')
xlabel('t/t_{crunch}')
ylabel('< V > / < k >','Rotation',0)
%title('The system potential energy over the average kintic energy ')