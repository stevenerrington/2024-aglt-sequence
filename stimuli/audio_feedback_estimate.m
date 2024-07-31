recObj = audiorecorder;




recDuration = 5;
recordblocking(recObj,recDuration);


y = getaudiodata(recObj);
plot(y);
