layout(matrix(1:2,ncol=1))
while(T){
    data=read.table("AutoSet_12-10-17_WTLiveFedRoti_148_007_tracks_1.csv",header=T);
    plot(data$frameN[(nrow(data)-600):nrow(data)],(data$EyeLDeg-data$EyeRDeg)[(nrow(data)-600):nrow(data)],type='l',ylim=c(-15,80));
    lines(data$frameN[(nrow(data)-600):nrow(data)],data$EyeLDeg[(nrow(data)-600):nrow(data)],col='green');
    lines(data$frameN[(nrow(data)-600):nrow(data)],data$EyeRDeg[(nrow(data)-600):nrow(data)],col='red');
    plot(data$frameN[(nrow(data)-600):nrow(data)],(data$AngleDeg)[(nrow(data)-600):nrow(data)],type='l',ylim=c(0,360));
}
