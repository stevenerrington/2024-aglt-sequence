figuren;
nsubplot(3,2,1,1)
plot(pc_out_auditory.window, pc_out_auditory.obs.pcs(:,1))
xlim([-100 2750])
ylabel('PC1')

nsubplot(3,2,2,1)
plot(pc_out_auditory.window, pc_out_auditory.obs.pcs(:,2))
xlim([-100 2750])
ylabel('PC2')

nsubplot(3,2,3,1)
plot(pc_out_auditory.window, pc_out_auditory.obs.pcs(:,3))
xlim([-100 2750])
xlabel('Time from sequence onset (ms)')
ylabel('PC3')

nsubplot(3,2,1,2)
plot(pc_out_auditory.window, pc_out_frontal.obs.pcs(:,1))
xlim([-100 2750])

nsubplot(3,2,2,2)
plot(pc_out_auditory.window, pc_out_frontal.obs.pcs(:,2))
xlim([-100 2750])

nsubplot(3,2,3,2)
plot(pc_out_auditory.window, pc_out_frontal.obs.pcs(:,3))
xlim([-100 2750])
