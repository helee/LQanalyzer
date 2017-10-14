sktree -q fastq -a diel -S DoubleEG -s SKTree_HNDiLepSkim -n 20
sktree -q fastq -a diel -S DoubleEG -fake True -s SKTree_HNDiLepSkim -n 20
sktree -a fastq -a diel -S DoubleEG -flip True -s SKTree_HNDiLepSkim -n 20
sktree -q fastq -a diel -list dilep_bkg -s SKTree_HNDiLepSkim -n 20
sktree -q fastq -a diel -list dilep_bkg2 -s SKTree_DiLepSkim -n 20
sktree -q fastq -a diel -list Tch2 -s SKTree_DiLepSkim -n 5 -SIG
#sktree -q fastq -a diel -list Sch2 -s SKTree_HNDiLepSkim -n 5 -SIG
