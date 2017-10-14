sktree -q fastq -a dimu -S DoubleMuon -s SKTree_HNDiLepSkim -n 20
sktree -q fastq -a dimu -S DoubleMuon -fake True -s SKTree_HNDiLepSkim -n 20
sktree -q fastq -a dimu -list dilep_bkg -s SKTree_HNDiLepSkim -n 20
sktree -q fastq -a dimu -list dilep_bkg2 -s SKTree_DiLepSkim -n 20
sktree -q fastq -a dimu -list Tch1 -s SKTree_DiLepSkim -n 5 -SIG
#sktree -q fastq -a dimu -list Sch1 -s SKTree_DiLepSkim -n 5 -SIG
