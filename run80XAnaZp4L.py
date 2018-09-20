import os, time

#dirMC = '/cms/data/store/user/dsperka/UFHZZAnalysisRun2/MC80X_M17_2l_Feb21/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amc-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2/170221_201244/0000/'
dirMC = 'root://cmsio5.rc.ufl.edu//store/user/dsperka/UFHZZAnalysisRun2/MC80X_M17_2l_Feb21/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amc-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2/170221_201244/0000/'

samplesMC  = [
        [
            dirMC,
            ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Run6MiniAODv2_%s"%i for i in range(1,279)],
        ]
]
dirZX = "/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2017/"
samplesZX = [
'SingleDoubleMuon_Run2017-17Nov2017-v1_NoDuplicates'
#'Data_Run2016-03Feb2017_4l'
]
for sample in samplesMC:
    print sample



#njobs = 6
#for job in range(1,njobs+1):

#  if (job>6): continue
  
#  for sample in samplesMC:
#    cmd = 'nohup ./ZZ4L_Ana.exe '+dirMC+'/'+sample+' Ntuples/'+sample+' 0 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
#    print cmd
#    os.system(cmd)

#  for sample in samplesZX:
#    cmd = 'nohup ./ZZ4L_Ana.exe '+dirZX+'/'+sample+' Ntuples/'+sample+' 1 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
#    print cmd
#    os.system(cmd)
