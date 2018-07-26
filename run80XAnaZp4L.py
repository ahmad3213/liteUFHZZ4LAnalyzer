import os, time

#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Feb2018/'
dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2017/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2016/'
#dirMC = '/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2016/official_signal_sample/'
samplesMC  = [
###2016######
#"ZpTomumu_M1_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M5_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M10_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M15_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M20_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M25_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M30_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M35_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M40_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M45_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M50_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M55_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M60_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M65_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M70_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M75_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M80_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M85_13TeV_MadGraph5_pythia8-v5",
#"ZpTomumu_M90_13TeV_MadGraph5_pythia8-v5",
#"ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2",
#"ZZTo4L_13TeV_powheg_pythia8_ext1_RunIISummer16MiniAODv2"
#"GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2"
#"ZZprimeTo4mu_M-45_13TeV-madgraph_RunIISummer16MiniAODv2",
#"ZZprimeTo4mu_M-15_13TeV-madgraph_RunIISummer16MiniAODv2"
###2017###
"ZpTomumu_M10_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M15_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M20_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M25_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M30_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M35_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M40_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M45_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M50_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M5_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M55_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M60_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M65_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZpTomumu_M70_13TeV_MadGraph5_pythia8-v5_muahmad-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8",
"ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2",
"ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1"
]
dirZX = "/raid/raid9/ahmad/RUN2_Analyzer/v2/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples_Input/2017/"
samplesZX = [
'SingleDoubleMuon_Run2017-17Nov2017-v1_NoDuplicates'
#'Data_Run2016-03Feb2017_4l'
]

njobs = 6
for job in range(1,njobs+1):

  if (job>6): continue
  
  for sample in samplesMC:
    cmd = 'nohup ./ZZ4L_Ana.exe '+dirMC+'/'+sample+' Ntuples/'+sample+' 0 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
    print cmd
    os.system(cmd)

#  for sample in samplesZX:
#    cmd = 'nohup ./ZZ4L_Ana.exe '+dirZX+'/'+sample+' Ntuples/'+sample+' 1 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
#    print cmd
#    os.system(cmd)

#'Data_Run2016-03Feb2017_hzz4l'
#'Data_Run2016-03Feb2017_hzz4l'
