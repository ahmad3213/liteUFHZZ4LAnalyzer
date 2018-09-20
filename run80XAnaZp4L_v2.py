import os, time

isMC = True

# ____________________________________________________________________________________________________________ ||
mcInfos = [
        [
            'root://cmsio5.rc.ufl.edu//store/user/dsperka/UFHZZAnalysisRun2/MC80X_M17_2l_Feb21/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amc-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2/170221_201244/0000/',
            ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Run6MiniAODv2_%s"%i for i in range(1,2)],
        ],
        ]

# ____________________________________________________________________________________________________________ ||
infos = dataInfos if not isMC else mcInfos
outputDir       = "/raid/raid9/ahmad/Tag_n_Probe/CMSSW_8_0_26_patch1/src/liteUFHZZ4LAnalyzer/Ntuples/"
if not os.path.exists(os.path.abspath(outputDir)):
    os.makedirs(outputDir)
njobs = 1
if njobs > 6: raise RuntimeError, "Too many resources required"
for inputDir,inputSamples in infos:
    for job in range(1,njobs+1):
        for sample in inputSamples:
            if isMC:
                cmd = 'nohup ./ZZ4L_Ana.exe '+inputDir+'/'+sample+' Ntuples/'+sample+' 0 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
            else:
                cmd = 'nohup ./ZZ4L_Ana.exe '+dirZX+'/'+sample+' Ntuples/'+sample+' 1 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &' 
            print cmd
            os.system(cmd)


