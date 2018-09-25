# liteUFHZZ4LAnalyzer

cmsrel CMSSW_8_0_26_patch1 <br/>
cd CMSSW_8_0_26_patch1/src <br/>
git clone https://github.com/ahmad3213/liteUFHZZ4LAnalyzer <br/>
cd liteUFHZZ4LAnalyzer/ <br/>
git checkout Tag_and_probe <br/>
// checking out recipe is finished <br/>

main code in file ZZ4L_Ana.cc <br/>
input tree variable/branches are decelared in file  include/AnalysisTree.h <br/>

// Every time when you change in the code, you have create compile and create executable fine again  by following command 
make <br/>
// to run over NTuple, following is the command  <br/>
python run80XAnaZp4L.py  <br/>
