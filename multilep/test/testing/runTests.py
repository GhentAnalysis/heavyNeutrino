#!/usr/bin/env python
import subprocess,ROOT,glob
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Ignore these branches (are based on random numbers)
ignoreBranches = ['_phRandomConeChargedIsolation']

# System command and retrieval of its output
def system(command):
  return subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)

# Saves the mean and RMS from the histos and branches in a dictionary
def extractFile(name):
  f = ROOT.TFile(name)
  f.cd('blackJackAndHookers')
  data = {}
  for key in ROOT.gDirectory.GetListOfKeys():
    object = key.ReadObj()
    try:
      for b in object.GetListOfBranches():
        if b.GetName() in ignoreBranches: continue
        object.Draw(b.GetName() + ' >> temp'+b.GetName())
        temp = ROOT.gDirectory.Get('temp'+b.GetName())
        data[b.GetName()] = (temp.GetMean(), temp.GetRMS())
    except:
      data[key.GetName()] = (object.GetMean(), object.GetRMS())
  return data

# Comparing the two ROOT files using the dictionaries
def compare(logger, name):
  newData = extractFile(name + '.root')
  refData = extractFile(name + '-ref.root')

  def compatible(tuple1, tuple2):
    return all([(abs(x-y)/x < 0.0001 if abs(x) > 0 else (abs(x-y)<0.0001)) for x,y in zip(tuple1, tuple2)])

  new     = sorted([i for i in newData if i not in refData])
  removed = sorted([i for i in refData if i not in newData])
  changed = sorted([i for i in refData if i in newData and not compatible(refData[i], newData[i])])
  if len(new):      logger.write('   New: '     + ','.join(new) + '\n')
  if len(removed):  logger.write('   Removed: ' + ','.join(removed) + '\n')
  if len(changed):  logger.write('   Changed:\n')
  for c in changed:
    logger.write('      %-50s mean: %-25s rms: %-25s\n' % (c, '%8.4f --> %8.4f' % (refData[c][0], newData[c][0]), '%8.4f --> %8.4f' % (refData[c][1], newData[c][1])))

# Compile
print system('eval `scram runtime -sh`;cd $CMSSW_BASE/src;scram b -j 10')

# Starting the test
with open('tests.log', 'w') as logFile:
  logFile.write(system("eval `scram runtime -sh`;git log -n 1;git diff -- . ':(exclude)*.log'"))

  def runTest(name, testFile):
    logFile.write('\n--------------------------------------------------------------------------------------------------\n\n')
    command = 'eval `scram runtime -sh`;cmsRun ../multilep.py inputFile=' + testFile + ' outputFile=noskim.root events=100 extraContent=storeLheParticles'
    logFile.write('Running test: ' + name)
    try:    
      system(command)
      system('mv noskim.root ' + name + '.root')
      logFile.write( ' --> OK\n')
      compare(logFile, name)
      system('mv ' + name + '.root ' + name + '-ref.root')
    except subprocess.CalledProcessError, e:
      logFile.write( ' --> FAILED\nOutput:')
      for line in e.output.splitlines():
        if '[arg' in line: continue
        logFile.write('   ' + line + '\n')

  # Tests files to run
  runTest('Run2018-17Sep2018',    'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/100000/42EFAC9D-DC91-DB47-B931-B6B816C60C21.root')
  runTest('Run2018-PromptReco',   'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v3/000/316/569/00000/0085320B-4E64-E811-A2D3-FA163E2A55D6.root')
  runTest('Run2017-31March2017',  'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2017C/DoubleEG/MINIAOD/31Mar2018-v1/00002/58841328-4838-E811-AC6D-00266CFFCAC0.root')
  runTest('Run2016-17Jul2018',    'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2016H/DoubleEG/MINIAOD/17Jul2018-v1/80000/B8EE4064-968D-E811-B077-0242AC1C0500.root')
  runTest('Autumn18MiniAOD',      'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIIAutumn18MiniAOD/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v2/110000/00707922-8E6F-3042-A709-2DD4DB9AEDED.root')
  runTest('Fall17MiniAODv2',      'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017RECOPF_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0A1754A2-256F-E811-AD07-6CC2173CAAE0.root')
  runTest('Summer16MiniAODv3',    'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-105To160_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/00000/2E242480-5C0D-E911-B9A6-90E2BACBAA90.root')

system('git add *ref.root')
system('git add runTests.py')
system('git add tests.log')
system('git commit -m"New test run: see tests.log"')
