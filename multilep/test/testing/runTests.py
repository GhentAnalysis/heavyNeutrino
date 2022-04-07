#!/usr/bin/env python
import subprocess,ROOT,glob,os
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
doCommit = True
with open('tests.log', 'w') as logFile:
  logFile.write(system("eval `scram runtime -sh`;git log -n 1;git diff -- . ':(exclude)*.log'"))

  def runTest(name, testFile):
    print 'runTest', name
    logFile.write('\n--------------------------------------------------------------------------------------------------\n\n')
    command = 'eval `scram runtime -sh`;cmsRun ../multilep.py inputFile=' + testFile + ' outputFile=noskim.root events=100 extraContent=storeLheParticles,storeParticleLevel,storeAllTauID,storeJecSources,storePrefireComponents,storeJetSubstructure'
    logFile.write('Running test: ' + name)
    try:
      system(command)
      system('mv noskim.root ' + name + '.root')
      logFile.write( ' --> OK\n')
      compare(logFile, name)
      system('mv ' + name + '.root ' + name + '-ref.root')
    except subprocess.CalledProcessError, e:
      global doCommit
      doCommit = False
      logFile.write( ' --> FAILED\nOutput:')
      for line in e.output.splitlines():
        if '[arg' in line: continue
        logFile.write('   ' + line + '\n')

  # Tests files to run
  if('lxplus' in os.environ['HOSTNAME']): testFileLocation = '/eos/user/t/lwezenbe/testFiles' 
  else:                                   testFileLocation = '/pnfs/iihe/cms/store/user/lwezenbe/heavyNeutrino/testFiles/'

  runTest('12Nov2019_UL2018',             'file://' + testFileLocation + '/store/data/Run2018C/SingleMuon/MINIAOD/12Nov2019_UL2018-v2/100000/08CD0000-EAC8-844F-96C6-A02E7F742007.root')
  runTest('09Aug2019_UL2017',             'file://' + testFileLocation + '/store/data/Run2017D/SingleElectron/MINIAOD/09Aug2019_UL2017-v1/260000/00A5C633-1806-0844-8D65-C31C779A57F6.root')
  runTest('21Feb2020_UL2016',             'file://' + testFileLocation + '/store/data/Run2016G/DoubleMuon/MINIAOD/21Feb2020_UL2016-v1/230000/0494EFE4-63FD-2448-8E9E-D7C2E5C1E1BE.root')
  runTest('21Feb2020_UL2016_HIPM',        'file://' + testFileLocation + '/store/data/Run2016D/DoubleMuon/MINIAOD/21Feb2020_UL2016_HIPM-v1/230000/0CA8DB06-DE75-4240-8DD1-BF576E8BCFD5.root')
  runTest('RunIISummer20UL18MiniAOD',     'file://' + testFileLocation + '/store/mc/RunIISummer20UL18MiniAOD/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/280000/02D2BCCD-2F6F-0240-96CB-4C9D2525A3C4.root')
  runTest('RunIISummer20UL17MiniAOD',     'file://' + testFileLocation + '/store/mc/RunIISummer20UL17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v1/280000/00E313D4-661C-D54E-86BD-25ABCE7639C5.root')
  runTest('RunIISummer20UL16MiniAOD',     'file://' + testFileLocation + '/store/mc/RunIISummer20UL16MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v13-v2/260000/003F1A76-9CDA-7644-A34E-923C4B1C0E5E.root')
  runTest('RunIISummer20UL16MiniAODAPV',  'file://' + testFileLocation + '/store/mc/RunIISummer20UL16MiniAODAPV/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v8-v1/270000/00974EBB-A8B3-4241-80C3-08C750C5838A.root')

if doCommit:
  system('git add *ref.root')
  system('git add runTests.py')
  system('git add tests.log')
  system('git commit -m"New test run: see tests.log"')
