from WMCore.Configuration import Configuration
import os

# Get parameters from submit script
productionLabel = os.environ['CRAB_PRODUCTIONLABEL']
dataset         = os.environ['CRAB_DATASET']
outputFile      = os.environ['CRAB_OUTPUTFILE']
lumiMask        = os.environ['CRAB_LUMIMASK']
extraContent    = os.environ['CRAB_EXTRACONTENT']

requestName = dataset.split('/')[2] + '_' + productionLabel
requestName = requestName.replace('RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3', 'MiniAOD2016v3')
requestName = requestName.replace('RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14', 'MiniAOD2017v2')
requestName = requestName.replace('RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14', 'MiniAOD2017v2NewPMX')
requestName = requestName.replace('RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15', 'MiniAOD2018')

# Crab configuration
config = Configuration()
config.section_('General')
config.General.transferLogs            = True
config.General.requestName             = requestName
config.General.workArea                = os.path.join('crab', productionLabel, dataset.split('/')[1])

config.section_('JobType')
config.JobType.psetName                = 'multilep.py'
config.JobType.pyCfgParams             = ['events=-1', 'outputFile='+outputFile, 'inputFile='+dataset] + (['extraContent='+extraContent] if extraContent else [])
config.JobType.pluginName              = 'analysis'
config.JobType.outputFiles             = [outputFile]
config.JobType.sendExternalFolder      = True
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
config.Data.inputDataset               = dataset
config.Data.unitsPerJob                = 1 if 'SIM' in dataset else 40
config.Data.splitting                  = 'FileBased' if 'SIM' in dataset else 'LumiBased'
config.Data.outLFNDirBase              = '/store/user/' + os.environ['USER'] + '/heavyNeutrino/'
config.Data.publication                = False
config.Data.lumiMask                   = lumiMask if not 'SIM' in dataset else None
config.Data.allowNonValidInputDataset  = True #allow unfinished samples (production) to be processed

config.section_('Site')
config.Site.storageSite                = 'T2_BE_IIHE'
