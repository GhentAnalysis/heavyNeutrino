from WMCore.Configuration import Configuration
import os

# Get parameters from submit script
productionLabel = os.environ['CRAB_PRODUCTIONLABEL']
dataset         = os.environ['CRAB_DATASET']
outputFile      = os.environ['CRAB_OUTPUTFILE']

requestName = dataset.split('/')[2] + '_' + productionLabel
requestName = requestName.replace('RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6','Moriond2017')
requestName = requestName.replace('RunIISummer16MiniAODv2-PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6','Moriond2017')
requestName = requestName.replace('RunIISummer16MiniAODv2-80X_mcRun2_asymptotic_2016_TrancheIV_v6', '80X')


# Crab configuration
config = Configuration()
config.section_('General')
config.General.transferLogs            = True
config.General.requestName             = requestName
config.General.workArea                = os.path.join('crab', productionLabel, dataset.split('/')[1])

config.section_('JobType')
config.JobType.psetName                = 'multilep.py'
config.JobType.pyCfgParams             = ['events=-1', 'outputFile='+outputFile, 'isData='+('False' if 'SIM' in dataset else 'True')]
config.JobType.pluginName              = 'analysis'
config.JobType.outputFiles             = [outputFile]
config.JobType.sendExternalFolder      = True
config.JobType.allowUndistributedCMSSW = True 

config.section_('Data')
config.Data.inputDataset               = dataset
config.Data.unitsPerJob                = 300 if 'SIM' in dataset else 50
config.Data.splitting                  = 'LumiBased'
config.Data.outLFNDirBase              = '/store/user/' + os.environ['USER'] + '/heavyNeutrino/'
config.Data.publication                = False

config.section_('Site')
config.Site.storageSite                = 'T2_BE_IIHE'
