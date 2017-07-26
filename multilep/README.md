# Structure of the heavyNeutrino/multilep module

This directory follows the general CMSSW module layout, see
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBuildFile#CmsswSrcCodeDir

In particular for the heavyNeutrino/multilep module following directories are used:

### data
Contains JSON and .txt configuration files needed to run.
Putting those in the data directory ensures they are shipped with crab.

### interface
Header-files for the classes which are implemented in ./src
In order to structure the code, this module uses sub-analyzers, each of them aimed at calculating and storing the variables for a given object type:
* GenAnalyzer
* JetAnalyzer
* LeptonAnalyzer
* LheAnalyzer
* PhotonAnalyzer
* TriggerAnalyzer

### src
Implementation of the sub-analyzers. Note that for the leptonAnalyzer, the identification and isolation functions are stored in separate .cc files in order to improve readability.

### plugins
Contains mulilep.h and multilep.cc, which form the main plugin of this module. Focusing on keeping track of all the tokens retrieved from the parameters the module is given, and organises the main order of how the sub-analyzers are run.
Note the LheAnalyzer should always be run before a skimming sub-analyzer, and that GenAnalyzer should be run before PhotonAnalyzer.

### python
Contains the jetSequence and egmSequence which group thogether jet and egm-related modules, with the goal of making the test/multilep.py file a bit lighter.
Also contains an updated version of PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties, a temporary workaround for a bug in this file

### test
Directory with the main python config file multilep.py as well as submission scripts
