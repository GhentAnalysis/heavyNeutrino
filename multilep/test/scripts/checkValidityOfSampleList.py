import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--sampleList',   action='store', default=False,  help='Input Sample List')

args = argParser.parse_args()

cms_das_names = [line.split('%')[0].strip() for line in open(args.sampleList, 'r')]
cms_das_names = [line.split('+')[0].strip() for line in cms_das_names if line]
cms_das_names = [line.split(':')[-1] for line in cms_das_names if line]

with open('sampleListCheck_{}'.format(args.sampleList.rsplit('/', 1)[-1]),"w") as f:
    f.write('INVALID DAS ENTRIES: \n \n')

    for cdn in cms_das_names:
        import subprocess, shlex
        command_line = '/cvmfs/cms.cern.ch/common/dasgoclient --query "dataset={0}"'.format(cdn)
        args = shlex.split(command_line)
        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        (out, err) = p.communicate()
        out_message = [x for x in out.decode().split('\n') if x]
        if len(out_message) != 1 or cdn not in out_message[0]:
            f.write(cdn +' \n')
        
    

