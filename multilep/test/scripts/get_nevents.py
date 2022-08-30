import os,sys,subprocess, glob, time
def system(command):
  try:
    return subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT).decode()
  except subprocess.CalledProcessError as e:
    print(e.output)
    return 'error'

search_string = sys.argv[1]
output = system('dasgoclient --query="dataset='+search_string+'"')
samples = output.split('\n')
samples = [sample for sample in samples if sample]

for sample in samples:
    output = system('dasgoclient -query="dataset='+sample+' | grep dataset.nevents,dataset.nfiles"').split('\n')
    output = [line.strip().split() for line in output if line and not '[' in line][0]
    output = [line for line in output if line.isalnum()]
    print('{}: {}'.format(sample, output[0]))

print('')
for sample in samples:
    print('{}'.format(sample))
