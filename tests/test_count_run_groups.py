#!/usr/bin/env python3
import collections,os,random,subprocess,tempfile,unittest
from pathlib import Path
CORE=Path(os.environ.get('MEGAHIT_TEST_CORE',Path(__file__).resolve().parents[1]/'build/megahit_core'))

class CountRunGroupsTest(unittest.TestCase):
 def test_exact_edges_contexts_and_mercy_candidates(self):
  rng=random.Random(1939)
  def dna(n):return ''.join(rng.choices('ACGT',k=n))
  def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
  with tempfile.TemporaryDirectory(prefix='megahit-run-groups-') as name:
   root=Path(name);sequences=[]
   for i in range(80):
    template=dna(220)
    for j in range(30):
     start=rng.randrange(50);length=rng.randrange(40,151)
     s=template[start:start+length]
     if j%3==0:s=rc(s)
     sequences.extend([s]*(1+j%9))
    sequences += [template[:100]]*140
    sequences += [template[:70]+'A'+template[71:100]]*2
   half=dna(20);sequences += [half+rc(half)]*4
   sequences += ['A'*150]*66000
   fasta=root/'reads.fa';fasta.write_text(''.join(f'>r{i}\n{s}\n' for i,s in enumerate(sequences)))
   config=root/'input.lib';config.write_text('test reads\nse '+str(fasta)+'\n')
   lib=root/'reads.lib'
   built=subprocess.run([str(CORE),'buildlib',str(config),str(lib),'4'],capture_output=True,text=True)
   self.assertEqual(built.returncode,0,built.stderr)
   for threshold in (1,2,3,127):
    with self.subTest(threshold=threshold):
     outputs={}
     for label in ('reference','grouped','filtered'):
      prefix=root/f'{threshold}-{label}'
      env=os.environ.copy();env['MEGAHIT_FORCE_HASH_COUNT']='1'
      env.pop('MEGAHIT_EXPERIMENTAL_COUNT_RUN_GROUPS',None)
      env.pop('MEGAHIT_EXPERIMENTAL_COUNT_SINGLETON_FILTER',None)
      if label!='reference':env['MEGAHIT_EXPERIMENTAL_COUNT_RUN_GROUPS']='1'
      if label=='filtered':env['MEGAHIT_EXPERIMENTAL_COUNT_SINGLETON_FILTER']='1'
      command=[str(CORE),'count','-k','39','-m',str(threshold),'--host_mem','2000000000','--mem_flag','1','--output_prefix',str(prefix),'--num_cpu_threads','4','--read_lib_file',str(lib)]
      run=subprocess.run(command,env=env,capture_output=True,text=True)
      self.assertEqual(run.returncode,0,run.stderr)
      self.assertIn('Exact hash count profile:',run.stderr)
      if label=='grouped':self.assertIn('Exact run groups:',run.stderr)
      records={}
      for suffix,width in (('edges',12),('endcand',13)):
       items=[]
       for file in root.glob(prefix.name+'.'+suffix+'.*'):
        if not file.suffix[1:].isdigit():continue
        data=file.read_bytes();self.assertEqual(len(data)%width,0)
        items.extend(data[i:i+width] for i in range(0,len(data),width))
       records[suffix]=collections.Counter(items)
      for suffix in ('counting','cand'):records[suffix]=Path(str(prefix)+'.'+suffix).read_bytes()
      outputs[label]=records
     for key in outputs['reference']:
      for mode in ('grouped','filtered'):self.assertEqual(outputs['reference'][key],outputs[mode][key],mode+' '+key)

if __name__=='__main__':unittest.main()
