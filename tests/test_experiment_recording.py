import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'automation'))
import record_experiment as record


class RecordingTests(unittest.TestCase):
    def fixture(self,root):
        run=root/'run_001';run.mkdir()
        (run/'comparison.csv').write_text('model,structure,pathway,response_rmse,absolute_rmse,API_KEY\nmlp,False,False,0.31,0.32,DO_NOT_EXPORT\nmlp,True,False,0.30,0.33,DO_NOT_EXPORT\n')
        (run/'config.json').write_text(json.dumps({'seed':42,'epochs':30,'API_KEY':'DO_NOT_EXPORT','prepared':'/private/path'}))
        return run

    def test_allowlisted_export_is_immutable_and_idempotent(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);run=self.fixture(root);repo=root/'repo'
            path,created=record.export_report(run,repo)
            self.assertTrue(created)
            combined='\n'.join(p.read_text() for p in path.iterdir())
            self.assertNotIn('DO_NOT_EXPORT',combined);self.assertNotIn('/private/path',combined)
            self.assertEqual(json.loads((path/'metadata.json').read_text())['experiment_config'],{'seed':42,'epochs':30})
            self.assertFalse(record.export_report(run,repo)[1])
            before=(path/'metrics.csv').read_bytes()
            source=run/'comparison.csv';source.write_text(source.read_text().replace('0.31','0.29'))
            with self.assertRaises(ValueError):record.export_report(run,repo)
            self.assertEqual(before,(path/'metrics.csv').read_bytes())

    def test_publisher_commits_only_report_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);repo=root/'repo';repo.mkdir();remote=root/'remote.git'
            def cmd(*args):subprocess.run(args,check=True,capture_output=True,text=True)
            cmd('git','init','--bare',str(remote));cmd('git','init','-b',record.BRANCH,str(repo))
            cmd('git','-C',str(repo),'config','user.name','Test');cmd('git','-C',str(repo),'config','user.email','test@localhost')
            cmd('git','-C',str(repo),'remote','add','origin',str(remote))
            (repo/'README.md').write_text('test\n');(repo/'.gitignore').write_text('*.csv\n');cmd('git','-C',str(repo),'add','README.md','.gitignore');cmd('git','-C',str(repo),'commit','-m','init')
            cmd('git','-C',str(repo),'push','-u','origin',record.BRANCH)
            run=self.fixture(root);path,_=record.export_report(run,repo);record.update_index(repo)
            (repo/'unrelated.txt').write_text('preserve me')
            self.assertEqual(record.publish_report(repo,path),'pushed')
            tracked=record.git(repo,'show','--pretty=format:','--name-only','HEAD').stdout.splitlines()
            self.assertTrue(all(p.startswith('reports/experiments/') for p in tracked if p))
            self.assertIn('reports/experiments/run_001/metrics.csv',tracked)
            self.assertTrue((repo/'unrelated.txt').exists())
            self.assertEqual(record.publish_report(repo,path),'unchanged')
