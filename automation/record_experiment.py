#!/usr/bin/env python3
"""Export aggregate IVIVE results; optionally commit/push a dedicated reports checkout.

No training, API calls, raw-expression copying, or automatic model selection.
"""
import argparse
import contextlib
import csv
from datetime import datetime
import fcntl
import hashlib
import html
import importlib.metadata
import io
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from zoneinfo import ZoneInfo

BRANCH='experiment-records'
REPOSITORY='kiseowon/practice'
LABELS=['model','family','structure','pathway','scope','group']
NUMBERS=['absolute_rmse','response_rmse','zero_response_rmse','ratio_to_no_change','pair_rmse','energy_score','mean_prediction_rmse','selected_epoch','parameters','discriminator_parameters','n_compounds','n_target_conditions','n_treated_conditions','generated_variance_per_probe','observed_variance_per_probe']
CONFIG_KEYS=['seed','model_seed','epochs','base_epochs','calibration_epochs','check_every','batch_size','fold','folds','threads','samples','steps','ddpm_steps','diffusion_steps','sampling_steps','latent','input_components','latent_components','graph_k','graph_threshold','sigma_floor','architecture','protocol','selection','input_mode','pathway_enabled','models','adversarial_weight']
PRIMARY=['absolute_rmse','response_rmse','energy_score','mean_prediction_rmse']


def git(repo,*args,check=True):
    result=subprocess.run(['git','-C',str(repo),*args],capture_output=True,text=True)
    if check and result.returncode:
        # Avoid forwarding remote URLs/credential-bearing stderr into reports.
        raise RuntimeError(f'git {args[0]} failed (exit {result.returncode}). Check Git authentication/connectivity in this dedicated checkout; no force push was attempted.')
    return result


def safe_id(value):
    result=re.sub(r'[^A-Za-z0-9_.-]+','_',value).strip('._')[:140]
    if not result:raise ValueError('Invalid experiment ID')
    return result


def parse_value(key,value):
    value=(value or '').strip()
    if key in ['structure','pathway'] and value.lower() in ['true','false']:return value.lower()=='true'
    if key in NUMBERS:
        if not value or value.lower() in ['nan','none','null']:return None
        number=float(value)
        if not math.isfinite(number):raise ValueError(f'Nonfinite aggregate metric: {key}')
        return int(number) if key in ['selected_epoch','parameters','discriminator_parameters','n_compounds','n_target_conditions','n_treated_conditions'] and number.is_integer() else number
    return value[:180]


def aggregate_metrics(source,family=None):
    with source.open(encoding='utf-8-sig',newline='') as handle:
        reader=csv.DictReader(handle)
        fields=[k for k in LABELS+NUMBERS if k in (reader.fieldnames or [])]
        if not any(k in fields for k in PRIMARY):raise ValueError('Expected comparison.csv or summary.csv with aggregate metrics')
        rows=[]
        for row in reader:
            if len(rows)>=250:raise ValueError('Too many rows: use aggregate results, not per-probe/sample data')
            clean={key:parse_value(key,row.get(key)) for key in fields}
            if not any(clean.get(k) is not None for k in PRIMARY):raise ValueError('A row has no usable performance metric')
            if family and not clean.get('family'):clean['family']=family
            rows.append(clean)
    if not rows:raise ValueError('Zero result rows')
    fields=[k for k in LABELS+NUMBERS if any(k in row for row in rows)]
    stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=fields,lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return fields,rows,stream.getvalue()


def reviewed_config(experiment):
    path=experiment/'config.json'
    if not path.is_file():return {},'not provided'
    raw=json.loads(path.read_text(encoding='utf-8'))
    allowed={}
    for key in CONFIG_KEYS:
        if key not in raw:continue
        value=raw[key]
        if isinstance(value,(str,int,float,bool)) or value is None:
            if isinstance(value,str) and len(value)>500:continue
            if isinstance(value,float) and not math.isfinite(value):continue
            allowed[key]=value
        elif isinstance(value,list) and all(isinstance(x,(str,int,float,bool)) for x in value):allowed[key]=value[:30]
    return allowed,'allowlisted fields only; filesystem/credential configuration omitted'


def md(value):
    if value is None:return '—'
    if isinstance(value,float):return f'{value:.6f}'
    return html.escape(str(value)).replace('|','\\|').replace('\n',' ').replace('\r',' ')


def observations(rows):
    groups={}
    for row in rows:groups.setdefault((row.get('scope',''),row.get('group','all')),[]).append(row)
    result=[]
    for (scope,group),part in groups.items():
        prefix='/'.join(x for x in [scope,group] if x)
        for metric in PRIMARY:
            candidates=[r for r in part if r.get(metric) is not None]
            if not candidates:continue
            winner=min(candidates,key=lambda r:r[metric])
            label=winner.get('family') or winner.get('model') or 'row'
            if 'structure' in winner:label+=f', structure={winner["structure"]}'
            if 'pathway' in winner:label+=f', pathway={winner["pathway"]}'
            result.append(f'- {md(prefix)} / {metric}: 최저 관측값 {winner[metric]:.6f} ({md(label)}).')
    return result


def build_report(run_id,fields,rows,metadata):
    selected=[k for k in ['family','model','structure','pathway','scope','group']+PRIMARY+['zero_response_rmse','ratio_to_no_change','selected_epoch','parameters'] if k in fields]
    lines=[f'# {md(run_id)}','',f'기록 시각: {metadata["recorded_at_seoul"]}','',
        f'출처: {md(metadata["provenance"])}. 입력 파일: {md(metadata["source_filename"])}.','',
        '## 결과','', '| '+' | '.join(selected)+' |','| '+' | '.join(['---']*len(selected))+' |']
    lines.extend('| '+' | '.join(md(row.get(k)) for k in selected)+' |' for row in rows)
    lines+=['','## 자동 확인 사항','']+observations(rows)
    lines+=['','## 설정','', '```json',json.dumps(metadata['experiment_config'],ensure_ascii=False,indent=2),'```','',
        '## 해석 범위','',
        '- 위 문장은 저장된 지표의 기술적 요약입니다. 유의성·재현성·원인이나 최적 모델을 확정하지 않습니다.',
        '- 다른 실행과 비교하기 전 분할·입력·전처리·평가 기준이 같은지 확인해야 합니다.',
        '- 실험 코드 버전은 별도 근거가 없으면 미확인으로 남깁니다. 기록 도구 버전을 학습 코드 버전으로 취급하지 않습니다.']
    if metadata.get('note'):lines+=['',f'기록 메모: {md(metadata["note"])}']
    return '\n'.join(lines)+'\n'


def export_report(experiment,repo,run_id=None,family=None,provenance='server aggregate file',note='',notebook=None):
    experiment=Path(experiment).resolve();repo=Path(repo).resolve()
    if not experiment.is_dir():raise ValueError('Experiment directory does not exist')
    source=next((experiment/name for name in ['comparison.csv','summary.csv'] if (experiment/name).is_file()),None)
    if source is None:raise ValueError('No comparison.csv or summary.csv found')
    fields,rows,metrics=aggregate_metrics(source,family)
    config,config_status=reviewed_config(experiment)
    run_id=safe_id(run_id or experiment.name)
    notebook_info=None
    if notebook:
        nb=Path(notebook)
        notebook_info={'filename':nb.name,'sha256':hashlib.sha256(nb.read_bytes()).hexdigest(),'status':'file observed at recording; not proof of execution'}
    fingerprint=hashlib.sha256(json.dumps({'rows':rows,'config':config,'provenance':provenance,'note':note,'notebook':notebook_info},sort_keys=True,ensure_ascii=False).encode()).hexdigest()
    destination=repo/'reports'/'experiments'/run_id
    if destination.exists():
        previous=json.loads((destination/'metadata.json').read_text())
        if previous.get('content_sha256')!=fingerprint:raise ValueError('Run ID already has different results; choose a new --run-id. Existing records were preserved.')
        return destination,False
    metadata=dict(schema_version=1,run_id=run_id,recorded_at_seoul=datetime.now(ZoneInfo('Asia/Seoul')).isoformat(timespec='seconds'),
        provenance=provenance,note=note,source_filename=source.name,content_sha256=fingerprint,
        source_file_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),experiment_config=config,config_status=config_status,
        training_code_revision=None,collection_runtime={'python':sys.version.split()[0]})
    for name in ['numpy','pandas','torch','scikit-learn']:
        try:metadata['collection_runtime'][name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:pass
    if notebook_info:metadata['notebook_fingerprint']=notebook_info
    destination.mkdir(parents=True,exist_ok=False)
    (destination/'metrics.csv').write_text(metrics,encoding='utf-8')
    (destination/'metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    (destination/'report.md').write_text(build_report(run_id,fields,rows,metadata),encoding='utf-8')
    return destination,True


def update_index(repo):
    root=repo/'reports'/'experiments';root.mkdir(parents=True,exist_ok=True)
    ignore=root/'.gitignore'
    if not ignore.exists():ignore.write_text('!**/metrics.csv\n',encoding='utf-8')
    lines=['# IVIVE experiment records','', '자동 생성된 실험 요약 목록입니다. 학습 실행은 연구자가 직접 수행합니다.','',
        '| Experiment | Recorded (Asia/Seoul) | Provenance |','| --- | --- | --- |']
    for path in sorted(root.iterdir(),reverse=True):
        if not path.is_dir() or not (path/'metadata.json').is_file():continue
        info=json.loads((path/'metadata.json').read_text())
        lines.append(f'| [{md(info["run_id"])}]({path.name}/report.md) | {md(info["recorded_at_seoul"])} | {md(info["provenance"])} |')
    (root/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def remote_repository(url):
    if url.startswith(('https://','http://','ssh://')):
        from urllib.parse import urlsplit
        parts=urlsplit(url);host=parts.hostname;path=parts.path
    elif url.startswith('git@github.com:'):host='github.com';path=url.split(':',1)[1]
    else:return None
    return path.strip('/').removesuffix('.git') if host=='github.com' else None


def validate_publisher(repo):
    if git(repo,'branch','--show-current').stdout.strip()!=BRANCH:raise ValueError(f'Use a separate checkout of branch {BRANCH}; this tool does not switch branches')
    remote=git(repo,'remote','get-url','origin').stdout.strip()
    if remote_repository(remote)!=REPOSITORY:raise ValueError('Unexpected origin repository')
    if git(repo,'diff','--cached','--name-only').stdout.strip():raise ValueError('Index contains staged changes; preserve them before publishing a report')
    # Do not bring unrelated local code changes into an automated report commit.
    pending=git(repo,'diff','--name-only',f'origin/{BRANCH}...HEAD',check=False).stdout.splitlines()
    if any(not p.startswith('reports/experiments/') for p in pending):raise ValueError('Unpublished commits contain non-report changes')
    changed=git(repo,'diff','--name-only').stdout.splitlines()
    if any(not p.startswith('reports/experiments/') for p in changed):raise ValueError('Dedicated reports checkout has unrelated tracked changes')


@contextlib.contextmanager
def recording_lock(repo,publish):
    repo=Path(repo)
    found=git(repo,'rev-parse','--absolute-git-dir',check=False)
    if found.returncode==0:
        directory=Path(found.stdout.strip())
    elif publish:
        raise ValueError('Publishing requires a dedicated Git checkout')
    else:
        directory=repo/'reports';directory.mkdir(parents=True,exist_ok=True)
    with (directory/'ivive-record.lock').open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX)
        try:yield
        finally:fcntl.flock(handle,fcntl.LOCK_UN)


def publish_report(repo,destination):
    paths=[str(destination.relative_to(repo)), 'reports/experiments/README.md','reports/experiments/.gitignore']
    git(repo,'add','--',*paths)
    if git(repo,'diff','--cached','--quiet',check=False).returncode==0:
        pending=git(repo,'rev-list','--count',f'origin/{BRANCH}..HEAD',check=False).stdout.strip()
        git(repo,'push','origin',f'HEAD:{BRANCH}')
        return 'pushed' if pending and pending!='0' else 'unchanged'
    if not git(repo,'config','user.name',check=False).stdout.strip():git(repo,'config','user.name','IVIVE automation')
    if not git(repo,'config','user.email',check=False).stdout.strip():git(repo,'config','user.email','ivive-automation@localhost')
    git(repo,'commit','-m',f'reports: record {destination.name}')
    git(repo,'push','origin',f'HEAD:{BRANCH}')
    return 'pushed'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--experiment',type=Path,required=True);p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--run-id');p.add_argument('--family');p.add_argument('--provenance',default='server aggregate file');p.add_argument('--note',default='');p.add_argument('--notebook',type=Path);p.add_argument('--publish',action='store_true')
    a=p.parse_args();repo=a.repo.resolve()
    with recording_lock(repo,a.publish):
        if a.publish:
            validate_publisher(repo)
            git(repo,'pull','--ff-only','origin',BRANCH)
        path,created=export_report(a.experiment,repo,a.run_id,a.family,a.provenance,a.note,a.notebook)
        update_index(repo)
        status=publish_report(repo,path) if a.publish else 'local-only'
    print(json.dumps({'report':str(path),'created':created,'status':status},ensure_ascii=False))


if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError,OSError) as error:
        print(f'Recording stopped: {error}',file=sys.stderr);sys.exit(1)
