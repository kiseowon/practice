# IVIVE 실험 기록

실험은 연구자가 직접 실행합니다. 기록 도구는 완료된 `comparison.csv` 또는 `summary.csv`에서 집계 지표와 허용된 설정 필드를 읽어 GitHub의 `experiment-records` 브랜치에 저장합니다. 학습·모델 선택·외부 AI API 호출은 수행하지 않습니다.

## 서버에서 한 번 준비

기존 실험 코드와 별도로 기록 전용 checkout을 만드세요. 아래 폴더가 이미 있으면 덮어쓰지 말고 확인하세요.

```bash
git clone --branch experiment-records --single-branch https://github.com/kiseowon/practice.git /home/seowon1211/data1/ivive_records
```

Git 읽기와 push 권한은 다릅니다. 서버의 기존 GitHub 인증을 사용합니다. 필요한 경우 서버에서 허용하는 Git credential manager나 GitHub CLI 로그인 절차를 사용하세요. 토큰을 노트북·명령 문자열·저장소에 넣지 마세요.

## 실험 완료 셀

`EXPERIMENT`만 해당 결과 폴더로 바꿉니다.

```python
from pathlib import Path
import subprocess, sys

RECORDS = Path('/home/seowon1211/data1/ivive_records')
EXPERIMENT = Path('/home/seowon1211/data1/conditioned_expression/pathway_structure_mlp_20261102_035406_327bbc')

subprocess.run([
    sys.executable, str(RECORDS/'automation/record_experiment.py'),
    '--experiment', str(EXPERIMENT),
    '--repo', str(RECORDS),
    '--publish',
], check=True)
```

이 셀은 다음을 수행합니다.

1. 기록 브랜치를 fast-forward로 갱신합니다.
2. `reports/experiments/<실행 폴더명>/`에 `metrics.csv`, `metadata.json`, `report.md`를 저장합니다.
3. 해당 기록과 자동 목록만 commit/push합니다.

같은 내용은 중복 등록하지 않습니다. 같은 run ID에서 지표나 설정이 달라지면 중단합니다. 수정된 결과는 새 `--run-id`로 구분하세요. 기존 파일은 보존합니다. `--notebook /path/to/file.ipynb`를 추가하면 노트북 파일명과 SHA256만 기록합니다.

`--publish`를 빼면 로컬 저장만 합니다. 인증·접속 실패 시 기존 실험은 그대로 남습니다. commit 이후 push가 실패하면 인증·접속을 복구한 뒤 같은 명령을 다시 실행할 수 있습니다. 분기된 Git 이력은 자동으로 merge/rebase/force-push하지 않고 중단합니다.

설정된 Git author가 없으면 기록용 checkout에 `IVIVE automation <ivive-automation@localhost>`를 설정합니다. 연구자의 개인 이메일로 추정하거나 다른 checkout의 설정을 변경하지 않습니다. 여러 노트북의 기록 명령은 같은 checkout의 lock으로 순차 처리합니다.

## 기록의 범위

- 집계 지표의 허용 열만 수집합니다. 전체 발현, 모델 가중치, 원본 로그는 전송하지 않습니다.
- seed·epoch·모델 구성 등의 허용 설정만 기록하고, 파일 경로·인증 설정은 복사하지 않습니다.
- 실험 코드 commit은 근거가 없으면 미확인으로 남깁니다. `collection_runtime`은 기록 시점의 환경이며 학습 환경의 증명이 아닙니다.
- 초기 4조건 MLP 기록은 대화에 붙여준 소수점 6자리 표에서 만들었습니다. 서버 원본 CSV/설정을 읽은 기록과 구분합니다. 실제 서버 파일은 별도 run ID로 추가할 수 있습니다.
- 자동 문장은 지표의 최솟값 등 사실 요약입니다. 통계적 유의성·원인·재현성·최적 모델을 판정하지 않습니다. 연구 해석은 사용자와 assistant가 함께 합니다.

상주 agent나 자동 학습 없이 완료 셀의 기록만 자동화합니다. 후속 연구 상담에서 이 브랜치의 기록을 근거로 참조할 수 있습니다.
