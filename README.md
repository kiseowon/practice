# IVIVE: structure and sparse pathway inputs

첫 연구미팅용 코드. 조건별 in vitro 발현과 약물 구조·실험 조건으로 in vivo 발현을 예측하고, KEGG pathway 표현의 효과를 MLP / cGAN / epsilon-DDPM에서 비교합니다.

## 실행

```bash
python -m pip install -r requirements.txt
jupyter lab
```

독립형 노트북은 저장소 scripts를 import하지 않아 서버에 노트북만 올려도 실행됩니다. 새 커널에서 Run All 하세요. 기본 실행은 CPU입니다.

1. `notebooks/pathway_pattern_mlp.ipynb`: 구조 입력 고정, pathway 없음 / 평균 / 실제 sparse MLP.
2. `notebooks/pathway_random_control.ipynb`: 기존 세 결과를 읽고 무작위 sparse MLP 하나만 추가.
3. `notebooks/sparse_pathway_gan_ddpm.ipynb`: MLP 결과 재사용, cGAN·DDPM × 없음/실제/무작위 6개 학습. 같은 OUT으로 재개하면 완료 모델은 건너뜁니다.

각 노트북 마지막 셀의 데이터·이전 결과 경로를 수정하세요. 노트북에는 당시 실제 서버 경로가 기본값으로 들어 있습니다. 실험 중인 서버 파일은 GitHub 정리와 별개입니다.

## 필요한 로컬 데이터

- prepared_controls: `vitro_expression.npy`, `vivo_expression.npy`, `vitro_conditions.csv`, `vivo_conditions.csv`, `pairs.csv.gz`, `probes.csv`.
- Morgan binary CSV: `compound`와 2048 fingerprint 컬럼. 누락 구조는 0과 presence flag로 처리.
- KEGG Rat snapshot: `projection.npy`, `pathway_audit.csv`, `provenance.json`. 기존 생성 snapshot을 사용합니다. 소속은 고정 외부 지식이며 probe 순서와 Rat 종을 검사합니다.
- 재사용 비교에는 이전 노트북의 comparison / config / compound_split / predictions 파일이 필요합니다.

원자료·annotation snapshot·체크포인트는 GitHub에 포함하지 않습니다. 따라서 새 환경에서 원자료 없이 실제 연구 결과를 재현할 수는 없습니다. 패키지 버전은 범위이며 당시 서버 환경의 lockfile은 확보하지 못했습니다.

## 읽는 순서

[실제 코드 해설](docs/code-walkthrough.md)을 먼저 보세요. `scripts/first_structure_mlp.py`는 전처리·분할·평가 공통 코드, `pathway_structure_mlp.py`는 기존 snapshot 로더와 평균 표현, `pathway_pattern_mlp.py`는 sparse encoder, `pathway_random_control.py`는 연결 대조군, `sparse_pathway_information_matrix.py`는 생성모델 학습입니다.

[실험 상세·SMILES 매핑·해석의 한계](docs/experimental-details.md)에는 사용자가 확인한 서버 매핑 절차, Morgan 생성 코드와 최종 입력의 미확인 연결, end-to-end 학습 방식 및 미실행 후속 계획을 정리했습니다. 최종 `morgan_binary.csv`의 radius 2 설정은 검증되지 않았습니다.

## 현재 결과와 범위

[9조건 결과](results/meeting_comparison.csv)는 사용자가 서버 출력으로 제공한 소수점 여섯 자리 요약입니다. 이 환경에서 실제 데이터로 다시 계산하지 않았습니다. [출처](results/provenance.json).

모든 조건에 구조를 제공했습니다. 실제 KEGG sparse MLP가 절대 RMSE 0.333559, 반응 RMSE 0.309485로 가장 좋았습니다. 변화 없음 반응 RMSE는 0.330116입니다. DDPM에서는 반응만 소폭 개선됐고 cGAN은 pathway 없음보다 개선되지 않았습니다. 한 seed·반복 관찰한 28개 테스트 화합물 결과이며 유의성이나 독립 검증을 주장하지 않습니다.

이는 단일 분할과 무작위 mask 1개의 탐색 결과입니다. MLP의 실제 KEGG는 이 분할과 단일 무작위 대조에서 오차가 낮았고, 생성모델의 pathway 효과는 제한적입니다. 9조건 모두 `structure=True`이므로 이 표만으로 구조 추가 효과를 주장할 수 없습니다. 검증으로 epoch를 골랐지만 같은 테스트 결과를 보며 다음 설계를 발전시켰으므로 완전히 독립적인 테스트는 아닙니다. 여러 분할·여러 무작위 mask·화합물 bootstrap CI·훈련 평균 반응 기준선은 아직 실행하지 않은 후속 계획입니다.

이는 AIVIVE/P-NET 원문 모델 재현이 아닙니다. Sparse 연결 원리 참고: [P-NET, Nature 2021](https://doi.org/10.1038/s41586-021-03922-4). 예측은 probe 수준이며 모든 probe가 고유 유전자에 일대일 대응하는 것은 아닙니다.

기존 main의 탐색 코드는 현재 파일 목록에서 제거했으며 Git 이력에 보존했습니다. 이전 보고 자동화는 `experiment-records` 브랜치에 남아 있습니다.
