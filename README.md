# IVIVE: structure and sparse pathway inputs

첫 연구미팅용 코드. 조건별 in vitro 발현과 약물 구조·실험 조건으로 in vivo 발현을 예측하고, KEGG pathway 표현의 효과를 MLP / cGAN / epsilon-DDPM에서 비교합니다.

## 실행

[현재 연구 상태](docs/research-state.md) · [원본 확보·검증 상태](docs/verification-status.md) · [메타데이터 보정·137개 선정 이력](docs/metadata-and-compound-history.md)

[추가 확인된 서버 원본](docs/server-source-update.md): [SMILES 매핑](notebooks/smiles_mapping_original.ipynb), [matched Morgan 생성](scripts/structure_features.py), [sparse 9조건 실행 기록](notebooks/sparse_pathway_server_record.ipynb)을 반영했습니다. [KEGG In89·In90 원본을 복구](docs/kegg-original-recovered.md)했습니다. 후속 첨부의 저장 KEGG 원문·annotation으로 snapshot의 321개 열과 audit가 정확히 일치함을 이 환경에서 확인했습니다. Morgan은 사용자 서버에서 저장 SMILES·config와 동일 RDKit 버전으로 재계산해 137개 binary 전부 일치했습니다. 이 작업 환경의 독립 재검증은 아닙니다.

매핑부터 확인할 때는 [매핑 준비·연구 재시작 안내](docs/mapping-and-restart.md)를 보세요. [매핑 재구성 코드](scripts/reconstruct_mapping.py)는 원본과 구분해 추가했고, [In[77]부터의 실험 원본 보존본](notebooks/research_restart_original.ipynb)은 필요한 In[76] 정의부터 포함합니다. 원본 보존본에는 학습 호출이 있으므로 실행 승인 없이 Run All하지 마세요.

미팅용 단계별 코드는 첨부 `연구미팅 1차 재현.ipynb`를 그대로 보존한 [연구미팅 재현 노트북](notebooks/research_meeting_reproduction.ipynb)입니다. 서버 데이터·기존 결과·무작위 mask가 필요하며, 함수 정의와 미학습 검증 예측을 포함합니다. 아래 기존 학습 노트북들과 용도를 구분하세요. 이번 갱신에서 추가 학습은 실행하지 않았습니다.

실제 시험 첨부 `디퓨전 2차 시험.ipynb`의 마지막 완료 12조건 결과는 [별도 CSV](results/diffusion_second_trial_mean_pathway_comparison.csv)에 옮겼습니다. 이 표는 **구조 on/off × 평균 pathway on/off**이며 아래 최종 9조건 sparse 표와 다른 실험입니다. [첨부 파일 반영 기록](docs/meeting-notebook-source.md)에 출처와 차이를 정리했습니다.

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

[실험 상세·SMILES 매핑·해석의 한계](docs/experimental-details.md)에는 사용자가 확인한 서버 매핑 절차, Morgan 생성 코드와 최종 입력의 미확인 연결, end-to-end 학습 방식 및 미실행 후속 계획을 정리했습니다. 추가 첨부의 matched config에는 radius 2가 기록되어 있으며, CSV checksum 및 독립 재생성 검증은 없습니다.

## 현재 결과와 범위

[9조건 결과](results/meeting_comparison.csv)는 사용자가 서버 출력으로 제공한 소수점 여섯 자리 요약입니다. 이 환경에서 실제 데이터로 다시 계산하지 않았습니다. [출처](results/provenance.json).

모든 조건에 구조를 제공했습니다. 실제 KEGG sparse MLP가 절대 RMSE 0.333559, 반응 RMSE 0.309485로 가장 좋았습니다. 변화 없음 반응 RMSE는 0.330116입니다. DDPM에서는 반응만 소폭 개선됐고 cGAN은 pathway 없음보다 개선되지 않았습니다. 한 seed·반복 관찰한 28개 테스트 화합물 결과이며 유의성이나 독립 검증을 주장하지 않습니다.

이는 단일 분할과 무작위 mask 1개의 탐색 결과입니다. MLP의 실제 KEGG는 이 분할과 단일 무작위 대조에서 오차가 낮았고, 생성모델의 pathway 효과는 제한적입니다. 9조건 모두 `structure=True`이므로 이 표만으로 구조 추가 효과를 주장할 수 없습니다. 검증으로 epoch를 골랐지만 같은 테스트 결과를 보며 다음 설계를 발전시켰으므로 완전히 독립적인 테스트는 아닙니다. 여러 분할·여러 무작위 mask·화합물 bootstrap CI·훈련 평균 반응 기준선은 아직 실행하지 않은 후속 계획입니다.

이는 AIVIVE/P-NET 원문 모델 재현이 아닙니다. Sparse 연결 원리 참고: [P-NET, Nature 2021](https://doi.org/10.1038/s41586-021-03922-4). 예측은 probe 수준이며 모든 probe가 고유 유전자에 일대일 대응하는 것은 아닙니다.

기존 main의 탐색 코드는 현재 파일 목록에서 제거했으며 Git 이력에 보존했습니다. 이전 보고 자동화는 `experiment-records` 브랜치에 남아 있습니다.
