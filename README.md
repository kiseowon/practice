# IVIVE: 구조·생물학적 관계를 학습에 반영하기

**연구 질문: 구조·pathway 정보가 미관측 화합물의 IVIVE 예측에 도움이 되는가?** Rat liver의 조건 평균 in vitro 발현으로 in vivo 발현을 예측합니다. 생물학적 정보를 예측 후 보정하는 대신 처음부터 모델 학습에 반영할 수 있는지 탐색합니다. 현재 encoder는 별도 사전학습 없이 최종 목적함수로 함께 end-to-end 학습합니다.

AIVIVE는 in vitro→in vivo 변환과 생물학적 신호 보완이라는 문제 설정의 출발점입니다. GEARS는 생물학적 관계를 모델에 활용하는 아이디어의 영감이며, 현재 구현이 GEARS의 그래프나 유전자 교란 예측 모델을 구현한 것은 아닙니다. **AIVIVE·GEARS를 그대로 재현한 연구가 아닙니다.** RAID는 구조 정보와 in vitro 발현을 이용한 예측의 선행연구이고, sparse 연결 구현은 P-NET 원리를 참고했습니다. 현재 KEGG 소속 행렬은 완성된 지식 그래프 기반 diffusion이나 별도 pretraining의 검증을 뜻하지 않습니다.

## 실험 발전과 현재 최종 비교

과거 24h·109 probe 탐색 → 조건 평균·3453 probe 연구 재시작 → pathway 없는 모델 비교 → **평균 pathway 12조건 → sparse pathway 9조건**으로 진행했습니다. 자세한 [연구 상태](docs/research-state.md)와 [원본·자료 이력](docs/metadata-and-compound-history.md)을 참고하세요.

| 단계 | 비교 범위 | 결과 |
| --- | --- | --- |
| 평균 pathway | MLP·cGAN·epsilon-DDPM × 구조 off/on × 평균 pathway off/on, 12조건 | [별도 결과 CSV](results/diffusion_second_trial_mean_pathway_comparison.csv) |
| 최종 sparse pathway | MLP·cGAN·epsilon-DDPM × pathway 없음/실제 KEGG sparse/무작위 sparse, 9조건. 구조는 모두 켬 | [6자리 요약](results/meeting_comparison.csv) · [서버 집계 CSV](results/server_sparse_comparison.csv) |

seed42 화합물 분할은 내부 학습89 / 검증23 / 테스트28입니다. 발현 3453 probe, vitro 1675·vivo 2224조건, pair 19,936개입니다. 검증으로 epoch를 고르고 112개 화합물로 다시 학습한 뒤 테스트를 평가하는 코드입니다. 이 표는 조건 평균의 예측 오차이며 생성 분포 품질 평가는 아닙니다.

| 모델 | pathway | 절대 RMSE ↓ | 반응 RMSE ↓ |
| --- | --- | ---: | ---: |
| MLP | 없음 | 0.337249 | 0.313829 |
| MLP | 실제 sparse | 0.333559 | 0.309485 |
| MLP | 무작위 sparse | 0.335970 | 0.313035 |
| cGAN | 없음 | 0.379863 | 0.332292 |
| cGAN | 실제 sparse | 0.388245 | 0.332922 |
| cGAN | 무작위 sparse | 0.410567 | 0.333861 |
| epsilon-DDPM | 없음 | 0.354855 | 0.330387 |
| epsilon-DDPM | 실제 sparse | 0.356288 | 0.328952 |
| epsilon-DDPM | 무작위 sparse | 0.357838 | 0.331531 |

이 분할과 단일 무작위 대조에서 실제 KEGG sparse MLP가 두 오차 모두 가장 낮았습니다. **DDPM에서도 실제 KEGG가 반응 RMSE를 소폭 낮췄고 무작위 연결보다 낮았습니다.** 다만 절대 오차는 pathway 없음보다 높았으며 cGAN에서는 개선되지 않았습니다. 변화 없음 반응 기준선은 0.330116입니다. 모든 9조건이 structure=True이므로 이 표만으로 구조 추가 효과를 주장하지 않습니다.

단일 분할·단일 mask·테스트28개이며 같은 테스트 결과를 보며 설계를 발전시킨 탐색 실험입니다. 독립적인 확증 검증이나 유의성을 주장하지 않습니다. 반복 분할·여러 mask·화합물 bootstrap CI·훈련 평균 반응 기준선은 미실행 후속 계획이고, 추가 학습·실험은 현재 승인되지 않았습니다.

## 파일 선택과 실행

[실행 안내](docs/execution-guide.md)에서 **최종 실행용 / 공부·미팅용 / 원본 보존용 / 재구성본**을 구분하고 입력·순서·출력을 확인하세요. 원본 보존 노트북에는 과거 학습·설치·외부 요청이 있으므로 일괄 Run All하지 않습니다. 기존 파일 경로는 유지했습니다.

최종 학습 경로는 pathway_pattern_mlp → pathway_random_control → sparse_pathway_gan_ddpm 노트북 순서입니다. prepared_controls 6파일, Morgan binary, 기존 KEGG snapshot, 단계별 참조 산출물이 필요합니다. 서버 절대경로를 사용하므로 새 환경에서는 각 실행 셀의 경로를 명시적으로 설정해야 합니다. 준비된 입력과 원시 발현·snapshot·체크포인트 전체는 공개 저장소에 없습니다.

requirements.txt는 의존성 범위이며 당시 서버 환경 lockfile이 아닙니다. 현재 클라우드에 설치된 버전을 과거 실험 버전으로 간주하지 않습니다.

## 반복 검증 코드

[반복 검증 실행기·설계](docs/repeated-validation.md)를 추가했습니다. 기본은 화합물 분할5개 × mask3개, 학습 seed 고정이며 계획만 출력합니다. 실제 반복 학습은 아직 실행하지 않았고, `--execute`는 별도 실행 승인 뒤에만 사용합니다. 기존 단일 실험 결과는 유지합니다.

## 검증 상태

[수행 환경·근거별 검증 상태](docs/verification-status.md)에 원본 확보와 산출물 검증을 구분했습니다.

- 사용자 서버: 원시 발현·현재 corrected 메타데이터로 준비 6파일 재생성, 전부 정확히 일치. NPY 최대 차이 0. 사용 코드 hash는 보존 원본과 일치.
- 사용자 서버: Morgan binary 137개 일치, 무작위 mask degree 보존·겹침률 검사 통과, 저장 예측의 9조건 지표 모두 일치 보고. 추가 모델 학습·예측 재생성 검증과 구분.
- 이 작업 환경: 첨부 파일 정합성과 KEGG 저장 원문·annotation·audit·projection의 정확한 일치 확인. 새 다운로드·원본 생성 코드 실행 없이 오프라인 메모리 대조.

전체 모델 재학습이나 새 환경의 end-to-end 재현을 완료한 것은 아닙니다. corrected 메타데이터의 생성 버전(In7/In10), reference 취득 출처, 당시 환경 전체, 과거 SVD64 생성 출처는 미확인입니다.

[코드 해설](docs/code-walkthrough.md) · [실험 상세](docs/experimental-details.md) · [KEGG 원본](docs/kegg-original-recovered.md) · [출처 기록](results/provenance.json)

현재 checkout은 과거 모든 실험의 완전한 보관소가 아닙니다. README의 이전 기록상 탐색 코드는 Git 이력에, 보고 자동화는 experiment-records 브랜치에 남아 있습니다. 해당 이력 전체를 감사한 것은 아닙니다.
