# 미팅용 코드와 실제 시험 결과의 출처

> 추가 첨부의 원본 코드·설정 기록을 반영한 최신 설명은 [서버 원본 추가 확인](server-source-update.md)을 보세요. 기존 미확인 사항 중 수동 CID·SMILES와 matched Morgan의 저장 설정은 새 근거로 확인되었습니다.

사용자가 지정한 `연구미팅 1차 재현.ipynb`는 [미팅용 노트북](../notebooks/research_meeting_reproduction.ipynb)으로 파일 내용을 그대로 복사했습니다. `디퓨전 2차 시험.ipynb`의 탐색·시험 코드는 GitHub에 추가하지 않았습니다. 기존 scripts·노트북·9조건 결과는 보존했습니다. 추가 학습·재평가는 실행하지 않았습니다.

## 미팅 노트북이 보여주는 내용

구조 입력 고정, MLP·cGAN·epsilon-DDPM의 sparse pathway 설명과 학습·평가 함수가 들어 있습니다. `RUN_TRAINING=False`이며 학습 함수 호출은 없지만, 데이터 로딩·모델 생성·미학습 모델의 검증 예측 등은 실행합니다. 따라서 단순한 Run All을 학습 완료나 실험 재현으로 해석하면 안 됩니다. 서버의 원자료, KEGG snapshot, 무작위 mask, 기존 결과 CSV·분할 파일이 필요합니다.

0부터 센 셀 2는 서버의 `sparse_pathway_gan_ddpm_20261102_062000_f4c529/comparison.csv`를 읽고, 셀 15는 기존 9조건의 반응 RMSE를 표시합니다. 그 저장 출력은 현재 저장소의 9조건 반응 RMSE와 일치하지만, 이 노트북에서 새로 학습한 결과는 아닙니다. 셀 14의 미학습 검증 점수는 최종 연구 결과에 포함하지 않습니다. 내부 조건 차원 16, refit 조건 차원 17이 기록되어 있어 단계별 차원 구분이 필요합니다.

## 실제 시험 첨부에서 확인한 결과

`디퓨전 2차 시험.ipynb`에는 여러 과거 탐색 실험이 섞여 있습니다. 마지막 완료 비교는 0부터 센 셀 92의 실행 로그와 셀 93의 12행 표이며 run은 `vanilla_information_matrix_20261102_040702_9210e0`입니다.

- MLP·cGAN·epsilon-DDPM × 구조 off/on × 평균 pathway off/on의 12조건입니다.
- 셀 90의 context는 `PathwayStructureMLP`와 pathway 평균 score를 사용합니다. 실제 sparse/무작위 sparse 9조건과 같은 표현이 아닙니다.
- 셀 93의 첫 HTML 표에 표시된 12행을 [결과 CSV](../results/diffusion_second_trial_mean_pathway_comparison.csv)로 추출했습니다. 표시된 소수점 여섯 자리 값이며 원시 예측으로 재계산하지 않았습니다.
- 반응 RMSE 최저는 MLP 구조 on·평균 pathway on의 0.309828, 절대 RMSE 최저는 MLP 구조 off·pathway off의 0.324497입니다. 지표별 최저 조건은 다릅니다.
- 첨부 시험 파일에는 최종 9조건 sparse 표가 없습니다. 따라서 기존 `meeting_comparison.csv`를 이 시험 파일의 값으로 덮어쓰거나, 이 시험 첨부를 9조건 전체의 검증 출처로 표시하지 않았습니다.

[별도 provenance](../results/diffusion_second_trial_provenance.json)에 파일 SHA-256, 셀·출력 위치, run과 한계를 기록했습니다. 두 실험을 합친 순위표나 반복 검증 결과로 해석하지 않습니다. 사용자는 후속 확인에서 12조건은 이전 실험이고 9조건은 별도 서버 완료 결과이며, 기존 `meeting_comparison.csv`는 그 9조건 출력의 요약이라고 명시했습니다. 제공한 반응 RMSE 9개는 기존 CSV 및 미팅 노트북의 저장 출력과 일치합니다. 서버 원본 `comparison.csv` 자체는 이번 첨부에 없으므로 원본 CSV 전체와 예측 파일을 직접 대조한 것은 아닙니다.
