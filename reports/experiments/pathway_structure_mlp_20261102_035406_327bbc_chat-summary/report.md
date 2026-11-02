# pathway_structure_mlp_20261102_035406_327bbc_chat-summary

기록 시각: 2026-11-02T13:24:21+09:00

출처: user-pasted-table. 입력 파일: comparison.csv.

## 결과

| family | model | structure | pathway | absolute_rmse | response_rmse | zero_response_rmse | ratio_to_no_change | selected_epoch | parameters |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MLP | expression | False | False | 0.324497 | 0.312818 | 0.330116 | 0.947600 | 10 | 839933 |
| MLP | expression+structure | True | False | 0.338225 | 0.311451 | 0.330116 | 0.943459 | 10 | 839933 |
| MLP | expression+pathway | False | True | 0.324664 | 0.316732 | 0.330116 | 0.959456 | 15 | 839933 |
| MLP | expression+structure+pathway | True | True | 0.333340 | 0.309828 | 0.330116 | 0.938541 | 10 | 839933 |

## 자동 확인 사항

- all / absolute_rmse: 최저 관측값 0.324497 (MLP, structure=False, pathway=False).
- all / response_rmse: 최저 관측값 0.309828 (MLP, structure=True, pathway=True).

## 설정

```json
{}
```

## 해석 범위

- 위 문장은 저장된 지표의 기술적 요약입니다. 유의성·재현성·원인이나 최적 모델을 확정하지 않습니다.
- 다른 실행과 비교하기 전 분할·입력·전처리·평가 기준이 같은지 확인해야 합니다.
- 실험 코드 버전은 별도 근거가 없으면 미확인으로 남깁니다. 기록 도구 버전을 학습 코드 버전으로 취급하지 않습니다.

기록 메모: 사용자가 대화에 공유한 소수점 6자리 결과표에서 기록했습니다. 서버 원본 CSV와 config.json은 읽지 않았으므로 실제 학습 설정과 코드 버전은 미확인입니다. 서버 원본은 원래 실행 폴더명으로 별도 등록할 수 있습니다.
