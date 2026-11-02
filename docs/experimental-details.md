# IVIVE 실험 상세·구조 입력의 출처·해석의 한계

> 추가 첨부의 원본 코드·설정 기록을 반영한 최신 설명은 [서버 원본 추가 확인](server-source-update.md)을 보세요. 기존 미확인 사항 중 수동 CID·SMILES와 matched Morgan의 저장 설정은 새 근거로 확인되었습니다.

## 대상과 기록의 범위

Rat liver의 in vitro 발현으로 in vivo 발현을 예측합니다. 사용자가 확인한 데이터 규모는 140개 화합물, 3453개 probe이며, 조건 평균 데이터는 vitro 1675행 / vivo 2224행입니다. 개별 replicate의 분포를 직접 학습하는 실험은 아닙니다. seed 42의 화합물 분할은 내부 학습 89 / 검증 23 / 테스트 28개입니다.

이 문서는 저장소의 최종 코드·노트북·[9조건 결과](../results/meeting_comparison.csv)와 사용자가 서버 원본에서 확인하여 전달한 기록을 구분해 정리합니다. 서버의 `Untitled.ipynb`, `make_structure_features.py`, 원자료 및 최종 구조 CSV를 이 클라우드에서 직접 열어 재검증한 것은 아닙니다. 이 문서 작업에서 추가 학습이나 실험은 실행하지 않았습니다.

## SMILES 매핑: 사용자가 확인한 서버 원본

서버 `Untitled.ipynb`에서 확인한 절차는 다음과 같습니다.

1. `pubchempy.get_compounds(name, "name")`으로 화합물 이름을 검색합니다.
2. 검색 결과의 첫 번째 compound를 선택합니다.
3. 입력 이름, CID, `c.connectivity_smiles`, status를 `compound_pubchem_smiles.csv`에 저장합니다.
4. 추가 첨부에서 확인한 원본 구현은 `range(3)`으로 최대 3번 시도합니다.
5. `acetamidofluorene`, `naphthyl isothiocyanate`는 CID와 SMILES를 수동 입력하고 `OK_MANUAL`로 기록합니다.

첫 검색 결과가 실험에 사용된 물질과 동일한지 검증하는 자동 단계는 확인되지 않았습니다. 이름의 모호성, 이성질체·염 형태 등으로 잘못된 구조가 연결될 가능성이 있으므로, 첫 검색 결과나 `OK_MANUAL` 표기만으로 물질 동일성이 검증됐다고 볼 수 없습니다. 수동 입력한 CID·SMILES의 구체적인 값은 [추가 확인 문서](server-source-update.md)에 기록했습니다. 물질 동일성의 별도 검증 근거는 확인하지 못했습니다.

위 내용은 원본에 대한 확인 기록을 설명한 것이며 실행 가능한 매핑 코드를 재구성해 제공하지 않습니다. 추후 코드를 재구성한다면 원본 그대로가 아닌 **재구성본**임을 표시해야 합니다.

## Morgan 생성 코드와 최종 학습 입력의 구분

사용자가 확인한 서버 `make_structure_features.py`는 다음 절차를 사용합니다.

- `Chem.MolFromSmiles`로 SMILES를 파싱합니다.
- `LargestFragmentChooser`로 가장 큰 fragment를 선택합니다.
- `AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)`로 fingerprint를 생성합니다.
- 파싱 실패 값은 열별 중앙값으로 채우고 상수 열을 제거하여 `struct_feats/struct_morgan2048.csv`에 저장합니다. 상수 열 제거로 출력 비트 수가 2048보다 줄어들 수 있습니다.

이 스크립트에서 확인된 설정은 radius 2, 생성 시 2048비트입니다. 그러나 **이 스크립트의 생성 과정이 최종 학습 파일의 생성 과정과 같은지는 미확인**입니다. 최종 학습 입력은 `conditioned_expression/morgan_matched/morgan_binary.csv`이며, 사용자가 확인한 크기는 `compound` 포함 2049열입니다. 최종 노트북의 구조 입력 경로도 이 파일을 가리킵니다. 추가 첨부의 `morgan_matched/config.json` 저장 출력에는 **radius 2, size 2048, chirality false, RDKit 2026.03.6**이 기록되어 있습니다. 이는 최종 경로의 설정 기록 확인이며 CSV checksum이나 독립 재생성 검증은 아닙니다.

추가 첨부에서 matched용 `structure_features.py`와 RDKit 버전 기록을 확인했습니다. 최종 CSV까지의 전체 변환·실행 이력과 입력 파일 checksum은 여전히 미확인입니다. 중간 생성 스크립트의 중앙값 보간을 최종 입력에도 적용했다고 단정할 수 없습니다. 현재 학습 로더는 구조 입력의 유한한 이진 값(0/1)과 이름 중복을 검사하고, 화합물 이름으로 매칭한 뒤 누락 구조를 0과 presence flag로 처리합니다. 이 검사 자체가 PubChem 물질 동일성이나 fingerprint 생성 설정을 검증하지는 않습니다.

## 학습 방식

발현·구조·pathway encoder는 별도 사전학습 없이 각 모델의 최종 학습 목적함수로 함께 end-to-end 학습합니다. Probe–gene 매핑과 KEGG의 gene–pathway 소속으로 구성한 probe–pathway 연결 위치는 고정하고 연결별 가중치는 학습합니다. 별도의 pathway 정답이나 pathway loss는 없습니다.

MLP는 가중 발현 MSE, cGAN은 조건부 적대적 BCE, epsilon-DDPM은 noise 예측 MSE를 사용합니다. cGAN의 generator와 discriminator는 각각 자신의 context encoder를 학습하며 별도 발현 재구성 loss는 없습니다. 자세한 연산은 [코드 해설](code-walkthrough.md)을 참고하세요.

학습 target은 각 probe의 학습 target 평균을 뺀 값이며, 예측 후 같은 평균을 더합니다. target 분산으로 나누는 표준화는 하지 않습니다. 내부 학습에서는 학습 89개 화합물로 전처리를 적합하고 검증 23개로 epoch를 선택합니다. 이후 112개 화합물에 전처리를 다시 적합하고 모델을 초기화해 선택한 epoch만큼 재학습한 뒤 테스트 28개를 평가합니다.

## 현재 결과에서 말할 수 있는 범위

현재 9조건은 MLP·cGAN·epsilon-DDPM × pathway 없음 / 실제 KEGG sparse / 무작위 sparse이며 **모두 `structure=True`**입니다. 이 표만으로 구조 입력을 추가한 효과를 주장할 수 없습니다. 무작위 대조는 probe별·pathway별 연결 수를 보존한 mask 1개입니다.

MLP의 실제 KEGG sparse가 9조건 중 절대 RMSE 0.333559, 반응 RMSE 0.309485로 가장 낮았습니다. **이 분할과 단일 무작위 대조에서 오차가 낮았다**고 해석합니다. 생성모델의 pathway 효과는 제한적입니다. DDPM의 실제 KEGG는 pathway 없음보다 반응 오차가 소폭 낮지만 절대 오차는 높았고, cGAN에서는 두 오차 모두 개선되지 않았습니다. 모델 전반에 pathway 효과가 일관되게 입증됐다고 주장하지 않습니다.

단일 분할, 테스트 화합물 28개, 무작위 mask 1개에 기반한 탐색 결과입니다. 코드 안에서는 검증 데이터로 epoch를 선택했지만, 같은 테스트 결과를 보며 다음 설계를 발전시켰습니다. 따라서 완전히 독립적인 테스트 또는 확증적 검증이라고 주장하지 않습니다. 여러 조건·pair가 있어도 독립적인 테스트 화합물 수가 늘어나는 것은 아닙니다. 유의성이나 새로운 분할·화합물 집단에 대한 일반화도 현재 결과만으로 확정할 수 없습니다.

예측은 probe 수준이며 probe와 고유 유전자가 일대일 대응한다고 가정할 수 없습니다. KEGG 연결은 고정 외부 지식이고 학습된 가중치가 인과적 활성·억제를 입증하지 않습니다. 생성모델은 생성 평균의 발현 오차를 평가했으므로 생성 분포의 품질을 입증하는 결과도 아닙니다.

결과 CSV는 사용자가 제공한 서버 출력의 소수점 여섯 자리 요약이며 이 환경에서 재계산하지 않았습니다. 원자료·annotation snapshot·체크포인트가 저장소에 없으므로 새 클라우드에서 접근 가능하거나 연구 결과를 즉시 재현할 수 있다고 가정하지 않습니다. [결과 출처](../results/provenance.json)를 함께 참고하세요.

## 후속 계획: 아직 실행하지 않음

- 여러 화합물 분할에서의 평가
- 여러 무작위 sparse mask를 사용한 대조
- 화합물을 재표집 단위로 하는 bootstrap 신뢰구간(CI)
- 훈련 평균 반응 기준선과의 비교

위 항목은 후속 계획이며 현재 결과에 포함되지 않습니다. 현재 표의 변화 없음 반응 기준선(`predicted response=0`, RMSE 0.330116)은 훈련 평균 반응 기준선과 다릅니다.
