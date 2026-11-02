# Python 코드 해설: 교수님께 설명하는 순서

## 1. 테이블에서 tensor로

`first_structure_mlp.Data`는 NPY 발현 배열과 CSV 조건·pair 테이블을 읽습니다. 배열 행 번호와 condition_id, probe 수가 맞는지 검사합니다. pair 한 행은 source vitro 조건 하나와 target vivo 조건 하나입니다. 같은 target에 source가 여러 개 연결될 수 있습니다.

```python
source = pack['source'][data.source_ids[rows]]
condition = pack['conditions'][rows]
structure = pack['structure'][data.compound_indices[rows]]
```

각 tensor의 첫 차원은 batch입니다. source는 처음 3453개 열이 표준화한 probe 발현이고, 그 뒤에 평균 pathway score가 붙습니다. sparse encoder는 앞의 probe 값만 사용하며 평균 score를 읽지 않습니다. 조건별 replicate 평균을 사용하며 개별 replicate 분포를 학습하는 실험은 아닙니다.

## 2. 학습 데이터에만 전처리를 적합

`prepare`에서 source StandardScaler는 고유 학습 source 조건으로 적합합니다. dose는 log1p 후 시간과 함께 표준화하고 단위·등급은 OneHotEncoder로 변환합니다. 서로 다른 단위의 dose를 물리적으로 환산하는 처리는 아닙니다. target은 학습 target 평균만 빼고 분산은 나누지 않습니다. raw probe 공간 오차를 유지하고 예측 후 평균을 더합니다. 추가 log/RMA 정규화는 이 코드에서 하지 않습니다.

하나의 target이 여러 pair로 반복되므로 학습 loss에 target별 pair 수의 역수를 가중치로 줍니다.

## 3. Sparse encoder는 실제로 어떤 연산인가

`pathway_pattern_mlp.PatternEncoder.__init__`:

```python
edges = np.argwhere(np.asarray(membership) > 0)
self.edge_weights = nn.Parameter(weights)  # [연결 수, 4]
```

3453×321 membership에서 실제 소속이 있는 위치만 추출합니다. 이번 snapshot에는 19,254개 연결이 있습니다. `register_buffer`로 probe/pathway index를 저장해 학습 parameter와 구분합니다. 연결별 가중치 4개는 optimizer가 갱신합니다. 실제 활성·억제 관계의 부호를 annotation에서 준 것은 아니며, 가중치 부호는 예측 손실로 배웁니다.

`pathway_units`:

```python
values = genes[:, self.probe_index]       # [batch, 연결 수]
weighted = values[:, :, None] * self.edge_weights[None, :, :]
aggregated = genes.new_zeros(batch, self.n_pathways, self.units)
aggregated = aggregated.index_add(1, self.pathway_index, weighted)
return F.silu(aggregated + self.bias)
```

`None`은 broadcasting용 차원을 추가합니다. 연결별 probe 발현에 서로 다른 가중치를 곱하고 `index_add`가 같은 pathway로 들어오는 항들을 더합니다. 노드 j,k는 `SiLU(sum_i x_i*w_ijk + b_jk)`입니다. pathway마다 4개 특징을 만들고 펼쳐 32차원으로 줄입니다. 고정 평균처럼 같은 가중치를 쓰지 않으므로 반대 방향 패턴을 구별할 수 있지만 압축이므로 완전한 정보 보존은 아닙니다. 같은 유전자에 여러 probe가 있으면 각각의 연결을 사용합니다.

`none`은 pathway branch 출력이 0입니다. `mean`은 개별 probe 대신 pathway 평균 score를 반복 입력합니다. 동일한 명목상 parameter 수라도 none branch는 비활성이고 mean의 유효 자유도도 sparse와 다릅니다.

## 4. 세 입력 표현을 합치는 방법

`PatternMLP.forward`는 발현 3453→128→64, Morgan+presence→64→32, pathway→32를 계산하고 실험 조건과 `torch.cat(..., dim=1)`로 이어 붙입니다. decoder는 결합 입력→64→3453을 예측합니다. 구조 유사도 그래프·타 약물 반응·사후 유전자 보정은 없습니다.

## 5. MLP / cGAN / DDPM 차이

MLP는 weighted MSE로 target 발현을 직접 예측합니다.

cGAN은 `Generator`가 context와 16차원 정규 노이즈를 받아 발현을 생성합니다. `Discriminator`는 같은 조건과 실제/생성 발현을 보고 logit을 출력합니다. BCEWithLogitsLoss로 판별기를 실제=1/생성=0으로 학습하고 생성기는 생성=1을 목표로 학습합니다. 재구성 MSE를 추가하지 않았습니다. Generator와 discriminator 각각의 context encoder도 학습합니다.

DDPM은 학습 중 t와 epsilon을 뽑아 `x_t=sqrt(alpha_bar_t)*x0+sqrt(1-alpha_bar_t)*epsilon`을 만들고 epsilon MSE를 최소화합니다. 기존 구현의 예측기는 `epsilon_hat = noisy + correction(context, noisy, time_embedding)`이며 마지막 correction layer를 0으로 초기화합니다. 따라서 원 논문의 U-Net 구조를 그대로 쓴 것은 아닙니다. 100단계에 맞춰 rescale한 linear beta와 전체 stochastic DDPM 역과정을 사용합니다. v 예측·DDIM·clipping·Gaussian 보정은 없습니다.

학습은 AdamW(lr=.001, weight_decay=.01), batch 128, gradient norm clipping 1, 최대 30 epoch입니다. 1,5,10,...,30 epoch에 내부 검증합니다. 생성모델은 pair별 4개 생성 평균으로 예측하므로 이 결과는 분포 품질 비교가 아닙니다.

## 6. epoch 선택·재학습·평가

seed 42로 compound를 나눠 내부 학습89 / 검증23 / 외부 테스트28을 구성합니다. 같은 화합물의 조건들이 분할을 넘어가지 않습니다. 내부 검증에서 source별 예측을 target별 평균한 절대 MSE로 epoch를 선택합니다. 전처리를 다시 112개 화합물에 적합하고 모델도 초기화해 선택 epoch만큼 재학습합니다. 테스트는 최종 결과 계산에 사용되지만 이 분할을 여러 설계에서 관찰했기 때문에 독립 confirmatory test는 아닙니다.

`predict`는 여러 source에서 같은 target을 예측한 값을 평균합니다. `pair_rmse`는 평균 전 pair 단위 진단이며 target 평균 RMSE와 다릅니다.

`evaluate`의 절대 RMSE는 target 발현 오차입니다. 반응 RMSE는 아래 변화량의 차이입니다.

```python
observed = actual_treated - actual_control
predicted = predicted_treated - predicted_control
```

실제 대조군은 정답 반응 계산에만 사용하며 예측 반응에 빼지 않습니다. target 조건을 동일하게 가중한 값이고 약물별 macro-average와 동일하다고 보장하지 않습니다. 변화 없음 기준은 predicted response=0이며, 절대 발현을 0이라고 예측한다는 뜻이 아닙니다.

## 7. 무작위 대조군

`shuffle_membership`는 (probe A, pathway B), (probe C, pathway D)를 (A,D),(C,B)로 교환합니다. 중복 연결은 허용하지 않아 probe별·pathway별 연결 수가 유지됩니다. 원래 edge 수×10번 성공한 교환 후 mask를 저장합니다. seed는 고정이며 성능으로 선택하지 않습니다. 원 KEGG edge가 약16.6% 남았고 균등 표본 추출을 보장하지 않습니다. 이는 논문 재현 항목이 아닌 우리가 추가한 대조 실험입니다.

## 발표 문장

“구조를 제공한 IVIVE 예측에서 pathway의 표현 방식을 비교했습니다. 고정 평균 대신 probe–pathway 소속 연결의 가중치를 학습하는 encoder가 MLP에서 절대 발현·반응 오차를 함께 줄였습니다. 연결 수를 보존한 무작위 대조군보다도 좋았으나, cGAN과 DDPM까지 일관된 개선은 아니었습니다.”
