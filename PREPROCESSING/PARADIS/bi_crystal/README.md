# ParaDiS Bicrystal Frank–Read Source Network Generator

이 MATLAB 스크립트는 ParaDiS에서 사용할 수 있는 **bicrystal 전위 네트워크 초기 조건**을 자동으로 생성합니다. 생성된 `.data`와 `.ctrl` 파일은 Frank–Read 소스 기반의 전위들을 포함하고 있으며, 두 결정립 사이에 `y = 0` 평면을 경계로 한 bicrystal 구조를 구성합니다.

---

## 📌 기본 특징

- 모든 기하학적 단위는 **Burgers 벡터 크기 `b`**를 기준으로 생성됨
- 결정립 회전을 위한 **Bunge Euler 각도** 입력 지원
- `y = 0` 평면에 수직한 grain boundary를 자동 삽입
- 설정된 전위 밀도에 도달할 때까지 **Frank–Read 소스를 무작위 배치**
- `.ctrl` 및 `.data` 파일 모두 자동 생성
- 시각화를 위한 3D 플로팅 기능 포함 (MATLAB 내에서 확인 가능)

---
## 🔧 주요 설정 항목 (`config`)

| 파라미터 | 설명 |
|----------|------|
| `config.maxstep` | `1`이면 DDM 사용, `>1`이면 ParaDiS 단독 실행 |
| `config.boxSize_b` | 전위 네트워크 도메인 크기 (단위: b) |
| `config.targetDensity` | 목표 전위 밀도 [m^-2] |
| `config.frs_length_alpha` | 평균 FRS 길이 조절 파라미터 |
| `config.crystalStructure` | 결정 구조 ('FCC', 'BCC') |
| `config.R` | **Grain 2의 결정 방향 회전을 정의하는 3×3 회전 행렬**입니다. 이 회전은 y>0 영역의 grain 의 방위를 설정합니다. 사용자는 Bunge Euler 각(`bunge_phi1`, `bunge_PHI`, `bunge_phi2`)을 DEGREE 단위로 지정한 후, `deg2rad()`로 라디안 변환하여 `bunge_euler_to_matrix()`를 통해 계산합니다. 회전 행렬 `R`은 grain 2의 Burgers 벡터와 면 법선 벡터를 기준 결정립 대비 회전시켜 bicrystal misorientation을 설정하는 데 사용됩니다. |
| `config.mobilityLaw` | ParaDiS mobility law 이름 |
| `config.burgMag` | 버거스 벡터 크기 (m) |
| `config.output_folder` | 생성된 파일 저장 위치 |

---

## 🖥️ 사용법

1. `main.m` 파일 (또는 포함된 스크립트)을 MATLAB에서 실행하세요.
2. 설정한 `config.output_folder` 경로에 `.data`와 `.ctrl` 파일이 생성됩니다.
3. ParaDiS 실행 시 이 파일들을 사용하여 bicrystal 시뮬레이션을 시작할 수 있습니다.

---

## 🎨 시각화

- `config.enable_plotting = true`로 설정하면 전위 네트워크가 3D로 시각화됩니다.
- 서로 다른 Burgers 벡터 방향은 색으로 구분됩니다.
- grain boundary(`y = 0`)는 반투명 평면으로 표시됩니다.

---

## 📂 생성 파일 예시

```
output/
├── Bicrystal_FCC_run_1.ctrl   ← ParaDiS 제어 파일
└── Bicrystal_FCC_run_1.data   ← 초기 전위 구조 파일
```

---

## ⚠️ 유의사항

- Frank–Read 소스의 길이는 자동 난수 생성되며, 도메인을 넘지 않도록 제한됩니다.
- 전위의 방향은 slip system에서 무작위로 선택됩니다.
- 전위가 grain boundary를 넘어가는 경우는 제거되어 다시 생성됩니다.

---

## 📞 문의

작성자: Gyujang Sim  
소속: MAMEL Lab, SNU  
문의: gyujang95@snu.ac.kr
